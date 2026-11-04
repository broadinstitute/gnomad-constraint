"""Test whether the histone genes' extra synonymous variants and DNMs copy a paralog.

Gene conversion, or reads from a paralog mapped to the wrong copy, puts variants where
a gene's paralogs differ from it, with the alt allele equal to a paralog's base. A
higher point-mutation rate has no reason to prefer those alleles.

For every histone transcript in ``syn_z_by_histone_group.tsv``:

1. Its coding sequence: the transcript's coding exons from Ensembl (GRCh38), read from
   hg38.2bit.
2. Its protein aligned to every other member of its family (H1, H2A, H2B, H3, H4), with
   codons mapped through the alignment, to get the paralog bases at each coding position.
3. Its possible synonymous SNVs, flagged where a paralog has a different base
   (``paralog_differs``), where the alt is any paralog's base (``alt_in_paralog``) and
   where it is the base of the most similar paralog (``alt_in_closest``).
4. Which are in gnomAD v4.1 exomes (observed = PASS with AC > 0; filtered calls kept
   for the filtering check) and which are synonymous DNMs (every cohort in the
   pipeline's DNM Table).

The nearby non-histone genes of ``syn_z_by_histone_group.tsv`` (no synonymous excess)
are the control genes: their share of possible variants seen in gnomAD, per
(trinucleotide context, alt), gives the histone variants' expected count.

Tests:

- gnomAD: logistic regression of being observed on ``paralog_differs`` and one paralog
  flag, with one term per (trinucleotide context, alt); and, at positions where a
  paralog differs, a Mantel-Haenszel odds ratio of the paralog's alt against the other
  synonymous alts of the same mutation type. Each is repeated with Roulette's own rate
  for the variant as a covariate, on the variants Roulette covers, in case paralog
  differences sit at sites more mutable than their trinucleotide suggests.
- DNMs: the share of DNMs whose alt is a paralog's base, against the share expected from
  each possible variant's mutation rate (Roulette's average for its trinucleotide
  context and alt, or its own rate for the variant).
- Observed against expected by allele class (copies a paralog or not): gnomAD against
  the control genes' rates, and DNMs against Roulette's rate for each variant scaled to
  the genome-wide DNM count.
- Frequency and filtering of the gnomAD variants that copy a paralog against those that
  don't. Reads mapped from a paralog tend to give common variants that fail gnomAD's
  filters (for example too many heterozygotes, or a low VQSR score).
"""

import argparse
import json
import logging
import os
import subprocess
import time
import urllib.error
import urllib.request
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from roulette_combo import BASE_CODE, LETTERS, read_chromosome
from scipy.stats import binomtest, chi2
from statsmodels.stats.contingency_tables import StratifiedTable

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

ENSEMBL = "https://rest.ensembl.org/lookup/id"
VCF = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/"
    "gnomad.exomes.v4.1.sites.{chrom}.vcf.bgz"
)
FAMILIES = [
    ("H1", r"^H1-"),
    ("H2A", r"^H2A"),
    ("H2B", r"^H2B"),
    ("H3", r"^H3"),
    ("H4", r"^H4"),
]
COMPLEMENT = str.maketrans("ACGT", "TGCA")
CONTROL = "non-histone, same chr6 region"
"""Group of ``syn_z_by_histone_group.tsv`` used as control genes: nearby, no excess."""


def revcomp(s: str) -> str:
    """Get the reverse complement of a DNA string."""
    return s.translate(COMPLEMENT)[::-1]


def ensembl_cds(transcripts: List[str]) -> Dict[str, dict]:
    """
    Get each transcript's chromosome, strand and coding segments from Ensembl.

    :param transcripts: Ensembl transcript IDs without versions.
    :return: {transcript: {"chrom", "strand", "segments": [(start, end), ...]}}, with
        1-based inclusive genomic segments in genomic order.
    """
    req = urllib.request.Request(
        ENSEMBL,
        data=json.dumps({"ids": transcripts, "expand": 1}).encode(),
        headers={"Content-Type": "application/json", "Accept": "application/json"},
    )
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req) as r:
                res = json.load(r)
            break
        except urllib.error.HTTPError as e:
            if attempt == 4 or e.code not in (429, 500, 502, 503, 504):
                raise
            logger.info("Ensembl returned %d; retrying.", e.code)
            time.sleep(10 * (attempt + 1))
    out = {}
    for t, v in res.items():
        if not v or "Translation" not in v:
            logger.warning("%s: no translation in Ensembl", t)
            continue
        lo, hi = v["Translation"]["start"], v["Translation"]["end"]
        segs = sorted(
            (max(e["start"], lo), min(e["end"], hi))
            for e in v["Exon"]
            if e["end"] >= lo and e["start"] <= hi
        )
        out[t] = {
            "chrom": f"chr{v['seq_region_name']}",
            "strand": v["strand"],
            "segments": segs,
        }
    return out


def coding_sequence(seq: np.ndarray, info: dict) -> tuple:
    """
    Read a transcript's coding sequence and the genomic position of each base.

    :param seq: Chromosome bases from ``read_chromosome``.
    :param info: Entry of ``ensembl_cds``.
    :return: (coding sequence 5'->3', genomic positions in the same order).
    """
    pos = np.concatenate([np.arange(s, e + 1) for s, e in info["segments"]])
    cds = seq[pos - 1].tobytes().decode()
    if info["strand"] < 0:
        return revcomp(cds), pos[::-1]
    return cds, pos


def codon_maps(prots: Dict[str, str]) -> Dict[tuple, Dict[int, int]]:
    """
    Align every ordered pair of proteins in a family and map their codons.

    :param prots: {transcript: protein without the stop}.
    :return: {(transcript, paralog): {codon index: paralog codon index}}.
    """
    aligner = Align.PairwiseAligner(
        mode="global",
        substitution_matrix=substitution_matrices.load("BLOSUM62"),
        open_gap_score=-10,
        extend_gap_score=-0.5,
    )
    out = {}
    for a, pa in prots.items():
        for b, pb in prots.items():
            if a == b:
                continue
            aln = aligner.align(pa, pb)[0]
            m = {}
            for (s1, e1), (s2, _) in zip(*aln.aligned):
                m.update({s1 + k: s2 + k for k in range(e1 - s1)})
            out[(a, b)] = m
    return out


def possible_synonymous(
    t: str,
    cds: str,
    pos: np.ndarray,
    info: dict,
    seq: np.ndarray,
    paralogs: Dict[str, str],
    maps: Dict[tuple, Dict[int, int]],
    closest: Optional[str],
) -> pd.DataFrame:
    """
    List a transcript's possible synonymous SNVs with the paralog flags.

    :param t: Transcript.
    :param cds: Its coding sequence.
    :param pos: Genomic position of each coding base.
    :param info: Its ``ensembl_cds`` entry.
    :param seq: Chromosome bases.
    :param paralogs: {paralog transcript: coding sequence} for the rest of its family.
    :param maps: Output of ``codon_maps``.
    :param closest: The paralog with the highest coding identity, or None for a gene
        without paralogs here.
    :return: One row per possible synonymous SNV.
    """
    rows = []
    n_codons = len(cds) // 3 - 1
    for i in range(n_codons):
        codon = cds[3 * i : 3 * i + 3]
        aa = str(Seq(codon).translate())
        for j in range(3):
            k = 3 * i + j
            par = [
                p_cds[3 * maps[(t, p)][i] + j]
                for p, p_cds in paralogs.items()
                if i in maps[(t, p)]
            ]
            near = (
                paralogs[closest][3 * maps[(t, closest)][i] + j]
                if closest and i in maps[(t, closest)]
                else None
            )
            for alt in "ACGT":
                if alt == codon[j]:
                    continue
                new = codon[:j] + alt + codon[j + 1 :]
                if str(Seq(new).translate()) != aa:
                    continue
                g = int(pos[k])
                tri = seq[g - 2 : g + 1].tobytes().decode()
                g_ref, g_alt = (
                    (codon[j], alt)
                    if info["strand"] > 0
                    else (codon[j].translate(COMPLEMENT), alt.translate(COMPLEMENT))
                )
                if g_ref in "GT":
                    ctx, c_alt = revcomp(tri), g_alt.translate(COMPLEMENT)
                else:
                    ctx, c_alt = tri, g_alt
                rows.append(
                    {
                        "transcript": t,
                        "chrom": info["chrom"],
                        "pos": g,
                        "ref": g_ref,
                        "alt": g_alt,
                        "cds_pos": k + 1,
                        "context": ctx,
                        "context_alt": c_alt,
                        "n_paralogs": len(par),
                        "paralog_differs": any(b != codon[j] for b in par),
                        "alt_in_paralog": alt in par,
                        "paralogs_with_alt": par.count(alt),
                        "alt_in_closest": alt == near,
                    }
                )
    return pd.DataFrame(rows)


def roulette_context_mu(path: str, twobit: str) -> pd.DataFrame:
    """
    Get Roulette's average rate per (strand-collapsed trinucleotide context, alt).

    :param path: Roulette ``all_hq_synonymous_variants.tsv.gz``.
    :param twobit: hg38.2bit path.
    :return: Rows with ``context``, ``context_alt`` and ``mu``.
    """
    d = pd.read_csv(
        path,
        sep="\t",
        usecols=["CHROM", "POS", "REF", "ALT", "mu"],
        dtype={"CHROM": str},
    )
    d = d[d.REF.str.len().eq(1) & d.ALT.str.len().eq(1)]
    out = []
    for chrom, g in d.groupby("CHROM", sort=False):
        seq = read_chromosome(twobit, f"chr{chrom}")
        p = g.POS.to_numpy()
        codes = BASE_CODE[seq[(p - 1)[:, None] + np.arange(-1, 2)]]
        alt = BASE_CODE[np.frombuffer("".join(g.ALT).encode(), dtype=np.uint8)]
        ok = (codes < 4).all(axis=1)
        codes, alt, mu = codes[ok], alt[ok], g.mu.to_numpy()[ok]
        flip = codes[:, 1] >= 2
        codes[flip] = (3 - codes[flip])[:, ::-1]
        alt[flip] = 3 - alt[flip]
        out.append(
            pd.DataFrame(
                {
                    "context": LETTERS[codes]
                    .astype(np.uint8)
                    .copy()
                    .view("S3")
                    .ravel(),
                    "context_alt": LETTERS[alt].astype(np.uint8).view("S1"),
                    "mu": mu,
                }
            )
        )
    r = pd.concat(out)
    r["context"] = r.context.str.decode("ascii")
    r["context_alt"] = r.context_alt.str.decode("ascii")
    return r.groupby(["context", "context_alt"], as_index=False).mu.mean()


def gnomad_variants(sites: pd.DataFrame, bcftools: str, tmp: str) -> pd.DataFrame:
    """
    Get the gnomAD v4.1 exome variants over the coding sites, with their filters.

    :param sites: Possible variants with ``chrom`` and ``pos``.
    :param bcftools: bcftools path.
    :param tmp: Path for the regions file.
    :return: Rows with ``chrom``, ``pos``, ``ref``, ``alt``, ``filter``, ``AC`` and
        ``AN``.
    """
    out = []
    for chrom, g in sites.groupby("chrom"):
        # Merge the coding positions into regions, so only exons are read remotely.
        p = np.unique(g.pos.to_numpy())
        brk = np.flatnonzero(np.diff(p) > 1000)
        regions = pd.DataFrame(
            {
                "c": chrom,
                "s": p[np.r_[0, brk + 1]],
                "e": p[np.r_[brk, len(p) - 1]],
            }
        )
        fmt = "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AC\t%INFO/AN\n"
        rows = []
        # Short remote reads in small batches, retried: long ones often drop mid-read.
        for b in range(0, len(regions), 20):
            regions.iloc[b : b + 20].to_csv(tmp, sep="\t", header=False, index=False)
            for attempt in range(5):
                # bcftools saves the remote index in its working directory.
                q = subprocess.run(
                    [bcftools, "query", "-R", tmp, "-f", fmt, VCF.format(chrom=chrom)],
                    capture_output=True,
                    cwd=os.path.dirname(tmp),
                )
                if q.returncode == 0:
                    break
                logger.info("bcftools failed on %s; retrying.", chrom)
                time.sleep(5)
            q.check_returncode()
            rows += [line.split("\t") for line in q.stdout.decode().splitlines()]
        out.append(
            pd.DataFrame(
                rows, columns=["chrom", "pos", "ref", "alt", "filter", "AC", "AN"]
            )
        )
        logger.info("%s: %d gnomAD variants", chrom, len(rows))
    v = pd.concat(out)
    v["pos"] = v.pos.astype(int)
    for c in ("AC", "AN"):
        v[c] = pd.to_numeric(v[c], errors="coerce")
    return v.drop_duplicates(["chrom", "pos", "ref", "alt"])


def roulette_site_mu(path: str, sites: pd.DataFrame, dnms: pd.DataFrame) -> tuple:
    """
    Get Roulette's own rate for each possible variant it covers, and the DNM scale.

    :param path: Roulette ``all_hq_synonymous_variants.tsv.gz``.
    :param sites: Possible variants with ``chrom`` and ``pos``.
    :param dnms: Synonymous DNMs with ``chrom``, ``pos``, ``ref`` and ``alt``.
    :return: (rows with ``chrom``, ``pos``, ``ref``, ``alt`` and ``site_mu``; DNMs per
        unit of Roulette rate over every synonymous variant Roulette covers).
    """
    keys = set(zip(sites.chrom.str.replace("chr", ""), sites.pos))
    dkeys = set(zip(dnms.chrom.str.replace("chr", ""), dnms.pos, dnms.ref, dnms.alt))
    out, total_mu, n_dnm = [], 0.0, 0
    for chunk in pd.read_csv(
        path,
        sep="\t",
        usecols=["CHROM", "POS", "REF", "ALT", "mu"],
        dtype={"CHROM": str},
        chunksize=2_000_000,
    ):
        total_mu += chunk.mu.sum()
        n_dnm += sum(
            k in dkeys for k in zip(chunk.CHROM, chunk.POS, chunk.REF, chunk.ALT)
        )
        out.append(chunk[[k in keys for k in zip(chunk.CHROM, chunk.POS)]])
    r = pd.concat(out)
    rates = pd.DataFrame(
        {
            "chrom": "chr" + r.CHROM,
            "pos": r.POS,
            "ref": r.REF,
            "alt": r.ALT,
            "site_mu": r.mu,
        }
    ).drop_duplicates(["chrom", "pos", "ref", "alt"])
    logger.info("%d synonymous DNMs at variants Roulette covers.", n_dnm)
    return rates, n_dnm / total_mu


def observed_vs_expected(sites: pd.DataFrame, dnm_per_mu: float) -> pd.DataFrame:
    """
    Compare histone variants with what their mutation type predicts, by allele class.

    gnomAD: expected is the share of possible variants of the same (trinucleotide
    context, alt) seen in the control genes, nearby non-histone genes with no excess.
    DNMs: expected is the variant's Roulette rate times the genome-wide DNMs per unit
    rate, over the variants Roulette covers.

    :param sites: Possible variants with ``histone``, the paralog flags, ``observed``,
        ``dnm``, ``stratum`` and ``site_mu``.
    :param dnm_per_mu: Output of ``roulette_site_mu``.
    :return: Observed, expected and their ratio with a 95% interval, per source and
        allele class.
    """
    ctrl = sites[~sites.histone]
    p_stratum = ctrl.groupby("stratum")["observed"].mean()
    h = sites[sites.histone].assign(
        exp_gnomad=lambda d: d.stratum.map(p_stratum).fillna(ctrl.observed.mean()),
        exp_dnm=lambda d: d.site_mu * dnm_per_mu,
    )
    classes = [
        ("all", np.ones(len(h), bool)),
        ("copies a paralog", h.alt_in_paralog.to_numpy()),
        ("copies no paralog", ~h.alt_in_paralog.to_numpy()),
        ("copies the closest paralog", h.alt_in_closest.to_numpy()),
    ]
    rows = []
    for name, m in classes:
        g = h[m]
        for source, obs, exp in (
            ("gnomAD v4.1 exomes", g.observed.sum(), g.exp_gnomad.sum()),
            (
                "DNMs",
                g[g.site_mu.notna()].dnm.sum(),
                g[g.site_mu.notna()].exp_dnm.sum(),
            ),
        ):
            # Exact (Garwood) Poisson interval.
            lo = chi2.ppf(0.025, 2 * obs) / 2 if obs else 0.0
            hi = chi2.ppf(0.975, 2 * obs + 2) / 2
            rows.append(
                {
                    "source": source,
                    "alleles": name,
                    "possible": int(m.sum()),
                    "observed": int(obs),
                    "expected": exp,
                    "ratio": obs / exp,
                    "ratio_lo": lo / exp,
                    "ratio_hi": hi / exp,
                }
            )
    return pd.DataFrame(rows)


def test_flag(d: pd.DataFrame, flag: str, site_mu: bool) -> dict:
    """
    Test whether variants carrying a paralog flag are seen more often than expected.

    :param d: Possible variants.
    :param flag: Flag column, e.g. ``alt_in_paralog``.
    :param site_mu: Whether to add log Roulette site rate as a covariate.
    :return: Odds ratios, the excess of flagged variants, the DNM test and the
        frequency and filtering comparisons.
    """
    d = d.assign(
        paralog_differs=d.paralog_differs.astype(int),
        flagged=d[flag].astype(int),
        observed=d.observed.astype(int),
    )
    formula = "observed ~ C(stratum) + paralog_differs + flagged"
    if site_mu:
        d = d.assign(log_mu=np.log(d.site_mu))
        formula += " + log_mu"
    # Strata where every variant is seen, or none is, carry no information on the flag.
    var = d.groupby("stratum")["observed"].transform("nunique") > 1
    fit = smf.glm(formula, d[var], family=sm.families.Binomial()).fit()
    v = d[d.paralog_differs == 1]
    tables = [
        pd.crosstab(g.flagged, g.observed)
        .reindex(index=[1, 0], columns=[1, 0], fill_value=0)
        .to_numpy()
        for _, g in v.groupby("stratum")
        if g.flagged.nunique() == 2
    ]
    mh = StratifiedTable(tables)
    # Observed flagged variants expected if the flag had no effect.
    dv = d[var]
    exp_flagged = (
        fit.predict(dv.assign(flagged=0))[dv.flagged == 1].sum()
        + d[~var & (d.flagged == 1)].observed.sum()
    )
    obs_v = d[d.observed == 1]
    n_dnm = int(d.dnm.sum())
    k_dnm = int((d.dnm & (d.flagged == 1)).sum())
    w = d.site_mu if site_mu else d.mu
    share = (w * d.flagged).sum() / w.sum()
    called = d[d.in_vcf]
    return {
        "transcripts": d.transcript.nunique(),
        "possible": len(d),
        "possible_flagged": int(d.flagged.sum()),
        "observed": int(d.observed.sum()),
        "observed_flagged": int(d.observed[d.flagged == 1].sum()),
        "expected_flagged_if_no_effect": exp_flagged,
        "or_paralog_differs": np.exp(fit.params["paralog_differs"]),
        "or_flagged": np.exp(fit.params["flagged"]),
        "or_flagged_lo": np.exp(fit.conf_int().loc["flagged", 0]),
        "or_flagged_hi": np.exp(fit.conf_int().loc["flagged", 1]),
        "p_flagged": fit.pvalues["flagged"],
        "mh_or_at_differing_positions": mh.oddsratio_pooled,
        "mh_or_lo": mh.oddsratio_pooled_confint()[0],
        "mh_or_hi": mh.oddsratio_pooled_confint()[1],
        "mh_p": mh.test_null_odds().pvalue,
        "dnms": n_dnm,
        "dnms_flagged": k_dnm,
        "dnm_expected_share_flagged": share,
        "dnm_p": binomtest(k_dnm, n_dnm, share).pvalue if n_dnm else np.nan,
        "median_ac_flagged": obs_v[obs_v.flagged == 1].AC.median(),
        "median_ac_other": obs_v[obs_v.flagged == 0].AC.median(),
        "share_ac10_flagged": (obs_v[obs_v.flagged == 1].AC >= 10).mean(),
        "share_ac10_other": (obs_v[obs_v.flagged == 0].AC >= 10).mean(),
        "share_filtered_flagged": (
            called[called.flagged == 1]["filter"] != "PASS"
        ).mean(),
        "share_filtered_other": (
            called[called.flagged == 0]["filter"] != "PASS"
        ).mean(),
    }


def main(args: argparse.Namespace) -> None:
    """Build the sites, run the tests and write the tables."""
    hz = pd.read_csv(args.histone_tsv, sep="\t")
    hz = hz[hz.group.str.startswith("histone") | (hz.group == CONTROL)].copy()
    hz["family"] = None
    for fam, pat in FAMILIES:
        hz.loc[
            hz.group.str.startswith("histone")
            & hz.gene.str.match(pat)
            & hz.family.isna(),
            "family",
        ] = fam
    hz["family"] = hz.family.fillna("control")
    info = ensembl_cds(hz.transcript.tolist())

    seqs, cds, pos = {}, {}, {}
    for t, v in info.items():
        if v["chrom"] not in seqs:
            seqs[v["chrom"]] = read_chromosome(args.twobit, v["chrom"])
        c, p = coding_sequence(seqs[v["chrom"]], v)
        prot = str(Seq(c).translate())
        if len(c) % 3 or not prot.endswith("*") or "*" in prot[:-1]:
            logger.warning("%s: coding sequence doesn't translate cleanly; dropped", t)
            continue
        cds[t], pos[t] = c, p
    hz = hz[hz.transcript.isin(cds)]

    rows, ident = [], []
    for t in hz[hz.family == "control"].transcript:
        rows.append(
            possible_synonymous(
                t, cds[t], pos[t], info[t], seqs[info[t]["chrom"]], {}, {}, None
            ).assign(family="control")
        )
    for fam, g in hz[hz.family != "control"].groupby("family"):
        prots = {t: str(Seq(cds[t]).translate())[:-1] for t in g.transcript}
        maps = codon_maps(prots)
        for t in g.transcript:
            par = {p: cds[p] for p in g.transcript if p != t}
            idt = {}
            for p in par:
                m = maps[(t, p)]
                idt[p] = np.mean(
                    [
                        cds[t][3 * i + j] == cds[p][3 * m[i] + j]
                        for i in m
                        for j in range(3)
                    ]
                )
                ident.append(
                    {
                        "family": fam,
                        "transcript": t,
                        "paralog": p,
                        "nt_identity": idt[p],
                    }
                )
            closest = max(idt, key=idt.get)
            sites = possible_synonymous(
                t, cds[t], pos[t], info[t], seqs[info[t]["chrom"]], par, maps, closest
            )
            rows.append(sites.assign(family=fam, closest_identity=idt[closest]))
    sites = pd.concat(rows, ignore_index=True)
    sites["histone"] = sites.family != "control"
    sites = sites.merge(hz[["transcript", "gene", "group"]], on="transcript")
    ident = pd.DataFrame(ident)
    logger.info(
        "%d possible synonymous SNVs in %d histone and %d control transcripts; nt"
        " identity to paralogs by family:\n%s",
        len(sites),
        sites[sites.histone].transcript.nunique(),
        sites[~sites.histone].transcript.nunique(),
        ident.groupby("family").nt_identity.describe().round(3).to_string(),
    )

    if args.roulette_mu:
        mu = pd.read_csv(args.roulette_mu, sep="\t")
    else:
        mu = roulette_context_mu(args.roulette, args.twobit)
        mu.to_csv(
            f"{args.out_prefix}.roulette_mu_by_context.tsv", sep="\t", index=False
        )
    sites = sites.merge(mu, on=["context", "context_alt"], how="left")

    if os.path.exists(args.gnomad_variants):
        gv = pd.read_csv(args.gnomad_variants, sep="\t")
    else:
        gv = gnomad_variants(
            sites,
            args.bcftools,
            os.path.join(
                os.path.dirname(os.path.abspath(args.gnomad_variants)), "regions.tsv"
            ),
        )
        gv.to_csv(args.gnomad_variants, sep="\t", index=False)
    sites = sites.merge(
        gv, on=["chrom", "pos", "ref", "alt"], how="left", validate="many_to_one"
    )
    sites["in_vcf"] = sites["filter"].notna()
    sites["observed"] = (sites["filter"] == "PASS") & (sites.AC > 0)
    dn = pd.read_csv(args.dnms, sep="\t", compression="gzip")
    dn["chrom"] = dn.locus.str.split(":").str[0]
    dn["pos"] = dn.locus.str.split(":").str[1].astype(int)
    al = dn.alleles.apply(json.loads)
    dn["ref"], dn["alt"] = al.str[0], al.str[1]
    dn = dn[["chrom", "pos", "ref", "alt"]].drop_duplicates().assign(dnm=True)
    sites = sites.merge(dn, on=["chrom", "pos", "ref", "alt"], how="left")
    sites["dnm"] = sites.dnm.notna()
    sites["stratum"] = sites.context + ">" + sites.context_alt
    sites[[c for c in sites.columns if c != "dnm"]].to_csv(
        f"{args.out_prefix}.sites.tsv.gz", sep="\t", index=False
    )

    rs, dnm_per_mu = roulette_site_mu(args.roulette, sites, dn)
    sites = sites.merge(rs, on=["chrom", "pos", "ref", "alt"], how="left")
    sites.drop(columns="dnm").to_csv(
        f"{args.out_prefix}.sites.tsv.gz", sep="\t", index=False
    )
    logger.info(
        "Roulette site rates for %.0f%% of the possible variants.",
        100 * sites.site_mu.notna().mean(),
    )
    excess = observed_vs_expected(sites, dnm_per_mu)
    excess.to_csv(f"{args.out_prefix}.excess.tsv", sep="\t", index=False)
    logger.info(
        "Observed against expected:\n%s", excess.round(3).to_string(index=False)
    )

    results = []
    hs = sites[sites.histone]
    for name, d in [
        ("all histone genes", hs),
        *hs.groupby("family"),
    ]:
        for flag in ("alt_in_paralog", "alt_in_closest"):
            results.append({"set": name, "flag": flag, **test_flag(d, flag, False)})
            results.append(
                {
                    "set": name,
                    "flag": f"{flag}, with Roulette site rate",
                    **test_flag(d[d.site_mu.notna()], flag, True),
                }
            )
    res = pd.DataFrame(results)
    res.to_csv(f"{args.out_prefix}.tests.tsv", sep="\t", index=False)
    ident.to_csv(f"{args.out_prefix}.paralog_identity.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    logger.info("Tests:\n%s", res.T.to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--histone-tsv", required=True)
    parser.add_argument("--twobit", required=True)
    parser.add_argument(
        "--roulette",
        required=True,
        help="Roulette all_hq_synonymous_variants.tsv.gz, for each variant's own rate"
        " (and the averages, without --roulette-mu).",
    )
    parser.add_argument(
        "--roulette-mu",
        help="Roulette's average rate per context and alt from an earlier run, instead"
        " of --roulette.",
    )
    parser.add_argument("--dnms", required=True, help="Pipeline syn_dnm.tsv.bgz.")
    parser.add_argument("--bcftools", required=True)
    parser.add_argument(
        "--gnomad-variants",
        required=True,
        help="TSV of the gnomAD variants over the coding sites: read if it exists,"
        " otherwise queried from the public VCFs and written here.",
    )
    parser.add_argument("--out-prefix", required=True)
    main(parser.parse_args())
