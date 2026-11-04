"""Redo the LOEUF and s_het decile summary with explained synonymous outliers removed.

Same transcripts and binning as ``remap_rr_to_shet.py --summarize`` (autosomal, with a
v4.1.2 LOEUF decile and an s_het, no release QC gene flag), built from the per-
transcript expected counts in ``syn_z_by_histone_group.tsv``
(``remap_rr_to_shet.py --syn-z``). Then transcripts that are synonymous outliers
(|z| > 8) under no correction or any synonymous rr, and whose excess has an explained
cause, are dropped:

- histone genes (``GENE_FAMILIES`` in remap_rr_to_shet.py)
- repeat, segdup or low-coverage genes (segdup fraction >= 0.5, AN90 < 0.5, or a MUC,
  OR or PSG gene)
- genes in a maternal C>G hotspot region (Jónsson et al. 2017)

Outliers without one of these causes stay. Gnocchi ``adj_r`` is not used to pick
outliers, since its outliers come from the correction itself.
"""

import argparse
import logging

import numpy as np
import pandas as pd
from outlier_drivers import MATERNAL_HOTSPOTS, REPEAT_FAMILIES
from remap_rr_to_shet import GENE_FAMILIES, read_release_flags, read_shet

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

CORRECTIONS = [
    "none",
    "adj_r",
    "syn_adj_r",
    "ruchit_17k",
    "ruchit_satterstrom",
    "recommended",
    "syn_adj_r_norm",
]
OUTLIER_CORRECTIONS = [
    "none",
    "syn_adj_r",
    "syn_adj_r_norm",
    "ruchit_17k",
    "ruchit_satterstrom",
    "recommended",
]


def explained_cause(df: pd.DataFrame) -> pd.Series:
    """
    Get the explained cause of each transcript's synonymous excess, if any.

    :param df: Transcripts with ``gene``, ``chromosome``, ``mid`` and gene quality metrics.
    :return: "histone gene", "repeat or segdup", "maternal hotspot region" or "".
    """
    cause = pd.Series("", index=df.index)
    for name, (c, s, e) in MATERNAL_HOTSPOTS.items():
        cause[(df.chromosome == c) & df.mid.between(s, e, inclusive="left")] = (
            "maternal hotspot region"
        )
    repeat = (
        (df["gene_quality_metrics.exome_prop_segdup"] >= 0.5)
        | (df["gene_quality_metrics.exome_prop_bp_AN90"] < 0.5)
        | df.gene.fillna("").str.match(REPEAT_FAMILIES)
    )
    cause[repeat] = "repeat or segdup"
    cause[df.gene.fillna("").str.match(GENE_FAMILIES["histone genes"])] = "histone gene"
    return cause


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sum syn obs and expected per s_het and LOEUF decile for each correction.

    :param df: Transcripts with decile columns, ``obs`` and ``exp_*``.
    :return: One row per binning, decile and correction.
    """
    rows = []
    for binning in ("shet_decile", "loeuf_decile"):
        g = df.groupby(binning)
        for corr in CORRECTIONS:
            s = g[["obs", f"exp_{corr}"]].sum()
            for dec, r in s.iterrows():
                rows.append(
                    {
                        "binning": binning.replace("_decile", ""),
                        "decile": dec,
                        "correction": corr,
                        "n": int(g.size()[dec]),
                        "obs": r.obs,
                        "exp": r[f"exp_{corr}"],
                        "oe": r.obs / r[f"exp_{corr}"],
                    }
                )
    return pd.DataFrame(rows)


def main(args: argparse.Namespace) -> None:
    """Write the decile summary with explained outliers removed."""
    zs = pd.read_csv(args.syn_z, sep="\t")
    tx = pd.read_csv(args.transcripts, sep="\t", compression="gzip", na_values=["NA"])
    tx = tx[tx.is_xy.astype(str).str.lower().eq("false") & tx.decile.notna()]
    rel = read_release_flags(args.release_tsv)
    rel = rel[rel.gene_flags.map(len) == 0][["transcript", "gene_id"]]
    tx = tx.merge(rel, on="transcript").merge(read_shet(args.shet), on="gene_id")
    tx["shet_decile"] = pd.qcut(
        tx.shet.rank(method="first", ascending=False), 10, labels=False
    )
    tx["loeuf_decile"] = tx.decile.astype(int)
    pos = pd.read_csv(
        args.release_tsv,
        sep="\t",
        compression="gzip",
        na_values=["NA"],
        usecols=[
            "transcript",
            "chromosome",
            "start_position",
            "end_position",
            "gene_quality_metrics.exome_prop_segdup",
            "gene_quality_metrics.exome_prop_bp_AN90",
        ],
    ).drop_duplicates("transcript")
    df = (
        tx[["transcript", "shet_decile", "loeuf_decile"]]
        .merge(zs, on="transcript")
        .merge(pos, on="transcript")
    )
    df["mid"] = (df.start_position + df.end_position) / 2
    logger.info("%d transcripts in the summary set.", len(df))
    if args.check_against:
        old = pd.read_csv(args.check_against, sep="\t")
        new = summarize(df).merge(
            old, on=["binning", "decile", "correction"], suffixes=("", "_old")
        )
        logger.info(
            "Full set vs existing summary: max |obs/exp diff| %.2g over %d rows.",
            (new.oe - new.oe_old).abs().max(),
            len(new),
        )

    outlier = np.zeros(len(df), dtype=bool)
    for c in OUTLIER_CORRECTIONS:
        outlier |= (df[f"z_{c}"].abs() > 8).to_numpy()
    df["cause"] = explained_cause(df)
    drop = outlier & (df.cause != "")
    logger.info(
        "%d outliers under any correction; %d removed as explained (%s); %d kept as unexplained.",
        outlier.sum(),
        drop.sum(),
        df[drop].cause.value_counts().to_dict(),
        (outlier & ~drop).sum(),
    )
    df[drop][["transcript", "gene", "cause"]].to_csv(
        args.removed_out, sep="\t", index=False
    )
    summarize(df[~drop]).to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--syn-z", required=True, help="syn_z_by_histone_group.tsv")
    parser.add_argument("--transcripts", required=True)
    parser.add_argument("--release-tsv", required=True)
    parser.add_argument(
        "--shet", required=True, help="Kar et al. 2026 Supplementary Table 3 xlsx."
    )
    parser.add_argument(
        "--check-against", help="Existing decile summary to check the full set against."
    )
    parser.add_argument("--out", required=True)
    parser.add_argument("--removed-out", required=True)
    main(parser.parse_args())
