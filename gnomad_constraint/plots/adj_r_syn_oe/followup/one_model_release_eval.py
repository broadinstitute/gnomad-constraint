"""Score the single cloglog models on pLoF and missense, as for the pipeline runs.

For the notes ``notes_one_model_for_expected.md``, linked from the rr report.

Input: ``one_model_fit.py --apply-cells`` output (``<prefix>.by_csq.per_transcript.tsv.gz``):
per Ensembl MANE Select transcript on the autosomes, the observed count and each
model's expected count for synonymous, missense and pLoF (LOFTEE HC, the release's
``lof``) variants, every model fit on the other chromosome half. For each model:

- LOEUF: the upper end of the 90% gamma interval of pLoF obs/exp (shape obs + 1, scale
  1 / exp), the pipeline's ``oe_ci_gamma``, checked against the release for the current
  stack;
- synonymous, missense and pLoF obs/exp in total, and by the release's LOEUF decile;
- against the current stack: LOEUF rank correlation, genes changing LOEUF decile (by
  rank among the scored transcripts), moving in or out of the most constrained decile,
  and crossing LOEUF 0.6; the AUC of LOEUF for ClinGen haploinsufficient genes;
- how many genes leave the most constrained LOEUF decile, and the share of them with
  the release's ``low_exome_mapping_quality`` flag;
- NDD and HI precision-recall with the flagship preprint's gene lists and rules
  (``disease_gene_pr.score``), on all genes and on genes with no release QC flag.

Writes ``data/one_model/one_model_release.{scores,deciles}.tsv``,
``data/one_model/one_model_pr.{auprc,curves}.tsv`` and
``data/one_model/one_model_pr_no_qc_flag.{auprc,curves}.tsv``.
"""

import argparse
import logging
import os

import numpy as np
import pandas as pd
from disease_gene_pr import read_list, score
from remap_rr_to_shet import read_release_flags
from scipy.stats import gamma, spearmanr
from score_release import auc, read_clingen

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

KEYS = {
    "current stack, no rr": "none",
    "current stack + ridge rr": "stack_ridge",
    "current stack x Roulette": "stack_roulette",
    "current stack x Roulette + ridge rr": "stack_roulette_ridge",
    "one model: groups + coverage": "om_groups_coverage",
    "one model: + window MQ": "om_mq",
    "one model: + DNM regional": "om_dnm",
    "one model: + DNM regional + window MQ": "om_dnm_mq",
    "one model: + gnomAD regional + window MQ": "om_gnomad_mq",
    "one model: + window MQ + Roulette": "om_mq_roulette",
    "one model: + DNM regional + window MQ + Roulette": "om_full",
}
"""Model names in ``one_model_fit.py`` and their short keys here."""


def loeuf(obs: pd.Series, exp: pd.Series) -> pd.Series:
    """Get the upper end of the 90% gamma interval of obs/exp."""
    return pd.Series(
        gamma.ppf(0.95, obs + 1, scale=1 / exp.where(exp > 0)), index=obs.index
    )


def main(args: argparse.Namespace) -> None:
    """Score every model."""
    d = pd.read_csv(args.per_transcript, sep="\t")
    models = [k for k in KEYS if k in d.columns]
    rel = pd.read_csv(
        args.release_tsv,
        sep="\t",
        compression="gzip",
        usecols=[
            "gene",
            "gene_id",
            "transcript",
            "mane_select",
            "lof.obs",
            "lof.exp",
            "lof.oe_ci.upper",
            "lof.oe_ci.upper_bin_decile",
        ],
    )
    rel = rel[rel.mane_select.astype(str).str.lower().eq("true")].drop_duplicates(
        "transcript"
    )
    rel = rel.set_index("transcript")

    wide = {
        csq: d[d.csq == csq].set_index("transcript") for csq in ("syn", "mis", "lof")
    }
    lof = wide["lof"]
    t = lof[lof["current stack, no rr"] > 0].join(rel, how="inner")
    # The current stack against the release.
    r = loeuf(t.observed, t["current stack, no rr"]) / t["lof.oe_ci.upper"]
    logger.info(
        "Current stack against the release, %d transcripts: pLoF obs equal %.4f,"
        " pLoF exp within 0.1%% %.4f, LOEUF within 0.1%% %.4f.",
        len(t),
        (t.observed == t["lof.obs"]).mean(),
        (abs(t["current stack, no rr"] / t["lof.exp"] - 1) < 0.001).mean(),
        (abs(r - 1) < 0.001).mean(),
    )

    g = t[["gene", "gene_id"]].copy()
    g["exp_none"] = t["current stack, no rr"]
    for m in models:
        g[KEYS[m]] = loeuf(t.observed, t[m])
    keys = [KEYS[m] for m in models]
    g = g[np.isfinite(g[keys]).all(axis=1)]
    flags = (
        read_release_flags(args.release_tsv)
        .drop_duplicates("transcript")
        .set_index("transcript")
        .gene_flags.reindex(g.index)
        .map(lambda f: f if isinstance(f, list) else [])
    )
    low_mq = flags.map(lambda f: "low_exome_mapping_quality" in f)
    no_qc_flag = flags.map(len).eq(0)
    clingen = read_clingen(args.clingen)
    hi_cg = g.gene.map(clingen).eq(3)
    n = len(g)
    dec = {k: (g[k].rank(method="first") - 1) * 10 // n for k in keys}

    rows, deciles = [], []
    for m in models:
        k = KEYS[m]
        lt, lt0 = g[k].lt(0.6), g["none"].lt(0.6)
        top, top0 = dec[k].eq(0), dec["none"].eq(0)
        row = {"model": m, "key": k, "n": n}
        for csq in ("syn", "mis", "lof"):
            w = wide[csq].reindex(g.index)
            row[f"{csq}_oe"] = float(w.observed.sum() / w[m].sum())
        row.update(
            {
                "lof_over_syn": row["lof_oe"] / row["syn_oe"],
                "mis_over_syn": row["mis_oe"] / row["syn_oe"],
                "loeuf_spearman_vs_none": spearmanr(g[k], g["none"]).statistic,
                "decile_changed": float(dec[k].ne(dec["none"]).mean()),
                "decile_changed_2plus": float(
                    (dec[k] - dec["none"]).abs().ge(2).mean()
                ),
                "into_top_decile": int((top & ~top0).sum()),
                "out_of_top_decile": int((~top & top0).sum()),
                "into_loeuf_lt_0.6": int((lt & ~lt0).sum()),
                "out_of_loeuf_lt_0.6": int((~lt & lt0).sum()),
                "hi_auc": auc(g[k][hi_cg], g[k][~hi_cg]),
                "left_top_decile_low_mq_share": (
                    float(low_mq[~top & top0].mean()) if (~top & top0).any() else np.nan
                ),
            }
        )
        rows.append(row)
        bins = t["lof.oe_ci.upper_bin_decile"].reindex(g.index)
        for csq in ("syn", "mis", "lof"):
            w = wide[csq].reindex(g.index).assign(bin=bins)
            s = w.groupby("bin")[["observed", m]].sum()
            for b, v in (s.observed / s[m]).items():
                deciles.append(
                    {"model": m, "key": k, "csq": csq, "decile": int(b), "oe": v}
                )
    out = pd.DataFrame(rows)
    out.to_csv(f"{args.out_prefix}_release.scores.tsv", sep="\t", index=False)
    pd.DataFrame(deciles).to_csv(
        f"{args.out_prefix}_release.deciles.tsv", sep="\t", index=False
    )

    # NDD and HI precision-recall, as disease_gene_pr.py does for the pipeline runs.
    lists = args.lists
    ndd = read_list(os.path.join(lists, "ndd.txt"), header=True)
    hi_severe = read_list(
        os.path.join(lists, "haploinsufficiency_severe_curated_2016.tsv")
    )
    hi_severe |= read_list(os.path.join(lists, "hi_added_genes.txt"))
    ces = set(pd.read_csv(os.path.join(lists, "ces_genes_list.tsv"), sep="\t").gene_id)
    hi = read_list(os.path.join(lists, "hi_genes_union.tsv"), header=True)
    for suffix, gs in (("", g), ("_no_qc_flag", g[no_qc_flag.reindex(g.index)])):
        gs = gs.drop_duplicates("gene")
        g_ndd = gs[~gs.gene.isin(hi_severe) & ~gs.gene_id.isin(ces)]
        rng = np.random.default_rng(0)
        pr_rows, curves = [], []
        for test, sub, pos in (
            ("NDD", g_ndd, g_ndd.gene.isin(ndd)),
            ("HI", gs, gs.gene.isin(hi)),
        ):
            logger.info(
                "%s%s: %d genes, %d positive.", test, suffix, len(sub), pos.sum()
            )
            a, c = score(sub, pos, test, rng, args.n_boot, keys)
            pr_rows += a
            curves += c
        pr = pd.DataFrame(pr_rows)
        pr.to_csv(f"{args.out_prefix}_pr{suffix}.auprc.tsv", sep="\t", index=False)
        pd.concat(curves).to_csv(
            f"{args.out_prefix}_pr{suffix}.curves.tsv", sep="\t", index=False
        )
    pd.set_option("display.width", 250)
    logger.info("\n%s", out.drop(columns="model").round(4).to_string(index=False))
    logger.info("\n%s", pr.round(4).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--per-transcript",
        default="data/one_model/one_model_full_v3.by_csq.per_transcript.tsv.gz",
    )
    parser.add_argument(
        "--release-tsv", default="../gnomad.v4.1.2.constraint_metrics.tsv.bgz"
    )
    parser.add_argument(
        "--clingen",
        default="data/ext_clingen/ClinGen_haploinsufficiency_gene_GRCh38.bed",
    )
    parser.add_argument("--lists", default="data/ext_flagship")
    parser.add_argument("--n-boot", type=int, default=1000)
    parser.add_argument("--out-prefix", default="data/one_model/one_model")
    main(parser.parse_args())
