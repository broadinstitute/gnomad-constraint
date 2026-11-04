"""Check the per-variant rate (Roulette) code paths on small synthetic data.

Runs locally (``hl.init`` with a local tmp dir) or on a cluster. Each check asserts
what the code should do and prints a line when it passes.
"""

import sys

import hail as hl
import numpy as np
import pandas as pd

from gnomad_constraint.resources.constants import (
    GENOMIC_FEATURE_COLUMNS,
    VARIANT_RATE_MULTIPLIER_FIELD,
)
from gnomad_constraint.utils.constraint import annotate_regional_rr_expected
from gnomad_constraint.utils.regional_rate import (
    RR_RATE_OFFSET_FIELD,
    compute_rr_normalization_factors,
    create_rr_training_ht,
    fit_rr_context_model,
    get_variant_rate_ht,
    rr_training_ht_to_pandas,
)

hl.init(tmp_dir=sys.argv[1] if len(sys.argv) > 1 else "/tmp/hail_tmp", quiet=True)
LOC = hl.tlocus("GRCh38")


def locus(pos: int) -> hl.Locus:
    """Get a chr1 locus."""
    return hl.Locus("chr1", pos, "GRCh38")


# --- get_variant_rate_ht -------------------------------------------------------
ctx_rows, rou_rows = [], []
for i in range(400):
    pos = 1_000_000 + i
    context = ["ACA", "ACG", "TAT"][i % 3]
    ref = "A" if context == "TAT" else "C"
    for alt, mu in (("T", 1.0 + (i % 5)), ("G", 0.5 + (i % 2))):
        ctx_rows.append(
            hl.Struct(
                locus=locus(pos),
                alleles=[ref, alt],
                context=context,
                ref=ref,
                alt=alt,
                methylation_level=0,
            )
        )
        filters = {"low"} if i % 50 == 0 else set()
        rou_rows.append(
            hl.Struct(
                locus=locus(pos),
                alleles=[ref, alt],
                roulette_mu=mu,
                roulette_filters=filters,
            )
        )
# A site outside the coding interval.
rou_rows.append(
    hl.Struct(
        locus=locus(5_000_000),
        alleles=["C", "T"],
        roulette_mu=9.0,
        roulette_filters=set(),
    )
)
context_ht = hl.Table.parallelize(
    ctx_rows,
    hl.tstruct(
        locus=LOC,
        alleles=hl.tarray(hl.tstr),
        context=hl.tstr,
        ref=hl.tstr,
        alt=hl.tstr,
        methylation_level=hl.tint32,
    ),
).key_by("locus", "alleles")
roulette_ht = hl.Table.parallelize(
    rou_rows,
    hl.tstruct(
        locus=LOC,
        alleles=hl.tarray(hl.tstr),
        roulette_mu=hl.tfloat64,
        roulette_filters=hl.tset(hl.tstr),
    ),
).key_by("locus", "alleles")
cds_ht = hl.Table.parallelize(
    [hl.Struct(interval=hl.Interval(locus(999_990), locus(1_000_500)))],
    hl.tstruct(interval=hl.tinterval(LOC)),
).key_by("interval")

rate_ht = get_variant_rate_ht(roulette_ht, context_ht, cds_ht)
df = rate_ht.to_pandas()
assert len(df) == len(ctx_rows), "the non-coding site should be dropped"
low = df.locus.astype(str).isin([str(locus(1_000_000 + i)) for i in range(0, 400, 50)])
assert df[low & (df.context != "ACG")][VARIANT_RATE_MULTIPLIER_FIELD].isna().all()
assert (df[df.context == "ACG"][VARIANT_RATE_MULTIPLIER_FIELD] == 1.0).all()
means = (
    df[~low & (df.context != "ACG")]
    .groupby(["context", "alt"])[VARIANT_RATE_MULTIPLIER_FIELD]
    .mean()
)
assert np.allclose(means, 1.0, atol=0.02), means
print(
    "PASS get_variant_rate_ht: non-coding dropped, low-quality missing, CpG 1, mean 1"
)

rate_cpg = get_variant_rate_ht(roulette_ht, context_ht, cds_ht, include_cpg=True)
cpg = rate_cpg.filter(rate_cpg.context == "ACG").to_pandas()
assert cpg[VARIANT_RATE_MULTIPLIER_FIELD].std() > 0.1
print("PASS get_variant_rate_ht(include_cpg=True): CpG multipliers vary")

# --- create_rr_training_ht with the offset -------------------------------------
syn_ht = context_ht.select("context", "ref", "alt", "methylation_level")
dnm_ht = syn_ht.filter(syn_ht.locus.position % 7 == 0).select("context")
exomes_ht = syn_ht.filter(syn_ht.locus.position % 5 == 0).select()
feat_ht = hl.Table.parallelize(
    [
        hl.Struct(
            interval=hl.Interval(locus(999_000), locus(1_001_000)),
            **{c: 1.0 for c in GENOMIC_FEATURE_COLUMNS},
        )
    ],
    hl.tstruct(
        interval=hl.tinterval(LOC), **{c: hl.tfloat64 for c in GENOMIC_FEATURE_COLUMNS}
    ),
).key_by("interval")
train_ht = create_rr_training_ht(
    syn_ht, dnm_ht, exomes_ht, feat_ht, 1.0, 0, variant_rate_ht=rate_ht
)
tr = train_ht.to_pandas()
m = (
    df.assign(key=df.locus.astype(str) + df.alleles.astype(str))
    .set_index("key")[VARIANT_RATE_MULTIPLIER_FIELD]
    .fillna(1.0)
)
dn = tr[tr.y == 1]
exp_dnm = np.log(m.loc[dn.locus.astype(str) + dn.alleles.astype(str)].to_numpy())
assert np.allclose(dn[RR_RATE_OFFSET_FIELD], exp_dnm)
bg = tr[tr.y == 0]
loc_mean = (
    df.assign(m=df[VARIANT_RATE_MULTIPLIER_FIELD].fillna(1.0))
    .groupby(df.locus.astype(str))
    .m.mean()
)
assert np.allclose(
    bg[RR_RATE_OFFSET_FIELD], np.log(loc_mean.loc[bg.locus.astype(str)].to_numpy())
)
assert RR_RATE_OFFSET_FIELD in rr_training_ht_to_pandas(train_ht).columns
print("PASS create_rr_training_ht: DNM offset = own variant, background = locus mean")

# --- fit_rr_context_model with the offset ---------------------------------------
rng = np.random.default_rng(1)
n = 200_000
z = rng.normal(size=(n, 2))
log_m = 0.5 * z[:, 0] + rng.normal(scale=0.3, size=n)
eta = -4 + log_m + 0.2 * z[:, 1]
y = (rng.random(n) < 1 / (1 + np.exp(-eta))).astype(int)
sim = pd.DataFrame(
    {"locus": "chr1:1", "f1": z[:, 0], "f2": z[:, 1], "y": y, "off": log_m}
)
with_off = fit_rr_context_model(sim, ["f1", "f2"], offset_col="off")
without = fit_rr_context_model(sim, ["f1", "f2"])
assert abs(with_off["coefs"][0]) < 0.03 and abs(with_off["coefs"][1] - 0.2) < 0.03
assert without["coefs"][0] > 0.3, "without the offset, f1 should absorb the rate"
assert with_off["rate_offset"] and not without["rate_offset"]
print(
    "PASS fit_rr_context_model: with offset coefs",
    np.round(with_off["coefs"], 3),
    "without",
    np.round(without["coefs"], 3),
)
zero = sim.assign(off=0.0)
a = fit_rr_context_model(zero, ["f1", "f2"], ridge_penalty=0.3, offset_col="off")
b = fit_rr_context_model(zero, ["f1", "f2"], ridge_penalty=0.3)
assert np.allclose(a["coefs"], b["coefs"], atol=1e-3), (a["coefs"], b["coefs"])
assert abs(a["intercept"] - b["intercept"]) < 1e-3
print("PASS ridge with a zero offset matches the sklearn ridge fit")

# --- normalization and application ------------------------------------------------
rr_ht = hl.Table.parallelize(
    [
        hl.Struct(
            interval=hl.Interval(locus(999_000), locus(1_001_000)),
            adj_r={"ACA": 2.0, "TAT": 1.0, "ACG": 1.0},
        )
    ],
    hl.tstruct(interval=hl.tinterval(LOC), adj_r=hl.tdict(hl.tstr, hl.tfloat64)),
).key_by("interval")
syn_exp = syn_ht.select("context", transcript="T1", obs=1, exp=1.0)
f_plain = compute_rr_normalization_factors(syn_exp, rr_ht)
syn_m = syn_exp.annotate(
    **{
        VARIANT_RATE_MULTIPLIER_FIELD: rate_ht[syn_exp.key][
            VARIANT_RATE_MULTIPLIER_FIELD
        ]
    }
)
f_rate = compute_rr_normalization_factors(syn_m, rr_ht)
aca = df[df.context == "ACA"][VARIANT_RATE_MULTIPLIER_FIELD].fillna(1.0).mean()
assert abs(f_plain["ACA"] - 2.0) < 1e-9 and abs(f_rate["ACA"] - 2.0 * aca) < 1e-9
print("PASS compute_rr_normalization_factors includes the multiplier")

pv = syn_ht.select("context", expected_variants=[1.0, 2.0])
out = annotate_regional_rr_expected(
    pv, rr_ht, "rr", variant_rate_ht=rate_ht, variant_rate_ht_path="rates"
)
res = out.to_pandas()
mm = res[VARIANT_RATE_MULTIPLIER_FIELD]
rr = res.context.map({"ACA": 2.0, "TAT": 1.0, "ACG": 1.0})
assert np.allclose(res.expected_variants.str[1], 2.0 * rr * mm)
assert (
    mm.notna().all()
    and hl.eval(out.regional_rr_globals.variant_rate_ht_path) == "rates"
)
plain = annotate_regional_rr_expected(pv, rr_ht, "rr").to_pandas()
assert np.allclose(plain.expected_variants.str[1], 2.0 * rr)
print("PASS annotate_regional_rr_expected multiplies by rr and the multiplier")
