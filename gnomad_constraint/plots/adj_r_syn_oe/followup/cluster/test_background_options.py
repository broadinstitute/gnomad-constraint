"""Check the rr background options on a small synthetic set of synonymous sites."""

import hail as hl

from gnomad_constraint.resources.constants import GENOMIC_FEATURE_COLUMNS
from gnomad_constraint.utils.regional_rate import create_rr_training_ht

hl.init(tmp_dir="gs://gnomad-tmp-4day")
rows = []
for i in range(6000):
    ctx = ["ACA", "ACA", "ACG", "TAT"][i % 4]
    for alt in ["A", "G"] if i % 2 else ["T"]:
        rows.append(
            hl.Struct(
                locus=hl.Locus("chr1", 1_000_000 + i, "GRCh38"),
                alleles=["C", alt],
                context=ctx,
            )
        )
schema = hl.tstruct(
    locus=hl.tlocus("GRCh38"), alleles=hl.tarray(hl.tstr), context=hl.tstr
)
syn = hl.Table.parallelize(rows, schema).key_by("locus", "alleles")
exomes = syn.filter(syn.locus.position % 3 == 0).select()
dnm = syn.filter(syn.locus.position % 97 == 0).select("context")
feat = hl.Table.parallelize(
    [
        hl.Struct(
            interval=hl.Interval(
                hl.Locus("chr1", 999_000, "GRCh38"),
                hl.Locus("chr1", 1_010_000, "GRCh38"),
            ),
            **{c: 1.0 for c in GENOMIC_FEATURE_COLUMNS},
        )
    ],
    hl.tstruct(
        interval=hl.tinterval(hl.tlocus("GRCh38")),
        **{c: hl.tfloat64 for c in GENOMIC_FEATURE_COLUMNS},
    ),
).key_by("interval")
n_variants = syn.count()
n_loci = syn.key_by("locus").distinct().count()
print(
    "variants",
    n_variants,
    "loci",
    n_loci,
    "observed",
    exomes.count(),
    "dnms",
    dnm.count(),
)
for kw in (
    {},
    {"background_sites": "all"},
    {"background_unit": "variant"},
    {"background_sites": "all", "background_unit": "variant"},
    {"context_sampling": "proportional"},
    {"context_sampling": "even"},
):
    ht = create_rr_training_ht(syn, dnm, exomes, feat, 0.5, 1, **kw)
    bg = ht.filter(ht.y == 0)
    print(
        kw,
        "background rows",
        bg.count(),
        "per context",
        dict(bg.aggregate(hl.agg.counter(bg.context))),
        "distinct loci",
        bg.select().distinct().count(),
        "DNM rows",
        ht.filter(ht.y == 1).count(),
    )
