"""Check remove_qc_flagged_genes on a few loci in and out of flagged genes."""

import hail as hl

from gnomad_constraint.resources.resource_utils import get_loeuf_decile_ht
from gnomad_constraint.utils.regional_rate import remove_qc_flagged_genes

hl.init(tmp_dir="gs://gnomad-tmp-4day")
metrics = get_loeuf_decile_ht().ht()
flagged = metrics.filter(
    metrics.mane_select
    & metrics.transcript.startswith("ENST")
    & (hl.len(metrics.gene_flags) > 0)
)
print("flagged MANE Ensembl transcripts:", flagged.count())
ex = flagged.take(2)
loci = [hl.Locus(r.chromosome, r.start_position + 10, "GRCh38") for r in ex] + [
    hl.Locus("chr1", 1_000_000, "GRCh38"),
    hl.Locus("chr7", 117_559_600, "GRCh38"),
]
ht = hl.Table.parallelize(
    [hl.Struct(locus=l, y=1) for l in loci],
    hl.tstruct(locus=hl.tlocus("GRCh38"), y=hl.tint32),
).key_by("locus")
out = remove_qc_flagged_genes(ht, metrics)
print("kept:", out.locus.collect(), "flagged genes used:", [r.gene for r in ex])
