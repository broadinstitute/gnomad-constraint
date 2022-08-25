# noqa: D100

import hail as hl

from gnomad.utils.constraint import (
    trimer_from_heptamer,
    annotate_variant_types,
    collapse_strand,
)

from gnomad.utils.vep import add_most_severe_consequence_to_consequence


def add_vep_context_annotations(ht: hl.Table, split_context_ht_path: str) -> hl.Table:
    """
    Add annotations from VEP context Table to gnomAD data.

    Function adds the following annotations:
        - context
        - methylation
        - coverage
        - gerp
        - pass_filters

    Function drops `a_index`, `was_split`, and`colocated_variants` annotations from gnomAD data.

    .. note::
        Function expects that multiallelic variants in the VEP context Table have been split.

    :param ht: gnomAD exomes or genomes public Hail Table.
    :param split_context_ht_path: Path to VEP context Table.
    :param output_ht_path: Path to output Table.
    :param overwrite: Whether to overwrite existing data. Defaults to False.
    """
    context_ht = hl.read_table(split_context_ht_path).drop("a_index", "was_split")
    context_ht = context_ht.annotate(vep=context_ht.vep.drop("colocated_variants"))
    return ht.annotate(**context_ht[ht.key], pass_filters=hl.len(ht.filters) == 0)


def prepare_ht_for_constraint_calculations(
    ht: hl.Table, trimers: bool = False, annotate_coverage: bool = True
):
    """
    Filter input Table and add annotations used in constraint calculations.

    Function filters to SNPs, removes rows with undefined contexts, collapses strands
    to deduplicate trimer or heptamer contexts, and annotates the input Table.

    The following annotations are added to the output Table:
        - ref
        - alt
        - methylation_level
        - exome_coverage
        - annotations added by `annotate_variant_types` and `collapse_strand`

    :param ht: Input Table to be annotated.
    :param trimers: Whether to use trimers or heptamers. Defaults to False.
    :param annotate_coverage: Whether to annotate the coverage of exome. Defaults to True.
    :return: Table with annotations.
    """
    if trimers:
        ht = trimer_from_heptamer(ht)

    ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1])
    ht = ht.filter(hl.is_snp(ht.ref, ht.alt))
    ht = annotate_variant_types(collapse_strand(ht), not trimers)

    annotation = {
        "methylation_level": hl.case()
        .when(ht.cpg & (ht.methylation.MEAN > 0.6), 2)
        .when(ht.cpg & (ht.methylation.MEAN > 0.2), 1)
        .default(0)
    }
    if annotate_coverage:
        annotation["exome_coverage"] = ht.coverage.exomes.median
    ht = ht.annotate(**annotation)
    ht = add_most_severe_csq_to_tc_within_ht(ht)
    return ht


def add_most_severe_csq_to_tc_within_ht(t):
    """
    Add most_severe_consequence annotation to 'transcript_consequences' within the vep annotation.

    :param t: Input Table or MatrixTable.
    :return: Input Table or MatrixTable with most_severe_consequence annotation added.
    """
    annotation = t.vep.annotate(
        transcript_consequences=t.vep.transcript_consequences.map(
            add_most_severe_consequence_to_consequence
        )
    )
    return (
        t.annotate_rows(vep=annotation)
        if isinstance(t, hl.MatrixTable)
        else t.annotate(vep=annotation)
    )