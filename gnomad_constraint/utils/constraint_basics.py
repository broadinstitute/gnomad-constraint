# noqa: D100

import hail as hl

from gnomad.utils.constraint import (
    trimer_from_heptamer,
    annotate_mutation_type,
    collapse_strand,
)
from gnomad.utils.vep import add_most_severe_csq_to_tc_within_vep_root


def add_vep_context_annotations(
    ht: hl.Table, annotated_context_ht_path: str
) -> hl.Table:
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
    :param annotated_context_ht_path: Path to VEP context Table.
    """
    context_ht = hl.read_table(annotated_context_ht_path).drop("a_index", "was_split")
    context_ht = context_ht.annotate(vep=context_ht.vep.drop("colocated_variants"))
    return ht.annotate(**context_ht[ht.key])


def prepare_ht_for_constraint_calculations(
    ht: hl.Table, annotate_coverage: bool = True
):
    """
    Filter input Table and add annotations used in constraint calculations.

    Function filters to SNPs, removes rows with undefined contexts, collapses strands
    to deduplicate trimer or heptamer contexts, and annotates the input Table.

    The following annotations are added to the output Table:
        - ref
        - alt
        - methylation_level
        - annotations added by `annotate_mutation_type()`, `collapse_strand()`, and `add_most_severe_csq_to_tc_within_vep_root()`
        - exome_coverage (optionally added)

    :param ht: Input Table to be annotated.
    :param annotate_coverage: Whether to annotate the coverage of exome. Defaults to True.
    :return: Table with annotations.
    """
    ht = trimer_from_heptamer(ht)

    if "filters" in ht.row_value.keys():
        ht = ht.annotate(pass_filters=hl.len(ht.filters) == 0)

    # Add annotations for 'ref' and 'alt'
    ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1])
    # Filter to SNPs and context fields where the bases are either A, T, C, or G
    ht = ht.filter(hl.is_snp(ht.ref, ht.alt) & ht.context.matches(f"[ATCG]{3}"))
    # Annotate mutation type (such as "CpG", "non-CpG transition", "transversion") and collapse strands to deduplicate the context
    ht = annotate_mutation_type(collapse_strand(ht))
    # Add annotation for the methylation level
    annotation = {
        "methylation_level": hl.case()
        .when(ht.cpg & (ht.methylation.MEAN > 0.6), 2)
        .when(ht.cpg & (ht.methylation.MEAN > 0.2), 1)
        .default(0)
    }
    # Add annotation for the median exome coverage if requested
    if annotate_coverage:
        annotation["exome_coverage"] = ht.coverage.exomes.median
    ht = ht.annotate(**annotation)
    # Add most_severe_consequence annotation to 'transcript_consequences' within the vep root annotation.
    ht = add_most_severe_csq_to_tc_within_vep_root(ht)
    return ht
