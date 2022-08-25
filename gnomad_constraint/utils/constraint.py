"""Script containing utility functions used in the constraint pipeline."""

import hail as hl
from gnomad.utils.constraint import (
    annotate_mutation_type,
    collapse_strand,
    trimer_from_heptamer,
)
from gnomad.utils.vep import add_most_severe_csq_to_tc_within_vep_root


def add_vep_context_annotations(
    ht: hl.Table, annotated_context_ht: hl.Table
) -> hl.Table:
    """
    Add annotations from VEP context Table to gnomAD data.

    Function adds the following annotations:
        - context
        - methylation
        - coverage
        - gerp

    Function drops `a_index`, `was_split`, and`colocated_variants` annotations from gnomAD data.

    .. note::
        Function expects that multiallelic variants in the VEP context Table have been split.

    :param ht: gnomAD exomes or genomes public Hail Table.
    :param annotated_context_ht: VEP context Table.
    :return: Table with annotations.
    """
    context_ht = annotated_context_ht.drop("a_index", "was_split")
    context_ht = context_ht.annotate(vep=context_ht.vep.drop("colocated_variants"))
    return ht.annotate(**context_ht[ht.key])


def prepare_ht_for_constraint_calculations(ht: hl.Table) -> hl.Table:
    """
    Filter input Table and add annotations used in constraint calculations.

    Function filters to SNPs, removes rows with undefined contexts, collapses strands
    to deduplicate trimer or heptamer contexts, and annotates the input Table.

    The following annotations are added to the output Table:
        - ref
        - alt
        - methylation_level
        - exome_coverage
        - pass_filters - Whether the variant passed all variant filters
        - annotations added by `annotate_mutation_type()`, `collapse_strand()`, and `add_most_severe_csq_to_tc_within_vep_root()`

    :param ht: Input Table to be annotated.
    :return: Table with annotations.
    """
    ht = trimer_from_heptamer(ht)

    if "filters" in ht.row_value.keys():
        ht = ht.annotate(pass_filters=hl.len(ht.filters) == 0)

    # Add annotations for 'ref' and 'alt'
    ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1])
    # Filter to SNPs and context fields where the bases are either A, T, C, or G
    ht = ht.filter(hl.is_snp(ht.ref, ht.alt) & ht.context.matches(f"[ATCG]{{{3}}}"))
    # Annotate mutation type (such as "CpG", "non-CpG transition", "transversion") and collapse strands to deduplicate the context
    ht = annotate_mutation_type(collapse_strand(ht))
    # Add annotation for the methylation level
    annotation = {
        "methylation_level": hl.case()
        .when(ht.cpg & (ht.methylation.MEAN > 0.6), 2)
        .when(ht.cpg & (ht.methylation.MEAN > 0.2), 1)
        .default(0)
    }
    # Add annotation for the median exome coverage
    annotation["exome_coverage"] = ht.coverage.exomes.median
    ht = ht.annotate(**annotation)

    # Add most_severe_consequence annotation to 'transcript_consequences' within the vep root annotation.
    ht = add_most_severe_csq_to_tc_within_vep_root(ht)

    # Filter out locus with undefined exome coverage
    ht = ht.filter(hl.is_defined(ht.exome_coverage))
    return ht

def get_proportion_observed_by_coverage(
    exome_ht: hl.Table,
    context_ht: hl.Table,
    mutation_ht: hl.Table,
    recompute_possible: bool = False,
    dataset: str = "gnomad",
    impose_high_af_cutoff_upfront: bool = True,
) -> hl.Table:
    """
    Count the observed variants and possible variants by exome coverage.

    :param exome_ht: Preprocessed exome Table.
    :param context_ht: Preprocessed context Table.
    :param mutation_ht: Preprocessed mutation rate Table.
    :param recompute_possible: Whether to use context Table to recompute the number of possible variants instead of using a precomputed intermediate Table if it exists. Defaults to False.
    :param dataset: Dataset to use when computing frequency index. Defaults to 'gnomad'.
    :param impose_high_af_cutoff_upfront: Whether to remove high frequency alleles. Defaults to True.
    :return: Table with observed variant and possible variant count.
    """
    exome_ht = add_most_severe_csq_to_tc_within_ht(exome_ht)
    context_ht = add_most_severe_csq_to_tc_within_ht(context_ht)

    context_ht = fast_filter_vep(context_ht).select(
        "context", "ref", "alt", "methylation_level", "exome_coverage"
    )
    context_ht = context_ht.filter(hl.is_defined(context_ht.exome_coverage))

    exome_ht = fast_filter_vep(exome_ht).select(
        "context",
        "ref",
        "alt",
        "methylation_level",
        "exome_coverage",
        "freq",
        "pass_filters",
    )

    grouping = ("exome_coverage",)
    af_cutoff = 0.001

    exome_join = exome_ht[context_ht.key]
    freq_index = exome_ht.freq_index_dict.collect()[0][dataset]

    def keep_criteria(ht):
        crit = (ht.freq[freq_index].AC > 0) & ht.pass_filters
        if impose_high_af_cutoff_upfront:
            crit &= ht.freq[freq_index].AF <= af_cutoff
        return crit

    context_ht = context_ht.filter(
        hl.is_missing(exome_join) | keep_criteria(exome_join)
    )

    exome_ht = exome_ht.filter(keep_criteria(exome_ht))

    possible_file = f"{root}/tmp/possible_coverage.ht"
    if recompute_possible:
        possible_ht = count_variants(
            context_ht, additional_grouping=grouping, force_grouping=True
        )
        possible_ht = annotate_with_mu(possible_ht, mutation_ht)
        possible_ht.transmute(possible_variants=possible_ht.variant_count).write(
            possible_file, True
        )

    possible_ht = hl.read_table(possible_file)
    ht = count_variants(
        exome_ht,
        additional_grouping=grouping,
        partition_hint=100,
        count_downsamplings=POPS,
        impose_high_af_cutoff_here=not impose_high_af_cutoff_upfront,
    )
    ht = ht.join(possible_ht, "outer")
    ht = annotate_variant_types(ht)
    return ht
