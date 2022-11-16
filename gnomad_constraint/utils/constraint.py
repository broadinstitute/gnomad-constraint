"""Script containing utility functions used in the constraint pipeline."""

from typing import Optional, Tuple

import hail as hl
from gnomad.utils.constraint import (
    annotate_mutation_type,
    annotate_with_mu,
    collapse_strand,
    count_variants_by_group,
    trimer_from_heptamer,
)
from gnomad.utils.filtering import (
    filter_by_frequency,
    filter_for_mu,
    filter_to_autosomes,
    remove_coverage_outliers,
)
from gnomad.utils.vep import (
    add_most_severe_csq_to_tc_within_vep_root,
    filter_vep_transcript_csqs,
)


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

    Function drops `a_index`, `was_split`, and`colocated_variants` annotations from
    gnomAD data.

    .. note::
        Function expects that multiallelic variants in the VEP context Table have been
        split.

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
        - annotations added by `annotate_mutation_type()`, `collapse_strand()`, and
          `add_most_severe_csq_to_tc_within_vep_root()`

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
    # Annotate mutation type (such as "CpG", "non-CpG transition", "transversion") and
    # collapse strands to deduplicate the context
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

    # Add most_severe_consequence annotation to 'transcript_consequences' within the
    # vep root annotation.
    ht = add_most_severe_csq_to_tc_within_vep_root(ht)

    # Filter out locus with undefined exome coverage
    ht = ht.filter(hl.is_defined(ht.exome_coverage))
    return ht


def create_constraint_training_dataset(
    exome_ht: hl.Table,
    context_ht: hl.Table,
    mutation_ht: hl.Table,
    max_af: float = 0.001,
    keep_annotations: Tuple[str] = (
        "context",
        "ref",
        "alt",
        "methylation_level",
    ),
    pops: Tuple[str] = (),
    grouping: Tuple[str] = ("exome_coverage",),
    partition_hint: int = 100,
) -> hl.Table:
    """
    Count the observed variants and possible variants by substitution, context, methylation level, and additional `grouping`.

    Prior to computing variant counts the following variants are removed:
        - Variants not observed by any samples in the dataset: `(freq_expr.AC > 0)`
        - Low-quality variants: `exome_ht.pass_filters`
        - Variants with allele frequency above `max_af` cutoff: `(freq_expr.AF <=
          max_af)`
        - Variants that are not synonymous or in the canonical transcript

    For each substitution, context, methylation level, and exome coverage, the rest of
    variants in `exome_ht` are counted and annotated as `observed_variants`, and the
    rest of variants in `context_ht` are counted and annotated as `possible_variants`.
    The final Table is the outer-join of the filtered `exome_ht` and `context_ht` with
    the `observed_variants` and `possible_variants` annotations.

    The returned Table includes the following annotations:
        - context - trinucleotide genomic context
        - ref - the reference allele
        - alt - the alternate base
        - methylation_level - methylation_level
        - observed_variants - observed variant counts in `exome_ht`
        - possible_variants - possible variant counts in `context_ht`
        - downsampling_counts_{pop} - variant counts in downsamplings for populations
          in `pops`
        - mu_snp - SNP mutation rate
        - annotations added by `annotate_mutation_type`

    :param exome_ht: Preprocessed exome Table.
    :param context_ht: Preprocessed context Table.
    :param mutation_ht: Preprocessed mutation rate Table.
    :param max_af: Maximum allele frequency for a variant to be included in returned
        counts. Default is 0.001.
    :param keep_annotations: Annotations to keep in the context Table.
    :param pops: List of populations to use for downsampling counts. Default is ().
    :param grouping: Annotations other than 'context', 'ref', 'alt', and
        `methylation_level` to group by when counting variants. Default is
        ('exome_coverage',).
    :param partition_hint: Target number of partitions for aggregation. Default is 100.
    :return: Table with observed variant and possible variant count.
    """
    # Allele frequency information for high-quality genotypes (GQ >= 20; DP >= 10; and
    # AB >= 0.2 for heterozygous calls) in all release samples in gnomAD.
    freq_expr = exome_ht.freq[0]

    # Set up the criteria to exclude variants not observed in the dataset, low-quality
    # variants, and variants with allele frequency above the `max_af` cutoff.
    keep_criteria = (
        (freq_expr.AC > 0) & exome_ht.pass_filters & (freq_expr.AF <= max_af)
    )
    keep_annotations += grouping

    # Keep variants that are synonymous, canonical, and satisfy the criteria above.
    # Count the observed variants in the entire Table and in each downsampling grouped
    # by `keep_annotations`.
    ht = count_variants_by_group(
        filter_vep_transcript_csqs(exome_ht.filter(keep_criteria)).select(
            *list(keep_annotations) + ["freq"]
        ),
        additional_grouping=grouping,
        partition_hint=partition_hint,
        count_downsamplings=pops,
        use_table_group_by=True,
    )
    ht = ht.transmute(observed_variants=ht.variant_count)

    # Keep variants that are synonymous and canonical in `context_ht`.
    context_ht = filter_vep_transcript_csqs(context_ht).select(*keep_annotations)
    # Filter the `exome_ht` to rows that don’t match the criteria above
    # Anti join the `context_ht` with filtered `exome_ht`, so that `context_ht` only
    # has rows that match the criteria above in the `exome_ht` or are never in
    # the `exome_ht`.
    context_ht = context_ht.anti_join(exome_ht.filter(keep_criteria, keep=False))

    # Count the possible variants in the context Table grouped by `keep_annotations`.
    possible_ht = count_variants_by_group(
        context_ht,
        additional_grouping=grouping,
        partition_hint=partition_hint,
        use_table_group_by=True,
    )
    possible_ht = annotate_with_mu(possible_ht, mutation_ht)
    possible_ht = possible_ht.transmute(possible_variants=possible_ht.variant_count)

    # Outer join the Tables with possible variant counts and observed variant counts.
    ht = ht.join(possible_ht, "outer")
    # Annotate the Table with cpg site and mutation_type (one of "CpG", "non-CpG
    # transition", or "transversion").
    ht = annotate_mutation_type(ht)
    # Annotate parameters used in the function as global annotations.
    ht = ht.annotate_globals(
        training_dataset_params=hl.struct(max_af=max_af, pops=pops)
    )
    return ht


def calculate_mu_by_downsampling(
    genome_ht: hl.Table,
    raw_context_ht: hl.MatrixTable,
    recalculate_all_possible_summary: bool = True,
    recalculate_all_possible_summary_unfiltered: bool = False,
    omit_methylation: bool = False,
    count_singletons: bool = False,
    remove_common_downsampled: bool = True,
    remove_common_ordinary: bool = False,
    grouping_variables: Optional[Tuple[str]] = (),
    summary_file: str = None,
) -> hl.Table:
    """
    _summary_.

    :param genome_ht: _description_
    :param raw_context_ht: _description_
    :param recalculate_all_possible_summary: _description_, defaults to True
    :param recalculate_all_possible_summary_unfiltered: _description_, defaults to False
    :param omit_methylation: _description_, defaults to False
    :param count_singletons: _description_, defaults to False
    :param remove_common_downsampled: _description_, defaults to True
    :param remove_common_ordinary: _description_, defaults to False
    :param grouping_variables: _description_, defaults to ()
    :param summary_file: _description_, defaults to None
    :return: _description_
    """
    context_ht = filter_to_autosomes(remove_coverage_outliers(raw_context_ht))
    genome_ht = filter_to_autosomes(remove_coverage_outliers(genome_ht))

    if grouping_variables:
        context_ht = context_ht.filter(
            hl.all(
                lambda x: x, [hl.is_defined(context_ht[x]) for x in grouping_variables]
            )
        )
        genome_ht = genome_ht.filter(
            hl.all(
                lambda x: x, [hl.is_defined(genome_ht[x]) for x in grouping_variables]
            )
        )
    else:
        context_ht = filter_for_mu(context_ht)
        genome_ht = filter_for_mu(genome_ht)

    context_ht = context_ht.select(
        "context", "ref", "alt", "methylation_level", *grouping_variables
    )
    genome_ht = genome_ht.select(
        "context",
        "ref",
        "alt",
        "methylation_level",
        "freq",
        "pass_filters",
        *grouping_variables,
    )

    genome_join = genome_ht[context_ht.key]
    if remove_common_downsampled:
        ac_cutoff = 5
        downsampling_level = 1000
        freq_index = hl.eval(
            hl.zip_with_index(genome_ht.freq_meta).find(
                lambda f: (f[1].get("downsampling") == str(downsampling_level))
                & (f[1].get("pop") == "global")
                & (f[1].get("group") == "adj")
                & (f[1].size() == 3)
            )
        )[0]
        # In the denominator, only keep variants not in the genome dataset, or with AC <= ac_cutoff and passing filters
        context_ht = context_ht.filter(
            hl.is_missing(genome_join)
            | (
                (genome_join.freq[freq_index].AC <= ac_cutoff)
                & genome_join.pass_filters
            )
        )
        # Keep only AC <= ac_cutoff in numerator
        genome_ht = genome_ht.filter(genome_ht.freq[freq_index].AC <= ac_cutoff)
    elif remove_common_ordinary:
        af_cutoff = 0.001
        context_ht = context_ht.filter(
            hl.is_missing(genome_join)
            | ((genome_join.freq[0].AF <= af_cutoff) & genome_join.pass_filters)
        )
        # Keep only AF < af_cutoff in numerator
        genome_ht = filter_by_frequency(
            genome_ht, "below", frequency=af_cutoff, adj=True
        )
    else:
        context_ht = context_ht.filter(
            hl.is_missing(genome_join) | genome_join.pass_filters
        )
    genome_ht = genome_ht.filter(genome_ht.pass_filters)

    if not summary_file:
        summary_file = all_possible_summary_pickle
    all_possible_dtype = count_variants_by_group(
        context_ht,
        omit_methylation=omit_methylation,
        additional_grouping=grouping_variables,
        return_type_only=True,
    )
    if recalculate_all_possible_summary:
        all_possible = count_variants_by_group(
            context_ht,
            omit_methylation=omit_methylation,
            additional_grouping=grouping_variables,
        ).variant_count
        with hl.hadoop_open(summary_file, "wb") as f:
            pickle.dump(all_possible, f)
    with hl.hadoop_open(summary_file, "rb") as f:
        all_possible = pickle.load(f)

    if recalculate_all_possible_summary_unfiltered:
        all_possible_unfiltered = count_variants_by_group(
            raw_context_ht.rows(), omit_methylation=omit_methylation
        ).variant_count
        all_possible_unfiltered = {
            x: y
            for x, y in list(all_possible_unfiltered.items())
            if x.context is not None
        }
        with hl.hadoop_open(all_possible_summary_unfiltered_pickle, "wb") as f:
            pickle.dump(all_possible_unfiltered, f)
    # Uncomment this line and next few commented lines to back-calculate total_mu from the old mutation rate dataset
    # all_possible_unfiltered = load_all_possible_summary(filtered=False)

    genome_ht = count_variants_by_group(
        genome_ht,
        count_downsamplings=["global", "nfe", "afr"],
        count_singletons=count_singletons,
        additional_grouping=grouping_variables,
        omit_methylation=omit_methylation,
    )

    old_mu_data = get_old_mu_data()
    ht = genome_ht.annotate(
        possible_variants=hl.literal(all_possible, dtype=all_possible_dtype)[
            genome_ht.key
        ],
        # possible_variants_unfiltered=hl.literal(all_possible_unfiltered, dtype=all_possible_dtype)[genome_ht.key],
        old_mu_snp=old_mu_data[
            hl.struct(context=genome_ht.context, ref=genome_ht.ref, alt=genome_ht.alt)
        ].mu_snp,
    )
    ht = ht.persist()
    total_bases = ht.aggregate(hl.agg.sum(ht.possible_variants)) // 3
    # total_bases_unfiltered = ht.aggregate(hl.agg.sum(ht.possible_variants_unfiltered)) // 3
    # total_mu = ht.aggregate(hl.agg.sum(ht.old_mu_snp * ht.possible_variants_unfiltered) / total_bases_unfiltered)
    total_mu = 1.2e-08

    correction_factors = ht.aggregate(
        total_mu / (hl.agg.array_sum(ht.downsampling_counts_global) / total_bases)
    )
    correction_factors_nfe = ht.aggregate(
        total_mu / (hl.agg.array_sum(ht.downsampling_counts_nfe) / total_bases)
    )
    correction_factors_afr = ht.aggregate(
        total_mu / (hl.agg.array_sum(ht.downsampling_counts_afr) / total_bases)
    )

    ht = annotate_mutation_type(
        ht.annotate(
            downsamplings_frac_observed=ht.downsampling_counts_global
            / ht.possible_variants,
            downsamplings_mu_snp=hl.literal(correction_factors)
            * ht.downsampling_counts_global
            / ht.possible_variants,
            downsamplings_mu_nfe=hl.literal(correction_factors_nfe)
            * ht.downsampling_counts_nfe
            / ht.possible_variants,
            downsamplings_mu_afr=hl.literal(correction_factors_afr)
            * ht.downsampling_counts_afr
            / ht.possible_variants,
        )
    )
    downsamplings = list(map(lambda x: x[1], get_downsamplings(ht)))
    index_1kg = downsamplings.index(1000)
    return ht.annotate(
        proportion_observed_1kg=ht.downsampling_counts_global[index_1kg]
        / ht.possible_variants,
        proportion_observed=ht.variant_count / ht.possible_variants,
        mu_snp=ht.downsamplings_mu_snp[index_1kg],
        mu_snp_nfe=ht.downsamplings_mu_nfe[index_1kg],
        mu_snp_afr=ht.downsamplings_mu_afr[index_1kg],
    )


# Data loads
def get_old_mu_data() -> hl.Table:
    """
    _summary_.

    :return: _description_
    """
    old_mu_data = hl.import_table(
        f"{root}/old_exac_data/fordist_1KG_mutation_rate_table.txt",
        delimiter=" ",
        impute=True,
    )
    return old_mu_data.transmute(
        context=old_mu_data["from"], ref=old_mu_data["from"][1], alt=old_mu_data.to[1]
    ).key_by("context", "ref", "alt")


def get_downsamplings(ht):
    """
    _summary_.

    :param ht: _description_
    :return: _description_
    """
    freq_meta = ht.freq_meta.collect()[0]
    downsamplings = [
        (i, int(x.get("downsampling")))
        for i, x in enumerate(freq_meta)
        if x.get("group") == "adj"
        and x.get("pop") == "global"
        and x.get("downsampling") is not None
    ]
    return downsamplings
