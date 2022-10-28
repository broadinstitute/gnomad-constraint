"""Script containing utility functions used in the constraint pipeline."""
# cSpell: disable
from typing import Tuple

import hail as hl
from gnomad.utils.constraint import (
    annotate_constraint_groupings,
    annotate_mutation_type,
    annotate_with_mu,
    collapse_strand,
    count_variants_by_group,
    trimer_from_heptamer,
)
from gnomad.utils.vep import (
    add_most_severe_csq_to_tc_within_vep_root,
    filter_vep_transcript_csqs,
)

from gnomad_constraint.resources.resource_utils import HIGH_COVERAGE_CUTOFF


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


def apply_models(
    exome_ht: hl.Table,
    context_ht: hl.Table,
    mutation_ht: hl.Table,
    plateau_models: hl.StructExpression,
    coverage_model: Tuple[float, float],
    max_af: float = 0.001,
    keep_annotations: Tuple[str] = (
        "context",
        "ref",
        "alt",
        "methylation_level",
    ),
    pops: Tuple[str] = (),
    partition_hint: int = 2000,
    custom_model: str = None,
    cov_cutoff: int = COVERAGE_CUTOFF,
) -> hl.Table:
    """
    Compute the expected number of variants and observed:expected ratio using plateau models and coverage model.

    This function sums the number of possible variants times the mutation rate for all variants, and applies the calibration
    model separately for CpG transitions and other sites. For sites with coverage lower than the coverage cutoff, the value obtained
    from the previous step is multiplied by the coverage correction factor. These values are summed across the set of variants
    of interest to obtain the expected number of variants.

    A brief view of how to get the expected number of variants:
        mu_agg = the number of possible variants * the mutation rate (all variants)
        predicted_proportion_observed = sum(plateau model slop * mu_agg + plateau model intercept) (separately for CpG transitions and other sites)
        if 0 < coverage < coverage cutoff:
            coverage_correction = coverage_model slope * log10(coverage) + coverage_model intercept
            expected_variants = sum(predicted_proportion_observed * coverage_correction)
        else:
            expected_variants = sum(predicted_proportion_observed)
        The expected_variants are summed across the set of variants of interest to obtain the final expected number of variants.

    Function adds the following annotations:
        - observed_variants - observed variant counts annotated by `count_variants` function grouped by groupings (output of `annotate_constraint_groupings`)
        - predicted_proportion_observed (including those for each population) - the sum of mutation rate adjusted by plateau models and possible variant counts grouped by groupings
        - possible_variants (including those for each population) - the sum of possible variant counts derived from the context Table grouped by groupings
        - expected_variants (including those for each population) - the sum of expected variant counts grouped by groupings
        - mu - sum(mu_snp * possible_variant * coverage_correction) grouped by groupings
        - obs_exp - observed:expected ratio
        - annotations annotated by `annotate_constraint_groupings`

    :param exome_ht: Exome site Table (output of `prepare_ht_for_constraint_calculations()`) filtered to autosomes
      and pseudoautosomal regions.
    :param context_ht: Context Table (output of `prepare_ht_for_constraint_calculations()`) filtered to autosomes and
      pseudoautosomal regions.
    :param mutation_ht: Mutation rate Table with 'mu_snp' field.
    :param plateau_models: Linear models (output of `build_models() in gnomad_methods`), with
      the values of the dictionary formatted as a StrucExpression of intercept and slope, that calibrates mutation rate to proportion observed for high coverage exome. It includes models for CpG site, non-CpG site, and each population in `POPS`.
    :param coverage_model: A linear model (output of `build_models()` in gnomad_methods), formatted
      as a Tuple of intercept and slope, that calibrates a given coverage level to observed:expected ratio. It's a correction factor for low coverage sites.
    :param max_af: Maximum allele frequency for a variant to be included in returned
      counts. Default is 0.001.
    :param keep_annotations: Annotations to keep in the context Table.
    :param pops: List of populations to use for downsampling counts. Default is ().
    :param partition_hint: Target number of partitions for aggregation. Default is 100.
    :param custom_model: The customized model (one of "standard" or "worst_csq" for
      now), defaults to None.
    :param cov_cutoff: Median coverage cutoff. Sites with coverage above this cutoff
      are considered well covered and will be used to build plateau models. Sites below
      this cutoff have low coverage and will be used to build coverage models. Defaults
      to `COVERAGE_CUTOFF`.
    :return: Table with `expected_variants` (expected variant counts) and `obs_exp`
      (observed:expected ratio) annotations.
    """
    # Add necessary constraint annotations for grouping
    if custom_model == "worst_csq":
        (
            context_ht,
            _,
        ) = annotate_constraint_groupings(context_ht, "worst_csq_by_gene")
        exome_ht, grouping = annotate_constraint_groupings(
            exome_ht, "worst_csq_by_gene"
        )

    else:
        context_ht, _ = annotate_constraint_groupings(
            context_ht, "transcript_consequences"
        )
        exome_ht, grouping = annotate_constraint_groupings(
            exome_ht, "transcript_consequences"
        )

    # Allele frequency information for high-quality genotypes (GQ >= 20; DP >= 10; and
    # AB >= 0.2 for heterozygous calls) in all release samples in gnomAD.
    freq_expr = exome_ht.freq[0]

    # Set up the criteria to exclude variants not observed in the dataset, low-quality
    # variants, and variants with allele frequency above the `max_af` cutoff.
    keep_criteria = (
        (freq_expr.AC > 0)
        & exome_ht.pass_filters
        & (freq_expr.AF <= max_af)
        & (exome_ht.coverage > 0)
    )
    keep_annotations += grouping

    # Keep variants that satisfy the criteria above.
    # Count the observed variants in the entire Table and in each downsampling grouped
    # by `keep_annotations`.
    observed_ht = count_variants_by_group(
        exome_ht.filter(keep_criteria).select(*list(keep_annotations) + ["freq"]),
        additional_grouping=grouping,
        partition_hint=partition_hint,
        count_downsamplings=pops,
        use_table_group_by=True,
    )
    observed_ht = observed_ht.transmute(observed_variants=observed_ht.variant_count)

    # Filter the `exome_ht` to rows that don’t match the criteria above.
    # Anti join the `context_ht` with filtered `exome_ht`, so that `context_ht` only
    # has rows that match the criteria above in the `exome_ht` or are never in
    # the `exome_ht`.
    context_ht = context_ht.select(*keep_annotations).anti_join(
        exome_ht.filter(keep_criteria, keep=False)
    )

    # Count the possible variants in the context Table grouped by `keep_annotations`.
    possible_ht = count_variants_by_group(
        context_ht,
        additional_grouping=grouping,
        partition_hint=partition_hint,
        use_table_group_by=True,
    )
    possible_ht = annotate_with_mu(possible_ht, mutation_ht)
    possible_ht = possible_ht.transmute(possible_variants=possible_ht.variant_count)

    # Apply plateau models and coverage model.
    possible_ht = annotate_mutation_type(
        possible_ht.annotate(mu_agg=possible_ht.mu_snp * possible_ht.possible_variants)
    )
    # Apply plateau and coverage models on the entire dataset.
    total_plateau_model = hl.literal(plateau_models.total)[possible_ht.cpg]
    ann_expr = {
        "predicted_proportion_observed": possible_ht.mu_agg * total_plateau_model[1]
        + total_plateau_model[0],
        "coverage_correction": hl.case()
        .when(possible_ht.coverage == 0, 0)
        .when(possible_ht.coverage >= cov_cutoff, 1)
        .default(
            coverage_model[1] * hl.log10(possible_ht.coverage) + coverage_model[0]
        ),
    }
    # Apply plateau models for specified populations.
    for pop in pops:
        pop_model = hl.literal(plateau_models[pop])
        slopes = hl.map(lambda f: f[possible_ht.cpg][1], pop_model)
        intercepts = hl.map(lambda f: f[possible_ht.cpg][0], pop_model)
        ann_expr[f"predicted_proportion_observed_{pop}"] = (
            possible_ht.mu_agg * slopes + intercepts
        )
    possible_ht = possible_ht.annotate(**ann_expr)
    # Compute expected variant counts on the entire dataset.
    ann_expr = {
        "expected_variants": possible_ht.predicted_proportion_observed
        * possible_ht.coverage_correction,
        "mu": possible_ht.mu_agg * possible_ht.coverage_correction,
    }
    # Compute expected variant counts for specified populations.
    for pop in pops:
        ann_expr[f"expected_variants_{pop}"] = (
            possible_ht[f"predicted_proportion_observed_{pop}"]
            * possible_ht.coverage_correction
        )
    possible_ht = possible_ht.annotate(**ann_expr)

    # Outer join the Tables with possible variant counts and observed variant counts.
    ht = observed_ht.join(possible_ht, "outer")

    grouping.remove("coverage")

    # Generate sum aggregators for 'observed_variants',
    # 'predicted_proportion_observed', 'possible_variants', 'expected_variants', and
    # 'mu' on the entire dataset.
    agg_expr = {
        "observed_variants": hl.agg.sum(ht.observed_variants),
        "predicted_proportion_observed": hl.agg.sum(ht.predicted_proportion_observed),
        "possible_variants": hl.agg.sum(ht.possible_variants),
        "expected_variants": hl.agg.sum(ht.expected_variants),
        "mu": hl.agg.sum(ht.mu),
    }
    # Generate sum aggregators for 'observed_variants',
    # 'predicted_proportion_observed', 'possible_variants', 'expected_variants', and
    # 'mu' for specified populations.
    for pop in pops:
        agg_expr[f"predicted_proportion_observed_{pop}"] = hl.agg.array_sum(
            ht[f"predicted_proportion_observed_{pop}"]
        )
        agg_expr[f"expected_variants_{pop}"] = hl.agg.array_sum(
            ht[f"expected_variants_{pop}"]
        )
        agg_expr[f"downsampling_counts_{pop}"] = hl.agg.array_sum(
            ht[f"downsampling_counts_{pop}"]
        )
    # Aggregate the sum aggregaters grouped by `grouping`
    ht = ht.group_by(*grouping).partition_hint(1000).aggregate(**agg_expr)

    # Annotate global annotations.
    ht = ht.annotate_globals(
        apply_model_params=hl.struct(
            max_af=max_af,
            pops=pops,
            plateau_models=plateau_models,
            coverage_model=coverage_model,
        )
    )
    # Compute the observed:expected ratio.
    return ht.annotate(obs_exp=ht.observed_variants / ht.expected_variants)
