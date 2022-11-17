"""Script containing utility functions used in the constraint pipeline."""
# cSpell: disable
import pickle
from typing import List, Optional, Tuple

import hail as hl
from gnomad.utils.constraint import (
    annotate_mutation_type,
    annotate_with_mu,
    collapse_strand,
    count_variants_by_group,
    get_downsamplings,
    trimer_from_heptamer,
)
from gnomad.utils.filtering import (
    filter_for_mu,
    filter_to_autosomes,
    remove_coverage_outliers,
)
from gnomad.utils.vep import (
    add_most_severe_csq_to_tc_within_vep_root,
    filter_vep_transcript_csqs,
)
from hail.utils.misc import new_temp_file

from gnomad_constraint.resources.resource_utils import mutation_rate_ht


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
    raw_context_ht: hl.Table,
    recalculate_all_possible_summary: bool = True,
    recalculate_all_possible_summary_unfiltered: bool = False,
    omit_methylation: bool = False,
    count_singletons: bool = False,
    summary_file: str = None,
    keep_annotations: Tuple[str] = (
        "context",
        "ref",
        "alt",
        "methylation_level",
    ),
    ac_cutoff: int = 5,
    downsampling_level: int = 1000,
    total_mu: float = 1.2e-08,
    pops: Tuple[str] = (),
) -> hl.Table:
    """
    Calculate mutation rate using the downsampling with size specified by `downsampling_level` in genome sites Table.

    Prior to computing mutation rate the only following variants are kept:
        - variants with the mean coverage in the gnomAD genomes was between 15X and 60X.
        - variants whoes the most severe annotation was intron_variant or
            intergenic_variant
        - variants with the GERP score between the 5th and 95th percentile of the
            genomewide distribution.
        - high-quality variants: `exome_ht.pass_filters`
        - Variants with allele count below `ac_cutoff`: `(freq_expr.AC <= ac_cutoff)`

    The returned Table includes the following annotations:
        - context - trinucleotide genomic context
        - ref - the reference allele
        - alt - the alternate base
        - methylation_level - methylation_level
        - downsampling_counts_{pop} - variant counts in downsamplings for populations
          in `pops`
        - mu_snp - SNP mutation rate
        - annotations added by `annotate_mutation_type`

    :param genome_ht: Genome sites Table.
    :param raw_context_ht: Context Table with locus that is on an autosome or in a
        pseudoautosomal region and sex chromosomes.
    :param recalculate_all_possible_summary: Whether to calculate possible variants
        using context Table with locus that is only on an autosome or in a
        pseudoautosomal region. Default is True.
    :param recalculate_all_possible_summary_unfiltered: Whether to calculate possible
        variants using raw context Table. Default is True.
    :param omit_methylation: Whether to omit 'methylation_level' from the grouping when
        counting variants. Default is False.
    :param count_singletons: Whether to count singletons. Default is False.
    :param summary_file: _description_, defaults to None
    :param keep_annotations: Annotations to keep in the context Table and genome sites
        Table.
    :param ac_cutoff: The cutoff of allele count when filtering context Table and genome sites Table.
    :param downsampling_level: The size of downsamplings will be used to count variants. Default is 1000.
    :param total_mu: The per-generation mutation rate. Default is 1.2e-08.
    :param pops: List of populations to use for downsampling counts. Default is ().
    :return: Mutation rate Table.
    """
    # keep only loci where the mean coverage in the gnomAD genomes was between 15X and
    # 60X.
    context_ht = filter_to_autosomes(remove_coverage_outliers(raw_context_ht))
    ## Is it necessary to call `filter_to_autosomes()` again here?
    genome_ht = filter_to_autosomes(remove_coverage_outliers(genome_ht))

    # filter the Table so that the most severe annotation was intron_variant or
    # intergenic_variant, and that the GERP score was between the 5th and 95th
    # percentile of the genomewide distribution.
    context_ht = filter_for_mu(context_ht)
    genome_ht = filter_for_mu(genome_ht)

    context_ht = context_ht.select(*keep_annotations)
    genome_ht = genome_ht.select(*list(keep_annotations) + ["freq", "pass_filters"])

    # downsampled the dataset to 1,000 genomes
    freq_index = hl.eval(
        hl.enumerate(genome_ht.freq_meta).find(
            lambda f: (f[1].get("downsampling") == str(downsampling_level))
            & (f[1].get("pop") == "global")
            & (f[1].get("group") == "adj")
            & (f[1].size() == 3)
        )
    )[0]
    freq_expr = genome_ht.freq[freq_index]
    # Set up the criteria to filtered out  low quality sites, and sites found in
    # greater than 5 copies in the downsampled set.
    keep_criteria = (freq_expr.AC <= ac_cutoff) & genome_ht.pass_filters

    # Only keep variants with AC <= ac_cutoff and passing filters in genome site
    # Table.
    genome_ht = genome_ht.filter(keep_criteria)

    # In the denominator, only keep variants not in the genome dataset, or with AC
    # <= ac_cutoff and passing filters.
    context_ht = context_ht.anti_join(genome_ht.filter(keep_criteria, keep=False))

    if not summary_file:
        ## Should we save possible variants in a random path?
        summary_file = new_temp_file(prefix="constraint", extension="he")

    ## Is it possible to get rid of the usage of dtype here?
    all_possible_dtype = count_variants_by_group(
        context_ht,
        omit_methylation=omit_methylation,
        return_type_only=True,
    )

    # Count possible variants in context Table.
    if recalculate_all_possible_summary:
        all_possible = count_variants_by_group(
            context_ht,
            omit_methylation=omit_methylation,
        ).variant_count
        # with hl.hadoop_open(summary_file, "wb") as f:
        #     pickle.dump(all_possible, f)
        hl.experimental.write_expression(all_possible, summary_file)
    # with hl.hadoop_open(summary_file, "rb") as f:
    #     all_possible = pickle.load(f)
    all_possible = hl.experimental.read_expression(summary_file)

    if recalculate_all_possible_summary_unfiltered:
        all_possible_unfiltered = count_variants_by_group(
            raw_context_ht.rows(), omit_methylation=omit_methylation
        ).variant_count
        all_possible_unfiltered = {
            x: y
            for x, y in list(all_possible_unfiltered.items())
            if x.context is not None
        }
        hl.experimental.write_expression(
            all_possible_unfiltered, new_temp_file(prefix="constraint", extension="he")
        )
        # with hl.hadoop_open(new_temp_file(prefix="constraint", extension="he"), "wb") as f:
        #     pickle.dump(all_possible_unfiltered, f)
    # Uncomment this line and next few commented lines to back-calculate total_mu from the old mutation rate dataset
    # all_possible_unfiltered = load_all_possible_summary(filtered=False)

    # Count the observed variants in the genome sites Table.
    genome_ht = count_variants_by_group(
        genome_ht,
        count_downsamplings=pops,
        count_singletons=count_singletons,
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

    for pop in pops:
        correction_factors = ht.aggregate(
            total_mu
            / (hl.agg.array_sum(ht[f"downsampling_counts_{pop}"]) / total_bases)
        )
        ht = ht.annotate(
            **{
                f"downsamplings_mu_{pop}": hl.literal(correction_factors)
                * ht[f"downsampling_counts_{pop}"]
                / ht.possible_variants
            }
        )

    ht = annotate_mutation_type(
        ht.annotate(
            downsamplings_frac_observed=ht.downsampling_counts_global
            / ht.possible_variants
        )
    )

    # Get the index of dowsamplings with size of 1000 genomes
    downsamplings = list(map(lambda x: x[1], get_downsamplings(ht.freq_meta)))
    index_1kg = downsamplings.index(1000)

    # Compute the absolute mutation rate
    for pop in pops:
        if pop == "global":
            ht = ht.annotate(**{"mu_snp": ht[f"downsamplings_mu_{pop}"][index_1kg]})
        else:
            ht = ht.annotate(
                **{f"mu_snp_{pop}": ht[f"downsamplings_mu_{pop}"][index_1kg]}
            )

    # Compute the proportion observed, which represents the relative mutability of each
    # variant class
    return ht.annotate(
        proportion_observed_1kg=ht.downsampling_counts_global[index_1kg]
        / ht.possible_variants,
        proportion_observed=ht.variant_count / ht.possible_variants,
    )


def get_old_mu_data(version) -> hl.Table:
    """
    Get the mutation rate Table of last version.

    :return: Mutation rate Table.
    """
    if version == "2.1.1":
        old_mu_data = hl.import_table(
            "gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/old_exac_data/fordist_1KG_mutation_rate_table.txt",
            delimiter=" ",
            impute=True,
        )
        return old_mu_data.transmute(
            context=old_mu_data["from"],
            ref=old_mu_data["from"][1],
            alt=old_mu_data.to[1],
        ).key_by("context", "ref", "alt")
    else:
        return mutation_rate_ht.ht()
