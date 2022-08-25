# noqa: D100
from typing import List, Tuple, Optional, Union, Any

import hail as hl

# cSpell: disable
def fast_filter_vep(
    t: hl.Table,
    vep_root: str = "vep",
    syn: bool = True,
    canonical: bool = True,
    filter_empty: bool = True,
) -> hl.Table:
    """
    Filter to variants where 'transcript_consequences' within the VEP annotation is not empty.

    Also filter to variants where 'most_severe_consequence' is 'synonymous' and the transcript is the
    canonical transcript, if 'syn' and 'canonical' parameter are set to True, respectively.

    :param t: Input Table or MatrixTable.
    :param vep_root: Name used for VEP annotation. Defaults to 'vep'.
    :param syn: Whether to filter to variants where the most severe consequence is "synonymous". Defaults to True.
    :param canonical: Whether to filter to only canonical transcripts. Defaults to True.
    :param filter_empty: Whether to filter out rows where 'transcript_consequences' is empty. Defaults to True.
    :return: Table or MatrixTable filtered to specified criteria.
    """
    transcript_csqs = t[vep_root].transcript_consequences
    criteria = [lambda csq: True]
    if syn:
        criteria.append(lambda csq: csq.most_severe_consequence == "synonymous_variant")
    if canonical:
        criteria.append(lambda csq: csq.canonical == 1)
    transcript_csqs = transcript_csqs.filter(lambda x: combine_functions(criteria, x))
    vep_data = t[vep_root].annotate(transcript_consequences=transcript_csqs)
    t = (
        t.annotate_rows(**{vep_root: vep_data})
        if isinstance(t, hl.MatrixTable)
        else t.annotate(**{vep_root: vep_data})
    )
    if not filter_empty:
        return t
    criteria = hl.is_defined(t.vep.transcript_consequences) & (
        hl.len(t.vep.transcript_consequences) > 0
    )
    return (
        t.filter_rows(criteria) if isinstance(t, hl.MatrixTable) else t.filter(criteria)
    )


def combine_functions(func_list, x, operator="and"):
    """Return booleans."""
    possible_operators = ("and", "or")
    if operator not in possible_operators:
        raise ValueError(
            f'combine_functions only allows operators: {", ".join(possible_operators)}'
        )
    cond = func_list[0](x)
    for c in func_list[1:]:
        if operator == "and":
            cond &= c(x)
        elif operator == "or":
            cond |= c(x)
    return cond


def annotate_with_mu(
    ht: hl.Table,
    mutation_ht: hl.Table,
    output_loc: str = "mu_snp",
    keys: Tuple[str] = ("context", "ref", "alt", "methylation_level"),
) -> hl.Table:
    """
    Annotate SNP mutation rate for the input Table.

    :param ht: Input Table.
    :param mutation_ht: Mutation rate Table.
    :param output_loc: Name for mutational rate annotation. Defaults to 'mu_snp'.
    :param keys: Common keys between mutation rate Table and input Table. Defaults to ('context', 'ref', 'alt', 'methylation_level').
    :return: Table with mutational rate annotation added (default name for annotation is 'mu_snp').
    """
    mu = hl.literal(
        mutation_ht.aggregate(
            hl.dict(
                hl.agg.collect(
                    (hl.struct(**{k: mutation_ht[k] for k in keys}), mutation_ht.mu_snp)
                )
            )
        )
    )
    mu = mu.get(hl.struct(**{k: ht[k] for k in keys}))
    return ht.annotate(
        **{output_loc: hl.case().when(hl.is_defined(mu), mu).or_error("Missing mu")}
    )


def count_variants(
    ht: hl.Table,
    count_singletons: bool = False,
    count_downsamplings: Optional[List[str]] = (),
    additional_grouping: Optional[List[str]] = (),
    partition_hint: int = 100,
    omit_methylation: bool = False,
    return_type_only: bool = False,
    force_grouping: bool = False,
    singleton_expression: hl.expr.BooleanExpression = None,
    impose_high_af_cutoff_here: bool = False,
    af_cutoff=0.001,
) -> Union[hl.Table, Any]:
    """
    Count number of observed variants by context, ref, alt, methylation_level.

    :param ht: Input Hail Table.
    :param count_singletons: Whether to count singletons. Defaults to False.
    :param count_downsamplings: List of populations to use for downsampling counts. Defaults to ().
    :param additional_grouping: Additional features to group by. i.e. exome_coverage. Defaults to ().
    :param partition_hint: Target number of partitions for aggregation. Defaults to 100.
    :param omit_methylation: Whether to omit 'methylation_level' from the grouping when counting variants. Defaults to False.
    :param return_type_only: Whether to only return the data type of 'variant_count'. Defaults to False.
    :param force_grouping: Whether to force grouping. Defaults to False.
    :param singleton_expression: Expression for defining a singleton. Defaults to None.
    :param impose_high_af_cutoff_here: Whether to filter to variants with an AF <= 0.001. Defaults to False.
    :param af_cutoff: Variants with AF above than AF cutoff will be removed, defaults to 0.001.
    :return: Table including 'variant_count' and downsampling counts if requested.
    """
    grouping = hl.struct(context=ht.context, ref=ht.ref, alt=ht.alt)
    if not omit_methylation:
        grouping = grouping.annotate(methylation_level=ht.methylation_level)
    for group in additional_grouping:
        grouping = grouping.annotate(**{group: ht[group]})

    if count_singletons:
        # singleton = hl.any(lambda f: (f.meta.size() == 1) & (f.meta.get('group') == 'adj') & (f.AC[1] == 1), ht.freq)
        if singleton_expression is None:
            singleton_expression = ht.freq[0].AC == 1

    if count_downsamplings or force_grouping:
        # Slower, but more flexible (allows for downsampling agg's)
        output = {
            "variant_count": hl.agg.count_where(ht.freq[0].AF <= af_cutoff)
            if impose_high_af_cutoff_here
            else hl.agg.count()
        }
        for pop in count_downsamplings:
            output[f"downsampling_counts_{pop}"] = downsampling_counts_expr(
                ht, pop, impose_high_af_cutoff=impose_high_af_cutoff_here
            )
        if count_singletons:
            output["singleton_count"] = hl.agg.count_where(singleton_expression)
            for pop in count_downsamplings:
                output[
                    f"singleton_downsampling_counts_{pop}"
                ] = downsampling_counts_expr(ht, pop, singleton=True)
        return (
            ht.group_by(**grouping).partition_hint(partition_hint).aggregate(**output)
        )
    else:
        agg = {"variant_count": hl.agg.counter(grouping)}
        if count_singletons:
            agg["singleton_count"] = hl.agg.counter(
                hl.agg.filter(singleton_expression, grouping)
            )

        if return_type_only:
            return agg["variant_count"].dtype
        else:
            return ht.aggregate(hl.struct(**agg))


def downsampling_counts_expr(
    ht: Union[hl.Table, hl.MatrixTable],
    pop: str = "global",
    variant_quality: str = "adj",
    singleton: bool = False,
    impose_high_af_cutoff: bool = False,
) -> hl.expr.ArrayExpression:
    """
    Downsample the variant count per given population.

    :param ht: Input Table.
    :param pop: Population. Defaults to 'global'.
    :param variant_quality: Variant quality for "group" key. Defaults to 'adj'.
    :param singleton: Whether to sum only alleles that are singletons. Defaults to False.
    :param impose_high_af_cutoff: Whether to sum only alleles with an allele frequency less than or equal to 0.001. Defaults to False.
    :return: Downsampling count for specified population.
    """
    indices = hl.zip_with_index(ht.freq_meta).filter(
        lambda f: (f[1].size() == 3)
        & (f[1].get("group") == variant_quality)
        & (f[1].get("pop") == pop)
        & f[1].contains("downsampling")
    )
    sorted_indices = hl.sorted(indices, key=lambda f: hl.int(f[1]["downsampling"])).map(
        lambda x: x[0]
    )

    def get_criteria(i):
        if singleton:
            return hl.int(ht.freq[i].AC == 1)
        elif impose_high_af_cutoff:
            return hl.int((ht.freq[i].AC > 0) & (ht.freq[i].AF <= 0.001))
        else:
            return hl.int(ht.freq[i].AC > 0)

    return hl.agg.array_sum(hl.map(get_criteria, sorted_indices))


def remove_unnecessary_variants(
    context_ht: hl.Table,
    exome_ht: hl.Table,
    dataset: str = "gnomad",
    af_cutoff=0.001,
    impose_high_af_cutoff_upfront: bool = True,
) -> Tuple(hl.Table, hl.Table):
    """
    Remove unnecessary variants in context Table and exome Table.

    Variants with no coverage, a low-quality variant, or a variant with AF above 0.1% will be removed.

    :param context_ht: The context Table.
    :param exome_ht: The exome Table.
    :param dataset: Dataset to use when computing frequency index. Defaults to 'gnomad'.
    :param af_cutoff: Variants with AF above than AF cutoff will be removed, defaults to 0.001.
    :param impose_high_af_cutoff_upfront: Whether to remove high frequency alleles. Defaults to True.
    :return: The context and exome Table with necessary variants.
    """
    exome_join = exome_ht[context_ht.key]
    freq_index = exome_ht.freq_index_dict.collect()[0][dataset]

    crit = (ht.freq[freq_index].AC > 0) & ht.pass_filters & (ht.coverage > 0)
    if impose_high_af_cutoff_upfront:
        crit &= ht.freq[freq_index].AF <= af_cutoff

    context_ht = context_ht.filter(hl.is_missing(exome_join) | crit)

    exome_ht = exome_ht.filter(crit)
    return context_ht, exome_ht
