# noqa: D100
from typing import List, Tuple, Optional, Union, Any

import hail as hl

# cSpell: disable


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
