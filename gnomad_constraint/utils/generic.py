# noqa: D100
from typing import List, Tuple, Optional, Union, Any

import hail as hl

# cSpell: disable


def get_all_pop_lengths(
    ht, pops: Tuple[str], prefix: str = "observed_", skip_assertion: bool = False
):
    """
    Get the minimum array length for specific per population annotations in `ht`.

    The annotations are specified by the combination of `prefix` and each population in `pops`.

    :param ht: Input Table.
    :param pops: List of populations. Defaults to `POPS`.
    :param prefix: Prefix of population variant count. Defaults to 'observed_'.
    :param skip_assertion: Whether to skip raising an AssertionError if all the arrays of variant counts within a population don't have the same length. Defaults to False.
    :return: A Dictionary with the minimum array length for each population.
    """
    ds_lengths = ht.aggregate(
        [hl.agg.min(hl.len(ht[f"{prefix}{pop}"])) for pop in pops]
    )
    # temp_ht = ht.take(1)[0]
    # ds_lengths = [len(temp_ht[f'{prefix}{pop}']) for pop in pops]
    pop_lengths = list(zip(ds_lengths, pops))
    print("Found: ", pop_lengths)
    if not skip_assertion:
        assert ht.all(
            hl.all(
                lambda f: f,
                [hl.len(ht[f"{prefix}{pop}"]) == length for length, pop in pop_lengths],
            )
        )
    return pop_lengths
