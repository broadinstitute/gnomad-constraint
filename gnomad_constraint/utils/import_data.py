"""
Functions to format the external data imported for the regional mutation rate (rr).

The data are imported by ``prepare_data.py --import-dnms`` and
``--import-genomic-features``:

    - De novo mutations (DNMs) from published trio cohorts.
    - The 1kb regional genomic feature matrix from Gnocchi
      (https://github.com/atgu/gnomad_nc_constraint).
"""

from typing import Any, Dict, List

import hail as hl

from gnomad_constraint.resources.constants import (
    GENOMIC_FEATURE_COLUMNS,
    GENOMIC_FEATURE_MISSING,
)


def get_dnm_import_args() -> Dict[str, Any]:
    """
    Return the ``hl.import_table`` arguments for the published DNM files.

    The files are gzip (not bgzip) compressed CSVs that quote fields with commas.

    :return: Keyword arguments for ``hl.import_table``.
    """
    return {"delimiter": ",", "quote": '"', "force": True}


def get_genomic_feature_import_args() -> Dict[str, Any]:
    """
    Return the ``hl.import_table`` arguments for the published genomic feature matrix.

    The feature columns are read as floats, and every spelling of a missing value in
    ``GENOMIC_FEATURE_MISSING`` as missing, so NaN never reaches the rr models.

    :return: Keyword arguments for ``hl.import_table``.
    """
    return {
        "types": {c: hl.tfloat64 for c in GENOMIC_FEATURE_COLUMNS},
        "missing": GENOMIC_FEATURE_MISSING,
    }


def union_dnm_hts(hts: List[hl.Table]) -> hl.Table:
    """
    Combine the imported DNM files into one Table of their ``Variant`` fields.

    :param hts: Imported DNM Tables, each with a ``Variant`` field.
    :return: Unkeyed Table with the ``Variant`` field of every row of ``hts``.
    """
    return hl.Table.union(*[ht.select("Variant") for ht in hts])


def format_dnm_ht(ht: hl.Table) -> hl.Table:
    """
    Key imported de novo mutations (DNMs) by variant, with one row per variant.

    :param ht: Imported DNMs with a GRCh38 ``Variant`` field ("chrom:pos:ref:alt").
    :return: Table keyed by ``locus`` and ``alleles`` with no other fields.
    """
    ht = ht.key_by(**hl.parse_variant(ht.Variant, reference_genome="GRCh38"))

    return ht.select().distinct()


def format_genomic_features(ht: hl.Table) -> hl.Table:
    """
    Key an imported 1kb regional genomic feature matrix by window interval.

    Each window's interval is the BED window in ``element_id`` ("chr-start-end") in
    1-based positions, ``[start + 1, end]``. (The Gnocchi rr Tables used
    ``[start, end)``, 1bp left of the BED window.)

    :param ht: Imported feature matrix with ``element_id`` and
        ``GENOMIC_FEATURE_COLUMNS``.
    :return: Feature Table keyed by ``interval`` with ``element_id`` and
        ``GENOMIC_FEATURE_COLUMNS``.
    """
    coords = ht.element_id.split("-")
    ht = ht.key_by(
        interval=hl.locus_interval(
            coords[0],
            hl.int(coords[1]) + 1,
            hl.int(coords[2]),
            includes_end=True,
            reference_genome="GRCh38",
        )
    )

    return ht.select("element_id", *GENOMIC_FEATURE_COLUMNS)


def check_unique_windows(ht: hl.Table) -> None:
    """
    Raise an error if a window appears more than once in the feature matrix.

    Counting runs the Table's whole plan, so checkpoint it first.

    :param ht: Feature Table keyed by ``interval``, as returned by
        ``format_genomic_features``.
    :return: None.
    """
    n_duplicated = ht.count() - ht.distinct().count()
    if n_duplicated:
        raise ValueError(f"{n_duplicated} windows are duplicated!")
