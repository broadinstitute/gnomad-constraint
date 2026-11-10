"""
Functions to format the external data imported for the regional mutation rate (rr).

The data are imported by ``import_data.py``:

    - De novo mutations (DNMs) from published trio cohorts.
    - The 1kb regional genomic feature matrix from Gnocchi
      (https://github.com/atgu/gnomad_nc_constraint).
"""

import hail as hl

from gnomad_constraint.resources.constants import GENOMIC_FEATURE_COLUMNS


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
