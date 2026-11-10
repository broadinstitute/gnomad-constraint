"""
Utility functions to format external data imported for the constraint pipeline.

The data are imported by ``import_data.py``:

    - De novo mutations (DNMs) from published trio cohorts, annotated with the
      constraint pipeline's context fields.
    - The 1kb regional genomic feature matrix from Gnocchi
      (https://github.com/atgu/gnomad_nc_constraint).
"""

import logging
from typing import Any, Dict

import hail as hl
from gnomad.utils.liftover import default_lift_data
from hail.utils.misc import new_temp_file

from gnomad_constraint.resources.constants import GENOMIC_FEATURE_COLUMNS

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("import_data_utils")
logger.setLevel(logging.INFO)


# ---------------------------------------------------------------------------
# De novo mutations
# ---------------------------------------------------------------------------


def format_dnm_cohort(cohort: str, ht: hl.Table) -> hl.Table:
    """
    Format one imported DNM cohort as a GRCh38 Table keyed by ``locus`` and ``alleles``.

    Supported cohorts (see
    :func:`gnomad_constraint.resources.resource_utils.get_dnm_cohort_resources`):

        - ``decode_psychencode``: deCODE + PsychENCODE DNMs used by Gnocchi
          (GRCh38). Each DNM is listed in both allele orientations; the wrong
          orientation is removed when joining to the context Table.
        - ``satterstrom_2026_asd``, ``satterstrom_2026_ddd``: Satterstrom 2026 ASD
          and DD DNMs (GRCh38, with a "chrom:pos:ref:alt" ``Variant`` column).
          Indels are removed when joining to the context Table.
        - ``gnomad_v4.1``: gnomAD v4.1 high-quality coding de novos (GRCh38 Hail
          Table).
        - ``chd_2017``, ``ee_2014``, ``cdh_2021``: GRCh37, lifted over to GRCh38.
        - ``ofc_2020``: GRCh38 positions.

    :param cohort: Cohort name.
    :param ht: Imported cohort Table.
    :return: Table keyed by ``locus`` and ``alleles`` with a ``cohort`` field.
    """
    if cohort == "gnomad_v4.1":
        ht = ht.key_by("locus", "alleles").select()
    elif cohort in {"satterstrom_2026_asd", "satterstrom_2026_ddd"}:
        # TODO: Unlike the gnomAD de novos, these calls aren't QC-filtered. The ASD
        #  file has Confidence and Role columns and the DDD file has LCR and
        #  gnomad_non_neuro_AF columns that a filter could use.
        variant = ht.Variant.split(":")
        ht = ht.select(
            locus=hl.locus(variant[0], hl.int(variant[1]), reference_genome="GRCh38"),
            alleles=variant[2:],
        ).key_by("locus", "alleles")
    elif cohort == "decode_psychencode":
        ht = ht.select(
            locus=hl.parse_locus(ht.locus, reference_genome="GRCh38"),
            alleles=hl.array([ht.ref, ht.alt]),
        ).key_by("locus", "alleles")
    elif cohort == "chd_2017":
        ht = ht.select(
            locus=hl.parse_locus(
                ht.CHROM + ":" + hl.str(ht.POS), reference_genome="GRCh37"
            ),
            alleles=hl.array([ht.REF, ht.ALT]),
        )
        ht = default_lift_data(ht.key_by("locus", "alleles")).select()
    elif cohort == "ofc_2020":
        ht = ht.select(
            locus=hl.parse_locus(
                "chr" + hl.str(ht.Chromosome) + ":" + hl.str(ht["Position (build 38)"]),
                reference_genome="GRCh38",
            ),
            alleles=hl.array([ht.Ref, ht.Alt]),
        ).key_by("locus", "alleles")
    elif cohort == "ee_2014":
        ht = ht.select(
            locus=hl.parse_locus(
                ht["hg19 coordinates (chr:position)"].replace("chr", ""),
                reference_genome="GRCh37",
            ),
            alleles=ht["Ref/Alt alleles"].split("/"),
        )
        ht = default_lift_data(ht.key_by("locus", "alleles")).select()
    elif cohort == "cdh_2021":
        ht = ht.select(
            locus=hl.parse_locus(
                hl.str(ht.chr) + ":" + hl.str(ht.position), reference_genome="GRCh37"
            ),
            alleles=hl.array([ht.ref, ht.alt]),
        )
        ht = default_lift_data(ht.key_by("locus", "alleles")).select()
    else:
        raise ValueError(f"Unknown DNM cohort: {cohort}!")

    return ht.annotate(cohort=cohort)


def combine_dnms(cohort_hts: Dict[str, hl.Table], context_ht: hl.Table) -> hl.Table:
    """
    Combine DNM cohorts into one Table of SNV DNMs with their context annotations.

    Only DNMs whose ``locus`` and ``alleles`` are in the annotated context Table are
    kept, which removes indels and alleles that don't match the reference. Each DNM
    site gets one row, with the number of DNMs reported at the site by each cohort,
    and is annotated with the strand-collapsed ``context``, ``ref``, ``alt``, and
    ``methylation_level`` and the VEP ``most_severe_consequence`` from the context
    Table.

    :param cohort_hts: Dictionary of cohort name to Table returned by
        :func:`format_dnm_cohort`.
    :param context_ht: Annotated context Table from the constraint pipeline
        (``--prepare-context-ht`` output).
    :return: DNM Table keyed by ``locus`` and ``alleles`` with ``cohort_counts``,
        ``context``, ``ref``, ``alt``, ``methylation_level``, and
        ``most_severe_consequence``.
    """
    hts = [ht.select("cohort") for ht in cohort_hts.values()]
    ht = hts[0].union(*hts[1:])
    ht = ht.checkpoint(new_temp_file("dnm.all_cohorts", "ht"))

    # The cohort Tables have few partitions (text files are imported as one), so
    # joining the context Table into them would scan most of the context Table in
    # each DNM partition. Instead, filter the context Table to the DNM sites in a
    # scan over its own partitions, checkpoint the small result, then coalesce it.
    context_ht = context_ht.semi_join(ht.select().distinct())
    context_ht = context_ht.select(
        "context",
        "ref",
        "alt",
        "methylation_level",
        most_severe_consequence=context_ht.vep.most_severe_consequence,
    ).select_globals()
    context_ht = context_ht.checkpoint(new_temp_file("dnm.context", "ht"))
    context_ht = context_ht.naive_coalesce(50)

    ht = ht.semi_join(context_ht)
    ht = ht.group_by("locus", "alleles").aggregate(
        cohort_counts=hl.agg.counter(ht.cohort)
    )

    return ht.annotate(**context_ht[ht.key])


# ---------------------------------------------------------------------------
# Regional genomic features
# ---------------------------------------------------------------------------


def format_genomic_features(ht: hl.Table) -> hl.Table:
    """
    Key an imported 1kb regional genomic feature matrix by window interval.

    The matrix has an ``element_id`` column ("chr-start-end", BED coordinates)
    followed by ``GENOMIC_FEATURE_COLUMNS``. NaN feature values become missing. Each
    window's interval is the BED window in 1-based positions, ``[start + 1, end]``.
    (The Gnocchi rr Tables used ``[start, end)``, 1bp left of the BED window.)

    :param ht: Imported feature matrix with ``element_id`` and the feature columns
        as floats.
    :return: Feature Table keyed by ``interval`` with ``element_id`` and the feature
        columns.
    """
    ht = ht.select("element_id", *GENOMIC_FEATURE_COLUMNS)
    ht = ht.annotate(
        **{
            c: hl.if_else(hl.is_nan(ht[c]), hl.missing(hl.tfloat64), ht[c])
            for c in GENOMIC_FEATURE_COLUMNS
        }
    )
    coords = ht.element_id.split("-")
    ht = ht.annotate(
        _contig=coords[0], _start=hl.int(coords[1]), _end=hl.int(coords[2])
    )
    ht = ht.key_by(
        interval=hl.locus_interval(
            ht._contig,
            ht._start + 1,
            ht._end,
            includes_end=True,
            reference_genome="GRCh38",
        )
    )

    return ht.drop("_contig", "_start", "_end")


def summarize_genomic_features(ht: hl.Table) -> Dict[str, Any]:
    """
    Summarize a genomic feature Table per contig.

    Returns the number of windows and the features that are entirely missing per
    contig, and the number of duplicated window intervals.

    :param ht: Feature Table returned by :func:`format_genomic_features`.
    :return: Dictionary with ``n_duplicated_windows`` and ``contigs`` (contig ->
        dict with ``n_windows`` and ``missing_features``).
    """
    n_dup = ht.count() - ht.select().distinct().count()

    contig_summary = ht.aggregate(
        hl.agg.group_by(
            ht.interval.start.contig,
            hl.struct(
                n_windows=hl.agg.count(),
                n_missing=hl.struct(
                    **{
                        c: hl.agg.count_where(hl.is_missing(ht[c]))
                        for c in GENOMIC_FEATURE_COLUMNS
                    }
                ),
            ),
        )
    )
    contigs = {}
    for contig, s in contig_summary.items():
        missing = sorted(
            {
                c.rsplit("_", 1)[0]
                for c in GENOMIC_FEATURE_COLUMNS
                if s.n_missing[c] == s.n_windows
            }
        )
        contigs[contig] = {"n_windows": s.n_windows, "missing_features": missing}

    return {"n_duplicated_windows": n_dup, "contigs": contigs}
