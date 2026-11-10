"""Tests for the import functions in ``gnomad_constraint.utils.import_data``."""

import math

import hail as hl
import pytest

from gnomad_constraint.resources.constants import GENOMIC_FEATURE_COLUMNS
from gnomad_constraint.utils.import_data import (
    combine_dnms,
    format_dnm_cohort,
    format_genomic_features,
)

VARIANT_KEY_TYPE = hl.tstruct(locus=hl.tlocus("GRCh38"), alleles=hl.tarray(hl.tstr))


def _cohort_ht(cohort: str, variants: list) -> hl.Table:
    """Return a formatted cohort Table of (position, alleles) on chr1."""
    ht = hl.Table.parallelize(
        [{"locus": hl.Locus("chr1", p, "GRCh38"), "alleles": a} for p, a in variants],
        VARIANT_KEY_TYPE,
        key=["locus", "alleles"],
    )
    return ht.annotate(cohort=cohort)


class TestFormatDnmCohort:
    """Test the format_dnm_cohort function on the GRCh38 cohorts."""

    @pytest.mark.parametrize(
        "cohort, row, schema",
        [
            (
                "satterstrom_2026_asd",
                {"Variant": "chr1:100:A:G", "Confidence": "high"},
                hl.tstruct(Variant=hl.tstr, Confidence=hl.tstr),
            ),
            (
                "satterstrom_2026_ddd",
                {"Variant": "chr1:100:A:G", "LCR": "False"},
                hl.tstruct(Variant=hl.tstr, LCR=hl.tstr),
            ),
            (
                "decode_psychencode",
                {"locus": "chr1:100", "ref": "A", "alt": "G"},
                hl.tstruct(locus=hl.tstr, ref=hl.tstr, alt=hl.tstr),
            ),
            (
                "ofc_2020",
                {"Chromosome": 1, "Position (build 38)": 100, "Ref": "A", "Alt": "G"},
                hl.tstruct(
                    **{
                        "Chromosome": hl.tint32,
                        "Position (build 38)": hl.tint32,
                        "Ref": hl.tstr,
                        "Alt": hl.tstr,
                    }
                ),
            ),
        ],
    )
    def test_grch38_cohorts(self, cohort: str, row: dict, schema: hl.tstruct) -> None:
        """Test that each cohort becomes a locus/alleles-keyed Table with its name."""
        ht = format_dnm_cohort(cohort, hl.Table.parallelize([row], schema))
        assert list(ht.key) == ["locus", "alleles"]
        assert ht.collect() == [
            hl.Struct(
                locus=hl.Locus("chr1", 100, "GRCh38"), alleles=["A", "G"], cohort=cohort
            )
        ]

    def test_gnomad_cohort(self) -> None:
        """Test that the gnomAD de novo Table keeps only its key and gets a cohort."""
        ht = hl.Table.parallelize(
            [{"locus": hl.Locus("chr1", 100, "GRCh38"), "alleles": ["A", "G"], "x": 1}],
            hl.tstruct(
                locus=hl.tlocus("GRCh38"), alleles=hl.tarray(hl.tstr), x=hl.tint32
            ),
        )
        ht = format_dnm_cohort("gnomad_v4.1", ht)
        assert ht.collect() == [
            hl.Struct(
                locus=hl.Locus("chr1", 100, "GRCh38"),
                alleles=["A", "G"],
                cohort="gnomad_v4.1",
            )
        ]

    @pytest.mark.parametrize(
        "cohort, row, schema",
        [
            (
                "chd_2017",
                {"CHROM": "1", "POS": 100, "REF": "A", "ALT": "G"},
                hl.tstruct(CHROM=hl.tstr, POS=hl.tint32, REF=hl.tstr, ALT=hl.tstr),
            ),
            (
                "ee_2014",
                {
                    "hg19 coordinates (chr:position)": "chr1:100",
                    "Ref/Alt alleles": "A/G",
                },
                hl.tstruct(
                    **{
                        "hg19 coordinates (chr:position)": hl.tstr,
                        "Ref/Alt alleles": hl.tstr,
                    }
                ),
            ),
            (
                "cdh_2021",
                {"chr": 1, "position": 100, "ref": "A", "alt": "G"},
                hl.tstruct(chr=hl.tint32, position=hl.tint32, ref=hl.tstr, alt=hl.tstr),
            ),
        ],
    )
    def test_grch37_cohorts(
        self,
        monkeypatch: pytest.MonkeyPatch,
        cohort: str,
        row: dict,
        schema: hl.tstruct,
    ) -> None:
        """
        Test the GRCh37 cohorts' locus and allele parsing before liftover.

        Liftover needs the chain file on GCS, so it is replaced by the identity: the
        Table it gets must be keyed by the GRCh37 locus and alleles.
        """
        monkeypatch.setattr(
            "gnomad_constraint.utils.import_data.default_lift_data", lambda ht: ht
        )
        ht = format_dnm_cohort(cohort, hl.Table.parallelize([row], schema))
        assert list(ht.key) == ["locus", "alleles"]
        assert ht.collect() == [
            hl.Struct(
                locus=hl.Locus("1", 100, "GRCh37"), alleles=["A", "G"], cohort=cohort
            )
        ]

    def test_unknown_cohort(self) -> None:
        """Test that an unknown cohort name raises an error."""
        ht = hl.Table.parallelize([], hl.tstruct(x=hl.tint32))
        with pytest.raises(ValueError, match="Unknown DNM cohort"):
            format_dnm_cohort("other", ht)


class TestCombineDnms:
    """Test the combine_dnms function."""

    def test_combine(self) -> None:
        """
        Test deduplication, cohort counts, the context join, and the annotations.

        DNMs not in the context Table (an indel, the wrong allele orientation, a site
        absent from it) are dropped. A site reported twice by one cohort counts 2.
        """
        cohort_hts = {
            "a": _cohort_ht(
                "a",
                [
                    (10, ["A", "G"]),
                    (10, ["A", "G"]),
                    (10, ["G", "A"]),
                    (20, ["C", "CT"]),
                    (99, ["A", "G"]),
                ],
            ),
            "b": _cohort_ht("b", [(10, ["A", "G"]), (30, ["C", "T"])]),
        }
        context_ht = hl.Table.parallelize(
            [
                {
                    "locus": hl.Locus("chr1", p, "GRCh38"),
                    "alleles": a,
                    "context": context,
                    "ref": a[0],
                    "alt": a[1],
                    "methylation_level": methylation,
                    "vep": {"most_severe_consequence": csq},
                }
                for p, a, context, methylation, csq in [
                    (10, ["A", "G"], "AAA", 0, "synonymous_variant"),
                    (30, ["C", "T"], "ACG", 2, "missense_variant"),
                    (50, ["A", "G"], "AAA", 0, "intron_variant"),
                ]
            ],
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                context=hl.tstr,
                ref=hl.tstr,
                alt=hl.tstr,
                methylation_level=hl.tint32,
                vep=hl.tstruct(most_severe_consequence=hl.tstr),
            ),
            key=["locus", "alleles"],
        )
        rows = combine_dnms(cohort_hts, context_ht).collect()
        assert [(r.locus.position, r.alleles) for r in rows] == [
            (10, ["A", "G"]),
            (30, ["C", "T"]),
        ]
        assert rows[0].cohort_counts == {"a": 2, "b": 1}
        assert rows[1].cohort_counts == {"b": 1}
        assert (rows[0].context, rows[0].most_severe_consequence) == (
            "AAA",
            "synonymous_variant",
        )
        assert (rows[1].context, rows[1].methylation_level) == ("ACG", 2)


class TestFormatGenomicFeatures:
    """Test the format_genomic_features function."""

    def test_windows_and_missing_values(self) -> None:
        """
        Test that each window is its BED window in 1-based positions, and NaN is missing.

        The BED window ``chr1-0-1000`` covers positions 1-1000 and ``chr1-1000-2000``
        positions 1001-2000.
        """
        values = dict.fromkeys(GENOMIC_FEATURE_COLUMNS, 1.0)
        ht = hl.Table.parallelize(
            [
                {"element_id": "chr1-0-1000", **values, "dist2telo_1k": math.nan},
                {"element_id": "chr1-1000-2000", **values},
            ],
            hl.tstruct(
                element_id=hl.tstr,
                **dict.fromkeys(GENOMIC_FEATURE_COLUMNS, hl.tfloat64),
            ),
        )
        ht = format_genomic_features(ht)
        loci_ht = hl.Table.parallelize(
            [{"pos": p} for p in [1, 1000, 1001, 2000, 2001]], hl.tstruct(pos=hl.tint32)
        )
        window = ht[hl.locus("chr1", loci_ht.pos, "GRCh38")]
        loci_ht = loci_ht.annotate(
            element_id=window.element_id, dist2telo_1k=window.dist2telo_1k
        )
        rows = {r.pos: (r.element_id, r.dist2telo_1k) for r in loci_ht.collect()}
        assert rows == {
            1: ("chr1-0-1000", None),
            1000: ("chr1-0-1000", None),
            1001: ("chr1-1000-2000", 1.0),
            2000: ("chr1-1000-2000", 1.0),
            2001: (None, None),
        }
