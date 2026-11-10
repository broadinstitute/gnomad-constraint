"""Tests for the regional mutation rate (rr) functions in ``gnomad_constraint.utils.regional_rate``."""

import hail as hl
import numpy as np
import pytest

from gnomad_constraint.resources.constants import (
    GENOMIC_FEATURE_COLUMNS,
    METHYLATION_CORRELATED_FEATURES,
    RR_CONTEXT_FEATURES,
)
from gnomad_constraint.utils.regional_rate import (
    create_rr_training_ht,
    remove_qc_flagged_genes,
    rr_training_ht_to_pandas,
)

INTERVAL_TYPE = hl.tinterval(hl.tlocus("GRCh38"))


def _locus(pos: int, contig: str = "chr1") -> hl.Locus:
    """Return a GRCh38 locus."""
    return hl.Locus(contig, pos, "GRCh38")


class TestRrContextFeatures:
    """Test the RR_CONTEXT_FEATURES constant."""

    def test_cpg_contexts_leave_out_methylation_features(self) -> None:
        """Test that only CpG contexts leave out the methylation-correlated features."""
        assert RR_CONTEXT_FEATURES["AAA"] == GENOMIC_FEATURE_COLUMNS
        assert RR_CONTEXT_FEATURES["ACG"] == [
            c
            for c in GENOMIC_FEATURE_COLUMNS
            if c.rsplit("_", 1)[0] not in METHYLATION_CORRELATED_FEATURES
        ]


class TestCreateRrTrainingHt:
    """Test the create_rr_training_ht function."""

    @staticmethod
    def _context_ht() -> hl.Table:
        """
        Return context Table rows at chr1 positions 10-70.

        10 and 20 are synonymous DNMs, 20 also in exomes. 30 has two synonymous
        alleles absent from exomes, 40 is in exomes and passes, 50 is in exomes as
        AC0, 60 is missense, and 70 is absent from exomes.
        """
        rows = [
            (10, ["A", "G"], "synonymous_variant", None),
            (20, ["C", "T"], "synonymous_variant", set()),
            (30, ["A", "G"], "synonymous_variant", None),
            (30, ["A", "C"], "synonymous_variant", None),
            (40, ["A", "G"], "synonymous_variant", set()),
            (50, ["A", "G"], "synonymous_variant", {"AC0"}),
            (60, ["C", "T"], "missense_variant", None),
            (70, ["A", "G"], "synonymous_variant", None),
        ]
        return hl.Table.parallelize(
            [
                {
                    "locus": _locus(pos),
                    "alleles": alleles,
                    "context": "AAA",
                    "vep": {"most_severe_consequence": csq},
                    "filters": {"exomes": filters, "genomes": None},
                }
                for pos, alleles, csq, filters in rows
            ],
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                context=hl.tstr,
                vep=hl.tstruct(most_severe_consequence=hl.tstr),
                filters=hl.tstruct(exomes=hl.tset(hl.tstr), genomes=hl.tset(hl.tstr)),
            ),
            key=["locus", "alleles"],
        )

    @staticmethod
    def _dnm_ht() -> hl.Table:
        """Return DNMs at positions 10, 20, and 60, and an indel at 70."""
        return hl.Table.parallelize(
            [
                {"locus": _locus(pos), "alleles": alleles}
                for pos, alleles in [
                    (10, ["A", "G"]),
                    (20, ["C", "T"]),
                    (60, ["C", "T"]),
                    (70, ["A", "AT"]),
                ]
            ],
            hl.tstruct(locus=hl.tlocus("GRCh38"), alleles=hl.tarray(hl.tstr)),
            key=["locus", "alleles"],
        )

    @staticmethod
    def _features_ht() -> hl.Table:
        """Return two windows, positions 1-25 (features 1.0) and 26-1000 (2.0)."""
        return hl.Table.parallelize(
            [
                {
                    "interval": hl.Interval(
                        _locus(start), _locus(end), includes_end=True
                    ),
                    **dict.fromkeys(GENOMIC_FEATURE_COLUMNS, value),
                }
                for start, end, value in [(1, 25, 1.0), (26, 1000, 2.0)]
            ],
            hl.tstruct(
                interval=INTERVAL_TYPE,
                **dict.fromkeys(GENOMIC_FEATURE_COLUMNS, hl.tfloat64),
            ),
            key="interval",
        )

    def _training_ht(self, background_fraction: float, seed: int = 0) -> hl.Table:
        """Return the training Table of the test rows."""
        return create_rr_training_ht(
            self._context_ht(),
            self._dnm_ht(),
            self._features_ht(),
            background_fraction,
            seed,
        )

    def test_training_sites(self) -> None:
        """
        Test the DNM and background sites of the training set.

        DNMs are the synonymous context sites in the DNM Table, even if in exomes
        (20); the missense DNM (60) and the indel (70) aren't. The background (all of
        it, with a fraction of 1) is the other synonymous sites with no exomes
        filters, one per locus (30, 70), so sites in exomes, passing (40) or AC0
        (50), aren't in it. Every site gets its window's features.
        """
        rows = self._training_ht(1.0).collect()
        assert sorted((r.locus.position, r.y, r.dist2telo_1k) for r in rows) == [
            (10, 1, 1.0),
            (20, 1, 1.0),
            (30, 0, 2.0),
            (70, 0, 2.0),
        ]

    def test_background_sample(self) -> None:
        """Test that a sampled background is a subset of it and repeats with the seed."""

        def _background(seed: int) -> list:
            ht = self._training_ht(0.5, seed)
            return sorted(r.locus.position for r in ht.collect() if r.y == 0)

        assert set(_background(0)) <= {30, 70}
        assert _background(0) == _background(0)


class TestRemoveQcFlaggedGenes:
    """Test the remove_qc_flagged_genes function."""

    @pytest.mark.parametrize("flagged", [True, False])
    def test_flagged_ensembl_mane_select_spans(self, flagged: bool) -> None:
        """
        Test that only the spans of flagged Ensembl MANE Select transcripts are removed.

        A span includes its start and end. A transcript with no span is skipped. With
        no flagged gene, every site is kept.
        """
        positions = [100, 149, 150, 200, 250, 251, 300]
        ht = hl.Table.parallelize(
            [{"locus": _locus(p)} for p in positions],
            hl.tstruct(locus=hl.tlocus("GRCh38")),
            key="locus",
        )
        flags = {"low_exome_coverage"} if flagged else set()
        metrics = [
            ("ENST1", True, flags, 150, 250),
            ("ENST2", True, set(), 50, 120),
            ("NM_1", True, {"low_exome_coverage"}, 280, 320),
            ("ENST3", False, {"low_exome_coverage"}, 90, 110),
            ("ENST4", True, {"low_exome_coverage"}, None, None),
        ]
        metrics_ht = hl.Table.parallelize(
            [
                {
                    "transcript": transcript,
                    "mane_select": mane_select,
                    "gene_flags": flags,
                    "chromosome": "chr1",
                    "start_position": start,
                    "end_position": end,
                }
                for transcript, mane_select, flags, start, end in metrics
            ],
            hl.tstruct(
                transcript=hl.tstr,
                mane_select=hl.tbool,
                gene_flags=hl.tset(hl.tstr),
                chromosome=hl.tstr,
                start_position=hl.tint32,
                end_position=hl.tint32,
            ),
        )
        out = remove_qc_flagged_genes(ht, metrics_ht).collect()
        expected = [100, 149, 251, 300] if flagged else positions
        assert [r.locus.position for r in out] == expected


class TestRrTrainingHtToPandas:
    """Test the rr_training_ht_to_pandas function."""

    def test_columns_types_and_folds(self) -> None:
        """
        Test the columns, the NumPy dtypes the fit relies on, and the folds.

        chr1 is in the odd fold and chr2 in the even one.
        """
        features = dict.fromkeys(GENOMIC_FEATURE_COLUMNS, 1.0)
        ht = hl.Table.parallelize(
            [
                {"locus": _locus(10), "context": "AAA", "y": 1, **features},
                {
                    "locus": _locus(10, "chr2"),
                    "context": "AAA",
                    "y": 0,
                    **features,
                    "dist2telo_1k": None,
                },
            ],
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                context=hl.tstr,
                y=hl.tint32,
                **dict.fromkeys(GENOMIC_FEATURE_COLUMNS, hl.tfloat64),
            ),
            key="locus",
        )
        df = rr_training_ht_to_pandas(ht)
        assert list(df.columns) == ["context", "y", *GENOMIC_FEATURE_COLUMNS, "fold"]
        assert df["y"].dtype == np.int64
        assert all(df[c].dtype == np.float64 for c in GENOMIC_FEATURE_COLUMNS)
        assert df["dist2telo_1k"].isna().sum() == 1
        assert list(df["fold"]) == ["odd", "even"]
