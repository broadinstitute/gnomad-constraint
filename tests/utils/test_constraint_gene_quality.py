"""Tests for the gene quality metrics functions in the constraint utility module."""

import math
from typing import Dict, List, Optional, Set, Tuple

import hail as hl
import pytest

from gnomad_constraint.utils.constraint import (
    _compute_coverage_metrics,
    _compute_region_flag_metrics,
    _compute_site_quality_metrics,
    compute_gene_quality_metrics,
)

# CDS positions per transcript. T1 and T2 share chr1:109 and chr1:110. T3 has no
# context coverage rows, and T4 has coverage but no exome SNVs.
CDS = {
    "T1": [f"chr1:{p}" for p in range(101, 111)],
    "T2": ["chr1:109", "chr1:110", "chr1:111"],
    "T3": ["chr2:101"],
    "T4": ["chr3:101", "chr3:102"],
}

# Exome coverage (AN percent) per locus. chr1:103 and chr1:107 sit exactly on the
# default threshold of 90, and chr1:109 has missing coverage.
COVERAGE = {
    "chr1:101": 100,
    "chr1:102": 95,
    "chr1:103": 90,
    "chr1:104": 89,
    "chr1:105": 0,
    "chr1:106": 50,
    "chr1:107": 90,
    "chr1:108": 91,
    "chr1:109": None,
    "chr1:110": 100,
    "chr1:111": 85,
    "chr3:101": 0,
    "chr3:102": 0,
}

# Exome sites: (locus, alleles, AS_MQ). The indel at chr1:102 and the SNV outside
# any CDS at chr1:150 must be ignored, and the missing AS_MQ at chr1:110 is skipped.
SITES = [
    ("chr1:101", ["C", "T"], 60.0),
    ("chr1:102", ["A", "G"], 40.0),
    ("chr1:102", ["A", "AT"], 0.0),
    ("chr1:110", ["G", "A"], None),
    ("chr1:111", ["T", "C"], 30.0),
    ("chr1:150", ["C", "G"], 10.0),
]

# Inclusive segdup and LCR intervals.
SEGDUP = [("chr1", 101, 103)]
LCR = [("chr1", 110, 111)]


def _locus(locus: str) -> hl.expr.LocusExpression:
    """Return a GRCh38 locus expression from a "contig:pos" string."""
    return hl.parse_locus(locus, reference_genome="GRCh38")


def _cds_ht(cds: Dict[str, List[str]] = CDS) -> hl.Table:
    """Return a GENCODE CDS positions Table keyed by locus with transcript IDs."""
    by_locus: Dict[str, List[str]] = {}
    for transcript, loci in cds.items():
        for locus in loci:
            by_locus.setdefault(locus, []).append(transcript)
    return hl.Table.parallelize(
        [{"locus": _locus(l), "transcript_id": t} for l, t in by_locus.items()],
        hl.tstruct(locus=hl.tlocus("GRCh38"), transcript_id=hl.tarray(hl.tstr)),
        key="locus",
    )


def _context_ht(coverage: Dict[str, Optional[int]] = COVERAGE) -> hl.Table:
    """Return a preprocessed-like context Table with three SNVs per locus."""
    rows = [
        {"locus": _locus(l), "alleles": ["A", alt], "exomes_coverage": c}
        for l, c in coverage.items()
        for alt in ["C", "G", "T"]
    ]
    return hl.Table.parallelize(
        rows,
        hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            alleles=hl.tarray(hl.tstr),
            exomes_coverage=hl.tint32,
        ),
        key=["locus", "alleles"],
    )


def _sites_ht(
    sites: List[Tuple[str, List[str], Optional[float]]] = SITES,
) -> hl.Table:
    """Return an exomes sites Table with ``info.AS_MQ``."""
    return hl.Table.parallelize(
        [
            {"locus": _locus(l), "alleles": a, "info": hl.Struct(AS_MQ=mq)}
            for l, a, mq in sites
        ],
        hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            alleles=hl.tarray(hl.tstr),
            info=hl.tstruct(AS_MQ=hl.tfloat64),
        ),
        key=["locus", "alleles"],
    )


def _intervals_ht(intervals: List[Tuple[str, int, int]]) -> hl.Table:
    """Return a Table keyed by inclusive GRCh38 locus intervals."""
    return hl.Table.parallelize(
        [
            {
                "interval": hl.locus_interval(
                    c, s, e, includes_end=True, reference_genome="GRCh38"
                )
            }
            for c, s, e in intervals
        ],
        hl.tstruct(interval=hl.tinterval(hl.tlocus("GRCh38"))),
        key="interval",
    )


def _by_transcript(ht: hl.Table) -> Dict[str, hl.Struct]:
    """Collect a transcript-keyed Table into a dict."""
    return {r.transcript: r for r in ht.collect()}


class TestComputeCoverageMetrics:
    """Test the _compute_coverage_metrics function."""

    @pytest.mark.parametrize(
        "threshold, expected",
        [
            # T1 has 9 covered positions (chr1:109 is missing); 6 are >= 90,
            # including the two at exactly 90. T2 has chr1:110 (100) and chr1:111
            # (85).
            (90, {"T1": 6 / 9, "T2": 1 / 2, "T4": 0.0}),
            (91, {"T1": 4 / 9, "T2": 1 / 2, "T4": 0.0}),
            (0, {"T1": 1.0, "T2": 1.0, "T4": 1.0}),
        ],
    )
    def test_prop_bp_an90(self, threshold: int, expected: Dict[str, float]) -> None:
        """Test the fraction of CDS positions at or above the coverage threshold."""
        out = _by_transcript(
            _compute_coverage_metrics(_context_ht(), _cds_ht(), threshold)
        )
        assert sorted(out) == sorted(expected)
        for transcript, prop in expected.items():
            assert out[transcript].prop_bp_AN90 == pytest.approx(prop)

    def test_one_row_per_locus(self) -> None:
        """Test that the three SNVs at a locus count as one CDS position."""
        ht = _compute_coverage_metrics(
            _context_ht({"chr1:101": 100, "chr1:102": 0}),
            _cds_ht({"T1": ["chr1:101", "chr1:102"]}),
        )
        assert _by_transcript(ht)["T1"].prop_bp_AN90 == pytest.approx(0.5)

    def test_positions_without_coverage_are_excluded(self) -> None:
        """Test that CDS positions with no coverage value leave the denominator.

        This pins current behavior: a position with missing ``exomes_coverage``, or
        no row in the context Table, is dropped rather than counted as not
        covered, and a transcript with no covered positions (T3) gets no row.
        """
        out = _by_transcript(_compute_coverage_metrics(_context_ht(), _cds_ht()))
        assert "T3" not in out
        assert out["T1"].prop_bp_AN90 == pytest.approx(6 / 9)


class TestComputeSiteQualityMetrics:
    """Test the _compute_site_quality_metrics function."""

    def test_mean_as_mq(self) -> None:
        """Test the mean AS_MQ over SNVs in each transcript's CDS.

        The indel at chr1:102 and the SNV outside the CDS are ignored, and the
        missing AS_MQ at chr1:110 is skipped by the mean.
        """
        out = _by_transcript(_compute_site_quality_metrics(_sites_ht(), _cds_ht()))
        assert sorted(out) == ["T1", "T2"]
        assert out["T1"].mean_AS_MQ == pytest.approx(50.0)
        assert out["T2"].mean_AS_MQ == pytest.approx(30.0)

    def test_all_missing_as_mq(self) -> None:
        """Test that a transcript whose only SNVs lack AS_MQ gets a NaN mean.

        This pins current behavior: ``hl.agg.mean`` over no defined values is NaN,
        not missing, and NaN is never flagged as low mapping quality.
        """
        ht = _compute_site_quality_metrics(
            _sites_ht([("chr1:101", ["C", "T"], None)]),
            _cds_ht({"T1": ["chr1:101"]}),
        )
        assert math.isnan(_by_transcript(ht)["T1"].mean_AS_MQ)


class TestComputeRegionFlagMetrics:
    """Test the _compute_region_flag_metrics function."""

    def test_proportions(self) -> None:
        """Test the fraction of CDS positions in segdups and LCRs per transcript.

        Region overlap uses every CDS position, including those with no coverage.
        """
        out = _by_transcript(
            _compute_region_flag_metrics(
                _cds_ht(), _intervals_ht(SEGDUP), _intervals_ht(LCR)
            )
        )
        expected = {
            "T1": (0.3, 0.1),
            "T2": (0.0, 2 / 3),
            "T3": (0.0, 0.0),
            "T4": (0.0, 0.0),
        }
        assert sorted(out) == sorted(expected)
        for transcript, (segdup, lcr) in expected.items():
            assert out[transcript].prop_segdup == pytest.approx(segdup)
            assert out[transcript].prop_LCR == pytest.approx(lcr)

    def test_no_intervals(self) -> None:
        """Test that empty interval Tables give proportions of 0."""
        out = _by_transcript(
            _compute_region_flag_metrics(
                _cds_ht(), _intervals_ht([]), _intervals_ht([])
            )
        )
        assert all(r.prop_segdup == 0 and r.prop_LCR == 0 for r in out.values())


# Flag boundary cases: (transcript, covered CDS positions out of 10, AS_MQ, flags).
# One covered position of 10 is a prop_bp_AN90 of exactly 0.1.
FLAG_CASES = [
    ("F0", 1, 50.0, set()),
    ("F1", 0, 50.0, {"low_exome_coverage"}),
    ("F2", 1, 49.9, {"low_exome_mapping_quality"}),
    ("F3", 0, 10.0, {"low_exome_coverage", "low_exome_mapping_quality"}),
]


def _run_gene_quality_metrics(
    context_ht: Optional[hl.Table] = None,
    sites_ht: Optional[hl.Table] = None,
    cds_ht: Optional[hl.Table] = None,
) -> hl.Table:
    """Run compute_gene_quality_metrics, on the default fixtures unless given."""
    return compute_gene_quality_metrics(
        _context_ht() if context_ht is None else context_ht,
        _sites_ht() if sites_ht is None else sites_ht,
        _cds_ht() if cds_ht is None else cds_ht,
        _intervals_ht(SEGDUP),
        _intervals_ht(LCR),
    )


@pytest.fixture(scope="module")
def gene_quality_ht() -> hl.Table:
    """Return compute_gene_quality_metrics on the default fixtures, checkpointed."""
    return _run_gene_quality_metrics().cache()


@pytest.fixture(scope="module")
def flag_cases_rows() -> Dict[str, hl.Struct]:
    """Return compute_gene_quality_metrics rows for the FLAG_CASES transcripts."""
    cds, coverage, sites = {}, {}, []
    for i, (transcript, n_covered, as_mq, _) in enumerate(FLAG_CASES):
        loci = [f"chr{4 + i}:{p}" for p in range(101, 111)]
        cds[transcript] = loci
        coverage.update({l: (100 if j < n_covered else 0) for j, l in enumerate(loci)})
        sites.append((loci[0], ["C", "T"], as_mq))
    return _by_transcript(
        _run_gene_quality_metrics(_context_ht(coverage), _sites_ht(sites), _cds_ht(cds))
    )


class TestComputeGeneQualityMetrics:
    """Test the compute_gene_quality_metrics function."""

    def test_schema(self, gene_quality_ht: hl.Table) -> None:
        """Test the key, row fields, and release field names."""
        ht = gene_quality_ht
        assert list(ht.key) == ["transcript"]
        assert list(ht.row) == ["transcript", "gene_quality_metrics", "gene_flags"]
        assert list(ht.gene_quality_metrics.dtype) == [
            "exome_prop_bp_AN90",
            "exome_mean_AS_MQ",
            "exome_prop_segdup",
            "exome_prop_LCR",
        ]

    def test_metrics(self, gene_quality_ht: hl.Table) -> None:
        """Test the combined metrics and flags.

        T1's mean AS_MQ is exactly 50, which is not low. T3 has no covered CDS
        positions, so it gets no row.
        """
        out = _by_transcript(gene_quality_ht)
        assert "T3" not in out

        t1 = out["T1"].gene_quality_metrics
        assert t1.exome_prop_bp_AN90 == pytest.approx(6 / 9)
        assert t1.exome_mean_AS_MQ == pytest.approx(50.0)
        assert t1.exome_prop_segdup == pytest.approx(0.3)
        assert t1.exome_prop_LCR == pytest.approx(0.1)
        assert out["T1"].gene_flags == set()

        t2 = out["T2"].gene_quality_metrics
        assert t2.exome_prop_bp_AN90 == pytest.approx(0.5)
        assert t2.exome_mean_AS_MQ == pytest.approx(30.0)
        assert t2.exome_prop_LCR == pytest.approx(2 / 3)
        assert out["T2"].gene_flags == {"low_exome_mapping_quality"}

    @pytest.mark.parametrize("transcript, n_covered, as_mq, expected", FLAG_CASES)
    def test_flag_thresholds(
        self,
        flag_cases_rows: Dict[str, hl.Struct],
        transcript: str,
        n_covered: int,
        as_mq: float,
        expected: Set[str],
    ) -> None:
        """Test the flag boundaries: AS_MQ < 50 and prop_bp_AN90 < 0.1."""
        row = flag_cases_rows[transcript]
        assert row.gene_quality_metrics.exome_prop_bp_AN90 == pytest.approx(
            n_covered / 10
        )
        assert row.gene_quality_metrics.exome_mean_AS_MQ == pytest.approx(as_mq)
        assert row.gene_flags == expected

    def test_missing_as_mq_keeps_coverage_flag(self, gene_quality_ht: hl.Table) -> None:
        """Test that a transcript with no exome SNVs still gets low_exome_coverage.

        T4 has CDS coverage of 0 and no exome SNVs, so its mean AS_MQ is missing.
        The missing mapping quality condition is skipped rather than making
        ``gene_flags`` missing.
        """
        t4 = _by_transcript(gene_quality_ht)["T4"]
        assert t4.gene_quality_metrics.exome_mean_AS_MQ is None
        assert t4.gene_flags == {"low_exome_coverage"}
