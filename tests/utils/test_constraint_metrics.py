"""Tests for the constraint metrics and release functions in ``gnomad_constraint.utils.constraint``."""

import os
from typing import Dict, List, Optional, Set, Tuple

import hail as hl
import numpy as np
import pytest
from scipy.stats import gamma

from gnomad_constraint.resources.constants import (
    CONSTRAINT_GRANULARITIES,
    RELEASE_GROUP_NAMES,
    RELEASE_KEY_ORDER,
    RELEASE_TOP_LEVEL_ANNOTATIONS,
)
from gnomad_constraint.utils.constraint import (
    _annotate_oe_ci_z,
    _compute_pli_scores,
    _compute_z_scores,
    _get_constraint_group_name,
    add_oe_upper_rank_and_bins,
    annotate_constraint_percentile_bins,
    compute_constraint_metrics,
    compute_constraint_percentile_bins,
    flatten_release_ht,
    get_transcript_filter_expr,
    lof_bin_thresholds_to_ht,
    prepare_release_ht,
    prepare_release_mutation_ht,
)

CONSTRAINT_GROUP_META = [
    {"csq_set": "syn"},
    {"csq_set": "mis"},
    {"lof": "classic"},
    {"lof": "hc_lc"},
    {"lof": "hc"},
]
"""Constraint group metadata in the order built by ``build_constraint_consequence_groups``."""

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

GEN_ANCS = ["global", "afr", "amr", "eas", "nfe", "sas"]

GROUP_NAMES = ["syn", "mis", "lof_classic", "lof_hc_lc", "lof_hc"]
"""Internal names of ``CONSTRAINT_GROUP_META``, in the same order."""

CI_TYPE = hl.tstruct(lower=hl.tfloat64, upper=hl.tfloat64)

RANK_TYPE = hl.tstruct(
    upper_rank=hl.tint64,
    upper_bin_percentile=hl.tint32,
    upper_bin_decile=hl.tint32,
    upper_bin_sextile=hl.tint32,
)

requires_qgamma = pytest.mark.skipif(
    not hasattr(hl, "qgamma"), reason="The gamma CI needs hl.qgamma (Hail >= 0.2.137)."
)


def _bin(value: Optional[float], thresholds: List[float]) -> Optional[int]:
    """Return the number of thresholds a value meets or exceeds."""
    return None if value is None else sum(value >= t for t in thresholds)


def _constraint_groups_ht(
    rows: List[Dict],
    oe_info_type: hl.tstruct,
    extra_cg_fields: Dict = None,
    meta: List[Dict[str, str]] = CONSTRAINT_GROUP_META,
) -> hl.Table:
    """Return a Table keyed by ``transcript`` with a ``constraint_groups`` array.

    :param rows: Dicts with ``transcript``, ``constraint_groups`` (one struct per
        group, each with an ``oe_info`` list) and any fields in the row type.
    :param oe_info_type: Type of each ``oe_info`` element.
    :param extra_cg_fields: Extra constraint group fields and their types.
    :param meta: Constraint group metadata. Default is ``CONSTRAINT_GROUP_META``.
    :return: Table with the constraint group globals.
    """
    cg_type = hl.tstruct(oe_info=hl.tarray(oe_info_type), **(extra_cg_fields or {}))
    row_fields = {
        k: v
        for k, v in {
            "gene_id": hl.tstr,
            "mane_select": hl.tbool,
            "canonical": hl.tbool,
            "no_variants": hl.tbool,
            "constraint_flags": hl.tset(hl.tstr),
        }.items()
        if k in rows[0]
    }
    ht = hl.Table.parallelize(
        rows,
        hl.tstruct(
            transcript=hl.tstr, **row_fields, constraint_groups=hl.tarray(cg_type)
        ),
        key="transcript",
    )
    return ht.annotate_globals(
        constraint_group_meta=meta, exomes_freq_meta=[{"group": "adj"}]
    )


class TestAnnotateOeCiZ:
    """Test the _annotate_oe_ci_z function."""

    Z_THRESHOLDS = {"syn": (-3.0, 3.0), "mis": (-3.0, None), "lof": (-3.0, None)}

    # (observed, expected) used for every group of the transcript.
    COUNTS = {
        "T_ok": (18, 20.0),
        "T_depleted": (0, 25.0),
        "T_excess": (50, 25.0),
        "T_no_exp": (0, 0.0),
        "T_missing_exp": (3, None),
    }

    @pytest.fixture(scope="class")
    def rows(self) -> Dict[str, hl.Struct]:
        """Return the annotated rows keyed by transcript."""
        oe_info_type = hl.tstruct(
            observed_variants=hl.tint64, expected_variants=hl.tfloat64
        )
        rows = [
            dict(
                transcript=t,
                constraint_groups=[
                    hl.Struct(
                        oe_info=[
                            hl.Struct(observed_variants=obs, expected_variants=exp)
                        ]
                    )
                ]
                * len(CONSTRAINT_GROUP_META),
            )
            for t, (obs, exp) in self.COUNTS.items()
        ]
        ht = _annotate_oe_ci_z(
            _constraint_groups_ht(rows, oe_info_type), self.Z_THRESHOLDS
        )
        return {r.transcript: r for r in ht.collect()}

    @requires_qgamma
    @pytest.mark.parametrize("transcript", ["T_ok", "T_depleted", "T_excess"])
    def test_oe_ci_and_z(self, rows: Dict[str, hl.Struct], transcript: str) -> None:
        """Test oe, both CIs, and z_raw against independent computations."""
        obs, exp = self.COUNTS[transcript]
        for cg in rows[transcript].constraint_groups:
            x = cg.oe_info[0]
            assert x.oe == pytest.approx(obs / exp)
            assert x.z_raw == pytest.approx(
                np.sqrt((obs - exp) ** 2 / exp) * (-1 if obs > exp else 1)
            )
            assert x.oe_ci_gamma.lower == pytest.approx(
                gamma.ppf(0.05, obs + 1, scale=1 / exp)
            )
            assert x.oe_ci_gamma.upper == pytest.approx(
                gamma.ppf(0.95, obs + 1, scale=1 / exp)
            )
            # The discretized Poisson CI only covers oe in [0, 2).
            ci = x.oe_ci_discretized_poisson
            if x.oe < 1.9:
                assert ci.lower <= x.oe <= ci.upper
            if obs == 0:
                assert ci.lower == 0

    @requires_qgamma
    @pytest.mark.parametrize("transcript", ["T_no_exp", "T_missing_exp"])
    def test_no_expected(self, rows: Dict[str, hl.Struct], transcript: str) -> None:
        """Test that a group with zero or missing expected has missing metrics."""
        for cg in rows[transcript].constraint_groups:
            x = cg.oe_info[0]
            assert x.oe is None
            assert x.z_raw is None
            assert x.oe_ci_gamma.lower is None and x.oe_ci_gamma.upper is None

    @requires_qgamma
    @pytest.mark.parametrize(
        "transcript, expected_flags",
        [
            ("T_ok", [set()] * 5),
            # z = 5: only syn has an upper threshold.
            ("T_depleted", [{"outlier_syn"}, set(), set(), set(), set()]),
            # z = -5: every group has a lower threshold, and LoF groups share one.
            (
                "T_excess",
                [{"outlier_syn"}, {"outlier_mis"}] + [{"outlier_lof"}] * 3,
            ),
            ("T_no_exp", [{"no_exp_syn"}, {"no_exp_mis"}] + [{"no_exp_lof"}] * 3),
            (
                "T_missing_exp",
                [{"no_exp_syn"}, {"no_exp_mis"}] + [{"no_exp_lof"}] * 3,
            ),
        ],
    )
    def test_flags(
        self,
        rows: Dict[str, hl.Struct],
        transcript: str,
        expected_flags: List[Set[str]],
    ) -> None:
        """Test the per-group no_exp and outlier flags."""
        assert [cg.flags for cg in rows[transcript].constraint_groups] == expected_flags


class TestComputeZScores:
    """Test the _compute_z_scores function."""

    # Raw z per transcript for (syn, mis, lof_classic, lof_hc_lc, lof_hc).
    Z = {
        "T1": [1.0, -1.0, -1.0, -1.0, -1.0],
        "T2": [-1.0, -3.0, -3.0, -3.0, -3.0],
        "T3": [2.0, 2.0, 2.0, 2.0, 2.0],
        "T4": [-2.0, 5.0, 5.0, 5.0, 5.0],
        "T5": [None, None, None, None, None],
        # Excluded from the SD: no variants, or flagged.
        "T_no_variants": [100.0, -100.0, -100.0, -100.0, -100.0],
        "T_flagged": [50.0, -50.0, -50.0, -50.0, -50.0],
    }

    FLAGS = {
        "T_flagged": [
            {"outlier_syn"},
            {"outlier_mis"},
            {"classic_only"},
            {"hc_lc_only"},
            {"outlier_lof"},
        ]
    }

    @pytest.fixture(scope="class")
    def result(self) -> Tuple[List[float], Dict[str, hl.Struct]]:
        """Return the sd_raw_z global and the annotated rows keyed by transcript."""
        rows = [
            dict(
                transcript=t,
                no_variants=t == "T_no_variants",
                constraint_groups=[
                    hl.Struct(
                        oe_info=[hl.Struct(z_raw=z)],
                        flags=self.FLAGS.get(t, [set()] * 5)[i],
                    )
                    for i, z in enumerate(zs)
                ],
            )
            for t, zs in self.Z.items()
        ]
        ht = _constraint_groups_ht(
            rows, hl.tstruct(z_raw=hl.tfloat64), {"flags": hl.tset(hl.tstr)}
        )
        ht = _compute_z_scores(ht)
        return hl.eval(ht.sd_raw_z), {r.transcript: r for r in ht.collect()}

    def test_sd_raw_z(self, result: Tuple[List[float], Dict[str, hl.Struct]]) -> None:
        """Test the SD of raw z, mirroring negative values for every group but syn."""
        sd, _ = result
        used = ["T1", "T2", "T3", "T4"]
        syn = [self.Z[t][0] for t in used]
        assert sd[0] == pytest.approx(np.std(syn))
        for i in range(1, 5):
            neg = [self.Z[t][i] for t in used if self.Z[t][i] < 0]
            assert sd[i] == pytest.approx(np.std(neg + [-z for z in neg]))

    def test_z_score(self, result: Tuple[List[float], Dict[str, hl.Struct]]) -> None:
        """Test that z_score is z_raw divided by the group's SD."""
        sd, rows = result
        for row in rows.values():
            for i, cg in enumerate(row.constraint_groups):
                z = self.Z[row.transcript][i]
                if z is None:
                    assert cg.z_score is None
                else:
                    assert cg.z_score == pytest.approx(z / sd[i])

    def test_constraint_flags(
        self, result: Tuple[List[float], Dict[str, hl.Struct]]
    ) -> None:
        """Test that constraint_flags is the union of the syn, mis, and HC LoF flags."""
        _, rows = result
        assert rows["T_flagged"].constraint_flags == {
            "outlier_syn",
            "outlier_mis",
            "outlier_lof",
        }
        assert rows["T1"].constraint_flags == set()


class TestGetTranscriptFilterExpr:
    """Test the get_transcript_filter_expr function."""

    # transcript: (gene_id, mane_select, canonical)
    TRANSCRIPTS = {
        "ENSTA1": ("A", True, True),
        "ENSTA2": ("A", False, True),
        "ENSTA3": ("A", False, False),
        "ENSTB1": ("B", False, True),
        "ENSTB2": ("B", False, False),
        "NM_C1": ("C", True, True),
        "ENSTC1": ("C", False, True),
        "ENSTD1": ("D", False, False),
    }

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            # MANE Select per gene, falling back to canonical when a gene has no
            # Ensembl MANE Select transcript (gene C's MANE Select is RefSeq).
            ({}, {"ENSTA1", "ENSTB1", "ENSTC1"}),
            (
                {"use_mane_select_over_canonical": False},
                {"ENSTA1", "ENSTA2", "ENSTB1", "ENSTC1"},
            ),
            ({"mane_select_only": True}, {"ENSTA1"}),
            (
                {"use_mane_select_over_canonical": False, "mane_select_only": True},
                {"ENSTA1"},
            ),
        ],
    )
    def test_selected_transcripts(self, kwargs: Dict, expected: Set[str]) -> None:
        """Test which transcripts each mode selects."""
        ht = hl.Table.parallelize(
            [
                dict(transcript=t, gene_id=g, mane_select=m, canonical=c)
                for t, (g, m, c) in self.TRANSCRIPTS.items()
            ],
            hl.tstruct(
                transcript=hl.tstr,
                gene_id=hl.tstr,
                mane_select=hl.tbool,
                canonical=hl.tbool,
            ),
            key="transcript",
        )
        ht = ht.annotate(selected=get_transcript_filter_expr(ht, **kwargs))
        assert {r.transcript for r in ht.collect() if r.selected} == expected


class TestAddOeUpperRankAndBins:
    """Test the add_oe_upper_rank_and_bins function."""

    # transcript: (gene_id, mane_select, gamma upper, poisson upper)
    TRANSCRIPTS = {
        "ENST1": ("G1", True, 0.8, 0.1),
        "ENST2": ("G2", True, 0.2, 0.3),
        "ENST3": ("G3", True, 0.5, 0.2),
        "ENST4": ("G4", True, None, 0.4),
        # Not selected: gene G1 already has a MANE Select transcript.
        "ENST5": ("G1", False, 0.01, 0.01),
    }

    @pytest.fixture(scope="class")
    def rows(self) -> Dict[str, hl.Struct]:
        """Return the ranked constraint group keyed by transcript."""
        oe_info_type = hl.tstruct(
            oe_ci_discretized_poisson=CI_TYPE, oe_ci_gamma=CI_TYPE
        )
        rows = [
            dict(
                transcript=t,
                gene_id=g,
                mane_select=m,
                canonical=m,
                constraint_groups=[
                    hl.Struct(
                        oe_info=[
                            hl.Struct(
                                oe_ci_discretized_poisson=hl.Struct(lower=0.0, upper=p),
                                oe_ci_gamma=hl.Struct(lower=0.0, upper=gu),
                            )
                        ]
                    )
                ],
            )
            for t, (g, m, gu, p) in self.TRANSCRIPTS.items()
        ]
        ht = add_oe_upper_rank_and_bins(_constraint_groups_ht(rows, oe_info_type))
        return {r.transcript: r.constraint_groups[0] for r in ht.collect()}

    def test_gamma_ranks(self, rows: Dict[str, hl.Struct]) -> None:
        """Test ranks and bins among selected transcripts with a defined value."""
        # Three ranked values (ENST4's is missing), so bins use n = 3.
        for t, rank in [("ENST2", 0), ("ENST3", 1), ("ENST1", 2)]:
            r = rows[t].oe_ci_gamma_rank
            assert r.upper_rank == rank
            assert r.upper_bin_percentile == int(rank * 100 / 3)
            assert r.upper_bin_decile == int(rank * 10 / 3)
            assert r.upper_bin_sextile == int(rank * 6 / 3)

        assert rows["ENST4"].oe_ci_gamma_rank is None

    def test_poisson_ranked_separately(self, rows: Dict[str, hl.Struct]) -> None:
        """Test that the Poisson CI upper bound gets its own ranking."""
        ranks = {
            t: cg.oe_ci_discretized_poisson_rank.upper_rank
            for t, cg in rows.items()
            if cg.oe_ci_discretized_poisson_rank is not None
        }
        assert ranks == {"ENST1": 0, "ENST3": 1, "ENST2": 2, "ENST4": 3}

    def test_unselected_transcript_not_ranked(self, rows: Dict[str, hl.Struct]) -> None:
        """Test that a transcript outside the ranking set gets missing ranks."""
        cg = rows["ENST5"]
        assert cg.oe_ci_gamma_rank is None
        assert cg.oe_ci_discretized_poisson_rank is None


def _thresholds(offset: float = 0.0) -> Dict[Tuple[str, str], List[float]]:
    """Return evenly spaced thresholds on [0, 1] for every granularity and metric."""
    out = {}
    for gran, bins in CONSTRAINT_GRANULARITIES.items():
        n_bins = len(bins) + 1
        for metric in ["syn", "mis", "lof"]:
            shift = offset if metric == "lof" else 0.0
            out[(gran, metric)] = [b / n_bins + shift for b in bins]
    return out


class TestAnnotateConstraintPercentileBins:
    """Test the annotate_constraint_percentile_bins function."""

    METRIC_GROUP_IDX = {"syn": 0, "mis": 1, "lof": 4}

    @pytest.mark.parametrize("value", [0.0, 0.5, 0.555, 0.99, 1.5, None])
    def test_bins(self, value: Optional[float]) -> None:
        """Test that each bin is the number of thresholds met, missing if no value."""
        rows = [
            dict(
                transcript="T1",
                constraint_groups=[
                    hl.Struct(
                        oe_info=[
                            hl.Struct(oe_ci_gamma=hl.Struct(lower=0.0, upper=value))
                        ]
                    )
                ]
                * 5,
            )
        ]
        ht = _constraint_groups_ht(rows, hl.tstruct(oe_ci_gamma=CI_TYPE))
        # LoF uses its own thresholds, shifted up by 0.2.
        thresholds = _thresholds(offset=0.2)
        ht = annotate_constraint_percentile_bins(ht, thresholds, self.METRIC_GROUP_IDX)
        bins = ht.collect()[0].constraint_bins
        for gran in CONSTRAINT_GRANULARITIES:
            for metric in ["syn", "mis", "lof"]:
                assert bins[gran][metric] == _bin(value, thresholds[(gran, metric)])


class TestComputeConstraintPercentileBins:
    """Test the compute_constraint_percentile_bins function."""

    @pytest.fixture(scope="class")
    def result(self) -> Tuple[hl.Struct, List[hl.Struct]]:
        """Return the thresholds global and the binned rows.

        The 30 MANE Select transcripts have values from 0.05 to 1.5. The outliers
        (flagged) and unselected transcripts (neither MANE Select nor canonical)
        have a value of 100, so including either in the thresholds would show. Only
        the syn, mis, and HC LoF groups are included, since each group is ranked
        separately.
        """
        transcripts = (
            [(f"ENSTM{i}", 0.05 * (i + 1), True, set()) for i in range(30)]
            + [(f"ENSTO{i}", 100.0, True, {"outlier_lof"}) for i in range(5)]
            + [(f"ENSTU{i}", 100.0, False, set()) for i in range(5)]
        )
        oe_info_type = hl.tstruct(
            oe_ci_discretized_poisson=CI_TYPE, oe_ci_gamma=CI_TYPE
        )
        rows = [
            dict(
                transcript=t,
                gene_id=f"G_{t}",
                mane_select=mane,
                canonical=mane,
                constraint_flags=flags,
                constraint_groups=[
                    hl.Struct(
                        oe_info=[
                            hl.Struct(
                                oe_ci_discretized_poisson=hl.Struct(lower=0.0, upper=v),
                                oe_ci_gamma=hl.Struct(lower=0.0, upper=v),
                            )
                        ]
                    )
                ]
                * 3,
            )
            for t, v, mane, flags in transcripts
        ]
        meta = [{"csq_set": "syn"}, {"csq_set": "mis"}, {"lof": "hc"}]
        ht = compute_constraint_percentile_bins(
            _constraint_groups_ht(rows, oe_info_type, meta=meta)
        )
        return hl.eval(ht.percentile_thresholds), ht.collect()

    def test_thresholds(self, result: Tuple[hl.Struct, List[hl.Struct]]) -> None:
        """Test the thresholds global: shape, order, and excluded transcripts."""
        thresholds, _ = result
        for metric in ["syn", "mis", "lof"]:
            for gran, bins in CONSTRAINT_GRANULARITIES.items():
                t = thresholds[metric][gran]
                assert len(t) == len(bins)
                assert t == sorted(t)
                # Outliers and unselected transcripts (value 100) are excluded.
                assert max(t) <= 1.5

    def test_bins(self, result: Tuple[hl.Struct, List[hl.Struct]]) -> None:
        """Test that every transcript, including excluded ones, is binned."""
        thresholds, rows = result
        for row in rows:
            value = row.constraint_groups[0].oe_info[0].oe_ci_gamma.upper
            for gran in CONSTRAINT_GRANULARITIES:
                for metric in ["syn", "mis", "lof"]:
                    assert row.constraint_bins[gran][metric] == _bin(
                        value, thresholds[metric][gran]
                    )
            if value == 100:
                assert row.constraint_bins.decile.lof == 9

    def test_ranks(self, result: Tuple[hl.Struct, List[hl.Struct]]) -> None:
        """Test that only selected transcripts are ranked."""
        _, rows = result
        ranked = {
            r.transcript: r.constraint_groups[2].oe_ci_gamma_rank.upper_rank
            for r in rows
            if r.constraint_groups[2].oe_ci_gamma_rank is not None
        }
        assert not any(t.startswith("ENSTU") for t in ranked)
        assert sorted(ranked.values()) == list(range(35))
        assert ranked["ENSTM0"] == 0


class TestGetConstraintGroupName:
    """Test the _get_constraint_group_name function."""

    @pytest.mark.parametrize(
        "meta, expected",
        [
            ({"csq_set": "syn"}, "syn"),
            ({"csq_set": "mis"}, "mis"),
            ({"lof": "hc"}, "lof_hc"),
            ({"lof": "hc_lc"}, "lof_hc_lc"),
            ({"lof": "hc", "csq_set": "mis"}, "lof_hc_mis"),
        ],
    )
    def test_names(self, meta: Dict[str, str], expected: str) -> None:
        """Test the names derived from constraint group metadata."""
        assert _get_constraint_group_name(meta) == expected


class TestComputePliScores:
    """Test the _compute_pli_scores function."""

    def _ht(self) -> hl.Table:
        """Return transcripts whose HC and HC+LC LoF counts differ."""
        cg_type = hl.tstruct(
            oe_info=hl.tarray(
                hl.tstruct(observed_variants=hl.tint64, expected_variants=hl.tfloat64)
            ),
            z_score=hl.tfloat64,
        )

        def _cg(obs: int, exp: float) -> hl.Struct:
            return hl.Struct(
                oe_info=[hl.Struct(observed_variants=obs, expected_variants=exp)],
                z_score=0.0,
            )

        # (HC obs, HC exp, HC+LC obs, HC+LC exp): constrained, unconstrained,
        # recessive-like, no LoF variants expected, and one with HC+LC much less
        # constrained than HC.
        counts = [
            (0, 20.0, 0, 21.0),
            (20, 20.0, 22, 23.0),
            (7, 10.0, 8, 11.0),
            (0, 0.0, 0, 0.0),
            (0, 15.0, 30, 35.0),
        ] * 4
        rows = [
            dict(
                transcript=f"T{i}",
                constraint_groups=[
                    _cg(100, 100.0),
                    _cg(50, 60.0),
                    _cg(hc_obs, hc_exp),
                    _cg(hc_lc_obs, hc_lc_exp),
                    _cg(hc_obs, hc_exp),
                ],
            )
            for i, (hc_obs, hc_exp, hc_lc_obs, hc_lc_exp) in enumerate(counts)
        ]
        ht = hl.Table.parallelize(
            rows,
            hl.tstruct(transcript=hl.tstr, constraint_groups=hl.tarray(cg_type)),
            key="transcript",
        )
        ht = ht.annotate_globals(
            constraint_group_meta=CONSTRAINT_GROUP_META,
            exomes_freq_meta=[{"group": "adj"}],
        )
        return _compute_pli_scores(ht)

    def test_scores_per_group(self) -> None:
        """Test that each LoF group gets scores from its own counts.

        Regression test: ``lof_hc_lc`` used to get a copy of the HC LoF scores.
        """
        rows = self._ht().collect()
        names = [_get_constraint_group_name(m) for m in CONSTRAINT_GROUP_META]
        for row in rows:
            groups = dict(zip(names, row.constraint_groups))
            for name in ["syn", "mis", "lof_classic"]:
                assert groups[name].pLI is None
                assert groups[name].pNull is None
                assert groups[name].pRec is None

            for name in ["lof_hc", "lof_hc_lc"]:
                g = groups[name]
                if g.oe_info[0].expected_variants == 0:
                    assert g.pLI is None
                else:
                    assert g.pLI + g.pNull + g.pRec == pytest.approx(1.0)

        differs = [
            r
            for r in rows
            if r.constraint_groups[3].pLI is not None
            and r.constraint_groups[3].pLI != r.constraint_groups[4].pLI
        ]
        assert differs

    def test_constrained_transcript_is_lof_intolerant(self) -> None:
        """Test that a transcript with no observed of 20 expected LoF has a high pLI."""
        ht = self._ht()
        row = ht.filter(ht.transcript == "T0").collect()[0]
        assert row.constraint_groups[4].pLI > 0.9
        assert row.constraint_groups[4].pNull < 0.01


@requires_qgamma
class TestComputeConstraintMetrics:
    """Test the compute_constraint_metrics function."""

    # transcript: (HC LoF obs, HC LoF exp). Some are above expected, so every group
    # has negative raw z-scores for the mirrored SD.
    LOF = {
        f"ENST{i}": c
        for i, c in enumerate(
            [
                (0, 10.0),
                (1, 12.0),
                (0, 8.0),
                (2, 15.0),
                (5, 4.0),
                (7, 5.0),
                (3, 3.5),
                (1, 9.0),
                (0, 20.0),
                (4, 6.0),
                (9, 6.0),
                (2, 2.5),
            ]
        )
    }

    def _gencode_ht(self) -> hl.Table:
        """Return GENCODE features for ENST0 and the PAR transcript ENST1."""

        def _feature(t, contig, start, end, feature, par_y=False):
            return dict(
                interval=hl.Interval(
                    hl.Locus(contig, start, "GRCh38"),
                    hl.Locus(contig, end, "GRCh38"),
                    includes_end=True,
                ),
                feature=feature,
                transcript_id=t,
                transcript_id_version=f"{t}.3" + ("_PAR_Y" if par_y else ""),
                gene_id_version="ENSG0.1",
                level=2,
                transcript_type="protein_coding",
            )

        rows = [
            _feature("ENST0", "chr1", 1000, 5000, "transcript"),
            _feature("ENST0", "chr1", 1100, 1199, "CDS"),
            _feature("ENST0", "chr1", 1300, 1349, "CDS"),
            _feature("ENST1", "chrX", 10000, 20000, "transcript"),
            _feature("ENST1", "chrX", 10100, 10199, "CDS"),
            _feature("ENST1", "chrY", 10000, 20000, "transcript", par_y=True),
            _feature("ENST1", "chrY", 10100, 10199, "CDS", par_y=True),
            _feature("ENST1", "chrY", 10300, 10399, "CDS", par_y=True),
        ]
        return hl.Table.parallelize(
            rows,
            hl.tstruct(
                interval=hl.tinterval(hl.tlocus("GRCh38")),
                feature=hl.tstr,
                transcript_id=hl.tstr,
                transcript_id_version=hl.tstr,
                gene_id_version=hl.tstr,
                level=hl.tint32,
                transcript_type=hl.tstr,
            ),
            key="interval",
        )

    @pytest.fixture(scope="class")
    def result(self) -> Tuple[List[float], Dict[str, hl.Struct]]:
        """Return the sd_raw_z global and the metrics rows keyed by transcript."""
        oe_info_type = hl.tstruct(
            observed_variants=hl.tint64,
            predicted_proportion_observed=hl.tfloat64,
            expected_variants=hl.tfloat64,
        )

        def _cg(obs, exp):
            return hl.Struct(
                mu_snp=1e-7,
                mu=1e-7,
                possible_variants=100,
                coverage_correction=1.0,
                oe_info=[
                    hl.Struct(
                        observed_variants=obs,
                        predicted_proportion_observed=0.5,
                        expected_variants=exp,
                    )
                ],
            )

        rows = [
            dict(
                transcript=t,
                no_variants=False,
                constraint_groups=[
                    _cg(90 + 2 * i, 100.0),
                    _cg(180 + 4 * i, 200.0),
                    _cg(obs, exp),
                    _cg(obs + 2, exp + 2.0),
                    _cg(obs, exp),
                ],
            )
            for i, (t, (obs, exp)) in enumerate(self.LOF.items())
        ]
        ht = _constraint_groups_ht(
            rows,
            oe_info_type,
            {
                "mu_snp": hl.tfloat64,
                "mu": hl.tfloat64,
                "possible_variants": hl.tint64,
                "coverage_correction": hl.tfloat64,
            },
        )
        gene_quality_ht = hl.Table.parallelize(
            [
                dict(
                    transcript="ENST0",
                    gene_quality_metrics=hl.Struct(exome_prop_bp_AN90=0.95),
                    gene_flags={"low_exome_coverage"},
                )
            ],
            hl.tstruct(
                transcript=hl.tstr,
                gene_quality_metrics=hl.tstruct(exome_prop_bp_AN90=hl.tfloat64),
                gene_flags=hl.tset(hl.tstr),
            ),
            key="transcript",
        )
        ht = compute_constraint_metrics(ht, self._gencode_ht(), gene_quality_ht)
        return hl.eval(ht.sd_raw_z), {r.transcript: r for r in ht.collect()}

    def test_metrics(self, result: Tuple[List[float], Dict[str, hl.Struct]]) -> None:
        """Test that each group gets oe, CIs, flags, z_score, and LoF pLI."""
        sd, rows = result
        assert set(rows) == set(self.LOF)
        for row in rows.values():
            for i, cg in enumerate(row.constraint_groups):
                x = cg.oe_info[0]
                assert x.oe == pytest.approx(x.observed_variants / x.expected_variants)
                assert x.oe_ci_gamma.lower < x.oe_ci_gamma.upper
                assert x.oe <= x.oe_ci_gamma.upper
                assert cg.z_score == pytest.approx(x.z_raw / sd[i])
                has_pli = GROUP_NAMES[i] in ("lof_hc", "lof_hc_lc")
                assert (cg.pLI is not None) == has_pli
            assert row.constraint_flags == set()

    def test_annotations(
        self, result: Tuple[List[float], Dict[str, hl.Struct]]
    ) -> None:
        """Test the gene quality and GENCODE annotations, dropping chrY PAR copies."""
        _, rows = result
        enst0, enst1, enst2 = rows["ENST0"], rows["ENST1"], rows["ENST2"]

        assert enst0.gene_quality_metrics.exome_prop_bp_AN90 == 0.95
        assert enst0.gene_flags == {"low_exome_coverage"}
        assert enst1.gene_quality_metrics is None

        assert enst0.chromosome == "chr1"
        assert enst0.transcript_id_version == "ENST0.3"
        assert (enst0.cds_length, enst0.num_coding_exons) == (150, 2)
        assert (enst1.chromosome, enst1.cds_length, enst1.num_coding_exons) == (
            "chrX",
            100,
            1,
        )
        assert enst2.chromosome is None


class TestPrepareReleaseHt:
    """Test the values in the release Table built by prepare_release_ht."""

    # Downsamplings deliberately out of order, to check the arrays follow the
    # downsamplings global.
    FREQ_META = [
        {"group": "adj"},
        {"group": "adj", "gen_anc": "global", "downsampling": "100"},
        {"group": "adj", "gen_anc": "global", "downsampling": "10"},
        {"group": "adj", "gen_anc": "afr", "downsampling": "10"},
    ]

    @pytest.fixture(scope="class")
    def internal_ht(self) -> hl.Table:
        """Return an internal metrics Table with distinct values for every group."""
        oe_info_type = hl.tstruct(
            observed_variants=hl.tint64,
            expected_variants=hl.tfloat64,
            oe=hl.tfloat64,
            oe_ci_discretized_poisson=CI_TYPE,
            oe_ci_gamma=CI_TYPE,
            z_raw=hl.tfloat64,
        )
        cg_type = hl.tstruct(
            mu_snp=hl.tfloat64,
            mu=hl.tfloat64,
            possible_variants=hl.tint64,
            oe_info=hl.tarray(oe_info_type),
            z_score=hl.tfloat64,
            oe_ci_gamma_rank=RANK_TYPE,
            pNull=hl.tfloat64,
            pRec=hl.tfloat64,
            pLI=hl.tfloat64,
        )

        def _oe_info(i, obs, exp):
            return hl.Struct(
                observed_variants=obs,
                expected_variants=exp,
                oe=obs / exp,
                oe_ci_discretized_poisson=hl.Struct(lower=9.0, upper=9.0),
                oe_ci_gamma=hl.Struct(lower=0.1 * (i + 1), upper=0.2 * (i + 1)),
                z_raw=float(i),
            )

        def _cg(i, possible):
            # oe_info: adj, then the downsamplings, where obs is the downsampling.
            ds = [int(m["downsampling"]) for m in self.FREQ_META[1:]]
            has_pli = GROUP_NAMES[i] in ("lof_hc", "lof_hc_lc")
            return hl.Struct(
                mu_snp=float(i + 1),
                mu=99.0,
                possible_variants=possible,
                oe_info=[_oe_info(i, 10 * (i + 1), 20.0 * (i + 1))]
                + [_oe_info(i, d, 1.5 * d) for d in ds],
                z_score=float(i) / 2,
                oe_ci_gamma_rank=hl.Struct(
                    upper_rank=i,
                    upper_bin_percentile=i,
                    upper_bin_decile=i,
                    upper_bin_sextile=i,
                ),
                pNull=0.5 if has_pli else None,
                pRec=0.5 - 0.01 * (i + 1) if has_pli else None,
                pLI=0.01 * (i + 1) if has_pli else None,
            )

        def _row(t, possible):
            return dict(
                gene="GENE",
                gene_id="ENSG1",
                transcript=t,
                canonical=True,
                mane_select=True,
                constraint_groups=[_cg(i, p) for i, p in enumerate(possible)],
                no_variants=False,
                transcript_id_version=f"{t}.2",
                gene_id_version="ENSG1.1",
                level=2,
                transcript_type="protein_coding",
                chromosome="chr1",
                start_position=100,
                end_position=900,
                cds_length=300,
                num_coding_exons=3,
                gene_quality_metrics=hl.Struct(exome_prop_bp_AN90=0.9),
                gene_flags=set(),
                constraint_flags={"outlier_mis"},
            )

        rows = [
            _row("ENST1", [100, 200, 300, 400, 500]),
            # Nothing possible in any group: dropped.
            _row("ENST2", [0, 0, 0, 0, 0]),
            # Possible variants in one group only: kept.
            _row("ENST3", [0, 7, 0, 0, 0]),
        ]
        ht = hl.Table.parallelize(
            rows,
            hl.tstruct(
                gene=hl.tstr,
                gene_id=hl.tstr,
                transcript=hl.tstr,
                canonical=hl.tbool,
                mane_select=hl.tbool,
                constraint_groups=hl.tarray(cg_type),
                no_variants=hl.tbool,
                transcript_id_version=hl.tstr,
                gene_id_version=hl.tstr,
                level=hl.tint32,
                transcript_type=hl.tstr,
                chromosome=hl.tstr,
                start_position=hl.tint32,
                end_position=hl.tint32,
                cds_length=hl.tint64,
                num_coding_exons=hl.tint64,
                gene_quality_metrics=hl.tstruct(exome_prop_bp_AN90=hl.tfloat64),
                gene_flags=hl.tset(hl.tstr),
                constraint_flags=hl.tset(hl.tstr),
            ),
            key=["gene", "transcript", "canonical"],
        )
        thresholds = hl.struct(
            percentile=[0.5] * 99, decile=[0.5] * 9, sextile=[0.5] * 5
        )
        return ht.annotate_globals(
            version="4.1.1",
            calculate_mu_globals=hl.struct(
                freq_meta=[{"group": "adj"}],
                genetic_ancestry_groups=["global"],
                downsampling_idx=1,
                ac_cutoff=5,
            ),
            build_models_globals=hl.struct(
                synonymous_transcript_filter_field="mane_select",
                skip_coverage_model=False,
                high_cov_cutoff=40,
            ),
            apply_models_globals=hl.struct(
                skip_coverage_model=False,
                groupings=["gene"],
                log10_coverage=True,
            ),
            exomes_freq_meta=self.FREQ_META,
            downsamplings=[10, 100],
            max_af=0.001,
            constraint_group_meta=CONSTRAINT_GROUP_META,
            sd_raw_z=[1.0, 2.0, 3.0, 4.0, 5.0],
            percentile_thresholds=hl.struct(
                syn=thresholds,
                mis=thresholds,
                lof=hl.struct(
                    percentile=[0.1] * 99, decile=[0.2] * 9, sextile=[0.3] * 5
                ),
            ),
        )

    @pytest.fixture(scope="class")
    def release(self, internal_ht: hl.Table) -> Tuple[hl.Table, Dict[str, hl.Struct]]:
        """Return the release Table and its rows keyed by transcript."""
        ht = prepare_release_ht(internal_ht)
        return ht, {r.transcript: r for r in ht.collect()}

    def test_rows_and_key(self, release: Tuple[hl.Table, Dict[str, hl.Struct]]) -> None:
        """Test the key, the row fields, and that rows with nothing possible drop."""
        ht, rows = release
        assert list(ht.key) == RELEASE_KEY_ORDER
        assert set(ht.row_value) == set(
            RELEASE_TOP_LEVEL_ANNOTATIONS + RELEASE_GROUP_NAMES
        )
        assert set(rows) == {"ENST1", "ENST3"}
        assert rows["ENST1"].transcript_version == "ENST1.2"
        assert rows["ENST1"].transcript_level == 2
        assert rows["ENST1"].constraint_flags == {"outlier_mis"}

    @pytest.mark.parametrize("name, idx", [("syn", 0), ("mis", 1), ("lof_hc_lc", 3)])
    def test_unranked_group_values(
        self, release: Tuple[hl.Table, Dict[str, hl.Struct]], name: str, idx: int
    ) -> None:
        """Test that unranked groups get their own values and only the gamma CI."""
        g = release[1]["ENST1"][name]
        assert g.mu == idx + 1
        assert g.possible == 100 * (idx + 1)
        assert (g.obs, g.exp) == (10 * (idx + 1), 20.0 * (idx + 1))
        assert (g.z_raw, g.z_score) == (float(idx), idx / 2)
        assert dict(g.oe_ci) == {
            "lower": pytest.approx(0.1 * (idx + 1)),
            "upper": pytest.approx(0.2 * (idx + 1)),
        }

    def test_lof_group(self, release: Tuple[hl.Table, Dict[str, hl.Struct]]) -> None:
        """Test that lof_hc becomes lof, with the gamma CI and its rank fields."""
        lof = release[1]["ENST1"].lof
        assert lof.mu == 5.0
        assert lof.oe_ci.lower == pytest.approx(0.5)
        assert lof.oe_ci.upper == pytest.approx(1.0)
        assert lof.oe_ci.upper_rank == 4
        assert lof.oe_ci.upper_bin_decile == 4

    def test_pli_fields(self, release: Tuple[hl.Table, Dict[str, hl.Struct]]) -> None:
        """Test that each LoF group has its own pLI and other groups have none."""
        ht, rows = release
        row = rows["ENST1"]
        assert row.lof.pLI == pytest.approx(0.05)
        assert row.lof_hc_lc.pLI == pytest.approx(0.04)
        assert "pLI" not in ht.syn.dtype and "pLI" not in ht.mis.dtype
        assert "upper_rank" not in ht.lof_hc_lc.oe_ci.dtype

    def test_downsampling_arrays(
        self, release: Tuple[hl.Table, Dict[str, hl.Struct]]
    ) -> None:
        """Test that gen_anc arrays are parallel to the downsamplings global."""
        ht, rows = release
        downsamplings = dict(hl.eval(ht.downsamplings))
        assert downsamplings == {"global": [100, 10], "afr": [10]}
        row = rows["ENST1"]
        for name in RELEASE_GROUP_NAMES:
            for gen_anc, levels in downsamplings.items():
                assert row[name].gen_anc_obs[gen_anc] == levels
                assert row[name].gen_anc_exp[gen_anc] == [1.5 * d for d in levels]

    @pytest.mark.parametrize(
        "release_version, expected", [(None, "4.1.1"), ("4.1.2", "4.1.2")]
    )
    def test_globals(
        self, internal_ht: hl.Table, release_version: Optional[str], expected: str
    ) -> None:
        """Test the release globals."""
        g = hl.eval(prepare_release_ht(internal_ht, release_version).globals)
        assert set(g) == {
            "version",
            "calculate_mu_params",
            "build_models_params",
            "apply_models_params",
            "downsamplings",
            "max_af",
            "sd_raw_z",
            "loeuf_percentile_thresholds",
        }
        assert g.version == expected
        assert dict(g.calculate_mu_params) == {"ac_cutoff": 5}
        assert dict(g.build_models_params) == {"high_cov_cutoff": 40}
        assert dict(g.apply_models_params) == {"log10_coverage": True}
        assert dict(g.sd_raw_z) == {
            "syn": 1.0,
            "mis": 2.0,
            "lof_hc_lc": 4.0,
            "lof": 5.0,
        }
        assert g.loeuf_percentile_thresholds.decile == [0.2] * 9
        assert g.max_af == 0.001


class TestReleaseColumns:
    """Test the release Table and TSV produced from the internal metrics Table."""

    def _internal_metrics_ht(self) -> hl.Table:
        """Return an empty internal metrics Table with the pipeline's schema and globals.

        The row schema is stored in ``tests/data/constraint_metrics_row_type.txt``.
        Only the schema and globals determine the release columns, so the Table has
        one row of missing values, which ``prepare_release_ht`` filters out.
        """
        with open(os.path.join(DATA_DIR, "constraint_metrics_row_type.txt")) as f:
            row_type = hl.dtype(f.read())

        ht = hl.utils.range_table(1)
        ht = ht.annotate(**{k: hl.missing(v) for k, v in row_type.items()})
        ht = ht.key_by(
            "gene", "gene_id", "transcript", "canonical", "mane_select"
        ).drop("idx")

        freq_meta = [{"group": "adj"}] + [
            {"group": "adj", "gen_anc": g, "downsampling": str(ds)}
            for g in GEN_ANCS
            for ds in (10, 100)
        ]
        thresholds = hl.struct(
            percentile=[0.5] * 99, decile=[0.5] * 9, sextile=[0.5] * 5
        )
        return ht.annotate_globals(
            calculate_mu_globals=hl.struct(
                freq_meta=[{"group": "adj"}],
                genetic_ancestry_groups=GEN_ANCS,
                downsampling_idx=1,
            ),
            build_models_globals=hl.struct(
                synonymous_transcript_filter_field="mane_select",
                skip_coverage_model=False,
            ),
            apply_models_globals=hl.struct(
                skip_coverage_model=False, groupings=hl.tuple(["gene"])
            ),
            exomes_freq_meta=freq_meta,
            downsamplings=[10, 100],
            max_af=0.001,
            constraint_group_meta=CONSTRAINT_GROUP_META,
            sd_raw_z=[1.0] * len(CONSTRAINT_GROUP_META),
            percentile_thresholds=hl.struct(
                syn=thresholds, mis=thresholds, lof=thresholds
            ),
        )

    def test_release_tsv_columns(self) -> None:
        """Test that the flattened release Table has the published TSV columns.

        ``tests/data/release_tsv_columns.txt`` is the header of the published v4.1.2
        release TSV, including the array-valued ``gen_anc_obs``/``gen_anc_exp``
        downsampling columns.
        """
        with open(os.path.join(DATA_DIR, "release_tsv_columns.txt")) as f:
            expected = f.read().split()

        release_ht = prepare_release_ht(self._internal_metrics_ht(), "test")
        assert list(flatten_release_ht(release_ht).row) == expected


class TestPrepareReleaseMutationHt:
    """Test the prepare_release_mutation_ht function."""

    def _ht(self, version: Optional[str]) -> hl.Table:
        """Return a mutation rate Table with internal fields and globals."""
        ht = hl.Table.parallelize(
            [
                dict(
                    context="ACG",
                    ref="C",
                    alt="T",
                    methylation_level=2,
                    mu_snp=1e-7,
                    possible_variants=10,
                    cpg=True,
                    transition=True,
                    mutation_type="CpG",
                )
            ],
            hl.tstruct(
                context=hl.tstr,
                ref=hl.tstr,
                alt=hl.tstr,
                methylation_level=hl.tint32,
                mu_snp=hl.tfloat64,
                possible_variants=hl.tint64,
                cpg=hl.tbool,
                transition=hl.tbool,
                mutation_type=hl.tstr,
            ),
            key=["context", "ref", "alt", "methylation_level"],
        )
        ht = ht.annotate_globals(
            calculate_mu_globals=hl.struct(
                freq_meta=[{"group": "adj"}],
                genetic_ancestry_groups=["global"],
                downsampling_idx=1,
                ac_cutoff=5,
            ),
            exomes_freq_meta=[{"group": "adj"}],
        )
        return ht.annotate_globals(version=version) if version else ht

    @pytest.mark.parametrize(
        "version, release_version, expected",
        [("4.1.1", None, "4.1.1"), ("4.1.1", "4.1.2", "4.1.2"), (None, None, None)],
    )
    def test_release_mutation_ht(
        self,
        version: Optional[str],
        release_version: Optional[str],
        expected: Optional[str],
    ) -> None:
        """Test the selected fields and the release globals."""
        ht = prepare_release_mutation_ht(self._ht(version), release_version)
        assert list(ht.row_value) == ["mu", "cpg", "transition", "mutation_type"]
        assert ht.collect()[0].mu == 1e-7

        g = hl.eval(ht.globals)
        assert dict(g.calculate_mu_params) == {"ac_cutoff": 5}
        assert g.get("version") == expected
        assert "exomes_freq_meta" not in g


class TestLofBinThresholdsToHt:
    """Test the lof_bin_thresholds_to_ht function."""

    def test_rows(self) -> None:
        """Test one row per granularity and 1-based bin, in global order."""
        thresholds = {
            gran: [i / 1000 for i in range(len(bins))]
            for gran, bins in CONSTRAINT_GRANULARITIES.items()
        }
        ht = hl.utils.range_table(1).annotate_globals(
            loeuf_percentile_thresholds=hl.struct(**thresholds)
        )
        rows = lof_bin_thresholds_to_ht(ht).collect()
        assert [(r.granularity, r.bin, r.threshold) for r in rows] == [
            (gran, i + 1, t)
            for gran, ts in thresholds.items()
            for i, t in enumerate(ts)
        ]
