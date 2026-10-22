"""Tests for the apply-models and aggregation functions in ``gnomad_constraint.utils.constraint``."""

import functools
import math
from typing import Dict, List, Optional, Tuple

import hail as hl
import pytest

from gnomad_constraint.resources.constants import AGGREGATE_SUM_FIELDS, MU_GROUPING
from gnomad_constraint.utils.constraint import (
    _annotate_apply_models_globals,
    _apply_constraint_models,
    _prepare_ht_for_apply_models,
    aggregate_by_constraint_groups,
    aggregate_per_variant_expected_ht,
    create_aggregated_expected_ht,
    create_per_variant_expected_ht,
)

CSQ_TYPE = hl.tstruct(
    most_severe_consequence=hl.tstr,
    lof=hl.tstr,
    polyphen_prediction=hl.tstr,
    gene_symbol=hl.tstr,
    gene_id=hl.tstr,
    transcript_id=hl.tstr,
    canonical=hl.tint32,
    mane_select=hl.tstr,
)

APPLY_MODEL_TYPE = hl.tstruct(
    high_or_low_coverage=hl.tstr,
    model_group=hl.tstruct(cpg=hl.tbool, genomic_region=hl.tstr),
)

PREPROCESSED_TYPE = hl.tstruct(
    locus=hl.tlocus("GRCh38"),
    alleles=hl.tarray(hl.tstr),
    context=hl.tstr,
    ref=hl.tstr,
    alt=hl.tstr,
    methylation_level=hl.tint32,
    cpg=hl.tbool,
    transition=hl.tbool,
    mutation_type=hl.tstr,
    mutation_type_model=hl.tstr,
    exomes_coverage=hl.tint32,
    vep=hl.tstruct(transcript_consequences=hl.tarray(CSQ_TYPE)),
    calibrate_mu=hl.tstruct(
        observed_variants=hl.tarray(hl.tint32),
        possible_variants=hl.tint32,
        apply_model=APPLY_MODEL_TYPE,
    ),
)

NON_CPG = {
    "context": "ATG",
    "ref": "T",
    "alt": "C",
    "methylation_level": 0,
    "cpg": False,
    "transition": True,
    "mutation_type": "non-CpG transition",
    "mutation_type_model": "non-CpG transition",
}
CPG = {
    "context": "ACG",
    "ref": "C",
    "alt": "T",
    "methylation_level": 2,
    "cpg": True,
    "transition": True,
    "mutation_type": "CpG",
    "mutation_type_model": "CpG",
}
MU_SNP = {False: 1e-8, True: 1e-7}
"""Mutation rate of the non-CpG and CpG substitutions in the fixture."""

PLATEAU_MODELS = {
    hl.Struct(cpg=False, genomic_region="autosome_or_par"): [[0.1, 2e7], [0.05, 1e7]],
    hl.Struct(cpg=True, genomic_region="autosome_or_par"): [[0.2, 1e6], [0.1, 5e5]],
}
"""One [intercept, slope] model per frequency group, so both give ppo [0.3, 0.15]."""

COVERAGE_MODEL = (0.1, 0.5)
"""[intercept, slope] of the coverage model on log10 coverage."""

HIGH_COV = 60
LOW_COV = 20
LOW_COV_CORRECTION = COVERAGE_MODEL[0] + COVERAGE_MODEL[1] * math.log10(LOW_COV)

GROUPINGS = [
    "annotation",
    "modifier",
    "gene",
    "gene_id",
    "transcript",
    "canonical",
    "mane_select",
]


def _csq(
    consequence: str,
    transcript: str = "T1",
    lof: Optional[str] = None,
    canonical: Optional[int] = 1,
    mane_select: Optional[str] = "NM_1",
) -> hl.Struct:
    """Return a VEP transcript consequence for gene G1."""
    return hl.Struct(
        most_severe_consequence=consequence,
        lof=lof,
        polyphen_prediction=None,
        gene_symbol="G1",
        gene_id="ENSG1",
        transcript_id=transcript,
        canonical=canonical,
        mane_select=mane_select,
    )


def _variant(
    pos: int,
    sub: Dict,
    coverage: Optional[int],
    observed: Optional[List[int]],
    possible: Optional[int],
    csqs: List[hl.Struct],
    has_apply_model: bool = True,
) -> Dict:
    """Return a preprocessed context row."""
    apply_model = None
    if has_apply_model and coverage is not None:
        apply_model = hl.Struct(
            high_or_low_coverage="high" if coverage >= 40 else "low",
            model_group=hl.Struct(cpg=sub["cpg"], genomic_region="autosome_or_par"),
        )
    return dict(
        locus=hl.Locus("chr1", pos, "GRCh38"),
        alleles=[sub["ref"], sub["alt"]],
        **sub,
        exomes_coverage=coverage,
        vep=hl.Struct(transcript_consequences=csqs),
        calibrate_mu=hl.Struct(
            observed_variants=observed,
            possible_variants=possible,
            apply_model=apply_model,
        ),
    )


def _preprocessed_ht() -> hl.Table:
    """Return a small Table shaped like ``prepare_ht_for_constraint_calculations``.

    One gene with a MANE Select transcript (T1) and a second transcript (T2).
    Synonymous variants are at high and low coverage, so the coverage model
    changes their expected count. The last four rows are dropped before the
    models are applied, or by the coding-consequence filter.
    """
    rows = [
        _variant(
            1,
            NON_CPG,
            HIGH_COV,
            [1, 0],
            1,
            [
                _csq("synonymous_variant"),
                _csq("missense_variant", "T2", canonical=None, mane_select=None),
            ],
        ),
        _variant(2, NON_CPG, LOW_COV, [0, 0], 1, [_csq("synonymous_variant")]),
        _variant(3, CPG, HIGH_COV, [1, 1], 1, [_csq("missense_variant")]),
        _variant(4, CPG, LOW_COV, [0, 0], 1, [_csq("stop_gained", lof="HC")]),
        _variant(5, NON_CPG, HIGH_COV, [1, 0], 1, [_csq("stop_gained", lof="LC")]),
        # No apply model: outside both coverage models.
        _variant(
            6,
            NON_CPG,
            30,
            [0, 0],
            1,
            [_csq("synonymous_variant")],
            has_apply_model=False,
        ),
        # Not possible: e.g. a common variant.
        _variant(7, NON_CPG, HIGH_COV, [0, 0], 0, [_csq("synonymous_variant")]),
        # Exome coverage undefined: observed and possible are missing.
        _variant(8, NON_CPG, None, None, None, [_csq("synonymous_variant")]),
        # Non-coding consequence.
        _variant(9, NON_CPG, HIGH_COV, [1, 0], 1, [_csq("intron_variant")]),
    ]
    return hl.Table.parallelize(rows, PREPROCESSED_TYPE, key=["locus", "alleles"])


def _mutation_ht() -> hl.Table:
    """Return a mutation rate Table for the two substitutions in the fixture."""
    return hl.Table.parallelize(
        [
            {**{k: sub[k] for k in MU_GROUPING}, "mu_snp": MU_SNP[sub["cpg"]]}
            for sub in [NON_CPG, CPG]
        ],
        hl.tstruct(
            context=hl.tstr,
            ref=hl.tstr,
            alt=hl.tstr,
            methylation_level=hl.tint32,
            mu_snp=hl.tfloat64,
        ),
        key=list(MU_GROUPING),
    )


def _coverage_model_expr(
    coverage_model: Optional[Tuple[float, float]],
) -> Optional[hl.expr.ArrayExpression]:
    """Return the coverage model as a Hail expression, or None to skip it."""
    return None if coverage_model is None else hl.literal(list(coverage_model))


def _rows_by_group(ht: hl.Table) -> Dict[Tuple[str, str, str], hl.Struct]:
    """Return rows keyed by (transcript, annotation, modifier)."""
    return {(r.transcript, r.annotation, r.modifier): r for r in ht.collect()}


@functools.lru_cache(maxsize=None)
def _prepared(use_mane_select: bool) -> Tuple[List[str], List[str], List[hl.Struct]]:
    """Return the groupings, row fields, and rows of ``_prepare_ht_for_apply_models``.

    Results are cached because each Hail run costs more than the assertions.
    """
    ht, groupings = _prepare_ht_for_apply_models(
        _preprocessed_ht(), use_mane_select=use_mane_select
    )
    return list(groupings), list(ht.row), ht.collect()


@functools.lru_cache(maxsize=None)
def _per_variant_expected(coverage_model: Optional[Tuple[float, float]]) -> Dict:
    """Return the fields, globals, and rows of ``create_per_variant_expected_ht``."""
    ht = create_per_variant_expected_ht(
        _preprocessed_ht(),
        _mutation_ht(),
        hl.literal(PLATEAU_MODELS),
        _coverage_model_expr(coverage_model),
        use_mane_select=True,
    )
    return {
        "fields": list(ht.row),
        "params": hl.eval(ht.apply_models_globals),
        "rows": {(r.locus.position, r.transcript): r for r in ht.collect()},
    }


@functools.lru_cache(maxsize=None)
def _path_sums(path: str, coverage_model: Optional[Tuple[float, float]]) -> Dict:
    """Return the key, row type, globals, and rows of one path's per-group sums.

    ``path`` is "per_variant" (apply the models per variant, then aggregate) or
    "aggregated" (aggregate the counts, then apply the models).
    """
    if path == "per_variant":
        ht = create_per_variant_expected_ht(
            _preprocessed_ht(),
            _mutation_ht(),
            hl.literal(PLATEAU_MODELS),
            _coverage_model_expr(coverage_model),
            use_mane_select=True,
        )
        ht = aggregate_per_variant_expected_ht(ht)
    else:
        ht = create_aggregated_expected_ht(
            _preprocessed_ht(),
            _mutation_ht(),
            hl.literal(PLATEAU_MODELS),
            _coverage_model_expr(coverage_model),
            use_mane_select=True,
        )
    return {
        "key": list(ht.key),
        "row_type": ht.row_value.dtype,
        "params": hl.eval(ht.apply_models_globals),
        "rows": _rows_by_group(ht),
    }


class TestPrepareHtForApplyModels:
    """Test the _prepare_ht_for_apply_models function."""

    @pytest.mark.parametrize(
        "use_mane_select, expected_groupings",
        [(True, GROUPINGS), (False, GROUPINGS[:-1])],
    )
    def test_groupings(
        self, use_mane_select: bool, expected_groupings: List[str]
    ) -> None:
        """Test the VEP grouping fields, with and without the MANE Select group."""
        groupings, fields, _ = _prepared(use_mane_select)
        assert groupings == expected_groupings
        assert all(g in fields for g in expected_groupings)
        assert "vep" not in fields

    def test_filters_and_explodes(self) -> None:
        """Test that unusable variants are dropped and consequences are exploded.

        Variants with no apply model, no possible count, or undefined exome
        coverage are dropped. A variant with two transcript consequences gives two
        rows.
        """
        _, _, rows = _prepared(True)
        assert sorted((r.locus.position, r.transcript) for r in rows) == [
            (1, "T1"),
            (1, "T2"),
            (2, "T1"),
            (3, "T1"),
            (4, "T1"),
            (5, "T1"),
            (9, "T1"),
        ]

    def test_promotes_calibrate_mu_fields(self) -> None:
        """Test that the calibrate_mu fields are available at the top level."""
        _, fields, _ = _prepared(True)
        for field in ["observed_variants", "possible_variants", "apply_model"]:
            assert field in fields

    def test_worst_csq_by_gene_with_mane_select(self) -> None:
        """Test that MANE Select grouping is rejected for worst_csq_by_gene."""
        with pytest.raises(ValueError, match="mane_select"):
            _prepare_ht_for_apply_models(
                _preprocessed_ht(),
                custom_vep_annotation="worst_csq_by_gene",
                use_mane_select=True,
            )


class TestApplyConstraintModels:
    """Test the _apply_constraint_models and _annotate_apply_models_globals functions."""

    PLATEAU_MODELS = {
        hl.Struct(cpg=False): [[0.1, 2e7]],
        hl.Struct(cpg=True): [[0.2, 1e7]],
    }

    def _ht(self, coverage_groups: List[str]) -> hl.Table:
        """Return sites at high coverage (60) or low coverage (20)."""
        rows = [
            dict(
                idx=i,
                mu_snp=1e-8,
                possible_variants=1,
                exomes_coverage=60 if g == "high" else 20,
                apply_model=hl.Struct(
                    high_or_low_coverage=g, model_group=hl.Struct(cpg=False)
                ),
            )
            for i, g in enumerate(coverage_groups)
        ]
        return hl.Table.parallelize(
            rows,
            hl.tstruct(
                idx=hl.tint32,
                mu_snp=hl.tfloat64,
                possible_variants=hl.tint64,
                exomes_coverage=hl.tint32,
                apply_model=hl.tstruct(
                    high_or_low_coverage=hl.tstr,
                    model_group=hl.tstruct(cpg=hl.tbool),
                ),
            ),
            key="idx",
        )

    def test_without_coverage_model(self) -> None:
        """Test that skipping the coverage model gives plateau-only expected counts.

        With ``--skip-coverage-model``, preprocessing assigns sites only to the high
        coverage model, so the correction is 1 and ``coverage_correction`` is still
        present for the downstream sums.
        """
        ht = _apply_constraint_models(
            self._ht(["high", "high"]), hl.literal(self.PLATEAU_MODELS), None
        )
        for row in ht.collect():
            assert row.coverage_correction == 1.0
            assert row.mu == pytest.approx(1e-8)
            assert row.expected_variants[0] == pytest.approx(0.3)

    def test_with_coverage_model(self) -> None:
        """Test that the coverage model corrects only low coverage sites."""
        ht = _apply_constraint_models(
            self._ht(["high", "low"]),
            hl.literal(self.PLATEAU_MODELS),
            hl.literal([0.1, 0.5]),
        )
        high, low = ht.order_by("idx").collect()
        assert high.coverage_correction == 1.0
        expected_correction = 0.1 + 0.5 * hl.eval(hl.log10(20))
        assert low.coverage_correction == pytest.approx(expected_correction)
        assert low.expected_variants[0] == pytest.approx(0.3 * expected_correction)

    @pytest.mark.parametrize(
        "coverage_model, expected", [(None, None), ([0.1, 0.5], [0.1, 0.5])]
    )
    def test_globals(
        self, coverage_model: Optional[List[float]], expected: Optional[List[float]]
    ) -> None:
        """Test that the coverage model global is a typed missing value when skipped."""
        ht = _annotate_apply_models_globals(
            self._ht(["high"]),
            hl.literal(self.PLATEAU_MODELS),
            hl.literal(coverage_model) if coverage_model else None,
            True,
            ["gene"],
        )
        assert ht.apply_models_globals.coverage_model.dtype == hl.tarray(hl.tfloat64)
        assert hl.eval(ht.apply_models_globals.coverage_model) == expected


class TestCreatePerVariantExpectedHt:
    """Test the create_per_variant_expected_ht function."""

    @pytest.mark.parametrize(
        "coverage_model, low_cov_correction",
        [(COVERAGE_MODEL, LOW_COV_CORRECTION), (None, 1.0)],
    )
    def test_expected_per_variant(
        self,
        coverage_model: Optional[Tuple[float, float]],
        low_cov_correction: float,
    ) -> None:
        """Test mu, expected, and coverage correction at high and low coverage.

        Without a coverage model, the low coverage site is not corrected.
        """
        rows = _per_variant_expected(coverage_model)["rows"]
        high, low = rows[(1, "T1")], rows[(2, "T1")]

        assert high.mu_snp == pytest.approx(MU_SNP[False])
        assert high.coverage_correction == pytest.approx(1.0)
        assert high.mu == pytest.approx(MU_SNP[False])
        assert high.predicted_proportion_observed == pytest.approx([0.3, 0.15])
        assert high.expected_variants == pytest.approx([0.3, 0.15])

        assert low.coverage_correction == pytest.approx(low_cov_correction)
        assert low.mu == pytest.approx(MU_SNP[False] * low_cov_correction)
        assert low.expected_variants == pytest.approx(
            [0.3 * low_cov_correction, 0.15 * low_cov_correction]
        )

    def test_drops_promoted_calibrate_mu_fields(self) -> None:
        """Test that only the nested calibrate_mu struct keeps the counts."""
        fields = _per_variant_expected(COVERAGE_MODEL)["fields"]
        assert "calibrate_mu" in fields
        for field in ["observed_variants", "possible_variants", "apply_model"]:
            assert field not in fields

    @pytest.mark.parametrize("coverage_model", [COVERAGE_MODEL, None])
    def test_globals(self, coverage_model: Optional[Tuple[float, float]]) -> None:
        """Test the model parameters and groupings stored in the globals."""
        params = _per_variant_expected(coverage_model)["params"]
        assert params.coverage_model == (
            None if coverage_model is None else list(coverage_model)
        )
        assert params.log10_coverage is True
        assert list(params.groupings) == GROUPINGS


class TestAggregatePerVariantExpectedHt:
    """Test the aggregate_per_variant_expected_ht function."""

    def test_sums_and_coding_filter(self) -> None:
        """Test the sums per transcript and consequence, and the coding filter."""
        rows = _path_sums("per_variant", COVERAGE_MODEL)["rows"]
        assert set(rows) == {
            ("T1", "synonymous_variant", "None"),
            ("T1", "missense_variant", "None"),
            ("T1", "stop_gained", "HC"),
            ("T1", "stop_gained", "LC"),
            ("T2", "missense_variant", "None"),
        }

        syn = rows[("T1", "synonymous_variant", "None")]
        assert syn.observed_variants == [1, 0]
        assert syn.possible_variants == 2
        assert syn.mu_snp == pytest.approx(2 * MU_SNP[False])
        assert syn.coverage_correction == pytest.approx(1 + LOW_COV_CORRECTION)
        assert syn.expected_variants == pytest.approx(
            [0.3 * (1 + LOW_COV_CORRECTION), 0.15 * (1 + LOW_COV_CORRECTION)]
        )

    def test_unobserved_group_sums_to_zero(self) -> None:
        """Test that a group whose only variant is unobserved sums to 0."""
        rows = _path_sums("per_variant", COVERAGE_MODEL)["rows"]
        assert rows[("T1", "stop_gained", "HC")].observed_variants == [0, 0]

    def test_include_mu_annotations_in_grouping(self) -> None:
        """Test that the mutation rate fields can be added to the grouping."""
        ht = create_per_variant_expected_ht(
            _preprocessed_ht(),
            _mutation_ht(),
            hl.literal(PLATEAU_MODELS),
            _coverage_model_expr(COVERAGE_MODEL),
            use_mane_select=True,
        )
        ht = aggregate_per_variant_expected_ht(
            ht, include_mu_annotations_in_grouping=True
        )
        assert list(ht.key) == [*MU_GROUPING, *GROUPINGS]
        assert ht.count() == 5


class TestCreateAggregatedExpectedHt:
    """Test the create_aggregated_expected_ht function."""

    def test_schema_matches_per_variant_path(self) -> None:
        """Test that both paths have the same key and summed fields."""
        aggregated = _path_sums("aggregated", COVERAGE_MODEL)
        per_variant = _path_sums("per_variant", COVERAGE_MODEL)
        assert aggregated["key"] == per_variant["key"] == GROUPINGS
        assert aggregated["row_type"] == per_variant["row_type"]
        assert set(aggregated["row_type"]) == set(AGGREGATE_SUM_FIELDS)

    @pytest.mark.parametrize("coverage_model", [COVERAGE_MODEL, None])
    def test_globals(self, coverage_model: Optional[Tuple[float, float]]) -> None:
        """Test the model parameters and groupings stored in the globals."""
        params = _path_sums("aggregated", coverage_model)["params"]
        assert params.coverage_model == (
            None if coverage_model is None else list(coverage_model)
        )
        assert list(params.groupings) == GROUPINGS


class TestPerVariantAndAggregatedPathsAgree:
    """Test that the per-variant and aggregated paths give the same sums.

    The two paths apply the models at different stages: per variant and then
    summed, or to summed counts. Every summed field should agree.
    """

    @pytest.mark.parametrize("coverage_model", [COVERAGE_MODEL, None])
    def test_same_sums(self, coverage_model: Optional[Tuple[float, float]]) -> None:
        """Test every summed field for every transcript and consequence."""
        per_variant = _path_sums("per_variant", coverage_model)["rows"]
        aggregated = _path_sums("aggregated", coverage_model)["rows"]
        assert set(per_variant) == set(aggregated)

        for group, pv_row in per_variant.items():
            for field in AGGREGATE_SUM_FIELDS:
                assert aggregated[group][field] == pytest.approx(pv_row[field]), (
                    group,
                    field,
                )


@functools.lru_cache(maxsize=None)
def _constraint_groups() -> Dict:
    """Return the globals, constraint group type, and rows of the grouped sums.

    The input has per-transcript, per-consequence sums. T1 has every consequence
    type. T2 has only unobserved variants. T3 has no observed and no expected
    variants. T4 has expected but no observed variants.
    """
    # (transcript, annotation, modifier, observed, expected)
    groups = [
        ("T1", "synonymous_variant", "None", 2, 3.0),
        ("T1", "missense_variant", "probably_damaging", 3, 4.0),
        ("T1", "missense_variant", "benign", 2, 2.0),
        ("T1", "stop_gained", "HC", 0, 2.5),
        ("T1", "stop_gained", "LC", 1, 0.5),
        ("T1", "splice_donor_variant", "HC", 1, 1.0),
        ("T1", "stop_gained", "None", 0, 0.7),
        ("T2", "synonymous_variant", "None", 0, 1.0),
        ("T2", "stop_gained", "HC", 0, 2.0),
        ("T3", "synonymous_variant", "None", 0, 0.0),
        ("T4", "synonymous_variant", "None", 0, 1.5),
    ]
    ht = hl.Table.parallelize(
        [
            dict(
                gene="G1",
                transcript=t,
                canonical=True,
                annotation=a,
                modifier=m,
                mu_snp=1e-8,
                mu=1e-8,
                observed_variants=[o, 0],
                possible_variants=10,
                predicted_proportion_observed=[e / 10, 0.0],
                coverage_correction=1.0,
                expected_variants=[e, 0.0],
            )
            for t, a, m, o, e in groups
        ],
        hl.tstruct(
            gene=hl.tstr,
            transcript=hl.tstr,
            canonical=hl.tbool,
            annotation=hl.tstr,
            modifier=hl.tstr,
            mu_snp=hl.tfloat64,
            mu=hl.tfloat64,
            observed_variants=hl.tarray(hl.tint64),
            possible_variants=hl.tint64,
            predicted_proportion_observed=hl.tarray(hl.tfloat64),
            coverage_correction=hl.tfloat64,
            expected_variants=hl.tarray(hl.tfloat64),
        ),
    )
    ht = aggregate_by_constraint_groups(ht)
    return {
        "meta": hl.eval(ht.constraint_group_meta),
        "cg_type": ht.constraint_groups.dtype.element_type,
        "rows": {r.transcript: r for r in ht.collect()},
    }


class TestAggregateByConstraintGroups:
    """Test the aggregate_by_constraint_groups function."""

    GROUP_NAMES = ["syn", "mis", "lof_classic", "lof_hc_lc", "lof_hc"]
    """Names of the groups in the order of the ``constraint_group_meta`` global."""

    def _group(self, transcript: str, group: str) -> hl.Struct:
        """Return one constraint group of one transcript."""
        row = _constraint_groups()["rows"][transcript]
        return row.constraint_groups[self.GROUP_NAMES.index(group)]

    def test_meta(self) -> None:
        """Test the constraint group metadata global."""
        assert _constraint_groups()["meta"] == [
            {"csq_set": "syn"},
            {"csq_set": "mis"},
            {"lof": "classic"},
            {"lof": "hc_lc"},
            {"lof": "hc"},
        ]

    @pytest.mark.parametrize(
        "group, observed, expected, n_variants",
        [
            ("syn", 2, 3.0, 1),
            ("mis", 5, 6.0, 2),
            ("lof_classic", 2, 4.7, 4),
            ("lof_hc_lc", 2, 4.0, 3),
            ("lof_hc", 1, 3.5, 2),
        ],
    )
    def test_group_sums(
        self, group: str, observed: int, expected: float, n_variants: int
    ) -> None:
        """Test the observed, expected, and possible sums of each group.

        Classic LoF is by consequence only, so it includes the stop_gained with no
        LOFTEE annotation. HC and HC+LC are by LOFTEE annotation.
        """
        cg = self._group("T1", group)
        assert cg.oe_info[0].observed_variants == observed
        assert cg.oe_info[0].expected_variants == pytest.approx(expected)
        assert cg.possible_variants == 10 * n_variants

    def test_oe_info_shape(self) -> None:
        """Test that the per-frequency arrays become one oe_info struct per entry."""
        cg_type = _constraint_groups()["cg_type"]
        assert set(cg_type) == {
            "mu_snp",
            "mu",
            "possible_variants",
            "coverage_correction",
            "oe_info",
        }
        assert set(cg_type["oe_info"].element_type) == {
            "observed_variants",
            "predicted_proportion_observed",
            "expected_variants",
        }
        assert len(self._group("T1", "syn").oe_info) == 2

    def test_unobserved_group_sums_to_zero(self) -> None:
        """Test that a group with only unobserved variants has 0 observed."""
        lof_hc = self._group("T2", "lof_hc")
        assert lof_hc.oe_info[0].observed_variants == 0
        assert lof_hc.oe_info[0].expected_variants == pytest.approx(2.0)

    def test_group_with_no_variants(self) -> None:
        """Test a group with no variants in the transcript.

        T2 has no missense variants, so its missense group sums no rows: possible
        is 0 and oe_info is missing.
        """
        mis = self._group("T2", "mis")
        assert mis.possible_variants == 0
        assert mis.oe_info is None

    def test_no_variants_and_filter(self) -> None:
        """Test the no_variants flag and the filter on observed or expected.

        T3 has no observed and no expected variants and is dropped. T4 has no
        observed variants but some expected, so it is kept with no_variants set.
        """
        rows = _constraint_groups()["rows"]
        assert set(rows) == {"T1", "T2", "T4"}
        assert rows["T1"].no_variants is False
        assert rows["T2"].no_variants is True
        assert rows["T4"].no_variants is True
