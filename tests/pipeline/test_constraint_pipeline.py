"""Tests for the constraint pipeline script."""

import logging

import hail as hl
import pytest

from gnomad_constraint.pipeline.constraint_pipeline import (
    check_constraint_group_input,
    check_skip_coverage_model,
)


class TestCheckSkipCoverageModel:
    """Test the check_skip_coverage_model function."""

    @pytest.mark.parametrize("skip_coverage_model", [False, True])
    def test_matching_setting(self, skip_coverage_model: bool) -> None:
        """Test that a setting matching preprocessing passes."""
        params_expr = hl.struct(skip_coverage_model=skip_coverage_model)
        check_skip_coverage_model(params_expr, skip_coverage_model)

    @pytest.mark.parametrize("skip_coverage_model", [False, True])
    def test_mismatched_setting(self, skip_coverage_model: bool) -> None:
        """Test that a setting different from preprocessing raises an error."""
        params_expr = hl.struct(skip_coverage_model=not skip_coverage_model)
        with pytest.raises(ValueError, match="Use the same setting"):
            check_skip_coverage_model(params_expr, skip_coverage_model)


class TestCheckConstraintGroupInput:
    """Test the check_constraint_group_input function."""

    @pytest.mark.parametrize(
        "aggregate_per_variant_expected, apply_models_aggregated,"
        " use_aggregated_expected, warned",
        [
            # The route that is read wasn't run, and the other one was: the step
            # reads a Table from an earlier run.
            (False, True, False, True),
            (True, False, True, True),
            # The route that is read was run in this command, or neither was.
            (True, False, False, False),
            (False, True, True, False),
            (True, True, False, False),
            (True, True, True, False),
            (False, False, False, False),
            (False, False, True, False),
        ],
    )
    def test_warning(
        self,
        caplog: pytest.LogCaptureFixture,
        aggregate_per_variant_expected: bool,
        apply_models_aggregated: bool,
        use_aggregated_expected: bool,
        warned: bool,
    ) -> None:
        """
        Test the warning when the grouping step reads the other route's output.

        For example, ``--apply-models-aggregated --aggregate-by-constraint-groups``
        without ``--use-aggregated-expected`` groups the
        ``--aggregate-per-variant-expected`` output, maybe from an older run,
        rather than the Table just written.
        """
        with caplog.at_level(logging.WARNING, logger="constraint_pipeline"):
            check_constraint_group_input(
                aggregate_per_variant_expected=aggregate_per_variant_expected,
                apply_models_aggregated=apply_models_aggregated,
                use_aggregated_expected=use_aggregated_expected,
            )
        assert bool(caplog.records) == warned
