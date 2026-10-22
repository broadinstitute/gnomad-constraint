"""Tests for the constraint pipeline script."""

import hail as hl
import pytest

from gnomad_constraint.pipeline.constraint_pipeline import check_skip_coverage_model


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
