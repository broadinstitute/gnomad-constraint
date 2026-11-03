"""Tests for the constraint pipeline script."""

import logging
import os
from types import SimpleNamespace
from unittest import mock

import hail as hl
import pytest
from gnomad.resources.resource_utils import TableResource

from gnomad_constraint.pipeline.constraint_pipeline import (
    check_apply_models_resources,
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


class TestCheckApplyModelsResources:
    """Test the check_apply_models_resources function."""

    @pytest.fixture(scope="class")
    @classmethod
    def preprocessed_path(cls, tmp_path_factory: pytest.TempPathFactory) -> str:
        """Return the path of a preprocessed Table made with --skip-coverage-model."""
        path = os.path.join(tmp_path_factory.mktemp("preprocessed"), "pre.ht")
        ht = hl.utils.range_table(1).annotate_globals(
            apply_models_globals=hl.struct(skip_coverage_model=True)
        )
        ht.write(path)
        return path

    @staticmethod
    def _res(path: str) -> SimpleNamespace:
        """Return an apply step's resources with a mock resource check."""
        return SimpleNamespace(
            temp_preprocess_data_ht=TableResource(path),
            check_resource_existence=mock.Mock(),
        )

    def test_mismatch_reported_before_resource_check(
        self, preprocessed_path: str
    ) -> None:
        """
        Test that a --skip-coverage-model mismatch is reported before the resource check.

        Regression test: the resource check ran first, so a mismatch was reported
        as the coverage model missing ("Please add --build-models to the command
        line") instead of as the mismatch.
        """
        res = self._res(preprocessed_path)
        with pytest.raises(ValueError, match="Use the same setting"):
            check_apply_models_resources(res, skip_coverage_model=False)
        res.check_resource_existence.assert_not_called()

    def test_matching_setting_checks_resources(self, preprocessed_path: str) -> None:
        """Test that a matching setting goes on to the resource check."""
        res = self._res(preprocessed_path)
        check_apply_models_resources(res, skip_coverage_model=True)
        res.check_resource_existence.assert_called_once()

    def test_missing_preprocessed_table_checks_resources(self, tmp_path: str) -> None:
        """Test that a missing preprocessed Table is left to the resource check."""
        res = self._res(os.path.join(tmp_path, "missing.ht"))
        check_apply_models_resources(res, skip_coverage_model=False)
        res.check_resource_existence.assert_called_once()


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
