"""Tests for the helpers of the constraint pipeline scripts."""

import argparse
import ast
import inspect
import logging
import os
import runpy
import sys
from collections import defaultdict
from types import SimpleNamespace
from typing import Any, Dict, List
from unittest import mock

import hail as hl
import pytest
from gnomad.resources.resource_utils import ExpressionResource, TableResource

import gnomad_constraint.pipeline as pipeline
import gnomad_constraint.resources.resource_utils as constraint_res
import gnomad_constraint.utils.constraint as constraint_utils
from gnomad_constraint.pipeline.constraint_metrics import (
    check_constraint_group_input,
    check_regional_rr_input,
    get_release_downsamplings,
)
from gnomad_constraint.pipeline.expected_counts import (
    check_apply_models_resources,
    get_coverage_model,
)
from gnomad_constraint.pipeline.expected_models import check_skip_coverage_model
from gnomad_constraint.pipeline.prepare_data import check_gerp_cutoffs_step
from gnomad_constraint.resources.constants import CURRENT_VERSION


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


class TestCheckGerpCutoffsStep:
    """Test the check_gerp_cutoffs_step function."""

    @pytest.mark.parametrize(
        "calculate_gerp_cutoffs, preprocess_data",
        [(True, True), (False, True), (False, False)],
    )
    def test_allowed(self, calculate_gerp_cutoffs: bool, preprocess_data: bool) -> None:
        """Test that the cutoffs can be calculated with preprocessing, or not at all."""
        check_gerp_cutoffs_step(calculate_gerp_cutoffs, preprocess_data)

    def test_without_preprocess_data(self) -> None:
        """
        Test that calculating the cutoffs without preprocessing raises an error.

        Regression test: the cutoffs are not saved and only ``--preprocess-data`` in
        the same command uses them, so on their own they were calculated and
        discarded, and a later ``--preprocess-data`` used the default cutoffs.
        """
        with pytest.raises(ValueError, match="--preprocess-data"):
            check_gerp_cutoffs_step(calculate_gerp_cutoffs=True, preprocess_data=False)


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
        with caplog.at_level(logging.WARNING, logger="constraint_metrics"):
            check_constraint_group_input(
                aggregate_per_variant_expected=aggregate_per_variant_expected,
                apply_models_aggregated=apply_models_aggregated,
                use_aggregated_expected=use_aggregated_expected,
            )
        assert bool(caplog.records) == warned


class TestCheckRegionalRrInput:
    """Test the check_regional_rr_input function."""

    MODELS = hl.Struct(coverage_model=[1.0, 0.5], log10_coverage=False)

    @classmethod
    def _per_variant_ht(cls, coverage_model: List[float]) -> hl.Table:
        """Return an empty per-variant Table with the given apply_models_globals."""
        ht = hl.utils.range_table(0)
        return ht.annotate_globals(
            apply_models_globals=cls.MODELS.annotate(coverage_model=coverage_model)
        )

    @classmethod
    def _rr_ht(cls, record: bool) -> hl.Table:
        """Return an empty rr Table, recording the models it was normalized on."""
        rr_globals = hl.struct(ridge_penalty=0.3)
        if record:
            rr_globals = rr_globals.annotate(
                per_variant_apply_models_globals=hl.literal(cls.MODELS)
            )
        return hl.utils.range_table(0).annotate_globals(rr_globals=rr_globals)

    @staticmethod
    def _no_rr_globals_ht() -> hl.Table:
        """Return an empty rr Table with no ``rr_globals``, like Gnocchi's."""
        return hl.utils.range_table(0)

    def test_matching_models(self) -> None:
        """Test that an rr Table normalized on the same models passes."""
        check_regional_rr_input(
            self._per_variant_ht(self.MODELS.coverage_model), self._rr_ht(True)
        )

    def test_changed_models(self) -> None:
        """
        Test that rr normalized on a per-variant Table with other models is rejected.

        This happens when ``--apply-models-per-variant --overwrite`` is rerun after
        ``--compute-rr``, which leaves the rr normalization factors stale.
        """
        with pytest.raises(ValueError, match="Rerun expected_counts.py"):
            check_regional_rr_input(self._per_variant_ht([1.0, 0.6]), self._rr_ht(True))

    def test_no_rr_globals(self, caplog: pytest.LogCaptureFixture) -> None:
        """Test that an rr Table with no ``rr_globals`` is only warned about."""
        with caplog.at_level(logging.WARNING, logger="constraint_metrics"):
            check_regional_rr_input(
                self._per_variant_ht(self.MODELS.coverage_model),
                self._no_rr_globals_ht(),
            )
        assert caplog.records

    def test_rr_globals_without_the_record(self) -> None:
        """
        Test that ``rr_globals`` without the recorded models is rejected.

        Regression test: such a Table, from rr code older than the record, passed
        with a warning, and then ``--prepare-release`` raised a KeyError when it
        dropped ``per_variant_apply_models_globals`` from the release globals.
        """
        with pytest.raises(ValueError, match="don't record"):
            check_regional_rr_input(
                self._per_variant_ht(self.MODELS.coverage_model), self._rr_ht(False)
            )


class TestGetReleaseDownsamplings:
    """Test the get_release_downsamplings function."""

    def test_downsamplings(self) -> None:
        """Test that the downsampling levels are returned."""
        ht = hl.utils.range_table(1).annotate_globals(
            downsamplings=hl.struct(nfe=[10, 100])
        )
        assert get_release_downsamplings(ht) == hl.Struct(nfe=[10, 100])

    def test_no_downsamplings(self) -> None:
        """
        Test that a release Table without downsamplings raises a clear error.

        Regression test: without ``--include-downsamplings`` at preprocessing, the
        release Table has no ``downsamplings`` global, and the downsampling export
        failed with an AttributeError.
        """
        with pytest.raises(ValueError, match="--include-downsamplings"):
            get_release_downsamplings(hl.utils.range_table(1))


class TestGetCoverageModel:
    """Test the get_coverage_model function."""

    COEFFICIENTS = [0.1, 0.02]

    # Production coverage models saved as only their coefficients, with the high
    # coverage cutoff and scale they were built with: v4.1.1 on linear AN percent,
    # and the v2.1.1 autosome/PAR model on log10 median depth.
    OLD_MODELS = [
        ([-0.14479061495071816, 0.014710539963350573], 90, False),
        ([0.2200905508803188, 0.48448910713160426], 40, True),
    ]

    @staticmethod
    def _res(tmp_path: str, coverage_model: hl.expr.Expression) -> SimpleNamespace:
        """Return an apply step's resources with the coverage model written to disk."""
        path = os.path.join(tmp_path, "coverage.he")
        hl.experimental.write_expression(coverage_model, path)
        return SimpleNamespace(model_coverage=ExpressionResource(path))

    def test_skip_coverage_model(self) -> None:
        """Test that no coverage model is read when it is skipped."""
        res = SimpleNamespace(model_coverage=None)
        assert get_coverage_model(res, True, False, 90) == (None, False)

    @pytest.mark.parametrize(
        "saved_log10_coverage, log10_coverage, warned",
        [(True, False, False), (True, True, False), (False, False, False)]
        + [(False, True, True)],
    )
    def test_saved_setting(
        self,
        tmp_path: str,
        caplog: pytest.LogCaptureFixture,
        saved_log10_coverage: bool,
        log10_coverage: bool,
        warned: bool,
    ) -> None:
        """
        Test that the coverage scale saved with the model is used.

        Regression test: only the coefficients were saved, and building and applying
        each read ``--use-logarithmic-coverage-model``, so leaving the flag off one of
        the two commands applied the model on the wrong coverage scale with no error.
        The flag can only ask for log10, so passing it for a linear model warns.
        """
        res = self._res(
            tmp_path,
            hl.struct(
                coefficients=self.COEFFICIENTS, log10_coverage=saved_log10_coverage
            ),
        )
        with caplog.at_level(logging.WARNING, logger="expected_counts"):
            coefficients, log10 = get_coverage_model(res, False, log10_coverage, 90)
        assert hl.eval(coefficients) == self.COEFFICIENTS
        assert log10 == saved_log10_coverage
        assert bool(caplog.records) == warned

    @pytest.mark.parametrize(
        "coefficients, high_cov_cutoff, log10_coverage", OLD_MODELS
    )
    def test_model_without_setting(
        self,
        tmp_path: str,
        caplog: pytest.LogCaptureFixture,
        coefficients: List[float],
        high_cov_cutoff: int,
        log10_coverage: bool,
    ) -> None:
        """Test that a model saved as only its coefficients uses the flag, and warns."""
        res = self._res(tmp_path, hl.literal(coefficients))
        with caplog.at_level(logging.WARNING, logger="expected_counts"):
            result, log10 = get_coverage_model(
                res, False, log10_coverage, high_cov_cutoff
            )
        assert hl.eval(result) == coefficients
        assert log10 == log10_coverage
        assert "saved without" in caplog.text

    @pytest.mark.parametrize(
        "coefficients, high_cov_cutoff, log10_coverage", OLD_MODELS
    )
    def test_model_without_setting_wrong_flag(
        self,
        tmp_path: str,
        coefficients: List[float],
        high_cov_cutoff: int,
        log10_coverage: bool,
    ) -> None:
        """
        Test that a flag contradicting an old model's coefficients raises an error.

        The coverage model is fit to o/e scaled so high coverage sites are at 1, so
        on the scale it was built on it gives a correction close to 1 at the high
        coverage cutoff. On the other scale it's far from 1: -0.12 for the v4.1.1
        model on log10 coverage, and 19.6 for the v2.1.1 model on linear coverage.
        """
        res = self._res(tmp_path, hl.literal(coefficients))
        with pytest.raises(ValueError, match="look like it was built on"):
            get_coverage_model(res, False, not log10_coverage, high_cov_cutoff)


SCRIPTS = ["prepare_data", "expected_models", "expected_counts", "constraint_metrics"]
"""The pipeline scripts, in the order they run."""

GET_CONSTRAINT_RESOURCES = constraint_res.get_constraint_resources
"""The real ``get_constraint_resources``, for its signature once a test replaces it."""

COMMON_FLAGS = [
    "--version",
    "4.1.1",
    "--overwrite",
    "--test",
    "--directory-post-fix",
    "dir",
    "--path-post-fix",
    "path",
]
"""Flags every script takes, set to values other than their defaults."""

COMMON_ARGS = {
    "version": "4.1.1",
    "overwrite": True,
    "test": True,
    "directory_post_fix": "dir",
    "path_post_fix": "path",
}
"""The ``get_constraint_resources`` arguments ``COMMON_FLAGS`` should give."""


class _StopRun(Exception):
    """Raised by a stand-in for a function a script calls, to end the script run."""


def _script_path(script: str) -> str:
    """
    Return the path of a pipeline script.

    :param script: Script name, without ``.py``.
    :return: Path of the script.
    """
    return os.path.join(os.path.dirname(pipeline.__file__), f"{script}.py")


def _run_script(monkeypatch: pytest.MonkeyPatch, script: str, flags: List[str]) -> None:
    """
    Run a pipeline script as ``__main__`` until a stand-in raises ``_StopRun``.

    :param monkeypatch: pytest monkeypatch fixture.
    :param script: Script name, without ``.py``.
    :param flags: Command line flags.
    :return: None.
    """
    monkeypatch.setattr(hl, "init", lambda *args, **kwargs: None)
    monkeypatch.setattr(hl, "copy_log", lambda *args, **kwargs: None)
    monkeypatch.setattr(sys, "argv", [f"{script}.py"] + flags)
    with pytest.raises(_StopRun):
        runpy.run_path(_script_path(script), run_name="__main__")


class TestScriptResources:
    """Test the resources each script asks for and the options it reads."""

    @pytest.mark.parametrize(
        "script, flags, expected",
        [
            *[(script, [], {}) for script in SCRIPTS],
            *[(script, COMMON_FLAGS, COMMON_ARGS) for script in SCRIPTS],
            *[(script, ["--test-gene-list"], {"test": True}) for script in SCRIPTS],
            # No prepare_data.py step reads the models, so it doesn't pass them.
            ("prepare_data", ["--skip-coverage-model"], {}),
            ("expected_models", ["--skip-coverage-model"], {"models": ["plateau"]}),
            ("expected_counts", ["--skip-coverage-model"], {"models": ["plateau"]}),
            (
                "constraint_metrics",
                ["--skip-pre-rank-metrics"],
                {"skip_pre_rank_metrics": True},
            ),
            (
                "constraint_metrics",
                ["--use-aggregated-expected"],
                {"use_aggregated_expected": True},
            ),
            (
                "constraint_metrics",
                ["--apply-regional-rr"],
                {"apply_regional_rr": True},
            ),
        ],
    )
    def test_get_constraint_resources_call(
        self,
        monkeypatch: pytest.MonkeyPatch,
        script: str,
        flags: List[str],
        expected: Dict[str, Any],
    ) -> None:
        """
        Test the whole ``get_constraint_resources`` call each script makes.

        The resources decide which Tables each step reads and checks, so a flag
        that doesn't reach the call makes a step check or read the wrong Tables.
        Regression test: ``constraint_metrics.py`` had no ``--skip-coverage-model``,
        and an earlier version of this test only checked ``models`` in three
        scripts, so it passed with ``use_aggregated_expected`` dropped from the
        ``constraint_metrics.py`` call.
        """
        calls = []

        def get_resources(*args: Any, **kwargs: Any) -> None:
            calls.append(
                inspect.signature(GET_CONSTRAINT_RESOURCES).bind(*args, **kwargs)
            )
            raise _StopRun

        monkeypatch.setattr(constraint_res, "get_constraint_resources", get_resources)
        _run_script(monkeypatch, script, flags)

        call = calls[0]
        call.apply_defaults()
        call.arguments["models"] = list(call.arguments["models"])
        assert call.arguments == {
            "version": CURRENT_VERSION,
            "overwrite": False,
            "test": False,
            "models": ["plateau", "coverage"],
            "directory_post_fix": None,
            "path_post_fix": None,
            "skip_pre_rank_metrics": False,
            "use_aggregated_expected": False,
            "apply_regional_rr": False,
            **expected,
        }

    @pytest.mark.parametrize("skip_coverage_model", [False, True])
    def test_preprocess_data_skip_coverage_model(
        self, monkeypatch: pytest.MonkeyPatch, skip_coverage_model: bool
    ) -> None:
        """
        Test that ``prepare_data.py --preprocess-data`` passes ``--skip-coverage-model`` on.

        ``prepare_data.py`` doesn't pass ``models`` to ``get_constraint_resources``,
        so the flag reaches the rest of the pipeline only through the preprocessed
        Table, which the later scripts check their own flag against.
        """
        calls = []

        def prepare(*args: Any, **kwargs: Any) -> None:
            calls.append(kwargs)
            raise _StopRun

        monkeypatch.setattr(
            constraint_res,
            "get_constraint_resources",
            lambda *args, **kwargs: mock.MagicMock(),
        )
        monkeypatch.setattr(
            constraint_utils, "prepare_ht_for_constraint_calculations", prepare
        )
        flags = ["--preprocess-data"]
        if skip_coverage_model:
            flags.append("--skip-coverage-model")
        _run_script(monkeypatch, "prepare_data", flags)

        assert calls[0]["skip_coverage_model"] == skip_coverage_model

    @pytest.mark.parametrize("script", SCRIPTS)
    def test_reads_only_its_own_options(
        self, monkeypatch: pytest.MonkeyPatch, script: str
    ) -> None:
        """
        Test that every ``args.<name>`` a script reads is one of its own options.

        After the split, a step block left reading an option that went to another
        script would only fail, with an AttributeError, when that step ran.
        """
        parsers = []

        def parse_args(parser: argparse.ArgumentParser, *args: Any) -> None:
            parsers.append(parser)
            raise _StopRun

        monkeypatch.setattr(argparse.ArgumentParser, "parse_args", parse_args)
        _run_script(monkeypatch, script, [])
        options = set(vars(parsers[0].parse_known_args([])[0]))

        with open(_script_path(script)) as f:
            tree = ast.parse(f.read())
        read = {
            node.attr
            for node in ast.walk(tree)
            if isinstance(node, ast.Attribute)
            and isinstance(node.value, ast.Name)
            and node.value.id == "args"
        }
        assert len(read) > 5
        assert sorted(read - options) == []


def _help_text(node: ast.expr) -> str:
    """
    Return the text of a help argument, with each f-string field replaced by "x".

    :param node: AST node of the ``help`` argument.
    :return: Help text.
    """
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.JoinedStr):
        return "".join(
            v.value if isinstance(v, ast.Constant) else "x" for v in node.values
        )
    raise ValueError(f"Unexpected help argument: {ast.dump(node)}")


class TestArgumentHelp:
    """Test the help text of the pipeline's arguments."""

    def test_help_formats(self) -> None:
        """
        Test that every flag's help, in every script, formats the way argparse does.

        argparse fills in each help string with ``%``, so a lone ``%`` must be
        written ``%%``. Regression test: "AN% of 20" in the
        ``--pipeline-low-coverage-filter`` help read as a ``%o`` conversion, so
        ``--help`` raised a TypeError.
        """
        helps = []
        for script in SCRIPTS:
            path = os.path.join(os.path.dirname(pipeline.__file__), f"{script}.py")
            with open(path) as f:
                tree = ast.parse(f.read())
            helps += [
                _help_text(kw.value)
                for node in ast.walk(tree)
                if isinstance(node, ast.Call)
                and getattr(node.func, "attr", None) == "add_argument"
                for kw in node.keywords
                if kw.arg == "help"
            ]
        assert len(helps) > 40
        for text in helps:
            text % defaultdict(str)
