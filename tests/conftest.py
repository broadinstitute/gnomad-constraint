"""Pytest configuration file for Hail tests."""

import os
import tempfile

# Force the local Spark backend before anything imports Hail. A developer's
# `~/.config/hail/config.ini` may set `[query] backend = batch`, and the
# environment variable overrides it, so tests never submit jobs to Hail Batch. It
# is assigned rather than defaulted so that an inherited value can't override it.
os.environ["HAIL_QUERY_BACKEND"] = "spark"

# Spark reads PYSPARK_SUBMIT_ARGS when Hail starts its JVM, so these must be set
# before Hail is imported:
#   - Bind the driver to the loopback address. On macOS the hostname often does not
#     resolve to a routable IP, and the driver then fails with a `BindException`.
#   - Give the driver a 4 GB heap. Spark's default of 1 GB runs out of memory
#     partway through the suite.
# Each argument is added only if it isn't already set, and each is checked on its
# own, so a developer's own setting for one doesn't drop the others.
_SPARK_SUBMIT_ARGS = {
    "spark.driver.bindAddress": "--conf spark.driver.bindAddress=127.0.0.1",
    "spark.driver.host": "--conf spark.driver.host=127.0.0.1",
    "--driver-memory": "--driver-memory 4g",
}
_existing_args = os.environ.get("PYSPARK_SUBMIT_ARGS", "pyspark-shell")
os.environ["PYSPARK_SUBMIT_ARGS"] = " ".join(
    [arg for key, arg in _SPARK_SUBMIT_ARGS.items() if key not in _existing_args]
    + [_existing_args]
)

import hail as hl  # noqa: E402
import pytest  # noqa: E402

_HAIL_LOG = os.path.join(
    tempfile.mkdtemp(prefix="gnomad_constraint_tests_"), "hail.log"
)
_STARTED_HAIL = False


def pytest_configure(config: pytest.Config) -> None:
    """
    Start Hail once for the test session, before tests are collected.

    Some test modules build Hail types such as ``hl.tlocus("GRCh38")`` when they are
    imported, which starts Hail with default settings if it isn't running yet. Hail
    is therefore started here rather than in a fixture, which would run after
    collection.

    :param config: Pytest config object.
    :return: None.
    """
    global _STARTED_HAIL
    if not hl.utils.java.Env.is_fully_initialized():
        hl.init(log=_HAIL_LOG, quiet=True)
        _STARTED_HAIL = True


def pytest_report_header(config: pytest.Config) -> str:
    """
    Show where the Hail log for this session is written.

    :param config: Pytest config object.
    :return: Header line with the Hail log path.
    """
    return f"hail log: {_HAIL_LOG}"


def pytest_unconfigure(config: pytest.Config) -> None:
    """
    Stop Hail at the end of the session if :func:`pytest_configure` started it.

    :param config: Pytest config object.
    :return: None.
    """
    if _STARTED_HAIL:
        hl.stop()
