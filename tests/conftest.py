"""Pytest configuration file for Hail tests."""

import os

# Force the local Spark backend before anything imports Hail. A developer's
# `~/.config/hail/config.ini` may set `[query] backend = batch`, and test
# modules import Hail at collection time, before any fixture runs. The
# environment variable overrides the config file, so tests never submit jobs to
# Hail Batch.
os.environ.setdefault("HAIL_QUERY_BACKEND", "spark")

# Bind Spark's driver to the loopback address before Hail/Spark are imported. On
# macOS the hostname often does not resolve to a routable IP, and the
# `sparkDriver` then fails to start with a `BindException`.
_SPARK_BIND_ARGS = (
    "--conf spark.driver.bindAddress=127.0.0.1 --conf spark.driver.host=127.0.0.1"
)
_existing_args = os.environ.get("PYSPARK_SUBMIT_ARGS", "")
if "spark.driver.bindAddress" not in _existing_args:
    if _existing_args:
        os.environ["PYSPARK_SUBMIT_ARGS"] = f"{_SPARK_BIND_ARGS} {_existing_args}"
    else:
        os.environ["PYSPARK_SUBMIT_ARGS"] = f"{_SPARK_BIND_ARGS} pyspark-shell"

import hail as hl  # noqa: E402
import pytest  # noqa: E402


@pytest.fixture(scope="session", autouse=True)
def setup_hail():
    """
    Set up Hail before running tests.

    Hail is initialized once per test session, because starting the JVM costs more
    than the tests themselves, and stopped after all tests are completed.
    """
    if not hl.utils.java.Env.is_fully_initialized():
        hl.init(log="/dev/null", quiet=True)

    yield
    hl.stop()
