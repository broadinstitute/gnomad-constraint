"""Install the rr pipeline's driver-side packages on the Dataproc master."""

import subprocess
import sys

subprocess.check_call(
    [
        sys.executable,
        "-m",
        "pip",
        "install",
        "statsmodels==0.14.4",
        "scikit-learn==1.5.2",
    ]
)
import sklearn  # noqa: E402
import statsmodels  # noqa: E402

print(
    "python",
    sys.executable,
    "statsmodels",
    statsmodels.__version__,
    "sklearn",
    sklearn.__version__,
)
