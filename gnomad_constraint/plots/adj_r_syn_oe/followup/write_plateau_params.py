"""Write the plateau fits that ``cluster/build_plateau_multipliers.py`` applies.

From ``plateau_models.fits.tsv`` (the stored lines and one line per substitution type)
and ``per_variant_expected.totals.tsv`` (the saturating and logistic curves, fit by
binomial likelihood), one ``[a, b]`` pair per model group, with the rate in units of
1e-8. Writes ``plateau_params.json``.
"""

import json
import os

import pandas as pd

OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def main() -> None:
    """Write the JSON."""
    fits = pd.read_csv(os.path.join(OUT_DIR, "plateau_models.fits.tsv"), sep="\t")
    tot = pd.read_csv(
        os.path.join(OUT_DIR, "per_variant_expected.totals.tsv"), sep="\t"
    )
    params = {"stored": {}, "saturating": {}, "logistic": {}, "by_substitution": {}}
    for name in ("stored", "by_substitution"):
        for _, r in fits[fits.fit == name].iterrows():
            params[name][r.line] = [r.intercept, r.slope * 1e-8]
    for name, fit in (
        ("saturating", "saturating, likelihood"),
        ("logistic", "logistic, likelihood"),
    ):
        for _, r in tot[(tot.fit == fit) & (tot.line != "both")].iterrows():
            params[name][r.line] = [r.a, r.b]
    with open(os.path.join(OUT_DIR, "plateau_params.json"), "w") as f:
        json.dump(params, f, indent=2)
        f.write("\n")


if __name__ == "__main__":
    main()
