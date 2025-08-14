# Gamma Distribution UDF for Hail

This directory contains Scala UDFs for Gamma distribution calculations in Hail, similar to how `hl.qchisqtail` works.

## Files

- `GammaUDF.scala`: Scala implementation of the UDFs
- `gamma_udf_registration.py`: Python module to register UDFs with Hail
- `README_GAMMA_UDF.md`: This documentation

## Usage

### Direct Usage

```python
import hail as hl
from gnomad_constraint.experimental.promis3d.gamma_udf_registration import register_gamma_udfs
from gnomad_constraint.experimental.promis3d.gamma_udf_registration import gamma_upper_ci

# Register UDFs first
register_gamma_udfs()

# Use like hl.qchisqtail
obs = 10
exp = 15.0
alpha = 0.05
upper_ci = gamma_upper_ci(obs, exp, alpha)
```

### In Hail Expressions

```python
ht = hl.utils.range_table(100)
ht = ht.annotate(
    obs=hl.rand_int(0, 20),
    exp=hl.rand_unif(5.0, 25.0)
)

# Calculate upper CI for each row
ht = ht.annotate(
    oe_upper=gamma_upper_ci(ht.obs, ht.exp, 0.05)
)
```

## Available Functions

### `gamma_upper_ci(obs, exp, alpha=0.05)`

Calculate the upper bound of the OE confidence interval using the Gamma distribution.

- `obs`: Observed count (int32)
- `exp`: Expected count (float64)
- `alpha`: Significance level (float64, default 0.05)
- Returns: Upper bound of confidence interval (float64)

### `gamma_ppf(p, shape, scale)`

Calculate the percent point function (inverse CDF) of the Gamma distribution.

- `p`: Probability value (0 ≤ p ≤ 1)
- `shape`: Shape parameter (alpha) of the Gamma distribution
- `scale`: Scale parameter (beta) of the Gamma distribution
- Returns: The value x such that P(X ≤ x) = p

## Dependencies

- `org.apache.commons.math3:commons-math3` (for GammaDistribution)
- Hail with Scala UDF support
- Python: hail, typing
