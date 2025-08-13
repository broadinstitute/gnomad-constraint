# Gamma Distribution UDFs for Hail

This directory contains Scala UDFs for Gamma distribution calculations in Hail, similar to how `hl.qchisqtail` works.

## Overview

The current implementation in `utils.py` uses Pandas UDFs for Gamma distribution calculations:

```python
def make_gamma_upper_ci(alpha: float = 0.05) -> hl.expr.ArrayExpression:
    @pandas_udf("double")
    def _udf(obs: pd.Series, exp: pd.Series) -> pd.Series:
        import numpy as np
        from scipy.stats import gamma as gamma_dist

        a = obs.to_numpy(dtype=np.float64) + 1.0
        b = exp.to_numpy(dtype=np.float64)

        return pd.Series(
            gamma_dist.ppf(1.0 - alpha, a=a, scale=1.0 / b), dtype="float64"
        )

    return _udf
```

This new implementation provides native Hail UDFs that:
1. Use Apache Commons Math's `GammaDistribution`
2. Integrate directly with Hail's IR system
3. Provide better performance by avoiding Python-JVM serialization
4. Work similarly to `hl.qchisqtail`

## Files

- `GammaUDF.scala`: Scala implementation of the UDFs
- `gamma_udf_registration.py`: Python module to register UDFs with Hail
- `gamma_udf_example.py`: Example usage and migration guide
- `README_GAMMA_UDF.md`: This documentation

## Setup

### 1. Compile the Scala UDFs

The Scala UDFs need to be compiled and included in your Hail JAR. You'll need to:

1. Add the Scala files to your Hail build
2. Ensure `org.apache.commons.math3` is in your dependencies
3. Compile and package with Hail

### 2. Register UDFs in Python

```python
from gnomad_constraint.experimental.promis3d.gamma_udf_registration import register_gamma_udfs

# Register the UDFs with Hail
register_gamma_udfs()
```

## Usage

### Direct Usage

```python
import hail as hl
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

### Replace Existing Implementation

Replace the current `make_gamma_upper_ci` function in `utils.py`:

```python
# Old implementation (Pandas UDF)
def make_gamma_upper_ci(alpha: float = 0.05) -> hl.expr.ArrayExpression:
    @pandas_udf("double")
    def _udf(obs: pd.Series, exp: pd.Series) -> pd.Series:
        # ... Pandas UDF implementation
    return _udf

# New implementation (Native Hail UDF)
def make_gamma_upper_ci(alpha: float = 0.05):
    def _udf(obs, exp):
        return gamma_upper_ci(obs, exp, alpha)
    return _udf
```

The function signature and usage pattern remain the same, so existing code doesn't need to change.

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

## Migration Guide

1. **Register UDFs**: Call `register_gamma_udfs()` at the start of your script
2. **Replace function**: Update `make_gamma_upper_ci` in `utils.py`
3. **Test**: Verify results match the Pandas UDF implementation
4. **Remove dependencies**: You can remove `scipy.stats` dependency if not used elsewhere

## Performance Benefits

Native Hail UDFs provide several advantages:

1. **No serialization overhead**: Avoid Python-JVM data transfer
2. **Better memory usage**: No intermediate Pandas DataFrames
3. **Faster execution**: Direct JVM computation
4. **Better integration**: Native Hail expression optimization

## Dependencies

- `org.apache.commons.math3:commons-math3` (for GammaDistribution)
- Hail with Scala UDF support
- Python: hail, typing

## Testing

Run the example script to test the implementation:

```bash
python gnomad_constraint/experimental/promis3d/gamma_udf_example.py
```

This will create test tables and demonstrate the usage patterns.
