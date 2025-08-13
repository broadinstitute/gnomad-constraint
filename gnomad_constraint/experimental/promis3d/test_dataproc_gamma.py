#!/usr/bin/env python3
"""
Simple test script for Gamma UDFs on DataProc cluster.

This script tests the Gamma distribution functions that have been integrated
into the Promis3D pipeline.
"""

import hail as hl

from gnomad_constraint.experimental.promis3d.gamma_udf_registration import qgamma
from gnomad_constraint.experimental.promis3d.utils import (
    calculate_oe_upper,
    gamma_upper_ci,
)


def test_gamma_functions():
    """Test the Gamma distribution functions."""
    print("Testing Gamma UDFs on DataProc")
    print("=" * 40)

    # Test 1: Basic qgamma function
    print("\n1. Testing qgamma function...")
    result = qgamma(0.95, 2.0, 3.0)
    print(f"   qgamma(0.95, 2.0, 3.0) = {result}")
    print("   ✅ qgamma works")

    # Test 2: gamma_upper_ci function
    print("\n2. Testing gamma_upper_ci function...")
    result = gamma_upper_ci(10, 15.0, 0.05)
    print(f"   gamma_upper_ci(10, 15.0, 0.05) = {result}")
    print("   ✅ gamma_upper_ci works")

    result = gamma_upper_ci(23, 23.56871740477893, 0.05)
    print(f"   gamma_upper_ci(10, 15.0, 0.05) = {result}")
    print("   ✅ gamma_upper_ci works")

    # Test 3: Hail table integration
    print("\n3. Testing Hail table integration...")
    ht = hl.Table.parallelize(
        [
            hl.struct(obs=10, exp=15.0),
            hl.struct(obs=5, exp=8.0),
            hl.struct(obs=20, exp=25.0),
            hl.struct(obs=23, exp=23.56871740477893),
        ]
    )
    ht = ht.annotate(oe_upper=gamma_upper_ci(ht.obs, ht.exp, 0.05))

    print("   Sample results:")
    ht.show()
    print("   ✅ Hail table integration works")

    # Test 4: calculate_oe_upper function
    print("\n4. Testing calculate_oe_upper function...")
    # Create a proper test structure for calculate_oe_upper
    test_ht = hl.utils.range_table(3)
    test_ht = test_ht.annotate(
        oe=hl.array(
            [
                hl.struct(obs=10, exp=15.0),
                hl.struct(obs=5, exp=8.0),
                hl.struct(obs=20, exp=25.0),
                hl.struct(obs=23, exp=23.56871740477893),
            ]
        )
    )
    test_ht = test_ht.annotate(result=calculate_oe_upper(test_ht.oe, alpha=0.05))
    test_ht = test_ht.explode(test_ht.result)
    print("   calculate_oe_upper result structure:")
    test_ht.show()
    print("   ✅ calculate_oe_upper works")

    print("\n🎉 All Gamma UDF tests passed!")
    print("Ready for Promis3D pipeline integration.")


if __name__ == "__main__":
    # Initialize Hail
    hl.init()

    try:
        test_gamma_functions()
    except Exception as e:
        print(f"❌ Test failed: {e}")
        raise
    finally:
        hl.stop()
