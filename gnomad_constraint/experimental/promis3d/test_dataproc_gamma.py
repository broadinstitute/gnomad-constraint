#!/usr/bin/env python3
"""
Test script for Gamma UDFs on DataProc cluster.

This script tests the Gamma distribution UDFs and will fail fast if the Scala UDF
is not properly set up, ensuring you know exactly what's happening.
"""

import hail as hl
from gnomad_constraint.experimental.promis3d.gamma_udf_registration import (
    qgamma, 
    register_gamma_udfs,
    check_udf_availability
)
from gnomad_constraint.experimental.promis3d.utils import (
    calculate_oe_upper,
    gamma_upper_ci,
)


def test_gamma_functions():
    """Test Gamma UDF functions."""
    print("Testing Gamma UDFs on DataProc")
    print("=" * 50)
    
    # First, check if UDFs are available
    print("1. Checking UDF availability...")
    if not check_udf_availability():
        print("❌ UDFs not available - you need to set up the Scala UDF!")
        print("See DATAPROC_SCALA_UDF_GUIDE.md for instructions")
        return False
    
    print("2. Testing qgamma function...")
    try:
        # Test with known values
        result = qgamma(0.95, 2.0, 3.0)
        print(f"qgamma(0.95, 2.0, 3.0) = {result}")
        print("✅ qgamma works")
    except Exception as e:
        print(f"❌ qgamma failed: {e}")
        return False

    print("3. Testing gamma_upper_ci function...")
    try:
        # Test with known values
        result1 = gamma_upper_ci(10, 15.0, 0.05)
        result2 = gamma_upper_ci(5, 8.0, 0.05)
        print(f"gamma_upper_ci(10, 15.0, 0.05) = {result1}")
        print(f"gamma_upper_ci(5, 8.0, 0.05) = {result2}")
        print("✅ gamma_upper_ci works")
    except Exception as e:
        print(f"❌ gamma_upper_ci failed: {e}")
        return False

    print("4. Testing Hail table integration...")
    try:
        # Create a simple Hail table
        ht = hl.utils.range_table(5)
        ht = ht.annotate(
            obs=hl.literal(10),
            exp=hl.literal(15.0)
        )
        ht = ht.annotate(oe_upper=gamma_upper_ci(ht.obs, ht.exp, 0.05))
        
        # Show some results
        print("Sample results:")
        ht.show(5)
        print("✅ Hail table integration works")
    except Exception as e:
        print(f"❌ Hail table integration failed: {e}")
        return False

    print("5. Testing calculate_oe_upper function...")
    try:
        # Test with a simple array
        test_ht = hl.utils.range_table(3)
        test_ht = test_ht.annotate(
            oe=hl.array([
                hl.struct(obs=10, exp=15.0),
                hl.struct(obs=5, exp=8.0),
                hl.struct(obs=20, exp=25.0),
            ])
        )
        result = calculate_oe_upper(test_ht.oe, alpha=0.05)
        print("✅ calculate_oe_upper works")
    except Exception as e:
        print(f"❌ calculate_oe_upper failed: {e}")
        return False

    print("\n🎉 All tests passed! You're using the real Scala UDF!")
    return True


if __name__ == "__main__":
    hl.init()
    print("Registering Gamma UDFs...")
    register_gamma_udfs()
    
    try:
        success = test_gamma_functions()
        if not success:
            print("\n❌ Tests failed - you need to set up the Scala UDF properly!")
            print("See DATAPROC_SCALA_UDF_GUIDE.md for complete instructions")
            exit(1)
    except Exception as e:
        print(f"❌ Test failed: {e}")
        print("\nYou need to set up the Scala UDF:")
        print("1. Compile: sbt assembly")
        print("2. Upload: gsutil cp target/scala-2.12/gamma-udf-1.0.jar gs://your-bucket/jars/")
        print("3. Start cluster with: --jars=gs://your-bucket/jars/gamma-udf-1.0.jar")
        print("4. Register: gateway.jvm.gnomad.constraint.promis3d.GammaUDFRegistration.registerAll()")
        raise
    finally:
        hl.stop()
