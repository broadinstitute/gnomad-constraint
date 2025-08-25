"""Module to register Gamma distribution UDFs with Hail."""

from typing import Optional

import hail as hl


def register_gamma_udfs():
    """
    Register Gamma distribution UDFs with Hail.
    
    This function attempts to register Python UDFs that provide accurate Gamma distribution
    calculations using scipy. However, for production use, the Scala UDF should be used.
    """
    try:
        from scipy.stats import gamma as gamma_dist
        
        # Create the UDF function
        def qgamma_udf(p, shape, scale):
            """Python UDF for Gamma quantile function."""
            return gamma_dist.ppf(p, shape, scale=scale)
        
        # Register with Hail (this is the key part!)
        hl.expr.functions.define_function(
            qgamma_udf, 
            hl.tfloat64, 
            hl.tfloat64, 
            hl.tfloat64, 
            hl.tfloat64
        )
        
        print("✅ Gamma UDFs registered successfully with Hail")
        return True
        
    except ImportError:
        print("⚠️  scipy not available")
        return False
    except Exception as e:
        print(f"❌ Failed to register UDFs: {e}")
        return False


def call_scala_gamma_udf(p, shape, scale):
    """
    Call the Scala Gamma UDF directly through the JVM gateway.
    
    This function calls the compiled Scala UDF for exact accuracy.
    """
    try:
        # Get the Spark context to access the JVM
        sc = hl.utils.java.Env.spark_session().sparkContext
        gateway = sc._gateway
        
        # Access the Scala UDF class via Py4J (this works!)
        gamma_udf = getattr(gateway.jvm.gnomad.constraint.promis3d, "GammaUDF$")
        
        # Call the static method
        result = gamma_udf.qgamma(float(p), float(shape), float(scale))
        return float(result)
        
    except Exception as e:
        raise RuntimeError(f"Failed to call Scala Gamma UDF: {e}") from e


def qgamma(
    p: hl.expr.Float64Expression,
    shape: hl.expr.Float64Expression,
    scale: hl.expr.Float64Expression,
) -> hl.expr.Float64Expression:
    """
    Calculate the quantile function of the Gamma distribution.

    This function requires the Scala UDF to be properly set up. If you're getting
    errors, you need to:
    
    1. Compile the Scala UDF: sbt assembly
    2. Upload JAR to GCS: gsutil cp target/scala-2.12/gamma-udf-1.0.jar gs://your-bucket/jars/
    3. Start DataProc with JAR: --jars=gs://your-bucket/jars/gamma-udf-1.0.jar
    4. Register UDF: gateway.jvm.gnomad.constraint.promis3d.GammaUDFRegistration.registerAll()

    :param p: Probability value (0 <= p <= 1)
    :param shape: Shape parameter (alpha) of the Gamma distribution
    :param scale: Scale parameter (beta) of the Gamma distribution
    :return: The value x such that P(X <= x) = p
    """
    # Try to use the registered Scala UDF first
    try:
        # Check if Scala UDF is registered by trying to call it
        return hl.expr.functions.call("qgamma", p, shape, scale)
    except Exception as e:
        # Try Python UDF as fallback
        try:
            return hl.expr.functions.call("qgamma_udf", p, shape, scale)
        except Exception as e2:
            # If neither UDF is available, try calling Scala UDF directly via Py4J
            try:
                # For Hail expressions, we need to use a different approach
                # This is a workaround for when the UDF isn't registered
                return hl.expr.functions.call("call_scala_gamma_udf", p, shape, scale)
            except Exception as e3:
                # If all else fails, fail with clear error message
                raise RuntimeError(
                    f"Gamma UDF not available! You need to set up the Scala UDF:\n"
                    f"1. Compile: sbt assembly\n"
                    f"2. Upload: gsutil cp target/scala-2.12/gamma-udf-1.0.jar gs://your-bucket/jars/\n"
                    f"3. Start cluster with: --jars=gs://your-bucket/jars/gamma-udf-1.0.jar\n"
                    f"4. Register: gateway.jvm.gnomad.constraint.promis3d.GammaUDFRegistration.registerAll()\n"
                    f"Original error: {e}\n"
                    f"Python UDF error: {e2}\n"
                    f"Direct call error: {e3}"
                ) from e


def gamma_ppf(
    p: hl.expr.Float64Expression,
    shape: hl.expr.Float64Expression,
    scale: hl.expr.Float64Expression,
) -> hl.expr.Float64Expression:
    """
    Calculate the percent point function (inverse CDF) of the Gamma distribution.

    :param p: Probability value (0 <= p <= 1)
    :param shape: Shape parameter (alpha) of the Gamma distribution
    :param scale: Scale parameter (beta) of the Gamma distribution
    :return: The value x such that P(X <= x) = p
    """
    return qgamma(p, shape, scale)


def check_udf_availability():
    """
    Check if the Gamma UDF is properly set up and available.
    
    Returns:
        bool: True if UDF is available, False otherwise
    """
    try:
        # Test if Scala UDF is available by calling it directly via Py4J
        sc = hl.utils.java.Env.spark_session().sparkContext
        gateway = sc._gateway
        
        # Access the Scala UDF class via Py4J (this works!)
        gamma_udf = getattr(gateway.jvm.gnomad.constraint.promis3d, "GammaUDF$")
        
        # Test the function
        result = gamma_udf.qgamma(0.5, 2.0, 1.0)
        print("✅ Scala Gamma UDF is available!")
        print(f"   Test result: qgamma(0.5, 2.0, 1.0) = {result}")
        return True
        
    except Exception as e:
        try:
            # Test if Python UDF is available
            test_expr = hl.expr.functions.call("qgamma_udf", 0.5, 2.0, 1.0)
            print("⚠️  Python Gamma UDF is available (use Scala UDF for production)")
            return True
        except Exception as e2:
            print("❌ No Gamma UDF available!")
            print("You need to set up the Scala UDF:")
            print("1. Compile: sbt assembly")
            print("2. Upload: gsutil cp target/scala-2.12/gamma-udf-1.0.jar gs://your-bucket/jars/")
            print("3. Start cluster with: --jars=gs://your-bucket/jars/gamma-udf-1.0.jar")
            print("4. Register: gateway.jvm.gnomad.constraint.promis3d.GammaUDFRegistration.registerAll()")
            return False
