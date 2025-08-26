#!/usr/bin/env python3
"""
Final success test using the working Py4J approach.
"""

import os, subprocess

import hail as hl
from pathlib import Path
from hail.ir.ir import register_function
from pyspark import SparkConf, SparkContext
from pyspark.sql import SparkSession

from hail.ir import Apply
from hail.expr.expressions import construct_expr, unify_all
from hail.expr.types import tfloat64, dtype

def _apply(name, ret_type, *args):
    idx, aggs = unify_all(*args)
    ir = Apply(name, ret_type, *(a._ir for a in args))  # Apply(function, return_type, *args)
    return construct_expr(ir, ret_type, idx, aggs)

def gamma_quantile(p, shape, scale):
    return _apply("gamma_quantile", tfloat64,
                  hl.float64(p), hl.float64(shape), hl.float64(scale))


def test_executor_jar_access_final():
    """
    Test if executors can access the JAR after multiple distribution attempts.
    """
    try:
        print("🔍 Testing executor JAR access...")
        sc = hl.spark_context()
        app_id = sc.getConf().get("spark.app.id")
        app_id = app_id.split("_")[-1]

        jvm = sc._jvm
        
        def check_jar_final(x):
            """Function to run on executors to check JAR access."""
            try:
                import os
                import glob
                
                # Check for our JAR in various locations
                jar_found = False
                jar_location = None
                
                # Check common Spark JAR locations
                spark_dirs = glob.glob("/tmp/spark-*")
                for spark_dir in spark_dirs:
                    # Check subdirectories
                    for subdir in ["jars", "files", "archives"]:
                        jar_path = os.path.join(spark_dir, subdir, "gamma-udf-1.0.jar")
                        if os.path.exists(jar_path):
                            jar_found = True
                            jar_location = jar_path
                            break
                    if jar_found:
                        break
                
                jar_path = "/tmp/gamma-udf-1.0.jar"
                if os.path.exists(jar_path):
                    jar_found = True
                    jar_location = jar_path

                jar_path = "gamma-udf-1.0.jar"
                if os.path.exists(jar_path):
                    jar_found = True
                    jar_location = jar_path
                
                # Also check if we can access the JAR through the classpath
                try:
                    print("🔍 Checking class access...")
                    
                    # Try to load our class
                    try:
                        import gnomad.constraint.promis3d.ir.functions.ExtraMathFunctions
                        class_accessible = "✅ Class loaded successfully"
                    except Exception as e:
                        class_accessible = f"❌ Error loading class 1: {e}"
                except Exception as e:
                    class_accessible = f"❌ Error loading class 2: {e}"
                
                return {
                    'partition_id': x,
                    'jar_found': jar_found,
                    'jar_location': jar_location,
                    'class_accessible': class_accessible
                }
                
            except Exception as e:
                return {
                    'partition_id': x,
                    'jar_found': False,
                    'jar_location': None,
                    'class_accessible': "❌ Error loading class 3",
                    'error': str(e)
                }
        
        print("  Running final JAR access check on executors...")
        
        # Create a simple RDD to test executor access
        rdd = sc.parallelize([1, 2])
        results = rdd.map(check_jar_final).collect()
        
        print("  ✅ Final executor JAR access results:")
        for result in results:
            print(f"    Partition {result['partition_id']}: JAR found={result['jar_found']}: class_accessible={result['class_accessible']}")
            if result['jar_found']:
                print(f"      Location: {result['jar_location']}")
            if 'error' in result:
                print(f"      Error: {result['error']}")
        
        # Check if any executor can access the class
        any_class_accessible = any(result['class_accessible'] for result in results)
        if any_class_accessible:
            print("  ✅ At least one executor can access the UDF class!")
        else:
            print("  ❌ No executor can access the UDF class")
        
        return any_class_accessible
        
    except Exception as e:
        print(f"❌ Final executor JAR access test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    hail_jar = str(Path(hl.__file__).with_name("backend") / "hail-all-spark.jar")
    
    # Copy gamma jar locally to remove GCS, but fall back to gs:// if copy not possible.
    user_jar_gs = "gs://gnomad-julia/jars/gamma-udf-1.0.jar"
    local_user_jar = "/tmp/gamma-udf-1.0.jar"
    subprocess.run(["gsutil", "-q", "cp", user_jar_gs, local_user_jar], check=True)
    user_jar = local_user_jar

    # ":" on Linux.
    cp_sep = os.pathsep  

    # Build a FRESH SparkConf: Kryo + Hail jar + gamma jar.
    conf = (
        SparkConf()
        .setMaster("yarn")
        .setAppName("hail+gamma")

        # Hail requirements.
        .set("spark.serializer", "org.apache.spark.serializer.KryoSerializer")
        .set("spark.kryo.registrator", "is.hail.kryo.HailKryoRegistrator")
        .set("spark.kryoserializer.buffer.max", "1g")

        # JAR shipping + classpath.
        .set("spark.jars", f"{hail_jar}:{user_jar_gs}")
        # Make ./hail-all-spark.jar visible on executors.
        .set("spark.yarn.dist.jars", f"file://{hail_jar},file://{local_user_jar}")
        .set("spark.repl.local.jars", f"file://{hail_jar},file://{local_user_jar}")
        # Driver sees Hail classes
        .set("spark.driver.extraClassPath", f"{hail_jar}{cp_sep}{user_jar_gs}")
        .set("spark.executor.extraClassPath", f"./hail-all-spark.jar{cp_sep}./gamma-udf-1.0.jar")

        .set("spark.driver.userClassPathFirst", "false")
        .set("spark.executor.userClassPathFirst", "false")

        # YARN client mode; aligns with working conf.
        .set("spark.submit.deployMode", "client")
        .set("spark.yarn.unmanagedAM.enabled", "true")
        .set("spark.dynamicAllocation.enabled", "true")
        .set("spark.yarn.secondary.jars", "hail-all-spark.jar,gamma-udf-1.0.jar")
    )

    # Start Spark with the new conf that includes Hail and the gamma jar.
    spark = SparkSession.builder.config(conf=conf).getOrCreate()
    hl.init(sc=spark.sparkContext)
    #hl.init()
    
    sc  = hl.spark_context()
    print("scala:", hl.spark_context()._jvm.scala.util.Properties.versionNumberString())

    print(sc.getConf().toDebugString())
    jvm = sc._jvm

    # Get Spark's user-class loader.
    loader = jvm.org.apache.spark.util.Utils.getContextOrSparkClassLoader()

    # Add the jar URL to THAT loader.
    loader.addURL(jvm.java.io.File("/tmp/gamma-udf-1.0.jar").toURI().toURL())
    
    # Force-load the Scala object; this runs its <clinit> which calls registerAll() once.
    jvm.java.lang.Class.forName("gnomad.constraint.promis3d.ir.functions.ExtraMathFunctions$", True, loader)

    register_function(
        "gamma_quantile",
        (
            dtype("float64"),
            dtype("float64"),
            dtype("float64"),
        ),
        dtype("float64"),
    )

    test_executor_jar_access_final()

    #result = hl.eval(gamma_quantile(0.5, 2.0, 1.0))
    #print(f"gamma_quantile(0.5, 2.0, 1.0) = {result}")
        
if __name__ == "__main__":
    main() 