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

def qgamma(p, shape, scale):
    return _apply("qgamma", tfloat64,
                  hl.float64(p), hl.float64(shape), hl.float64(scale))

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
        .set("spark.jars", f"{hail_jar},{user_jar}")
        # Make ./hail-all-spark.jar visible on executors.
        .set("spark.yarn.dist.jars", f"file://{hail_jar},file://{user_jar}")
        # Driver sees Hail classes
        .set("spark.driver.extraClassPath", f"{hail_jar}{cp_sep}{user_jar}")
        .set("spark.executor.extraClassPath", f"./hail-all-spark.jar{cp_sep}./gamma-udf-1.0.jar")

        .set("spark.driver.userClassPathFirst", "true")
        .set("spark.executor.userClassPathFirst", "true")

        # YARN client mode; aligns with working conf.
        .set("spark.submit.deployMode", "client")
        .set("spark.yarn.unmanagedAM.enabled", "true")
        .set("spark.dynamicAllocation.enabled", "true")
    )

    # Start Spark with the new conf that includes Hail and the gamma jar.
    spark = SparkSession.builder.config(conf=conf).getOrCreate()
    hl.init(sc=spark.sparkContext)
    
    sc  = hl.spark_context()
    jvm = sc._jvm

    # Get Spark's user-class loader.
    loader = jvm.org.apache.spark.util.Utils.getContextOrSparkClassLoader()

    # Add the jar URL to THAT loader.
    loader.addURL(jvm.java.io.File("/tmp/gamma-udf-1.0.jar").toURI().toURL())
    
    # Force-load the Scala object; this runs its <clinit> which calls registerAll() once.
    jvm.java.lang.Class.forName("gnomad.constraint.promis3d.ir.functions.ExtraMathFunctions$", True, loader)

    register_function(
        "qgamma",
        (
            dtype("float64"),
            dtype("float64"),
            dtype("float64"),
        ),
        dtype("float64"),
    )

    hl._set_flags(
        use_ssa_logs="1",
        no_whole_stage_codegen="1",
        no_ir_logging=None,
        print_ir_on_worker="1",
    )
    result = hl.eval(qgamma(0.5, 2.0, 1.0))
    print(f"qgamma(0.5, 2.0, 1.0) = {result}")
        
if __name__ == "__main__":
    main() 