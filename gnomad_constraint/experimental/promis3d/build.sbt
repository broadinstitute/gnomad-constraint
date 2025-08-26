name := "gamma-udf"

version := "1.0"

scalaVersion := "2.12.20"

// Add local Hail JAR
unmanagedJars in Compile += file("/Users/jgoodric/miniconda3/envs/gnomad_constraint/lib/python3.10/site-packages/hail/backend/hail-all-spark.jar")

// Dependencies
libraryDependencies ++= Seq(
  "org.apache.commons" % "commons-math3" % "3.6.1",
  "net.sourceforge.jdistlib" % "jdistlib" % "0.4.5"
)

// Assembly plugin for creating fat JAR
assembly / assemblyMergeStrategy := {
  case PathList("META-INF", xs @ _*) => MergeStrategy.discard
  case x => MergeStrategy.first
}

// JAR name
assembly / assemblyJarName := "gamma-udf-1.0.jar" 