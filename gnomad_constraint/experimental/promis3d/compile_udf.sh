#!/bin/bash

# Compile Scala UDF to JAR
echo "Compiling Gamma UDF to JAR..."

# Install sbt if not available
if ! command -v sbt &> /dev/null; then
    echo "sbt not found. Please install sbt first:"
    echo "  macOS: brew install sbt"
    echo "  Ubuntu: sudo apt-get install sbt"
    echo "  Or download from: https://www.scala-sbt.org/download.html"
    exit 1
fi

# Compile and create JAR
echo "Running sbt assembly..."
sbt assembly

if [ $? -eq 0 ]; then
    echo "✅ JAR compiled successfully!"
    echo "JAR location: target/scala-2.12/gamma-udf-1.0.jar"
    
    # Upload to GCS bucket for DataProc
    echo ""
    echo "To use with DataProc, upload the JAR to GCS:"
    echo "gsutil cp target/scala-2.12/gamma-udf-1.0.jar gs://your-bucket/jars/"
else
    echo "❌ Compilation failed"
    exit 1
fi 