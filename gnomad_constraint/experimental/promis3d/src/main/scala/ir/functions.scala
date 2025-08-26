package gnomad.constraint.promis3d.ir.functions

import is.hail.expr.ir.functions._
import is.hail.types.virtual._

import java.util.concurrent.atomic.AtomicBoolean

object ExtraMathFunctions extends RegistryFunctions {
      private val didRegister = new AtomicBoolean(false)

    def registerAll(): Unit = {
        println("ExtraMathFunctions.registerAll() called!")
        if (didRegister.compareAndSet(false, true)) {
            try {
                val statsPackageClass = Class.forName("gnomad.constraint.promis3d.stats.package$")

                registerScalaFunction(
                    "gamma_quantile",  // Changed from "qgamma" to avoid conflicts
                    Array(TFloat64, TFloat64, TFloat64),
                    TFloat64,
                    null
                )(
                    statsPackageClass,
                    "qgamma"
                )
                println("✅ Successfully registered gamma_quantile function!") 
            } catch {
                case e: ClassNotFoundException =>
                    println(s"❌ ClassNotFoundException: $e")
                case e: Exception =>
                    println(s"❌ Exception in registerAll: $e")
                    e.printStackTrace()
            }
        }
    }
    
    // Call registerAll when this object is loaded
    registerAll()
} 