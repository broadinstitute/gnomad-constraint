package gnomad.constraint.promis3d

import is.hail.expr.ir.functions.IRFunction
import is.hail.expr.ir.{EmitCodeBuilder, IEmitCode, EmitCode}
import is.hail.types.physical.{PFloat64, PInt32}
import is.hail.types.virtual.{TFloat64, TInt32}
import org.apache.commons.math3.distribution.GammaDistribution

/**
 * Gamma distribution UDFs for Hail
 *
 * These UDFs provide Gamma distribution functionality similar to hl.qchisqtail
 * but using the Apache Commons Math GammaDistribution implementation.
 */

object GammaPPF extends IRFunction {
  val name = "gamma_ppf"
  val argTypes = Seq(TFloat64, TFloat64, TFloat64) // p, shape, scale
  val returnType = TFloat64

  def apply(cb: EmitCodeBuilder, args: IndexedSeq[IEmitCode]): IEmitCode = {
    val p = args(0)
    val shape = args(1)
    val scale = args(2)

    p.flatMap(cb) { pCode =>
      shape.flatMap(cb) { shapeCode =>
        scale.flatMap(cb) { scaleCode =>
          cb.emb.newEmitCode(
            returnType,
            cb => {
              val pValue = cb.memoize(pCode.asFloat64.value)
              val shapeValue = cb.memoize(shapeCode.asFloat64.value)
              val scaleValue = cb.memoize(scaleCode.asFloat64.value)

              // Create Gamma distribution with shape and scale parameters
              val gammaDist = cb.memoize(
                cb.emb.newInstance[GammaDistribution](shapeValue, scaleValue)
              )

              // Calculate PPF (inverse CDF)
              val result = cb.memoize(
                gammaDist.inverseCumulativeProbability(pValue)
              )

              cb.assign(cb.emb.newRVariable("result"), result)
              cb.emb.getCode(result)
            }
          )
        }
      }
    }
  }
}

object QGamma extends IRFunction {
  val name = "qgamma"
  val argTypes = Seq(TFloat64, TFloat64, TFloat64) // p, shape, scale
  val returnType = TFloat64

  def apply(cb: EmitCodeBuilder, args: IndexedSeq[IEmitCode]): IEmitCode = {
    val p = args(0)
    val shape = args(1)
    val scale = args(2)

    p.flatMap(cb) { pCode =>
      shape.flatMap(cb) { shapeCode =>
        scale.flatMap(cb) { scaleCode =>
          cb.emb.newEmitCode(
            returnType,
            cb => {
              val pValue = cb.memoize(pCode.asFloat64.value)
              val shapeValue = cb.memoize(shapeCode.asFloat64.value)
              val scaleValue = cb.memoize(scaleCode.asFloat64.value)

              // Create Gamma distribution with shape and scale parameters
              val gammaDist = cb.memoize(
                cb.emb.newInstance[GammaDistribution](shapeValue, scaleValue)
              )

              // Calculate PPF (inverse CDF)
              val result = cb.memoize(
                gammaDist.inverseCumulativeProbability(pValue)
              )

              cb.assign(cb.emb.newRVariable("result"), result)
              cb.emb.getCode(result)
            }
          )
        }
      }
    }
  }
}

/**
 * Companion object to register the UDFs with Hail
 */
object GammaUDFRegistration {

  /**
   * Register all Gamma distribution UDFs with Hail
   */
  def registerAll(): Unit = {
    // Register the UDFs with Hail's function registry
    // This would typically be called during Hail initialization
    GammaPPF.register()
    QGamma.register()
  }
}
