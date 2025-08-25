package gnomad.constraint.promis3d

import is.hail.expr.ir.functions._
import org.apache.commons.math3.distribution.GammaDistribution

package object stats {
  def qgamma(p: Double, shape: Double, scale: Double): Double = {
    val gammaDist = new GammaDistribution(shape, scale)
    gammaDist.inverseCumulativeProbability(p)
  }
}