/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */

package uk.ac.sussex.gdsc.smlm.function;

import uk.ac.sussex.gdsc.core.utils.MathUtils;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

/**
 * Implements the probability density function for a combined Poisson distribution used to model a
 * Poisson-Gaussian mixture by modelling the Gaussian as an independent Poisson process.
 *
 * <p>This is a simple implementation of the LikelihoodFunction interface.
 *
 * <p>The likelihood function is designed to model on-chip amplification of a EMCCD/CCD/sCMOS camera
 * which captures a Poisson process of emitted light, converted to electrons on the camera chip,
 * amplified by a gain and then read.
 *
 * <p>This uses the MLE-sCMOS formula from Huang, et al (2013), Supplementary Notes Eq 3.3:<br>
 * P_sCMOS (x=[(Di-oi)/gi + vari/gi^2]|ui,vari,gi,oi) = e^-(ui+vari/gi^2) (ui+vari/gi^2)^x /
 * gamma(x+1) <br> Where:<br> i = the pixel index <br> vari = the variance of the pixel <br> gi =
 * the gain of the pixel <br> oi = the offset of the pixel <br> ui = the function value (expected
 * number of photons) <br> Di = the observed value at the pixel x = the observed random variable
 * (observed number of photons adjusted by a pixel dependent constant) <br>
 *
 * <p>The log-likelihood function is: <br> LL(P_sCMOS (x=[(Di-oi)/gi + vari/gi^2]|ui,vari,gi,oi))
 * <br> = -(ui+vari/gi^2) + x * ln(ui+vari/gi^2) - ln(gamma(x+1)) <br>
 *
 * <p>The negative log-likelihood (and gradient) can be evaluated over the entire set of observed
 * values or for a chosen observed value.
 *
 * <p>To allow a likelihood to be computed: (a) when the function predicts negative photon count
 * data the function prediction is set to zero; (b) if the observed random variable (x) is negative
 * it is also set to zero. This occurs when true signal readout from the sCMOS camera is low enough
 * to be negated by readout noise. In this case the noise can be ignored.
 *
 * @see "Hunag, et al (2013) Video-rate nanoscopy using sCMOS camera–specific single-molecule localization algorithms. Nature Methods 10, 653–658."
 */
public class PoissonPoissonFunction implements LikelihoodFunction, LogLikelihoodFunction {
  /**
   * The inverse of the on-chip gain multiplication factor.
   */
  final double alpha;

  /**
   * The log of the inverse on-chip gain multiplication factor.
   */
  final double logAlpha;

  /** The variance divided by the gain squared (i.e. variance in electron units) */
  final double var_g2;

  /**
   * Instantiates a new poisson poisson function.
   *
   * @param alpha The inverse of the on-chip gain multiplication factor
   * @param variance the variance
   */
  private PoissonPoissonFunction(double alpha, double variance) {
    if (variance <= 0) {
      throw new IllegalArgumentException("Gaussian variance must be strictly positive");
    }
    this.alpha = Math.abs(alpha);
    logAlpha = Math.log(this.alpha);
    var_g2 = variance * MathUtils.pow2(this.alpha);
  }

  /**
   * Creates the with standard deviation.
   *
   * @param alpha The inverse of the on-chip gain multiplication factor
   * @param s The standard deviation of the Gaussian distribution at readout
   * @return the poisson poisson function
   * @throws IllegalArgumentException if the variance is zero or below
   */
  public static PoissonPoissonFunction createWithStandardDeviation(final double alpha,
      final double s) {
    return new PoissonPoissonFunction(alpha, s * s);
  }

  /**
   * Creates the with variance.
   *
   * @param alpha The inverse of the on-chip gain multiplication factor
   * @param var The variance of the Gaussian distribution at readout (must be positive)
   * @return the poisson poisson function
   * @throws IllegalArgumentException if the variance is zero or below
   */
  public static PoissonPoissonFunction createWithVariance(final double alpha, final double var) {
    return new PoissonPoissonFunction(alpha, var);
  }

  /** {@inheritDoc} */
  @Override
  public double likelihood(double o, double e) {
    // convert to photons
    e = e + var_g2;
    o = o * alpha + var_g2;
    if (o < 0 || e <= 0) {
      return 0;
    }

    // Allow non-integer observed value using the gamma function to provide a
    // factorial for non-integer values. Standard Poisson:
    // PMF(l,k) = C * e^-l * l^k / gamma(k+1)
    // log(PMF) = -l + k * log(l) - logGamma(k+1)

    // P_sCMOS (x=[(Di-oi)/gi + vari/gi^2]|ui,vari,gi,oi)
    // = e^-(ui+vari/gi^2) (ui+vari/gi^2)^x / gamma(x+1) <br>

    // LL(P_sCMOS (x=[(Di-oi)/gi + vari/gi^2]|ui,vari,gi,oi)) <br>
    // = -(ui+vari/gi^2) + x * ln(ui+vari/gi^2) - ln(gamma(x+1)) <br>

    double ll = -e;
    if (o > 0) {
      ll += o * Math.log(e) - logFactorial(o);
    }

    // Scale the output so the cumulative probability is 1
    return FastMath.exp(ll) * alpha;
  }

  /**
   * Return the log of the factorial for the given real number, using the gamma function.
   *
   * @param k the number
   * @return the log factorial
   */
  public static double logFactorial(double k) {
    if (k <= 1) {
      return 0;
    }
    return Gamma.logGamma(k + 1);
  }

  /**
   * Return the factorial for the given real number, using the gamma function.
   *
   * @param k the number
   * @return the factorial
   */
  public static double factorial(double k) {
    if (k <= 1) {
      return 1;
    }
    return Gamma.gamma(k + 1);
  }

  /** {@inheritDoc} */
  @Override
  public double logLikelihood(double o, double e) {
    // convert to photons
    e = e + var_g2;
    o = o * alpha + var_g2;
    if (o < 0 || e <= 0) {
      return Double.NEGATIVE_INFINITY;
    }

    // Allow non-integer observed value using the gamma function to provide a
    // factorial for non-integer values. Standard Poisson:
    // PMF(l,k) = C * e^-l * l^k / gamma(k+1)
    // log(PMF) = -l + k * log(l) - logGamma(k+1)

    // P_sCMOS (x=[(Di-oi)/gi + vari/gi^2]|ui,vari,gi,oi)
    // = e^-(ui+vari/gi^2) (ui+vari/gi^2)^x / gamma(x+1) <br>

    // LL(P_sCMOS (x=[(Di-oi)/gi + vari/gi^2]|ui,vari,gi,oi)) <br>
    // = -(ui+vari/gi^2) + x * ln(ui+vari/gi^2) - ln(gamma(x+1)) <br>

    double ll = -e;
    if (o > 0) {
      ll += o * Math.log(e) - logFactorial(o);
    }

    // Scale the output so the cumulative probability is 1
    return ll + logAlpha;
  }
}
