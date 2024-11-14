/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import uk.ac.sussex.gdsc.smlm.math3.analysis.integration.CustomSimpsonIntegrator;
import uk.ac.sussex.gdsc.smlm.utils.GaussianKernel;
import uk.ac.sussex.gdsc.smlm.utils.StdMath;

/**
 * This is a function to compute the likelihood assuming a Poisson-Gamma-Gaussian distribution.
 *
 * <p>For each observed value the log-likelihood is computed from the Poisson-Gamma-Gaussian
 * distribution (a Poisson convolved with a Gamma distribution convolved with a Gaussian). The
 * Poisson-Gamma distribution is derived analytically in the paper Ulbrich &amp; Isacoff (2007).
 * Nature Methods 4, 319-321 to explain the probability distribution of ADUs given a fixed photon
 * level per pixel and set gain in an EM-CCD camera (The Poisson distribution models the photon
 * count and the Gamma distribution models the EM-gain). This is then numerically convolved with a
 * Gaussian distribution to model the read noise of the camera.
 *
 * <p>The distribution of Ulbrich &amp; Isacoff has no analytical solution to the convolution with a
 * Gaussian. However the convolution with a Gaussian has the most effect when the counts are low.
 * The Poisson-Gamma-Gaussian can be approximated using the Poisson-Gamma and a partial convolution
 * with a Gaussian at low counts. This method is provided as Python source code within the
 * supplementary information of the paper Mortensen, et al (2010) Nature Methods 7, 377-383. This
 * Java implementation is based on the Python code.
 *
 * <p>The mean of the Poisson distribution is set using the expected value. The scale (EM-gain) for
 * the Gamma distribution and standard deviation of the Gaussian is fixed and set in the
 * constructor. The mean of the Gaussian is assumed to be zero.
 *
 * <p>The likelihood function is designed to model on-chip amplification of a EMCCD camera which
 * captures a Poisson process of emitted light, converted to electrons on the camera chip, amplified
 * by gain modelled by the Gamma distribution and then read with Gaussian noise.
 */
public class PoissonGammaGaussianFunction implements LikelihoodFunction, LogLikelihoodFunction {
  /**
   * The minimum read-noise. Values below this will be set to zerro and read-noise is not modelled.
   */
  public static final double MIN_READ_NOISE = 1e-3;

  /** Math.sqrt(2 * Math.PI). */
  private static final double SQRT_2PI = Math.sqrt(2 * Math.PI);

  /** The convolution mode. */
  private ConvolutionMode convolutionMode = ConvolutionMode.APPROXIMATION;

  /** The pmf mode flag for the convolution of the Dirac delta function at c=0. */
  private boolean pmfMode = true;

  /**
   * The scale of the Gamma distribution (e.g. the on-chip gain multiplication factor)
   */
  private final double gain;
  /**
   * The standard deviation of the Gaussian (e.g. Width of the noise distribution in the EMCCD
   * output)
   */
  private final double sigma;

  /** The two sigma 2. */
  private final double twoSigma2;

  /** The sqrt 2 sigma 2. */
  private final double sqrt2sigma2;

  /** The sqrt 2 pi sigma 2. */
  private final double sqrt2piSigma2;

  /** The minimum probability. */
  private double minimumProbability = Double.MIN_VALUE;

  /** The integrator. */
  private UnivariateIntegrator integrator;

  /** The kernel. */
  private double[] kernel;

  /**
   * Define the convolution mode for combining the Poisson-Gamma PMF with the Gaussian PDF.
   */
  public enum ConvolutionMode {
    /**
     * Use the approximation described in the Python source code within the supplementary
     * information of the paper Mortensen, et al (2010) Nature Methods 7, 377-383.
     */
    APPROXIMATION,

    /**
     * Convolve the Poisson-Gamma PDF with the Gaussian PDF using integer intervals. This method is
     * accurate when the read noise is far above 1. If the noise is below 1 then the
     * {@link #DISCRETE_PMF} is used. This computes a PDF(X=x).
     */
    DISCRETE_PDF,

    /**
     * Convolve the Poisson-Gamma PMF on discrete (integer) intervals with the Gaussian cumulative
     * probability density function. This method is accurate for all read noise. This computes a
     * PMF(X=x) on the integer scale.
     */
    DISCRETE_PMF,

    /**
     * Convolve the Poisson-Gamma PDF with the Gaussian PDF using Simpson's Rule. If integration
     * fails then the approximation will be used.
     */
    SIMPSON_PDF,

    /**
     * Convolve the Poisson-Gamma with the Gaussian PDF using a <a
     * href="http://mathworld.wolfram.com/Legendre-GaussQuadrature.html"> Legendre-Gauss</a>
     * quadrature. If integration fails then the approximation will be used.
     */
    LEGENDRE_GAUSS_PDF
  }

  /**
   * Provide the Poisson-Gamma-Gaussian function as a {@link UnivariateFunction}.
   */
  private class PggFunction implements UnivariateFunction {
    int counter;
    final double obs;
    final double exp;
    final double gain;

    PggFunction(double obs, double exp, double gain) {
      this.obs = obs;
      this.exp = exp;
      this.gain = gain;
    }

    @Override
    public double value(double count) {
      counter++;
      final double pg = PoissonGammaFunction.poissonGammaN(count, exp, gain);
      return (pg == 0) ? 0 : pg * gaussianPdf(count - obs);
    }
  }

  /**
   * Instantiates a new poisson gamma gaussian function.
   *
   * @param alpha Inverse gain of the EMCCD chip
   * @param sd The Gaussian standard deviation at readout
   * @throws IllegalArgumentException If the gain is below 1 (the Gamma distribution is not modelled
   *         for a scale below 1).
   */
  public PoissonGammaGaussianFunction(double alpha, double sd) {
    if (!(alpha > 0 && alpha <= 1)) {
      throw new IllegalArgumentException("Gain must be above 1");
    }
    this.gain = 1.0 / alpha;
    sd = Math.abs(sd);
    // Ignore tiny read noise
    if (sd < MIN_READ_NOISE) {
      sigma = twoSigma2 = sqrt2sigma2 = sqrt2piSigma2 = 0;
    } else {
      sigma = sd;
      twoSigma2 = 2 * sd * sd;
      sqrt2sigma2 = Math.sqrt(2 * sd * sd);
      sqrt2piSigma2 = Math.sqrt(2 * Math.PI * sd * sd);
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>This code is adapted from the Python source code within the supplementary information of the
   * paper Mortensen, et al (2010) Nature Methods 7, 377-383.
   *
   * <p>The output is a PMF. Ideally the input x should be discrete but this is not a requirement.
   */
  @Override
  public double likelihood(final double obs, final double exp) {
    // This did not speed up MLE fitting so has been commented out.
    // // When the observed ADUs and expected ADUs are much higher than the sigma then
    // // there is no point in convolving with a Gaussian
    //
    // double mySigma = sigma;
    // final double sLimit = sigma * 10;
    // if (o > sLimit && e > sLimit)
    // {
    // //System.out.println("Skipping convolution");
    // //mySigma = sigma;
    // }
    //
    // if (mySigma == 0)

    if (sigma == 0) {
      // No convolution with a Gaussian. Simply evaluate for a Poisson-Gamma distribution.
      // This can handle e<=0.
      return checkMinProbability(PoissonGammaFunction.poissonGamma(obs, exp, gain));
    }

    // If no Poisson mean then just use the Gaussian (Poisson-Gamma p=1 at x=0, p=0 otherwise)
    if (exp <= 0) {
      // Output as PMF
      return checkMinProbability(gaussianCdf(obs - 0.5, obs + 0.5));
    }

    if (convolutionMode == ConvolutionMode.APPROXIMATION) {
      return mortensenApproximation(obs, exp);
    }

    // Integrate to infinity is not necessary. The convolution of the function with the
    // Gaussian should be adequately sampled using a nxSD around the function value.
    // Find a bracket around the value.
    final double range = 5 * sigma;
    final double upperu = obs + range;
    if (upperu < 0) {
      return 0;
    }
    ConvolutionMode mode = convolutionMode;
    double loweru = obs - range;
    if (loweru < 0) {
      loweru = 0;
      // Edge case
      if (upperu == 0) {
        mode = ConvolutionMode.DISCRETE_PMF;
      }
    }

    double pvalue = 0;

    // PDF method does not work for small sigma
    if (mode == ConvolutionMode.DISCRETE_PDF && sigma < 1) {
      mode = ConvolutionMode.DISCRETE_PMF;
    }

    if (mode == ConvolutionMode.DISCRETE_PDF) {
      // Check
      if (loweru == upperu) {
        throw new IllegalStateException();
      }

      // Use a cached kernel using a range of 5
      final double[] g = createKernel(5);

      // Perform a simple convolution at point o with the kernel
      for (int i = 0, j = -g.length / 2; i < g.length; i++, j++) {
        final double u = obs - j;
        if (u >= 0.5) {
          // This approximates [u-0.5 to u+0.5] as a single point
          pvalue += PoissonGammaFunction.poissonGammaN(u, exp, gain) * g[i];
        } else if (u >= 0) {
          // Only count [0 to u+0.5] as a single point as the Poisson-Gamma
          // function is a step at x=0
          pvalue += PoissonGammaFunction.poissonGammaN(u, exp, gain) * g[i] * (0.5 + u);
        }
      }
    } else if (mode == ConvolutionMode.DISCRETE_PMF) {
      // Integrate the Poisson-Gamma to a PMF.
      // Convolve with the Gaussian CDF over the range (i.e. a discrete Gaussian PMF).

      // Use a simple integration by adding the points in the range.
      // Use the error function to obtain the integral of the Gaussian
      final int upper = (int) Math.ceil(upperu);
      int lower = (int) Math.floor(loweru);

      // Do an integration of the Poisson-Gamma PMF.
      // Trapezoid integration underestimates the total probability when the
      // Poisson-Gamma curve has no roots on the gradient (e.g. p<1).
      // (Since the trapezoid lines always miss part of the curve as it decays to zero.)
      // Simpson integration could be used to improve this.
      // Make this an option. For now just set to true as this mode should not be used
      // anyway. The Simpson integrator should be faster.
      final boolean doSimpson = true;

      double pg;
      // This is the CDF of the Gaussian
      double gauss;

      // If at zero then the Poisson-Gamma PMF approximation for u=0
      // is the integral from 0 to 0.5
      if (lower == 0) {
        // Lower = -0.5
        final double prevPg = PoissonGammaFunction.poissonGammaN(0, exp, gain);
        final double prevGauss = gaussianErf(-0.5 - obs);
        // Upper = 0.5
        pg = PoissonGammaFunction.poissonGammaN(0.5, pvalue, gain);
        gauss = gaussianErf(0.5 - obs);
        final double sum = (doSimpson)
            ? (prevPg + pg + 4 * PoissonGammaFunction.poissonGammaN(0.25, exp, gain)) / 12
            : (prevPg + pg) / 4;
        pvalue += sum * (gauss - prevGauss) / 2;
        lower++;
      } else {
        pg = PoissonGammaFunction.poissonGammaN(lower - 0.5, exp, gain);
        gauss = gaussianErf(lower - 0.5 - obs);
      }

      // For the rest of the range the Poisson-Gamma PMF approximation for u
      // is the integral from u-0.5 to u+0.5
      for (int u = lower; u <= upper; u++) {
        final double prevPg = pg;
        final double prevGauss = gauss;
        pg = PoissonGammaFunction.poissonGammaN(u + 0.5, exp, gain);
        gauss = gaussianErf(u + 0.5 - obs);
        final double sum =
            (doSimpson) ? (prevPg + pg + 4 * PoissonGammaFunction.poissonGammaN(u, exp, gain)) / 6
                : (prevPg + pg) / 2;
        pvalue += sum * (gauss - prevGauss) / 2;
      }
    } else {
      // Use integrator

      // Note that the Poisson-Gamma function has a delta function at u=0.
      // This prevents integration close to the zero boundary as the function
      // is not smooth.
      // So integration is done using the non-delta function version.

      // The integrator may be the fastest method when the range (upper-lower)
      // is large as it uses fewer points.

      // Specify the function to integrate.
      final PggFunction f = new PggFunction(obs, exp, gain);

      try {
        pvalue += createIntegrator().integrate(2000, f, loweru, upperu);
      } catch (final TooManyEvaluationsException ex) {
        Logger.getLogger(getClass().getName()).log(Level.WARNING,
            () -> String.format("Integration failed: o=%g, e=%g, eval=%d", obs, exp, f.counter));
        return mortensenApproximation(obs, exp);
      }
    }

    // Special case:
    // Due to the Poisson-Gamma function delta function at u=0
    // The probability at u=0 may be very large compared to the rest of
    // the Poisson-Gamma when e is low. To compensate always compute the
    // at u=0.

    // If this function is to be used as a PMF (with discrete integer observed values)
    // then use the Gaussian CDF convolution. This is the best option for
    // EM-CCD data fitting.
    // If this function is to be used as a PDF (e.g. for integration routines)
    // then use the Gaussian PDF convolution.

    if (pmfMode) {
      final double erf2 = gaussianErf(-obs + 0.5);
      if (erf2 != -1) {
        // Assume u==0
        pvalue += PoissonGammaFunction.dirac(exp) * (erf2 - gaussianErf(-obs - 0.5)) * 0.5;
      }
    } else {
      pvalue += PoissonGammaFunction.dirac(exp) * gaussianPdf(-obs);
    }

    return checkMinProbability(pvalue);
  }

  private double checkMinProbability(double pvalue) {
    return (pvalue > minimumProbability) ? pvalue : minimumProbability;
  }

  /**
   * Mortensen approximation.
   *
   * @param cij the cij
   * @param eta the eta
   * @return the double
   */
  private double mortensenApproximation(final double cij, final double eta) {
    // This code is adapted from the Python source code within the supplementary information of
    // the paper Mortensen, et al (2010) Nature Methods 7, 377-383.

    // The implementation of the approximation is not documented.
    // This is meant to be convolving a PMF of a Poisson-Gamma mixture with the PDF of a Gaussian.

    // See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
    // G n>0 (c) = sum n { (1 / n!) p^n e^-p (1 / ((n-1!)m^n)) c^n-1 e^-c/m }
    // G n>0 (c) = sqrt(p/(c*m)) * exp(-c/m - p) * Bessel.I1(2 * sqrt(c*p/m))
    // G n>0 (c) = sqrt(p/(c*m)) * exp(-c/m - p) * exp(2 * sqrt(c*p/m)) / sqrt(2*pi*sqrt(c*p/m))
    // G n=0 (c) = exp(-p)

    // p = eta
    // m = 1/alpha
    // c = cij

    // This is the value of the Poisson-Gamma at c=0:
    // PoissonGammaFunction.poissonGammaN(0, eta, m);
    final double exp_eta = StdMath.exp(-eta);
    final double f0 = exp_eta * eta / gain;

    // ?
    final double fp0 = f0 * 0.5 * (eta - 2) / gain;

    // The cumulative normal distribution of the read noise
    // at the observed count
    final double conv0 = gaussianCdf(cij);

    // [Noise * Gaussian PMF at observed count] +
    // [observed count * cumulative distribution of read noise at observed count]
    // [sigma*StdMath.exp(-cij**2/(twoSigma2))/Math.sqrt(2*pi)] + [cij*conv0]
    final double conv1 = sigma * StdMath.exp(-(cij * cij) / twoSigma2) / SQRT_2PI + cij * conv0;

    // ?
    double temp = f0 * conv0 + fp0 * conv1 + exp_eta * gaussianPdf(cij);

    // // TESTING
    // // Simple method:
    // temp = StdMath.exp(-eta) * gauss(cij); // G(c==0) * Gaussian;
    // if (cij <= 0)
    // return temp;
    // // Reset. The remaining will be the Poisson-Gamma and no convolution
    // f0 = fp0 = 0;
    //
    // // Q. How to normalise so that at low cij there is a mixture and at high cij there is no
    // mixture and the result is the Poisson-Gamma. Perhaps this is what the above code is doing.

    if (cij > 0) {
      temp += PoissonGammaFunction.poissonGammaN(cij, eta, gain) - f0 - fp0 * cij;
    }

    // XXX : Debugging: Store the smallest likelihood we ever see.
    // This can be used to set a limit for the likelihood
    // if (pMinObserved > temp && temp > 0)
    // {
    // pMinObserved = temp;
    // }

    return checkMinProbability(temp);
  }

  /**
   * {@inheritDoc}
   *
   * <p>This computes the log of {@link #likelihood(double, double)}.
   *
   * <p>The output is a PMF. Ideally the input x should be discrete but this is not a requirement.
   *
   * @param obs The observed count
   * @param exp The expected count
   * @return The log-likelihood
   *
   * @see #likelihood(double, double)
   */
  @Override
  public double logLikelihood(final double obs, final double exp) {
    return Math.log(likelihood(obs, exp));
  }

  /**
   * Gaussian PDF.
   *
   * @param x the x
   * @return the density
   */
  double gaussianPdf(final double x) {
    return StdMath.exp(-(x * x) / twoSigma2) / sqrt2piSigma2;
  }

  /**
   * Gaussian CDF.
   *
   * @param x the x
   * @return the cumulative density
   */
  double gaussianCdf(final double x) {
    // return 0.5 * org.apache.commons.numbers.gamma.Erfc.value(-x / sqrt2sigma2);
    // This may not be precise enough.
    // Absolute error is <3e-7. Not sure what relative error is.
    // The standard Erf is much slower.
    return 0.5 * (1 + Erf.erf(x / sqrt2sigma2));
  }

  /**
   * Gaussian CDF.
   *
   * @param x the x
   * @param x2 the x 2
   * @return the cumulative density
   */
  double gaussianCdf(final double x, final double x2) {
    // return 0.5 * org.apache.commons.numbers.gamma.ErfDifference.value(x / sqrt2sigma2, x2 / sqrt2sigma2);
    // This may not be precise enough.
    // Absolute error is <3e-7. Not sure what relative error is.
    // The standard Erf is much slower.
    return 0.5 * Erf.erf(x / sqrt2sigma2, x2 / sqrt2sigma2);
  }

  /**
   * Gaussian CDF.
   *
   * @param x the x
   * @return the cumulative density
   */
  double gaussianErf(final double x) {
    // return org.apache.commons.numbers.gamma.Erf.value(x / sqrt2sigma2);
    // This may not be precise enough.
    // Absolute error is <3e-7. Not sure what relative error is.
    // The standard Erf is much slower.
    return Erf.erf(x / sqrt2sigma2);
  }

  /**
   * Gets the alpha.
   *
   * @return the alpha
   */
  public double getAlpha() {
    return 1 / gain;
  }

  /**
   * Gets the sigma.
   *
   * @return the sigma
   */
  public double getSigma() {
    return sigma;
  }

  /**
   * Gets the minimum probability that will ever be returned. Setting this above zero allows the use
   * of Math.log() on the likelihood value.
   *
   * @return the minimum probability
   */
  public double getMinimumProbability() {
    return minimumProbability;
  }

  /**
   * Sets the minimum probability that will ever be returned. Setting this above zero allows the use
   * of Math.log() on the likelihood value.
   *
   * @param probability the new minimum probability
   */
  public void setMinimumProbability(double probability) {
    this.minimumProbability = probability;
  }

  /**
   * Gets the convolution mode.
   *
   * @return the convolution mode
   */
  public ConvolutionMode getConvolutionMode() {
    return convolutionMode;
  }

  /**
   * Sets the convolution mode.
   *
   * @param convolutionMode the new convolution mode
   */
  public void setConvolutionMode(ConvolutionMode convolutionMode) {
    this.convolutionMode = convolutionMode;
    integrator = null;
  }

  /**
   * Checks if is PMF mode. This is a hint on how to convolve the Dirac delta contribution at
   * observed count zero (c=0) with the Gaussian. If using the function for full integration then it
   * should be set to false. The default is true for modelling discrete count data from an EM-CCD.
   *
   * @return true, if PMF mode
   * @see #likelihood(double, double)
   * @see #getConvolutionMode()
   */
  public boolean isPmfMode() {
    return pmfMode;
  }

  /**
   * Sets the PMF mode flag. This is a hint on how to convolve the Dirac delta contribution at
   * observed count zero (c=0) with the Gaussian. If using the function for full integration then it
   * should be set to false. The default is true for modelling discrete count data from an EM-CCD.
   *
   * @param pmfMode the new PMF mode
   * @see #likelihood(double, double)
   * @see #setConvolutionMode(ConvolutionMode)
   */
  public void setPmfMode(boolean pmfMode) {
    this.pmfMode = pmfMode;
  }

  private UnivariateIntegrator createIntegrator() {
    UnivariateIntegrator in = integrator;
    if (in == null) {
      // This is the integrator for the Poisson-Gamma when observed count x>=1
      // i.e. not at the boundary x=0.

      final double relativeAccuracy = 1e-4;
      final double absoluteAccuracy = 1e-16;
      int minimalIterationCount;

      switch (convolutionMode) {
        case SIMPSON_PDF:
          // This is a CustomSimpsonIntegrator that computes 1 refinement
          // on the first iteration.
          // Number of function evaluations = 2^(iteration+1) + 1
          // => 5 for 1 iterations
          // => 9 for 2 iterations
          minimalIterationCount = 1;
          in = new CustomSimpsonIntegrator(relativeAccuracy, absoluteAccuracy,
              minimalIterationCount, CustomSimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
          break;

        case LEGENDRE_GAUSS_PDF:
          // Not sure how to configure this.
          // The integration points are used for each sub-interval.
          // Function evaluations = integrationpoints * intervals.
          // The intervals start at 1,2 and increase by at least 4 at each stage after that.
          // At least 1 stage is done thus 3 * integrationpoints functions evaluations
          // will be done for minimalIterationCount=1.
          minimalIterationCount = 1;
          final int maximalIterationCount = 32;
          final int integrationpoints = 8;
          in = new IterativeLegendreGaussIntegrator(integrationpoints, relativeAccuracy,
              absoluteAccuracy, minimalIterationCount, maximalIterationCount);
          break;

        default:
          throw new IllegalStateException();
      }

      integrator = in;
    }
    return in;
  }

  private double[] createKernel(int range) {
    double[] knl = kernel;
    if (knl == null) {
      kernel = knl = GaussianKernel.makeGaussianKernel(sigma, range, true);
    }
    return knl;
  }
}
