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

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.util.FastMath;

import uk.ac.sussex.gdsc.smlm.math3.analysis.integration.CustomSimpsonIntegrator;
import uk.ac.sussex.gdsc.smlm.utils.GaussianKernel;

/**
 * This is a function to compute the likelihood assuming a Poisson-Gamma-Gaussian distribution. <p>
 * For each observed value the log-likelihood is computed from the Poisson-Gamma-Gaussian
 * distribution (a Poisson convolved with a Gamma distribution convolved with a Gaussian). The
 * Poisson-Gamma distribution is derived analytically in the paper Ulbrich &amp; Isacoff (2007).
 * Nature Methods 4, 319-321 to explain the probability distribution of ADUs given a fixed photon
 * level per pixel and set gain in an EM-CCD camera (The Poisson distribution models the photon
 * count and the Gamma distribution models the EM-gain). This is then numerically convolved with a
 * Gaussian distribution to model the read noise of the camera. <p> The distribution of Ulbrich
 * &amp; Isacoff has no analytical solution to the convolution with a Gaussian. However the
 * convolution with a Gaussian has the most effect when the counts are low. The
 * Poisson-Gamma-Gaussian can be approximated using the Poisson-Gamma and a partial convolution with
 * a Gaussian at low counts. This method is provided as Python source code within the supplementary
 * information of the paper Mortensen, et al (2010) Nature Methods 7, 377-383. This Java
 * implementation is based on the Python code. <P> The mean of the Poisson distribution is set using
 * the expected value. The scale (EM-gain) for the Gamma distribution and standard deviation of the
 * Gaussian is fixed and set in the constructor. The mean of the Gaussian is assumed to be zero. <p>
 * The likelihood function is designed to model on-chip amplification of a EMCCD camera which
 * captures a Poisson process of emitted light, converted to electrons on the camera chip, amplified
 * by gain modelled by the Gamma distribution and then read with Gaussian noise.
 */
public class PoissonGammaGaussianFunction implements LikelihoodFunction, LogLikelihoodFunction {
  /**
   * The minimum read-noise. Values below this will be set to zerro and read-noise is not modelled.
   */
  public static final double MIN_READ_NOISE = 1e-3;

  /** Math.sqrt(2 * Math.PI). */
  private static final double sqrt2pi = Math.sqrt(2 * Math.PI);

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

  private class PGGFunction implements UnivariateFunction {
    int i = 0;
    final double o, e;

    public PGGFunction(double o, double e) {
      this.o = o;
      this.e = e;
    }

    @Override
    public double value(double u) {
      i++;
      final double pg = PoissonGammaFunction.poissonGammaN(u, e, m);
      return (pg == 0) ? 0 : pg * gaussianPDF(u - o);
    }
  }

  /** The convolution mode. */
  private ConvolutionMode convolutionMode = ConvolutionMode.APPROXIMATION;

  /** The pmf mode flag for the convolution of the Dirac delta function at c=0. */
  private boolean pmfMode = true;

  /**
   * The scale of the Gamma distribution (e.g. the on-chip gain multiplication factor)
   */
  private final double m;
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
  private UnivariateIntegrator integrator = null;

  /** The kernel. */
  private double[] kernel = null;

  /**
   * Instantiates a new poisson gamma gaussian function.
   *
   * @param alpha Inverse gain of the EMCCD chip
   * @param s The Gaussian standard deviation at readout
   * @throws IllegalArgumentException If the gain is below 1 (the Gamma distribution is not modelled
   *         for a scale below 1).
   */
  public PoissonGammaGaussianFunction(double alpha, double s) throws IllegalArgumentException {
    if (!(alpha > 0 && alpha <= 1)) {
      throw new IllegalArgumentException("Gain must be above 1");
    }
    this.m = 1.0 / alpha;
    s = Math.abs(s);
    // Ignore tiny read noise
    if (s < MIN_READ_NOISE) {
      sigma = twoSigma2 = sqrt2sigma2 = sqrt2piSigma2 = 0;
    } else {
      sigma = s;
      twoSigma2 = 2 * s * s;
      sqrt2sigma2 = Math.sqrt(2 * s * s);
      sqrt2piSigma2 = Math.sqrt(2 * Math.PI * s * s);
    }
  }

  /**
   * {@inheritDoc} <p> This code is adapted from the Python source code within the supplementary
   * information of the paper Mortensen, et al (2010) Nature Methods 7, 377-383. <p> The output is a
   * PMF. Ideally the input x should be discrete but this is not a requirement.
   *
   * @see uk.ac.sussex.gdsc.smlm.function.LikelihoodFunction#likelihood(double, double)
   */
  @Override
  public double likelihood(final double o, final double e) {
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
      return checkMinProbability(PoissonGammaFunction.poissonGamma(o, e, m));
    }

    // If no Poisson mean then just use the Gaussian (Poisson-Gamma p=1 at x=0, p=0 otherwise)
    if (e <= 0) {
      // Output as PMF
      return checkMinProbability(gaussianCDF(o - 0.5, o + 0.5));
    }

    if (convolutionMode == ConvolutionMode.APPROXIMATION) {
      return mortensenApproximation(o, e);
    }

    ConvolutionMode mode = convolutionMode;

    // Integrate to infinity is not necessary. The convolution of the function with the
    // Gaussian should be adequately sampled using a nxSD around the function value.
    // Find a bracket around the value.
    final double range = 5 * sigma;
    final double upperu = o + range;
    if (upperu < 0) {
      return 0;
    }
    double loweru = o - range;
    if (loweru < 0) {
      loweru = 0;
      // Edge case
      if (upperu == 0) {
        mode = ConvolutionMode.DISCRETE_PMF;
      }
    }

    double p = 0;

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
        final double u = o - j;
        if (u >= 0.5) {
          // This approximates [u-0.5 to u+0.5] as a single point
          p += PoissonGammaFunction.poissonGammaN(u, e, m) * g[i];
        } else if (u >= 0) {
          // Only count [0 to u+0.5] as a single point as the Poisson-Gamma
          // function is a step at x=0
          p += PoissonGammaFunction.poissonGammaN(u, e, m) * g[i] * (0.5 + u);
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
      double g;

      // If at zero then the Poisson-Gamma PMF approximation for u=0
      // is the integral from 0 to 0.5
      if (lower == 0) {
        // Lower = -0.5
        final double prevPG = PoissonGammaFunction.poissonGammaN(0, e, m);
        final double prevG = gaussianErf(-0.5 - o);
        // Upper = 0.5
        pg = PoissonGammaFunction.poissonGammaN(0.5, p, m);
        g = gaussianErf(0.5 - o);
        final double sum =
            (doSimpson) ? (prevPG + pg + 4 * PoissonGammaFunction.poissonGammaN(0.25, e, m)) / 12
                : (prevPG + pg) / 4;
        p += sum * (g - prevG) / 2;
        lower++;
      } else {
        pg = PoissonGammaFunction.poissonGammaN(lower - 0.5, e, m);
        g = gaussianErf(lower - 0.5 - o);
      }

      // For the rest of the range the Poisson-Gamma PMF approximation for u
      // is the integral from u-0.5 to u+0.5
      for (int u = lower; u <= upper; u++) {
        final double prevPG = pg;
        final double prevG = g;
        pg = PoissonGammaFunction.poissonGammaN(u + 0.5, e, m);
        g = gaussianErf(u + 0.5 - o);
        final double sum =
            (doSimpson) ? (prevPG + pg + 4 * PoissonGammaFunction.poissonGammaN(u, e, m)) / 6
                : (prevPG + pg) / 2;
        p += sum * (g - prevG) / 2;
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
      final PGGFunction f = new PGGFunction(o, e);

      try {
        p += createIntegrator().integrate(2000, f, loweru, upperu);
      } catch (final TooManyEvaluationsException ex) {
        System.out.printf("Integration failed: o=%g, e=%g, eval=%d\n", o, e, f.i);
        return mortensenApproximation(o, e);
      }
      // System.out.printf("Integration eval=%d\n", f.i);
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
      final double erf2 = gaussianErf(-o + 0.5);
      if (erf2 != -1) {
        // Assume u==0
        p += PoissonGammaFunction.dirac(e) * (erf2 - gaussianErf(-o - 0.5)) * 0.5;
      }
    } else {
      p += PoissonGammaFunction.dirac(e) * gaussianPDF(-o);
    }

    return checkMinProbability(p);
  }

  private double checkMinProbability(double p) {
    return (p > minimumProbability) ? p : minimumProbability;
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
    final double exp_eta = FastMath.exp(-eta);
    final double f0 = exp_eta * eta / m;

    // ?
    final double fp0 = f0 * 0.5 * (eta - 2) / m;

    // The cumulative normal distribution of the read noise
    // at the observed count
    final double conv0 = gaussianCDF(cij);

    // [Noise * Gaussian PMF at observed count] +
    // [observed count * cumulative distribution of read noise at observed count]
    // [sigma*FastMath.exp(-cij**2/(twoSigma2))/Math.sqrt(2*pi)] + [cij*conv0]
    final double conv1 = sigma * FastMath.exp(-(cij * cij) / twoSigma2) / sqrt2pi + cij * conv0;

    // ?
    double temp = f0 * conv0 + fp0 * conv1 + exp_eta * gaussianPDF(cij);

    // // TESTING
    // // Simple method:
    // temp = FastMath.exp(-eta) * gauss(cij); // G(c==0) * Gaussian;
    // if (cij <= 0)
    // return temp;
    // // Reset. The remaining will be the Poisson-Gamma and no convolution
    // f0 = fp0 = 0;
    //
    // // Q. How to normalise so that at low cij there is a mixture and at high cij there is no
    // mixture
    // // and the result is the Poisson-Gamma. Perhaps this is what the above code is doing.

    if (cij > 0.0) {
      temp += PoissonGammaFunction.poissonGammaN(cij, eta, m) - f0 - fp0 * cij;
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
   * {@inheritDoc} <p> This computes the log of {@link #likelihood(double, double)}. <p> The output
   * is a PMF. Ideally the input x should be discrete but this is not a requirement.
   *
   * @param o The observed count
   * @param e The expected count
   * @return The log-likelihood
   *
   * @see #likelihood(double, double)
   * @see uk.ac.sussex.gdsc.smlm.function.LogLikelihoodFunction#logLikelihood(double, double)
   */
  @Override
  public double logLikelihood(final double o, final double e) {
    return Math.log(likelihood(o, e));
  }

  /**
   * Gaussian PDF.
   *
   * @param x the x
   * @return the density
   */
  double gaussianPDF(final double x) {
    return FastMath.exp(-(x * x) / twoSigma2) / sqrt2piSigma2;
  }

  /**
   * Gaussian CDF.
   *
   * @param x the x
   * @return the cumulative density
   */
  double gaussianCDF(final double x) {
    // return 0.5 * (1 + org.apache.commons.math3.special.Erf.erf(x / sqrt2sigma2));
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
  double gaussianCDF(final double x, final double x2) {
    // return 0.5 * (org.apache.commons.math3.special.Erf.erf(x / sqrt2sigma2, x2 / sqrt2sigma2));
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
    // return org.apache.commons.math3.special.Erf.erf(x / sqrt2sigma2);
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
    return 1 / m;
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
   * @param p the new minimum probability
   */
  public void setMinimumProbability(double p) {
    this.minimumProbability = p;
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
    UnivariateIntegrator i = integrator;
    if (i == null) {
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
          i = new CustomSimpsonIntegrator(relativeAccuracy, absoluteAccuracy, minimalIterationCount,
              CustomSimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
          break;

        case LEGENDRE_GAUSS_PDF:
          // Not sure how to configure this.
          // The integration points are used for each sub-interval.
          // Function evaluations = integrationPoints * intervals.
          // The intervals start at 1,2 and increase by at least 4 at each stage after that.
          // At least 1 stage is done thus 3 * integrationPoints functions evaluations
          // will be done for minimalIterationCount=1.
          minimalIterationCount = 1;
          final int maximalIterationCount = 32;
          final int integrationPoints = 8;
          i = new IterativeLegendreGaussIntegrator(integrationPoints, relativeAccuracy,
              absoluteAccuracy, minimalIterationCount, maximalIterationCount);
          break;

        default:
          throw new IllegalStateException();
      }

      integrator = i;
    }
    return i;
  }

  private double[] createKernel(int range) {
    double[] g = kernel;
    if (g == null) {
      kernel = g = GaussianKernel.makeGaussianKernel(sigma, range, true);
    }
    return g;
  }
}
