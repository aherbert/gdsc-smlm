/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleUnaryOperator;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optim.AbstractConvergenceChecker;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.util.FastMath;
import uk.ac.sussex.gdsc.core.math.NumberUtils;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.utils.Convolution;
import uk.ac.sussex.gdsc.smlm.utils.Convolution.DoubleConvolutionValueProcedure;
import uk.ac.sussex.gdsc.smlm.utils.GaussianKernel;
import uk.ac.sussex.gdsc.smlm.utils.StdMath;

/**
 * Calculate the Fisher information for a Poisson-Gamma-Gaussian distribution.
 *
 * <p>Uses a modified form of the equation of Chao, et al (2013) Nature Methods, 10, 335-338, SI Eq
 * S8.
 *
 * <p>Performs a convolution with a finite Gaussian kernel.
 */
public class PoissonGammaGaussianFisherInformation extends BasePoissonFisherInformation {
  /**
   * Class to scale the lists.
   */
  private static class ScalingFunction implements DoubleUnaryOperator {
    final double scale;

    ScalingFunction(double scale) {
      this.scale = scale;
    }

    @Override
    public double apply(double value) {
      return value * scale;
    }
  }

  private static final int UPPER_LIMIT = Integer.MAX_VALUE - 1;

  /** The default minimum range for the Gaussian kernel (in units of SD). */
  public static final int DEFAULT_MIN_RANGE = 6;

  /** The default maximum range for the Gaussian kernel (in units of SD). */
  public static final int DEFAULT_MAX_RANGE = 38;

  /** The maximum range for the Gaussian kernel (in units of SD). */
  public static final int MAX_RANGE = 38;

  /** The log2 of the maximum scale for the Gaussian kernel (in units of SD). */
  public static final int LOG_2_MAX_SCALE = 16;

  /** The maximum scale for the Gaussian kernel (in units of SD). */
  public static final int MAX_SCALE = 1 << LOG_2_MAX_SCALE;

  /** The default threshold for the relative probability. */
  public static final double DEFAULT_RELATIVE_PROBABILITY_THRESHOLD = 1e-5;

  /** The default threshold for the cumulative probability. */
  public static final double DEFAULT_CUMULATIVE_PROBABILITY = 1 - 1e-6;

  /**
   * The default sampling of the Gaussian kernel. The kernel will be sampled at s/sampling, i.e.
   * this is the number of samples to take per standard deviation unit.
   */
  public static final int DEFAULT_SAMPLING = 2;

  /** The default relative accuracy for convergence. */
  public static final double DEFAULT_RELATIVE_ACCURACY = 1e-6;

  /**
   * The default number of maximum iterations. The Gaussian sampling will be doubled during each
   * iteration.
   */
  public static final int DEFAULT_MAX_ITERATIONS = 6;

  /**
   * Minimum relative tolerance.
   */
  private static final double MIN_RELATIVE_TOLERANCE = 2 * FastMath.ulp(1d);

  /** The gain multiplication factor. */
  public final double gain;

  /** The standard deviation of the Gaussian. */
  public final double sd;

  /** The minimum range of the Gaussian kernel (in SD units). */
  private int minRange = DEFAULT_MIN_RANGE;

  /** The maximum range of the Gaussian kernel (in SD units). */
  private int maxRange = DEFAULT_MAX_RANGE;

  /** The default scale of the kernel. */
  private final int defaultScale;

  /** Working space to store the probabilities. */
  private final DoubleArrayList listP = new DoubleArrayList();
  /** Working space to store the gradient. */
  private final DoubleArrayList listA = new DoubleArrayList();

  /**
   * The relative probability threshold of the Poisson-Gamma distribution that is used. Any
   * probability less than this is ignored. This prevents using values that do not contribute to the
   * sum.
   */
  private double relativeProbabilityThreshold = DEFAULT_RELATIVE_PROBABILITY_THRESHOLD;

  /** The upper mean threshold for the switch to half the Poisson Fisher information. */
  private double upperMeanThreshold = 200;

  /** Set to true to use Simpson's 3/8 rule for cubic interpolation of the integral. */
  private boolean use38 = true;

  /** Flag to indicate no possible convolution with the gaussian. */
  private final boolean noGaussian;

  /** The gaussian kernel. */
  private GaussianKernel gaussianKernel;

  /** The relative accuracy for convergence. */
  private double relativeAccuracy = DEFAULT_RELATIVE_ACCURACY;

  /** The max iterations. */
  private int maxIterations = DEFAULT_MAX_ITERATIONS;

  /** The last Gaussian kernel. */
  private double[] lastG;

  /** The scale of the last Gaussian kernel. */
  private int lastScale;

  /** The mean of the last Fisher information. */
  private double lastT;

  /**
   * Instantiates a new poisson gamma gaussian fisher information.
   *
   * @param gain the gain multiplication factor
   * @param sd the standard deviation of the Gaussian
   * @throws IllegalArgumentException If the standard deviation is not strictly positive
   */
  public PoissonGammaGaussianFisherInformation(double gain, double sd) {
    this(gain, sd, DEFAULT_SAMPLING);
  }

  /**
   * Instantiates a new poisson gamma gaussian fisher information.
   *
   * @param gain the gain multiplication factor
   * @param sd the standard deviation of the Gaussian
   * @param sampling The number of Gaussian samples to take per standard deviation
   * @throws IllegalArgumentException If the gain is not strictly positive
   * @throws IllegalArgumentException If the standard deviation is below 1
   * @throws IllegalArgumentException If the sampling is below 1
   * @throws IllegalArgumentException If the maximum kernel size after scaling is too large
   */
  public PoissonGammaGaussianFisherInformation(double gain, double sd, double sampling) {
    if (!(gain > 0 && gain <= Double.MAX_VALUE)) {
      throw new IllegalArgumentException("Gain multiplication factor must be strictly positive");
    }
    if (!(sd >= 0 && sd <= Double.MAX_VALUE)) {
      throw new IllegalArgumentException("Gaussian standard deviation must be positive");
    }
    if (!(sampling >= 1 && sampling <= Double.MAX_VALUE)) {
      throw new IllegalArgumentException("Gaussian sampling must at least 1");
    }

    this.gain = gain;
    this.sd = sd;
    noGaussian = (sd * MAX_RANGE < 1);
    if (noGaussian) {
      // This is OK. Just return the information for the Poisson-Gamma.
      defaultScale = 0;
    } else {
      // This is set to work for reasonable values of the Gaussian kernel and sampling
      // e.g. s [0.5:20], sampling from [1:8].

      defaultScale = getPow2Scale(sampling / sd);

      // Don'theta support excess scaling caused by small kernels
      if (defaultScale * sd * MAX_RANGE > 1000000000) {
        throw new IllegalArgumentException(
            "Maximum Gaussian kernel size too large: " + defaultScale * sd * MAX_RANGE);
      }

      gaussianKernel = new GaussianKernel(sd);
    }
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected PoissonGammaGaussianFisherInformation(PoissonGammaGaussianFisherInformation source) {
    gain = source.gain;
    sd = source.sd;
    minRange = source.minRange;
    maxRange = source.maxRange;
    defaultScale = source.defaultScale;
    relativeProbabilityThreshold = source.relativeProbabilityThreshold;
    upperMeanThreshold = source.upperMeanThreshold;
    use38 = source.use38;
    noGaussian = source.noGaussian;
    if (source.gaussianKernel != null) {
      gaussianKernel = source.gaussianKernel.copy();
    }
    relativeAccuracy = source.relativeAccuracy;
    maxIterations = source.maxIterations;
    lastG = source.lastG;
    lastScale = source.lastScale;
    lastT = source.lastT;
  }

  @Override
  public PoissonGammaGaussianFisherInformation copy() {
    return new PoissonGammaGaussianFisherInformation(this);
  }

  /**
   * Gets a value using the next integer power of 2. This is limited to a size of 2^16.
   *
   * @param scaleValue the value
   * @return the next power of 2
   */
  protected static int getPow2Scale(double scaleValue) {
    final double scale = Math.ceil(scaleValue);
    if (scale > MAX_SCALE) {
      return MAX_SCALE;
    }
    return MathUtils.nextPow2((int) scale);
  }

  /**
   * {@inheritDoc}
   *
   * <p>The input parameter refers to the mean of the Poisson distribution.
   *
   * <p>The Fisher information is computed using the equation of Chao, et al (2013) Nature Methods,
   * 10, 335-338, SI Eq S8. Note that that equation computes the noise coefficient relative to a
   * Poisson, this computes the Fisher information. To get the noise coefficient multiply by the
   * input parameter.
   */
  @Override
  public double getFisherInformation(double theta) {
    final double I = getPoissonGammaGaussianI(theta);

    // Check limits.

    // It should be worse than the Poisson Fisher information (upper limit).
    // Note a low Fisher information is worse as this is the amount of information
    // carried about the parameter.
    final double upper = 1.0 / theta; // PoissonFisherInformation.getPoissonI(theta)
    // When the mean is high then the lower bound should be half the
    // Poisson Fisher information. Otherwise check against 0.
    // This is not true when the gain is low so leave out for now.
    // final double lower = (theta > 10) ? 0.5 * upper : 0;
    final double lower = 0;
    return MathUtils.clip(lower, upper, I);
  }

  @Override
  public double getAlpha(double theta) {
    // Don'theta check against MIN_MEAN as this is done later if it is necessary
    // to compute the Fisher information.

    if (theta > upperMeanThreshold) {
      // Use an approximation as half the Poisson Fisher information
      return 0.5;
    }

    return theta * getFisherInformation(theta);
  }

  /**
   * Gets the Poisson-Gamma-Gaussian Fisher information.
   *
   * <p>The input parameter refers to the mean of the Poisson distribution.
   *
   * <p>The Fisher information is computed using the equation of Chao, et al (2013) Nature Methods,
   * 10, 335-338, SI Eq S8. Note that that equation computes the noise coefficient relative to a
   * Poisson, this computes the Fisher information. To get the noise coefficient multiply by the
   * input parameter.
   *
   * <p>Note: This uses a convolution of an infinite integral over a finite range. It may
   * under-estimate the information when the mean is large.
   *
   * @param theta the Poisson mean
   * @return the Poisson Gaussian Fisher information
   * @throws IllegalArgumentException the illegal argument exception
   */
  public double getPoissonGammaGaussianI(double theta) {
    // Reset
    listP.clear();
    listA.clear();
    lastG = null;
    lastT = theta;

    if (theta <= 0) {
      throw new IllegalArgumentException("Poisson mean must be positive");
    }
    if (theta < MIN_MEAN) {
      return Double.POSITIVE_INFINITY;
    }

    if (theta > upperMeanThreshold) {
      return largeApproximation(theta);
    }

    final double dirac = PoissonGammaFunction.dirac(theta);

    // Special case where the start of the range will have no probability.
    // This occurs when mean > limit of exp(-p) => 746.
    if (1 / dirac == Double.POSITIVE_INFINITY) {
      return largeApproximation(theta);
    }

    // This computes the convolution of a Poisson-Gamma PDF and a Gaussian PDF.
    // The value of this is p(z).

    // The Poisson-Gamma-Gaussian must be differentiated to compute the Fisher information:
    // Expected [ (d ln(p(z)) dv)^2 ]
    // = integral [ (1/p(z) . d p(z) dv)^2 p(z) dz ]
    // = integral [ 1/p(z) . (d p(z) dv)^2 dz ]

    // Gaussian standard deviation = s

    // Note: (fg)' => fg' + f'g
    // e^-v => -e^-v
    // (e^-v * g(v))' => e^-v * g'(v) - e^-v * g(v)
    // So for a probability distribution with a e^-v component:
    // p(z|v) = e^-v * g(v)
    // p'(z|v) = e^-v * g'(v) - p(z|v)

    // -=-=-=-

    // This Fisher information is based on Chao, et al (2013) Nature Methods, 10, 335-338, SI.

    // Chao et al, Eq 4:
    // (Poisson Geometric Gaussian convolution)
    // P(z) = e^-v / (sqrt(2pi)*s) * [ e^-0.5*(z/s)^2 +
    // sum_l=1_inf [ e^-0.5*((z-l)/s)^2 *
    // sum_j=0_l-1 [ (l-1)!*(1-1/g)^(l-1-j)*(v/g)^(j+1) / ( j!*(l-1-j)!*(j+1)! ) ] ]
    // ]
    // Note: dividing by (g/v)^(j+1) has been changed to multiplying by (v/g)^(j+1)

    // Gradient of the Poisson Geometric Gaussian is:
    // P'(z|v) = e^-v / (sqrt(2pi)*s) * [
    // sum_l=1_inf [ e^-0.5*((z-l)/s)^2 *
    // sum_j=0_l-1 [ (l-1)!*(1-1/g)^(l-1-j)*(v/g)^j / ( j!*(l-1-j)!*j!*g ) ] ]
    // ]
    // - P(z)
    // This is because d((v/g)^(j+1)) dv = (j+1)*(v/g)^j*(1/g)

    // This component can be seen in Chao et al, Eq S8:
    // FI = E e^-2v / P(z) * [
    // sum_l=1_inf [ e^-0.5*((z-l)/s)^2 / (sqrt(2pi)*s*g) *
    // sum_j=0_l-1 [ (l-1)!*(1-1/g)^(l-1-j)*(v/g)^j / ( j!*(l-1-j)!*j! ) ] ]
    // ]^2
    // - 1

    // This equation is:
    // integral [ 1/p(z) . a(z|v)^2 dz ] - 1
    // Where a(z|v) is the partial gradient of the Poisson-Geometric-Gaussian with respect to v
    // without the subtraction of P(z).

    // Note that the partial gradient of the Poisson-Geometric-Gaussian is just the
    // partial gradient of the Poisson-Geometric convolved with a Gaussian (since the Gaussian does
    // not have a component v).

    // The Poisson-Geometric is computationally intense. Thus this can be substituted using
    // the Poisson-Gamma convolution which can be computed efficiently using Bessel functions.

    // Poisson-Gamma-Gaussian:
    // P(z) = e^-v / (sqrt(2pi)*s) * [ e^-0.5*(z/s)^2 +
    // sum_l=0_inf [ e^-0.5*((z-l)/s)^2 * e^-l/m *
    // sum_n=1_inf [ (1 / n!(n-1)!) * v^n * l^(n-1) / m^n ) ] ]
    // ]

    // P'(z) = e^-v / (sqrt(2pi)*s) * [
    // sum_l=0_inf [ e^-0.5*((z-l)/s)^2 * e^-l/m *
    // 1/m * sum_n=0_inf [ (1 / (n!)^2) * v^n * l^n / m^n ) ] ]
    // ]
    // - P(z)

    // Note:
    // sum_n=0_inf [ (1 / (n!)^2) * v^n * l^n / m^n ) ] = sum_n=0_inf [ (1 / (n!)^2) * (v*l/m)^n) ]
    // Bessel I0 = sum_n=0_inf [ (1 / (n!)^2) * (x^2 / 4)^n
    // set x = 2 * sqrt(v*l/m)
    // = I0(x)

    // P'(z) = e^-v / (sqrt(2pi)*s) * [
    // sum_l=0_inf [ e^-0.5*((z-l)/s)^2 * e^-l/m * 1/m * I0(x) ]
    // ]
    // - P(z)

    // Fisher information is:
    // FI = E e^-2v / (P(z)) * [ sum_l=0_inf [ e^-0.5*((z-l)/s)^2 / (sqrt(2pi)*s) *
    // e^-l/m * 1/m * I0(x) ] ]^2
    // - 1
    // (where the Gaussian normalisation has been moved inside the power 2 bracket)

    // Note the Poisson-Gamma can be computed using the approximation defined in
    // Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
    // PG(l|v) = sqrt(v/(l*m)) * e^(-l/m -v) * I1(x)
    // x = 2*sqrt(v*l/m)
    // The partial gradient as defined above is:
    // PG'(l|v) = e^(-l/m -v) * I0(x) / m - PG(l|v)

    // We need the convolution of:
    // - the Poisson-Gamma with the Gaussian (P):
    // - the Poisson-Gamma partial gradient component with the Gaussian (A):
    // This is combined to a function A^2/P which is integrated.

    // Note that A & P should sum to 1 over the full range used for integration.
    // In practice a relative height threshold of the function can be used
    // to determine the range where it is significant and the sums normalised to 1.

    // Fisher information is:
    // FI = E [ 1 / (P(z)) * [ A(z) ]^2 ] - 1

    // -=-=-=-

    // Note: e^-v is a scaling factor used in P and A.
    // Compute P and A without the scaling factor as Pu and Au.
    // This stabilises the computation at high mean (v).

    // Fisher information is:
    // FI = [e^-v * E [ 1 / (Pu(z)) * [ Au(z) ]^2 ] ] - 1

    // Pu and Au will sum to 1/e^v.

    // -=-=-=-

    // Note that when the Poisson mean is low then the contribution from the
    // Dirac delta function at 0 will be large. As the mean increases
    // the Dirac will have little significance, the function is very smooth
    // and may have a large range. Thus convolution with an extended Gaussian kernel
    // may not be necessary over the full range.
    // So the integral can be computed as [integral range1] + [integral range2]

    // Determine the extent of the function: A^2/P
    final double[] maximum = findMaximum(theta, 1e-4);
    final double[] upper = findUpperLimit(theta, maximum, relativeProbabilityThreshold);
    double upperLimit = upper[0];

    if (upperLimit > UPPER_LIMIT) {
      throw new IllegalStateException("Unsupported upper limit: " + upperLimit);
    }

    // Note: This method works by using the Poisson-Gamma PDF as a substitute
    // for a Poisson-Geometric PMF. For simplicity the PDF is integrated over
    // integer intervals to create a discrete PMF. This is a fair approximation
    // of how an EM-CCD camera would work as input electrons should produce
    // discrete output electrons.

    // The exponent provides a rough idea of the size of the mean (i.e. log2(mean))
    final int exp = NumberUtils.getSignedExponent(theta);
    // As the mean reduces the Poisson distribution is more skewed
    // and the extent of the kernel must change. Just increase the range
    // for the kernel for each power of 2 the number is below 1.
    int range1 = minRange;
    for (int e = exp; range1 < maxRange && e <= 0; e++, range1++) {
      // All work in the loop update clause
    }

    if (noGaussian || sd * range1 < 1) {
      // Special case. Compute E [ A^2/P ] using only the Poisson-Gamma PMF.
      if (!computePmf(dirac, theta, maximum, (int) Math.ceil(upperLimit))) {
        return largeApproximation(theta);
      }

      // No interpolated integration necessary as this is discrete
      double sum = 0;
      final double[] p = listP.elements();
      final double[] a = listA.elements();
      for (int i = 0, size = listP.size(); i < size; i++) {
        sum += getF(p[i], a[i]);
      }
      if (sum == Double.POSITIVE_INFINITY) {
        return extremeLimit(theta);
      }
      // Remember to rescale to the correct range
      return dirac * sum - 1;
    }

    // Extend the upper limit using the configured kernel range
    upperLimit += sd * range1;
    if (upperLimit > UPPER_LIMIT) {
      throw new IllegalStateException("Unsupported upper limit: " + upperLimit);
    }

    // // TODO:
    // // Set range1 using the range of the Gaussian.
    // // This is the range around c=0 where the Dirac has a large effect.
    // // Set range2 using the rest of the range.
    // // Iteratively integrate range 1 up to a maximum kernel size.
    // // Iteratively Integrate range 2.
    // // Combine the sums.
    //
    // // Range 1:
    // // The range around the Dirac delta function is computed using a full convolution
    // // with the Gaussian. This is iterated until convergence using more Gaussian samples.
    //
    // // Range 2:
    // // The remaining range of the Poisson-Gamma is integrated using an integrator until the
    // // the Gaussian kernel is non-trivial.
    //
    // // For integration to work efficiently the Gaussian scaling should be in powers of 2.
    // // Adjust the ranges appropriately so that the ranges start and end on powers of 1.
    // int iEndRange1 = Maths.nextPow2((int) Math.ceil(s * range1));
    // int iEndRange2 = iEndRange1;
    // if (upperLimit > iEndRange1)
    // iEndRange2 += Maths.nextPow2((int) Math.ceil(upperLimit - iEndRange1));

    // Compute the Poisson-Gamma over the range. This is a PDF but for
    // simplicity it is integrated to a discrete PMF.
    final int iEndRange = (int) Math.ceil(upperLimit);
    if (!computePmf(dirac, theta, maximum, iEndRange)) {
      return largeApproximation(theta);
    }

    final double[] p = listP.toDoubleArray();
    final double[] a = listA.toDoubleArray();

    // Find the minimum value to determine if A^2/P can be unchecked
    // for divide by zero error.
    // The min is always at the end.
    final double minP = Math.min(Math.min(p[0], p[1]), p[p.length - 1]);

    // int iMinP = SimpleArrayUtils.findMinIndex(p);
    // int iMinA = SimpleArrayUtils.findMinIndex(a);
    // int iMaxP = SimpleArrayUtils.findMaxIndex(p);
    // int iMaxA = SimpleArrayUtils.findMaxIndex(a);
    // minP = p[iMinP];
    // double minA = a[iMinA];
    // double maxP = p[iMaxP];
    // double maxA = a[iMaxA];
    // System.out.printf("m=%g s=%s p=%g n=%d MinP=%s @ %d, MaxP=%s @ %d, MinA=%s @ %d, MaxA=%s @
    // %d\n", m, s,
    // theta, p.length, minP, iMinP, maxP, iMaxP, minA, iMinA, maxA, iMaxA);

    // Initial sum
    int scale = defaultScale;
    double sum = compute(scale, range1, p, a, minP);

    if (sum == Double.POSITIVE_INFINITY) {
      return extremeLimit(theta);
    }

    // Iterate
    for (int iteration = 1; iteration <= maxIterations && scale < MAX_SCALE; iteration++) {
      scale *= 2;
      if (iteration == maxIterations || scale == MAX_SCALE) {
        // System.out.printf("s=%g theta=%g final iteration=%d scale=%d\n", s, theta, iteration,
        // scale);
      }
      final double oldSum = sum;
      try {
        sum = compute(scale, range1, p, a, minP);
        if (sum == Double.POSITIVE_INFINITY) {
          sum = oldSum;
          break;
        }
      } catch (final IllegalArgumentException ex) {
        // Occurs when the convolution has grown too big
        break;
      }
      final double delta = Math.abs(sum - oldSum);
      // System.out.printf("s=%g theta=%g Iteration=%d sum=%s oldSum=%s change=%s\n", s, theta,
      // iteration,
      // sum, oldSum,
      // delta / (Math.abs(oldSum) + Math.abs(sum)) * 0.5);
      // final double rLimit = getRelativeAccuracy() * (Math.abs(oldSum) + Math.abs(sum)) *
      // 0.5;
      // Avoid problems with very large sums by scaling each individually
      final double rLimit =
          (getRelativeAccuracy() * Math.abs(oldSum) + getRelativeAccuracy() * Math.abs(sum))
              * 0.5;
      if (delta <= rLimit) {
        break;
      }
    }

    // Remember to rescale to the correct range
    return dirac * sum - 1;
  }

  /**
   * Use an approximation as half the Poisson Fisher information when the mean is large.
   *
   * @param theta the mean
   * @return The approximate Fisher information
   */
  private static double largeApproximation(double theta) {
    return 1.0 / (2 * theta);
  }

  /**
   * Return the extreme limit for the given mean.
   *
   * <p>This can be called when an infinite sum has occurred on the unscaled A^2/P integral.
   *
   * @param theta the mean
   * @return the double
   */
  private static double extremeLimit(double theta) {
    // Infinite sum occurs when at the upper limit or lower limit
    return (theta > 1) ? largeApproximation(theta) : Double.POSITIVE_INFINITY;
  }

  private boolean computePmf(double dirac, double theta, double[] max, int endX) {
    double value;
    final double[] gradient = new double[1];

    // Initialise at x=0
    value = PoissonGammaFunction.unscaledPoissonGammaPartial(0, theta, gain, gradient);

    // At x=0 integrate from 0 to 0.5 without the Dirac contribution
    // (i.e. make it a a smooth function for integration).
    // Note the Dirac contribution for the unscaled function is 1.
    final double diracContribution = 1;
    value -= diracContribution;
    value = add(theta, value, gradient, 0, 0.25, 0.5);
    listP.elements()[0] += diracContribution;

    // The next x value after the max
    final int maxx = (int) Math.ceil(max[0]);
    // Threshold for switch to single point sampling
    final double threshold = max[1] * 0.5;

    // For the start of the range integrate from x-0.5 to x+0.5
    int x = 1;
    for (; x <= endX; x++) {
      value = add(theta, value, gradient, x - 0.5, x, x + 0.5);
      if (value == 0) {
        listP.removeDouble(x);
        listA.removeDouble(x);
        break;
      } else if (x >= maxx) {
        // Compute f
        final double f = getF(listP.getDouble(x), listA.getDouble(x));
        // If this is less than half the max then break and switch to a single
        // point evaluation.
        if (f < threshold) {
          break;
        }
      }
    }
    if (value != 0) {
      // For rest of the range just sample at point x
      for (x = listP.size(); x <= endX; x++) {
        value = PoissonGammaFunction.unscaledPoissonGammaPartial(x, theta, gain, gradient);
        if (value == 0) {
          break;
        }
        listP.add(value);
        listA.add(gradient[0]);
      }
    }

    // The PoissonGammaFunction can switch to an approximation
    // of the Bessel functions and the totals may not be correct.
    // Ensure the the sum of P and A are both 1 as this is required
    // for the integration to work.
    // This will be a factor when the mean is large. If the sums are nearly 1
    // then this will have little effect so do it anyway.

    // The sum of the unscaled P and A should be 1/Dirac
    final double expected = 1 / dirac;

    // At large mean the summation of the unscaled function underestimates.
    // This may be due to not using enough of the range. It may be due
    // to the Bessel approximation failing.

    final double sumP = listP.doubleStream().sum();
    if (sumP == Double.POSITIVE_INFINITY || sumP / expected < 0.75) {
      return false;
    }
    final double sumA = listA.doubleStream().sum();
    if (sumA == Double.POSITIVE_INFINITY || sumA / expected < 0.75) {
      return false;
    }

    // System.out.printf("theta=%g 0=%d Expected=%s SumP=%s (%s) SumA=%s (%s)\n", theta, endX,
    // expected, sumP, sumP / expected, sumA, sumA / expected);
    SimpleArrayUtils.apply(listP.elements(), 0, listP.size(), new ScalingFunction(expected / sumP));
    SimpleArrayUtils.apply(listA.elements(), 0, listA.size(), new ScalingFunction(expected / sumA));

    return true;
  }

  /**
   * Adds the integral of the Poisson-Gamma PDF to the working lists for P and A. It is assumed that
   * the input is the value and gradient of the lower bound. When finished return the value and
   * gradient of the upper bound.
   *
   * @param value the previous computed point value
   * @param gradient the previous computed point gradient
   * @param lower the lower bound
   * @param mid the mid point
   * @param upper the upper bound
   * @return the double
   */
  private double add(double theta, double value, double[] gradient, double lower, double mid,
      double upper) {
    // TODO - does this need refinement?

    // Do a single Simpson sum (f0 + 4 * fxi + fn) * deltaX / 3.
    double sumG = value;
    double sumA = gradient[0];
    value = PoissonGammaFunction.unscaledPoissonGammaPartial(mid, theta, gain, gradient);
    sumG += 4 * value;
    sumA += 4 * gradient[0];
    value = PoissonGammaFunction.unscaledPoissonGammaPartial(upper, theta, gain, gradient);
    final double deltaXover3 = (upper - lower) / 6;
    listP.add((sumG + value) * deltaXover3);
    listA.add((sumA + gradient[0]) * deltaXover3);
    return value;
  }

  private abstract static class IntegrationProcedure implements DoubleConvolutionValueProcedure {
    int index;

    protected double getF(double pz, double az) {
      // Compute with respect to the ultimate limit

      // p is always below 1.
      // - (az / pz) * az can tend towards infinity when pz is small
      // - (az * az) / pz can tend towards zero when az is small
      // Note: Computing infinity will result in an infinite fisher information
      // which is only possible for extremely low mean. Check for this.
      final double f = (az / pz) * az;
      if (f == Double.POSITIVE_INFINITY) {
        // System.out.printf("p=%s a=%s f=%s\n", pz, az, az * az / pz);
        return (az * az) / pz;
      }
      return f;
    }

    public abstract double getSum();
  }

  private static class SimpsonIntegrationProcedure extends IntegrationProcedure {
    double sum2;
    double sum4;

    @Override
    public boolean execute(double pz, double az) {
      ++index;
      if (pz > 0) {
        // Simpson's rule.
        // This computes the sum as:
        // h/3 * [ f(x0) + 4f(x1) + 2f(x2) + 4f(x3) + 2f(x4) ... + 4f(xn-1) + f(xn) ]
        if (index % 2 == 0) {
          sum2 += super.getF(pz, az);
        } else {
          sum4 += super.getF(pz, az);
        }
      }
      return true;
    }

    @Override
    public double getSum() {
      // Assume the end function values are zero
      // return (sum4 * 4 + sum2 * 2) / 3;
      // Stabilise for high sums by dividing first
      sum4 /= 3;
      sum2 /= 3;
      return sum4 * 4 + sum2 * 2;
    }
  }

  private static class Simpson38IntegrationProcedure extends IntegrationProcedure {
    double sum2;
    double sum3;

    @Override
    public boolean execute(double pz, double az) {
      ++index;
      if (pz > 0) {
        // Simpson's 3/8 rule based on cubic interpolation has a lower error.
        // This computes the sum as:
        // 3h/8 * [ f(x0) + 3f(x1) + 3f(x2) + 2f(x3) + 3f(x4) + 3f(x5) + 2f(x6) + ... + f(xn) ]
        if (index % 3 == 0) {
          sum2 += super.getF(pz, az);
        } else {
          sum3 += super.getF(pz, az);
        }
      }
      return true;
    }

    @Override
    public double getSum() {
      // Assume the end function values are zero
      // return (3.0 / 8) * (sum3 * 3 + sum2 * 2);
      // Stabilise for high sums by dividing first
      sum3 /= 8;
      sum2 /= 8;
      return sum3 * 9 + sum2 * 6;
    }
  }

  private static class UncheckedSimpsonIntegrationProcedure extends SimpsonIntegrationProcedure {
    @Override
    public boolean execute(double pz, double az) {
      // Simpson's rule.
      // This computes the sum as:
      // h/3 * [ f(x0) + 4f(x1) + 2f(x2) + 4f(x3) + 2f(x4) ... + 4f(xn-1) + f(xn) ]
      if (++index % 2 == 0) {
        sum2 += super.getF(pz, az);
      } else {
        sum4 += super.getF(pz, az);
      }
      return true;
    }
  }

  private static class UncheckedSimpson38IntegrationProcedure
      extends Simpson38IntegrationProcedure {
    @Override
    public boolean execute(double pz, double az) {
      // Simpson's 3/8 rule based on cubic interpolation has a lower error.
      // This computes the sum as:
      // 3h/8 * [ f(x0) + 3f(x1) + 3f(x2) + 2f(x3) + 3f(x4) + 3f(x5) + 2f(x6) + ... + f(xn) ]
      if (++index % 3 == 0) {
        sum2 += super.getF(pz, az);
      } else {
        sum3 += super.getF(pz, az);
      }
      return true;
    }
  }

  /**
   * Compute the integral.
   *
   * @param scale the scale of the Gaussian kernel
   * @param range the range of the Gaussian kernel
   * @param pd the Poisson distribution
   * @param a the a
   * @param minP the minimum P-value in the Poisson distribution
   * @return the integral
   * @throws IllegalArgumentException If the convolution will be too large
   */
  private double compute(int scale, int range, double[] pd, double[] a, double minP) {
    this.lastG = gaussianKernel.getGaussianKernel(scale, range, true);
    this.lastScale = scale;

    // If the product of the minimum kernel value and the min value is non-zero
    // then the convolution will never be zero.
    final boolean unchecked = lastG[0] * minP > 0;

    // If zeros can occur then the convolution could be reduced but the ends of
    // the distribution p that can be removed are unknown without convolution.
    // So the resulting convolution output is checked for zero.

    final IntegrationProcedure ip = createProcedure(unchecked);

    Convolution.convolve(lastG, pd, a, scale, ip);

    // Convolution should not create excess P
    // if (ip.sumP > 1.0001)
    // System.out.printf("theta=%g m=%g s=%g Bad convolution? = %s\n", lastT, m, s, ip.sumP);

    return ip.getSum();
  }

  private IntegrationProcedure createProcedure(final boolean unchecked) {
    if (use38) {
      return (unchecked) ? new UncheckedSimpson38IntegrationProcedure()
          : new Simpson38IntegrationProcedure();
    }
    return (unchecked) ? new UncheckedSimpsonIntegrationProcedure()
        : new SimpsonIntegrationProcedure();
  }

  private static double getF(double valueP, double gradientA) {
    if (valueP == 0) {
      return 0;
    }
    return MathUtils.pow2(gradientA) / valueP;
  }

  /**
   * Gets the fisher information function. This is the function that is integrated to get the Fisher
   * information. This can be called after a called to {@link #getFisherInformation(double)}. It
   * will return the last value of the function that was computed (before/after convolution with the
   * Gaussian kernel).
   *
   * <p>The function should be smooth with zero at the upper end. The lower end may not be smooth
   * due to the Dirac delta function at c=0.
   *
   * @param withGaussian Set to true to return the function after convolution with the gaussian.
   *        This re-performs the convolution.
   * @return the fisher information function
   * @throws IllegalStateException if no previous Fisher information computation is available
   */
  public double[][] getFisherInformationFunction(boolean withGaussian) {
    if (listP.isEmpty()) {
      throw new IllegalStateException("No previous Fisher information computation is available");
    }

    final double[] p = listP.toDoubleArray();
    final double[] a = listA.toDoubleArray();

    double[] gradientA;
    double[] valueP;
    double dx; // step size
    double offset;
    double scaleFactor;
    if (withGaussian && lastG != null) {
      // Do the convolution
      final double[][] result = Convolution.convolve(lastG, p, a, lastScale);
      valueP = result[0];
      gradientA = result[1];
      dx = 1.0 / lastScale;
      // Requires integer division. The kernel should be an odd size.
      final int radius = lastG.length / 2;
      offset = radius * -dx;
      scaleFactor = gaussianKernel.getConversionFactor(lastG);
    } else {
      valueP = p;
      gradientA = a;
      dx = 1;
      offset = 0;
      scaleFactor = 1;
    }

    // Correct the unscaled function to the scaled function
    final double dirac = StdMath.exp(-lastT);
    scaleFactor *= dirac;

    final double[] f = new double[gradientA.length];
    final double[] x = new double[gradientA.length];
    for (int i = 0; i < gradientA.length; i++) {
      x[i] = i * dx + offset;
      f[i] = getF(valueP[i], gradientA[i]) * scaleFactor;
    }
    return new double[][] {x, f};
  }

  /**
   * Gets the min range of the Gaussian kernel in SD units. The Gaussian kernel is scaled
   * dynamically depending on the shape of the Poisson distribution.
   *
   * @return the min range
   */
  public int getMinRange() {
    return minRange;
  }

  /**
   * Sets the min range of the Gaussian kernel in SD units. The Gaussian kernel is scaled
   * dynamically depending on the shape of the Poisson distribution.
   *
   * @param minRange the new min range
   */
  public void setMinRange(int minRange) {
    this.minRange = checkRange(minRange);
  }

  private static int checkRange(int range) {
    // Gaussian = Math.exp(-0.5 * x^2)
    // StdMath.exp(-746) == 0
    // => range for the Gaussian is sqrt(2*746) = 38.6
    return MathUtils.clip(1, MAX_RANGE, range);
  }

  /**
   * Gets the max range of the Gaussian kernel in SD units. The Gaussian kernel is scaled
   * dynamically depending on the shape of the Poisson distribution.
   *
   * @return the max range
   */
  public int getMaxRange() {
    return maxRange;
  }

  /**
   * Sets the max range of the Gaussian kernel in SD units. The Gaussian kernel is scaled
   * dynamically depending on the shape of the Poisson distribution.
   *
   * @param maxRange the new max range
   */
  public void setMaxRange(int maxRange) {
    this.maxRange = checkRange(maxRange);
  }

  /**
   * Gets the mean threshold for the switch to half the Poisson Fisher information.
   *
   * <p>When the mean is high then there is no Dirac delta contribution to the probability function.
   * This occurs when exp(-mean) is zero, which is approximately 746. In practice there is not much
   * difference when the mean is above 100.
   *
   * @return the mean threshold
   */
  public double getMeanThreshold() {
    return upperMeanThreshold;
  }

  /**
   * Sets the mean threshold for the switch to half the Poisson Fisher information.
   *
   * <p>When the mean is high then there is no Dirac delta contribution to the probability function.
   * This occurs when exp(-mean) is zero, which is approximately 746. In practice there is not much
   * difference when the mean is above 100.
   *
   * @param meanThreshold the new mean threshold
   */
  public void setMeanThreshold(double meanThreshold) {
    this.upperMeanThreshold = meanThreshold;
  }

  /**
   * Gets the relative probability threshold of the Poisson-Gamma distribution that is used. Any
   * probability less than this is ignored. This prevents using values that do not contribute to the
   * sum.
   *
   * @return the relative probability threshold
   */
  public double getRelativeProbabilityThreshold() {
    return relativeProbabilityThreshold;
  }

  /**
   * Sets the relative probability threshold of the Poisson-Gamma distribution that is used. Any
   * probability less than this is ignored. This prevents using values that do not contribute to the
   * sum.
   *
   * <p>This must be in the range 0 to 0.5. Larger values would prevent using any slope of the
   * Poisson-Gamma distribution. A value of zero requires all probability values to be used and may
   * result in long computation times.
   *
   * @param relativeProbabilityThreshold the new relative probability threshold
   */
  public void setRelativeProbabilityThreshold(double relativeProbabilityThreshold) {
    if (!(relativeProbabilityThreshold >= 0 && relativeProbabilityThreshold <= 0.5)) {
      throw new IllegalArgumentException("P must be in the range 0-0.5");
    }
    this.relativeProbabilityThreshold = relativeProbabilityThreshold;
  }

  /**
   * If true, use Simpson's 3/8 rule for cubic interpolation of the integral. False uses Simpson's
   * rule for quadratic interpolation.
   *
   * @return the use 38
   */
  public boolean getUse38() {
    return use38;
  }

  /**
   * Set to true to use Simpson's 3/8 rule for cubic interpolation of the integral. False uses
   * Simpson's rule for quadratic interpolation.
   *
   * @param use38 the new use 38
   */
  public void setUse38(boolean use38) {
    this.use38 = use38;
  }

  /**
   * Gets the relative accuracy for convergence during iteration.
   *
   * @return the relative accuracy
   */
  public double getRelativeAccuracy() {
    return relativeAccuracy;
  }

  /**
   * Sets the relative accuracy for convergence during iteration.
   *
   * <p>Set to below zero to prevent convergence check. This results in a fixed number of
   * iterations.
   *
   * @param relativeAccuracy the new relative accuracy
   */
  public void setRelativeAccuracy(double relativeAccuracy) {
    this.relativeAccuracy = relativeAccuracy;
  }

  /**
   * Gets the max iterations for iteration.
   *
   * @return the max iterations
   */
  public int getMaxIterations() {
    return maxIterations;
  }

  /**
   * Sets the max iterations for iteration. Set to zero to prevent iteration.
   *
   * @param maxIterations the new max iterations
   */
  public void setMaxIterations(int maxIterations) {
    this.maxIterations = maxIterations;
  }

  private static class QuickConvergenceChecker
      extends AbstractConvergenceChecker<UnivariatePointValuePair> {
    int fixedIterations;
    UnivariatePointValuePair previous;
    UnivariatePointValuePair current;

    public QuickConvergenceChecker(int fixedIterations) {
      super(0, 0);
      this.fixedIterations = fixedIterations;
    }

    @Override
    public boolean converged(int iteration, UnivariatePointValuePair previous,
        UnivariatePointValuePair current) {
      this.previous = previous;
      this.current = current;
      return iteration > fixedIterations;
    }

    UnivariatePointValuePair getBest() {
      if (previous == null || previous.getValue() < current.getValue()) {
        return current;
      }
      return previous;
    }
  }

  /**
   * Find the maximum of the integrated function A^2/P. P is the Poisson-Gamma convolution, A is the
   * partial gradient.
   *
   * <p>When both A and P are convolved with a Gaussian kernel, the integral of this function - 1 is
   * the Fisher information.
   *
   * <p>This method is used to determine the maximum of the function. The integral of A should equal
   * 1. The range can be determined using a fraction of the maximum, or integrating until the sum is
   * 1.
   *
   * @param theta the Poisson mean
   * @param rel Relative threshold
   * @return [max,max value,evaluations]
   */
  public double[] findMaximum(final double theta, double rel) {
    if (theta < MIN_MEAN) {
      throw new IllegalArgumentException("Poisson mean must be positive");
    }

    final double dirac = PoissonGammaFunction.dirac(theta);

    final UnivariateFunction f = new UnivariateFunction() {
      double[] gradient = new double[1];

      @Override
      public double value(double x) {
        // Return F without the dirac. This is so that the maximum of the function
        // is known for relative height analysis.
        double value;
        if (x == 0) {
          gradient[0] = dirac / gain;
          value = gradient[0] * theta;
        } else {
          value = PoissonGammaFunction.unscaledPoissonGammaPartial(x, theta, gain, gradient);
        }
        return getF(value, gradient[0]);
      }
    };

    // Determine the extent of the function: A^2/P
    // Without convolution this will have a maximum that is above, and
    // converges to [Poisson mean * amplification].
    final double mean = theta * gain;

    final int iterations = 10;
    final QuickConvergenceChecker checker = new QuickConvergenceChecker(iterations);
    final BrentOptimizer opt = new BrentOptimizer(rel, Double.MIN_VALUE, checker);
    UnivariatePointValuePair pair;
    try {
      final double start = (theta < 1) ? 0 : mean;
      pair = opt.optimize(new SearchInterval(0, mean * 1.5, start), GoalType.MAXIMIZE,
          new MaxEval(50000), new UnivariateObjectiveFunction(f));
    } catch (final TooManyEvaluationsException ex) {
      // This should not occur with the quick checker fixed iterations
      // - just get the current best
      pair = checker.getBest();
    }

    final double[] max = new double[3];
    max[0] = pair.getPoint();
    max[1] = pair.getValue();
    max[2] = opt.getEvaluations();
    return max;
  }

  /**
   * Find the upper limit of the integrated function A^2/P. P is the Poisson-Gamma convolution, A is
   * the partial gradient.
   *
   * <p>When both A and P are convolved with a Gaussian kernel, the integral of this function - 1 is
   * the Fisher information.
   *
   * <p>This method is used to determine the upper limit of the function using a binary search.
   *
   * @param theta the Poisson mean
   * @param max the max of the function (returned from {@link #findMaximum(double, double)})
   * @param rel Relative threshold
   * @return [upper,upper value,evaluations]
   */
  public double[] findUpperLimit(final double theta, double[] max, double rel) {
    if (rel < MIN_RELATIVE_TOLERANCE) {
      throw new IllegalArgumentException("Relative tolerance too small: " + rel);
    }

    final UnivariateFunction f = new UnivariateFunction() {
      double[] dgDp = new double[1];

      @Override
      public double value(double x) {
        final double G = PoissonGammaFunction.unscaledPoissonGammaPartial(x, theta, gain, dgDp);
        return getF(G, dgDp[0]);
      }
    };

    int eval = 0;

    // Increase from the max until the tolerance is achieved.

    // Use the mean to get a rough initial step size
    final double mean = theta * gain;
    double step = Math.max(mean, 1);

    double upper = max[0];
    double upperValue = max[1];
    final double threshold = upperValue * rel;

    double lower = upper;
    while (upperValue > threshold) {
      lower = upper;
      upper += step;
      step *= 2;
      eval++;
      upperValue = f.value(upper);
    }

    // Binary search the bracket between lower and upper
    while (lower + 1 < upper) {
      final double mid = (lower + upper) * 0.5;
      eval++;
      final double midg = f.value(mid);
      if (midg > threshold) {
        lower = mid;
      } else {
        upper = mid;
        upperValue = midg;
      }
    }

    return new double[] {upper, upperValue, eval};
  }
}
