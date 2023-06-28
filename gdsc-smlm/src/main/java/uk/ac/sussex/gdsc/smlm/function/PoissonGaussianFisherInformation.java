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
/*
 *
 */

package uk.ac.sussex.gdsc.smlm.function;

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import java.util.Arrays;
import uk.ac.sussex.gdsc.core.math.NumberUtils;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.smlm.math3.distribution.PoissonDistribution;
import uk.ac.sussex.gdsc.smlm.utils.Convolution;
import uk.ac.sussex.gdsc.smlm.utils.Convolution.ConvolutionValueProcedure;
import uk.ac.sussex.gdsc.smlm.utils.GaussianKernel;

/**
 * Calculate the Fisher information for a Poisson-Gaussian distribution.
 *
 * <p>Uses the equation of Chao, et al (2013) Nature Methods, 10, 335-338, SI Eq S7.
 *
 * <p>Performs a convolution with a finite Gaussian kernel.
 *
 * <p>An optimisation is used when the mean of the Poisson is above a threshold. In this case the
 * Poisson can be approximated as a Gaussian and the Fisher information is returned for the
 * Gaussian-Gaussian convolution.
 */
public class PoissonGaussianFisherInformation extends BasePoissonFisherInformation {
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

  /** The default cumulative probability for the Poisson distribution. */
  public static final double DEFAULT_CUMULATIVE_PROBABILITY = 1 - 1e-10;

  /**
   * The default sampling of the Gaussian kernel. The kernel will be sampled at s/sampling, i.e.
   * this is the number of samples to take per standard deviation unit.
   */
  public static final int DEFAULT_SAMPLING = 4;

  /** The default relative accuracy for convergence. */
  public static final double DEFAULT_RELATIVE_ACCURACY = 1e-6;

  /**
   * The default number of maximum iterations. The Gaussian sampling will be doubled during each
   * iteration.
   */
  public static final int DEFAULT_MAX_ITERATIONS = 4;

  /**
   * Store the limit of the Poisson distribution for small mean for the default cumulative
   * probability.
   */
  private static final int[] DEFAULT_LIMITS;
  /**
   * Store the limit of the Poisson distribution for tiny mean for the default cumulative
   * probability.
   */
  private static final int[] DEFAULT_TINY_LIMITS;

  static {
    final PoissonDistribution pd = new PoissonDistribution(1);
    DEFAULT_LIMITS = new int[101];
    for (int i = 1; i < DEFAULT_LIMITS.length; i++) {
      DEFAULT_LIMITS[i] = computeLimit(pd, i, DEFAULT_CUMULATIVE_PROBABILITY);
      // System.out.printf("[%d] = %d scale=%d\n", i, defaultLimits[i], getScale(Math.sqrt(i)));
    }

    // Use exponent of the mean down to -20
    DEFAULT_TINY_LIMITS = new int[21];
    for (int i = 1; i < DEFAULT_TINY_LIMITS.length; i++) {
      DEFAULT_TINY_LIMITS[i] = computeTinyLimit(pd, -i, DEFAULT_CUMULATIVE_PROBABILITY);
      // System.out.printf("[%d] = %d : p(0) = %g : cumul = %s : next = %g\n", i,
      // defaultTinyLimits[i],
      // pd.probability(0), pd.cumulativeProbability(defaultTinyLimits[i]),
      // pd.probability(defaultTinyLimits[i] + 1));
    }
  }

  /** The standard deviation of the Gaussian. */
  public final double sd;

  /** The minimum range of the Gaussian kernel (in SD units). */
  private int minRange = DEFAULT_MIN_RANGE;

  /** The maximum range of the Gaussian kernel (in SD units). */
  private int maxRange = DEFAULT_MAX_RANGE;

  /** The scale of the kernel. This is a power of 2. */
  private final int defaultScale;

  /** The poisson distribution used to generate the Poisson probabilities. */
  private final PoissonDistribution pd = new PoissonDistribution(1);

  /** Working space to store the Poisson probabilities. */
  private final DoubleArrayList list = new DoubleArrayList();

  /** The mean threshold for the switch to a Gaussian-Gaussian convolution. */
  private double meanThreshold = 100;

  /** The cumulative probability of the Poisson distribution that is used. */
  private double cumulativeProbability = DEFAULT_CUMULATIVE_PROBABILITY;

  /** Set to true to use Simpson's 3/8 rule for cubic interpolation of the integral. */
  private boolean use38 = true;

  /** Store the limit of the Poisson distribution for small mean for the cumulative probability. */
  private int[] limits = DEFAULT_LIMITS;

  /** Store the limit of the Poisson distribution for tiny mean for the cumulative probability. */
  private int[] tinyLimits = DEFAULT_TINY_LIMITS;

  /** Flag to indicate no possible convolution with the gaussian. */
  private final boolean noGaussian;

  /** The gaussian kernel. */
  private GaussianKernel gaussianKernel;

  /** The relative accuracy for convergence. */
  private double relativeAccuracy = DEFAULT_RELATIVE_ACCURACY;

  /** The max iterations. */
  private int maxIterations = DEFAULT_MAX_ITERATIONS;

  /**
   * Instantiates a new poisson gaussian fisher information.
   *
   * @param sd the standard deviation of the Gaussian
   * @throws IllegalArgumentException If the standard deviation is not strictly positive
   */
  public PoissonGaussianFisherInformation(double sd) {
    this(sd, DEFAULT_SAMPLING);
  }

  /**
   * Instantiates a new poisson gaussian fisher information.
   *
   * @param sd the standard deviation of the Gaussian
   * @param sampling The number of Gaussian samples to take per standard deviation
   * @throws IllegalArgumentException If the standard deviation is not strictly positive
   * @throws IllegalArgumentException If the sampling is below 1
   * @throws IllegalArgumentException If the maximum kernel size after scaling is too large
   */
  public PoissonGaussianFisherInformation(double sd, double sampling) {
    if (!(sd >= 0 && sd <= Double.MAX_VALUE)) {
      throw new IllegalArgumentException("Gaussian standard deviation must be positive");
    }
    if (!(sampling >= 1 && sampling <= Double.MAX_VALUE)) {
      throw new IllegalArgumentException("Gaussian sampling must at least 1");
    }

    this.sd = sd;
    noGaussian = (sd * MAX_RANGE < 1);
    if (noGaussian) {
      // This is OK. Just return the information for a Poisson.
      this.defaultScale = 0;
    } else {
      // This is set to work for reasonable values of the Gaussian kernel and sampling
      // e.g. s [0.5:20], sampling from [1:8].

      this.defaultScale = getPow2Scale(sampling / sd);

      // Don't support excess scaling caused by small kernels
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
  protected PoissonGaussianFisherInformation(PoissonGaussianFisherInformation source) {
    sd = source.sd;
    minRange = source.minRange;
    maxRange = source.maxRange;
    defaultScale = source.defaultScale;
    meanThreshold = source.meanThreshold;
    cumulativeProbability = source.cumulativeProbability;
    use38 = source.use38;
    limits = source.limits;
    tinyLimits = source.tinyLimits;
    noGaussian = source.noGaussian;
    gaussianKernel = source.gaussianKernel.copy();
    relativeAccuracy = source.relativeAccuracy;
    maxIterations = source.maxIterations;
  }

  @Override
  public PoissonGaussianFisherInformation copy() {
    return new PoissonGaussianFisherInformation(this);
  }

  /**
   * Compute the limit of a usable probability above 0.
   *
   * <p>Find the point where probability will return above 0. Using StdMath.exp this is -746.
   * However sub-normal output occurs at -709. This is a good limit for computation.
   *
   * @param mean the mean
   * @return the limit
   */
  public static int computeLimit(double mean) {
    final PoissonDistribution pd = new PoissonDistribution(mean);
    int x = (int) mean;
    while (pd.logProbability(x + 1) > -709) {
      x++;
    }
    return x;
  }

  /**
   * Compute the limit of a usable probability above 0.
   *
   * @param pd the pd
   * @param mean the mean
   * @param cumulativeProbability the cumulative probability
   * @return the limit
   */
  private static int computeLimit(PoissonDistribution pd, double mean,
      double cumulativeProbability) {
    pd.setMeanUnsafe(mean);
    return pd.inverseCumulativeProbability(cumulativeProbability);
  }

  /**
   * Compute the limit of a usable probability above 0.
   *
   * @param pd the pd
   * @param exp the exponent of the mean (in base 2)
   * @param cumulativeProbability the cumulative probability
   * @return the limit
   */
  private static int computeTinyLimit(PoissonDistribution pd, int exp,
      double cumulativeProbability) {
    // Fill all bits of the mantissa
    final long bits = 0xffffffffffffffL;
    final double mean = Double.longBitsToDouble(bits | (long) (exp + 1023) << 52);
    pd.setMeanUnsafe(mean);
    return pd.inverseCumulativeProbability(cumulativeProbability);
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
   * 10, 335-338, SI Eq S7. Note that that equation computes the noise coefficient relative to a
   * Poisson, this computes the Fisher information. To get the noise coefficient multiply by the
   * input parameter.
   */
  @Override
  public double getFisherInformation(double theta) {
    final double I = getPoissonGaussianI(theta);

    // Check limits.

    // It should be worse than the Poisson Fisher information (upper limit) but
    // better than the Poisson-Gaussian Fisher information (lower limit).
    // Note a low Fisher information is worse as this is the amount of information
    // carried about the parameter.
    final double lower = 1.0 / (theta + sd * sd); // getPoissonGaussianApproximationI(theta)
    final double upper = 1.0 / theta; // PoissonFisherInformation.getPoissonI(theta)
    return MathUtils.clip(lower, upper, I);
  }

  @Override
  public double getAlpha(double theta) {
    // Simple implementation
    return (noGaussian) ? 1 : theta * getFisherInformation(theta);
  }

  /**
   * Gets the Poisson-Gaussian Fisher information.
   *
   * <p>The input parameter refers to the mean of the Poisson distribution.
   *
   * <p>The Fisher information is computed using the equation of Chao, et al (2013) Nature Methods,
   * 10, 335-338, SI Eq S7. Note that that equation computes the noise coefficient relative to a
   * Poisson, this computes the Fisher information. To get the noise coefficient multiply by the
   * input parameter.
   *
   * <p>Note: This uses a convolution of an infinite integral over a finite range. It may
   * under-estimate the information when the mean is large. Use
   * {@link #getFisherInformation(double)} for a checked return value, clipped to the expected range
   * for a Poisson and the Poisson-Gaussian approximation.
   *
   * @param theta the Poisson mean
   * @return the Poisson Gaussian Fisher information
   * @throws IllegalArgumentException the illegal argument exception
   */
  public double getPoissonGaussianI(double theta) {
    if (theta <= 0) {
      throw new IllegalArgumentException("Poisson mean must be positive");
    }
    if (theta < MIN_MEAN) {
      return Double.POSITIVE_INFINITY;
    }

    if (noGaussian) {
      // No Gaussian convolution
      // Get the Fisher information for a Poisson.
      return 1.0 / theta;
    }

    if (theta > meanThreshold) {
      // Use an approximation when the Poisson mean is large
      return 1.0 / (theta + sd * sd); // getPoissonGaussianApproximationI(theta);
    }

    // This computes the convolution of a Poisson PMF and a Gaussian PDF.
    // The value of this is p(z).

    // The Poisson-Gaussian must be differentiated to compute the Fisher information:
    // Expected [ (d ln(p(z)) dv)^2 ]
    // = integral [ (1/p(z) . d p(z) dv)^2 p(z) dz ]
    // = integral [ 1/p(z) . (d p(z) dv)^2 dz ]

    // Gaussian standard deviation = s

    // Chao et al, S5:
    // p(z) = 1/sqrt(2pi)s sum_j=0:Inf e^-v . v^j / j! . e^-1/2((z-j)/s)^2

    // This is the sum over j of the probability of Poisson(j) * probability of Gaussian(z-j)

    // Note: (fg)' => f'g + fg'
    // e^-v => -e^-v
    // v^j => j.v^(j-1)
    // e^-v v^j / j! => e^-v v^(j-1) / (j-1)! - e^-v v^j / j!

    // d p(z) dv = 1/sqrt(2pi)s sum_j=1:Inf e^-v . v^(j-1) / (j-1)! . e^-1/2((z-j)/s)^2 -
    // sum_j=0:Inf e^-v . v^j / j! . e^-1/2((z-j)/s)^2
    // Note: j=0 differentiates to -e^v since v^j / j! = 1. This removes j=0 from the first sum
    // but not the second.
    // Over the sum the second term adds up to become p(z) so:
    // d p(z) dv = (1/sqrt(2pi)s sum_j=1:Inf e^-v . v^(j-1) / (j-1)! . e^-1/2((z-j)/s)^2 ) - p(z)

    // Set the first term to A, the second to P:
    // d p(z) dv = A - P

    // E = integral [ (1/p(z) . d p(z) dv)^2 p(z) dz ]
    // = integral [ (1/P . (A - P))^2 * P ]
    // = integral [ (1/P^2 . (A^2 - 2AP + p^2) * P ]
    // = integral [ (A^2/P^2 - 2A/P + 1) * P ]
    // = integral [ A^2/P - 2A + P ]
    // = integral [A^2/P] - integral [2A] + integral [P]

    // Note that the integral of P==1.
    // Since the integral of A is just P offset by j-1, integral of A==1

    // E = integral [A^2/P] - 1

    // P(z) = 1/sqrt(2pi)s sum_j=0:Inf e^-v . v^j / j! . e^-1/2((z-j)/s)^2
    // A(z) = 1/sqrt(2pi)s sum_j=1:Inf e^-v . v^(j-1) / (j-1)! . e^-1/2((z-j)/s)^2
    // A(z) = 1/sqrt(2pi)s sum_j=0:Inf e^-v . v^j / j! . e^-1/2((z-(j+1))/s)^2
    // A(z) = P(z-1)

    // We need the convolution of the Poisson with the Gaussian

    // The Poisson is a PMF. The Gaussian is a PDF.
    // The expected value is integrated over all real z, -Inf:Inf,
    // (the full range of the Gaussian).

    // Sample the values of the full range and compute a sum using Simpson integration.

    // Only use part of the cumulative distribution.
    // Start at zero for a 1-tailed truncation of the cumulative distribution.
    // This code will be for low mean values where the contribution at zero is
    // significant.
    int minx = 0;
    // Find the limit. These can be cached (or may be the defaults).
    // TODO - Determine if the Poisson can be truncated. We may have to use more of
    // the values (for example those returned by computeLimit(...).
    int maxx;
    // The exponent provides a rough idea of the size of the mean
    final int exp = NumberUtils.getSignedExponent(theta);
    if (theta < 1) {
      int ex = -exp;
      if (ex >= tinyLimits.length) {
        ex = tinyLimits.length - 1;
      }
      if (tinyLimits[ex] == 0) {
        tinyLimits[ex] = computeTinyLimit(pd, -ex, cumulativeProbability);
      }
      maxx = tinyLimits[ex];
    } else {
      final int x = (int) Math.ceil(theta);
      if (x < limits.length) {
        if (limits[x] == 0) {
          limits[x] = computeLimit(pd, x, cumulativeProbability);
        }
        maxx = limits[x];
      } else {
        // For large mean the distribution will be far from zero.
        // In this case use a 2-tailed limit.
        final double lower = (1 - cumulativeProbability) / 2;
        minx = computeLimit(pd, x, lower);
        maxx = computeLimit(pd, x, 1 - lower);
      }
    }

    // Build the Poisson distribution.
    pd.setMeanUnsafe(theta);
    list.clear();

    // XXX - check this is needed
    // For small theta the tail of the distribution is important
    if (maxx < 10) {
      maxx = 10;
    }

    for (int x = minx; x <= maxx; x++) {
      final double pp = pd.probability(x);
      if (pp == 0) {
        break;
      }
      list.add(pp);
    }

    if (list.size() < 2) {
      // Extreme case where there is no Poisson for convolution.
      // Assume a Gaussian distribution. Return the Fisher information
      // for the Gaussian with mean 0. This will happen when the cumulative
      // probability has been altered from the default.
      return getGaussianI();
    }

    // Unscaled Poisson
    final double[] p = list.toDoubleArray();

    // Convolve with the Gaussian kernel.
    // As the mean reduces the Poisson distribution is more skewed
    // and the extent of the kernel must change. Just increase the range
    // for the kernel for each power of 2 the number is below 1.
    int range = minRange;
    for (int e = exp; range < maxRange && e <= 0; e++, range++) {
      // The counters are incremented in the for loop increment clause
    }
    // Ensure the kernel range covers multiple values of the Poisson distribution.
    // Only applicable to small kernels
    while (range < maxRange && range * sd < 1) {
      range++;
    }

    // In order for A(z) = P(z-1) to work sum A(z) must be 1
    double sum;
    // sum = list.sum();
    // System.out.printf("Normalisation (theta=%g) = %s\n", theta, sum);
    // for (int i = 0; i < p.length; i++)
    // p[i] /= sum;

    // Initial sum
    int scale = defaultScale;
    sum = compute(scale, range, p);

    // Iterate
    for (int iteration = 1; iteration <= maxIterations && scale < MAX_SCALE; iteration++) {
      scale *= 2;
      final double oldSum = sum;
      try {
        sum = compute(scale, range, p);
      } catch (final IllegalArgumentException ex) {
        // Occurs when the convolution has grown too big
        return sum;
      }
      final double delta = Math.abs(sum - oldSum);
      // System.out.printf("s=%g theta=%g Iteration=%d sum=%s oldSum=%s change=%s\n", s, theta,
      // iteration, sum, oldSum,
      // delta / (Math.abs(oldSum) + Math.abs(sum)) * 0.5);
      final double rLimit = getRelativeAccuracy() * (Math.abs(oldSum) + Math.abs(sum)) * 0.5;
      if (delta <= rLimit) {
        break;
      }
    }

    return sum;
  }

  /**
   * Define an integration procedure for a convolution result.
   */
  private abstract static class IntegrationProcedure implements ConvolutionValueProcedure {
    /**
     * The scale mask to compute a modulus using the counter where the scale is assumed to be a
     * power of 2.
     */
    private final int scaleMask;
    final double[] pz1;
    int counter;

    /**
     * Create an instance.
     *
     * @param scale the scale (assumed to be a power of 2)
     */
    IntegrationProcedure(int scale) {
      scaleMask = scale - 1;

      // E = integral [A^2/P] - 1
      // P(z) = 1/sqrt(2pi)s sum_j=0:Inf e^-v . v^j / j! . e^-1/2((z-j)/s)^2
      // A(z) = 1/sqrt(2pi)s sum_j=1:Inf e^-v . v^(j-1) / (j-1)! . e^-1/2((z-j)/s)^2
      // A(z) = 1/sqrt(2pi)s sum_j=0:Inf e^-v . v^j / j! . e^-1/2((z-(j+1))/s)^2
      // A(z) = P(z-1)

      // Store p(z-1). This is initialised to 0.
      pz1 = new double[scale];
    }

    @Override
    public boolean execute(double pz) {
      final int index = counter & scaleMask;
      counter++;
      final double az = pz1[index];
      pz1[index] = pz;
      if (pz > 0) {
        // Compute with respect to the ultimate limit.
        // Both az and pz should be < 1
        // if az < 1 : az^2 -> 0
        // if pz < 1 then dividing first will reduce the chance of computing zero.
        sum((az / pz) * az);
      }
      return true;
    }

    /**
     * Sum the value.
     *
     * @param value the value
     */
    abstract void sum(double value);

    /**
     * Gets the sum.
     *
     * @return the sum
     */
    abstract double getSum();
  }

  /**
   * Perform Simpson 1/3 integration.
   *
   * @see <a href="https://en.wikipedia.org/wiki/Simpson%27s_rule">Simpsons rule</a>
   */
  private static class SimpsonIntegrationProcedure extends IntegrationProcedure {
    double sum2;
    double sum4;

    SimpsonIntegrationProcedure(int scale) {
      super(scale);
    }

    @Override
    void sum(double value) {
      // Simpson's rule.
      // This computes the sum as:
      // h/3 * [ f(x0) + 4f(x1) + 2f(x2) + 4f(x3) + 2f(x4) ... + 4f(xn-1) + f(xn) ]
      if ((counter & 0x1) == 0) {
        sum2 += value;
      } else {
        sum4 += value;
      }
    }

    @Override
    double getSum() {
      // Assume the end function values are zero
      // return (sum4 * 4 + sum2 * 2) / 3;
      // Stabilise for high sums by dividing first
      sum4 /= 3;
      sum2 /= 3;
      return sum4 * 4 + sum2 * 2;
    }
  }

  /**
   * Perform Simpson 3/8 integration.
   *
   * @see <a href="https://en.wikipedia.org/wiki/Simpson%27s_rule">Simpsons rule</a>
   */
  private static class Simpson38IntegrationProcedure extends IntegrationProcedure {
    double sum2;
    double sum3;

    Simpson38IntegrationProcedure(int scale) {
      super(scale);
    }

    @Override
    void sum(double value) {
      // Simpson's 3/8 rule based on cubic interpolation has a lower error.
      // This computes the sum as:
      // 3h/8 * [ f(x0) + 3f(x1) + 3f(x2) + 2f(x3) + 3f(x4) + 3f(x5) + 2f(x6) + ... + f(xn) ]
      if (counter % 3 == 0) {
        sum2 += value;
      } else {
        sum3 += value;
      }
    }

    @Override
    double getSum() {
      // Assume the end function values are zero
      // return (3.0 / 8) * (sum3 * 3 + sum2 * 2);
      // Stabilise for high sums by dividing first
      sum3 /= 8;
      sum2 /= 8;
      return sum3 * 9 + sum2 * 6;
    }
  }

  /**
   * Compute the integral.
   *
   * @param scale the scale of the Gaussian kernel
   * @param range the range of the Gaussian kernel
   * @param pd the poisson distribution
   * @return the integral
   * @throws IllegalArgumentException If the convolution will be too large
   */
  private double compute(int scale, int range, double[] pd) {
    final double[] g = gaussianKernel.getGaussianKernel(scale, range, true);

    final IntegrationProcedure ip =
        (use38) ? new Simpson38IntegrationProcedure(scale) : new SimpsonIntegrationProcedure(scale);

    Convolution.convolve(g, pd, scale, ip);

    // Subtract the final 1
    return ip.getSum() - 1;
  }

  /**
   * Gets the approximate Poisson-Gaussian Fisher information.
   *
   * <p>Approximate the Poisson as a Gaussian with {@code u=theta} and {@code var=theta}.
   * Gaussian-Gaussian convolution: s<sub>a</sub> * s<sub>b</sub> =&gt; s<sub>c</sub> =
   * sqrt(s<sub>a</sub><sup>2</sup>+s<sub>b</sub><sup>2</sup>). Fisher information of Gaussian mean
   * is 1/variance.
   *
   * <p>The returned value is: {@code 1.0 / (theta + s * s)} with {@code theta} the Poisson mean and
   * {@code s} the Gaussian standard deviation.
   *
   * @param theta the poisson mean
   * @return the Poisson Gaussian Fisher information
   */
  public double getPoissonGaussianApproximationI(double theta) {
    if (theta <= 0) {
      throw new IllegalArgumentException("Poisson mean must be positive");
    }
    return 1.0 / (theta + sd * sd);
  }

  /**
   * Gets the approximate Poisson-Gaussian Fisher information.
   *
   * <p>Approximate the Poisson as a Gaussian with {@code u=theta} and {@code var=theta}.
   * Gaussian-Gaussian convolution: s<sub>a</sub> * s<sub>b</sub> =&gt; s<sub>c</sub> =
   * sqrt(s<sub>a</sub><sup>2</sup>+s<sub>b</sub><sup>2</sup>). Fisher information of Gaussian mean
   * is 1/variance.
   *
   * <p>The returned value is: {@code 1.0 / (theta + s * s)} with {@code theta} the Poisson mean and
   * {@code s} the Gaussian standard deviation.
   *
   * @param theta the poisson mean
   * @param sd the Gaussian standard deviation (no check made for negatives as this is squared)
   * @return the Poisson Gaussian Fisher information
   */
  public static double getPoissonGaussianApproximationI(double theta, double sd) {
    if (theta <= 0) {
      throw new IllegalArgumentException("Poisson mean must be positive");
    }
    return 1.0 / (theta + sd * sd);
  }

  /**
   * Gets the Gaussian Fisher information for mean 0. Fisher information of Gaussian mean is
   * 1/variance.
   *
   * @return the Gaussian Fisher information
   */
  public double getGaussianI() {
    return 1.0 / (sd * sd);
  }

  /**
   * Gets the Gaussian Fisher information for mean 0. Fisher information of Gaussian mean is
   * 1/variance.
   *
   * @param sd the Gaussian standard deviation (no check made for negatives as this is squared)
   * @return the Gaussian Fisher information
   */
  public static double getGaussianI(double sd) {
    return 1.0 / (sd * sd);
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
   * Gets the mean threshold for the switch to a Gaussian-Gaussian convolution.
   *
   * @return the mean threshold
   */
  public double getMeanThreshold() {
    return meanThreshold;
  }

  /**
   * Sets the mean threshold for the switch to a Gaussian-Gaussian convolution.
   *
   * @param meanThreshold the new mean threshold
   */
  public void setMeanThreshold(double meanThreshold) {
    this.meanThreshold = meanThreshold;
  }

  /**
   * Gets the cumulative probability of the Poisson distribution that is used.
   *
   * @return the cumulative probability
   */
  public double getCumulativeProbability() {
    return cumulativeProbability;
  }

  /**
   * Sets the cumulative probability of the Poisson distribution that is used.
   *
   * @param cumulativeProbability the new cumulative probability
   */
  public void setCumulativeProbability(double cumulativeProbability) {
    if (!(cumulativeProbability > 0 && cumulativeProbability <= 1)) {
      throw new IllegalArgumentException("P must be in the range 0-1");
    }
    if (this.cumulativeProbability != cumulativeProbability) {
      this.cumulativeProbability = cumulativeProbability;
      if (limits == DEFAULT_LIMITS) {
        limits = new int[DEFAULT_LIMITS.length];
        tinyLimits = new int[DEFAULT_TINY_LIMITS.length];
      } else {
        // Reset
        Arrays.fill(limits, 0);
        Arrays.fill(tinyLimits, 0);
      }
    }
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
}
