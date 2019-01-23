/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.utils;

import uk.ac.sussex.gdsc.core.utils.MathUtils;

import gnu.trove.list.array.TDoubleArrayList;

import org.apache.commons.math3.util.FastMath;

import java.util.Arrays;

/**
 * Store a Gaussian kernel for use in convolution.
 */
public class GaussianKernel {
  /** The maximum size of half the Gaussian kernel. */
  public static final int HALF_WIDTH_LIMIT = 1 << 28;

  private static final double ONE_OVER_ROOT_2_PI = 1.0 / Math.sqrt(2 * Math.PI);

  /** The standard deviation. */
  public final double sd;

  private final double var2;

  private int currentScale;

  private TDoubleArrayList halfKernel;

  /**
   * Instantiates a new gaussian kernel.
   *
   * @param sd The standard deviation
   */
  public GaussianKernel(double sd) {
    // Check not infinite or NaN
    if (!(sd > 0 && sd < Double.MAX_VALUE)) {
      throw new IllegalArgumentException("Standard deviation must be positive");
    }
    this.sd = sd;
    currentScale = 1;
    halfKernel = new TDoubleArrayList();
    // Initialise the first value exp(-0)
    halfKernel.add(1);
    // Precompute exponential scaling factor
    var2 = -(sd * sd * 2);
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected GaussianKernel(GaussianKernel source) {
    this.sd = source.sd;
    this.var2 = source.var2;
    currentScale = source.currentScale;
    halfKernel = new TDoubleArrayList(source.halfKernel);
  }

  /**
   * Get a copy.
   *
   * @return the copy
   */
  public GaussianKernel copy() {
    return new GaussianKernel(this);
  }

  /**
   * Gets the normalisation for the Gaussian:
   *
   * <pre>
   * 1 / sqrt(2 * pi * s ^ 2)
   * </pre>
   *
   * <p>This is the normalisation component that should be applied to the Gaussian exponential
   * component:
   *
   * <pre>
   *  exp^(-x^2 / 2s^2).
   * </pre>
   *
   * <p>Note that the kernel computed is normalised to 1 irrespective of the number of samples. This
   * is so that convolution using the kernel does not alter the total signal in the input function.
   * However the raw Gaussian can be deduced using the difference in the magnitude at the centre of
   * the kernel and this normalisation factor.
   *
   * @return the normalisation
   */
  public double getNormalisation() {
    return ONE_OVER_ROOT_2_PI / sd;
  }

  /**
   * Gets the conversion factor. Thus is the factor that should be applied to the kernel values to
   * convert them to the true Gaussian value:
   *
   * <pre>
   * 1/sqrt(2*pi*s^2) * exp^(-x^2 / 2s^2)
   * </pre>
   *
   * @param kernel the kernel
   * @return the conversion factor
   */
  public double getConversionFactor(double[] kernel) {
    return getNormalisation() / kernel[kernel.length / 2];
  }

  /**
   * Create a 1-dimensional normalized Gaussian kernel with standard deviation s*scale. To avoid a
   * step due to the cutoff at a finite value, the near-edge values are replaced by a 2nd-order
   * polynomial with its minimum=0 at the first out-of-kernel pixel. Thus, the kernel function has a
   * smooth 1st derivative in spite of finite length.
   *
   * @param scale the scale (must be a power of 2)
   * @param range the range (in units of standard deviation)
   * @param edgeCorrection Set to true to perform the edge correction
   * @return The kernel, decaying towards zero, which would be reached at the first out of kernel
   *         index
   * @throws IllegalArgumentException If the input scale is not a power of 2
   */
  public double[] getGaussianKernel(int scale, double range, boolean edgeCorrection) {
    if (!MathUtils.isPow2(scale)) {
      throw new IllegalArgumentException("Scale must be a power of 2: " + scale);
    }

    // Limit range for the Gaussian
    if (range < 1) {
      range = 1;
    } else if (range > 38) {
      range = 38;
    }

    final int kradius = getGaussianHalfWidth(sd * scale, range) + 1;
    increaseScale(scale);
    increaseKernel(kradius);

    // Create kernel
    // Note: The stored values in the halfKernel are always non-zero.
    final double[] kernel = new double[2 * kradius - 1];
    kernel[0] = 1;
    if (currentScale == scale) {
      for (int i = 1, size = Math.min(kradius, halfKernel.size()); i < size; i++) {
        kernel[i] = halfKernel.getQuick(i);
      }
    } else {
      final double step = 1.0 / scale;
      final int sample = currentScale / scale;
      for (int i = 1, j = sample; i < kradius; i++, j += sample) {
        // In case sampling requires a different end point in the kernel
        // check the size
        if (j < halfKernel.size()) {
          kernel[i] = halfKernel.getQuick(j);
        } else {
          kernel[i] = FastMath.exp(MathUtils.pow2(i * step) / var2);
          // Check if zero
          if (kernel[i] == 0) {
            break;
          }
        }
      }
    }

    return buildKernel(kernel, kradius, edgeCorrection);
  }

  /**
   * Create a 1-dimensional normalized Gaussian kernel with standard deviation s/scale. To avoid a
   * step due to the cutoff at a finite value, the near-edge values are replaced by a 2nd-order
   * polynomial with its minimum=0 at the first out-of-kernel pixel. Thus, the kernel function has a
   * smooth 1st derivative in spite of finite length.
   *
   * <p>Note: This can lead to inefficiency if the kernel has been upscaled (i.e. the current scale
   * is above 1) and the cached range is below the input range since the current kernel must be
   * expanded. It is recommended to only call this method with the same range that has been used for
   * {@link #getGaussianKernel(int, double, boolean)}. The result is just a sparsely sampled kernel.
   *
   * <p>This method is not recommended when the scale is above the standard deviation as the kernel
   * will be too sparsely sampled.
   *
   * @param scale the scale
   * @param range the range (in units of standard deviation)
   * @param edgeCorrection Set to true to perform the edge correction
   * @return The kernel, decaying towards zero, which would be reached at the first out of kernel
   *         index
   * @throws IllegalArgumentException If the input scale is above zero
   */
  public double[] getDownscaleGaussianKernel(int scale, double range, boolean edgeCorrection) {
    if (scale < 1) {
      throw new IllegalArgumentException("Scale must be strictly positive: " + scale);
    }

    // Limit range for the Gaussian
    if (range < 1) {
      range = 1;
    } else if (range > 38) {
      range = 38;
    }

    // Check if the current half-kernel would be too large if expanded to the range
    if ((double) currentScale * scale * range > HALF_WIDTH_LIMIT) {
      return makeErfGaussianKernel(sd / scale, range);
    }

    final int sample = currentScale * scale;

    // Expand the kernel to cover the range at the current scale.
    // Note: This can lead to inefficiency if the kernel has been upscaled
    // (i.e. the current scale is above 1).
    int kradius = getGaussianHalfWidth(sample, range) + 1;
    increaseKernel(kradius);

    // Now get the radius of the downscaled kernel
    kradius = getGaussianHalfWidth(sd / scale, range) + 1;

    // Create kernel
    // Note: The stored values in the halfKernel are always non-zero.
    final double[] kernel = new double[2 * kradius - 1];
    kernel[0] = 1;
    final double step = scale;
    for (int i = 1, j = sample; i < kradius; i++, j += sample) {
      // In case sampling requires a different end point in the kernel
      // check the size
      if (j < halfKernel.size()) {
        kernel[i] = halfKernel.getQuick(j);
      } else {
        kernel[i] = FastMath.exp(MathUtils.pow2(i * step) / var2);
        // Check if zero
        if (kernel[i] == 0) {
          break;
        }
      }
    }

    return buildKernel(kernel, kradius, edgeCorrection);
  }

  /**
   * Get half the width of the region smoothed by a Gaussian filter for the specified standard
   * deviation. The full region size is 2N + 1
   *
   * @param sigma the sigma
   * @param range the range
   * @return The half width
   */
  public static int getGaussianHalfWidth(double sigma, double range) {
    final double limit = Math.ceil(sigma * range);
    // Ensure the kernel is clipped to the size of an array
    return (limit < HALF_WIDTH_LIMIT) ? (int) limit : HALF_WIDTH_LIMIT;
  }

  private void increaseScale(int scale) {
    if (currentScale < scale) {
      // Up sample the current kernel
      final int upsample = scale / currentScale;
      currentScale = scale;

      final double[] g = halfKernel.toArray();
      final int kradius = g.length;
      final double step = 1.0 / currentScale;
      halfKernel.resetQuick();
      halfKernel.add(1);

      for (int i = 1, j = 0; i < kradius; i++, j += upsample) {
        for (int k = 1; k < upsample; k++) {
          halfKernel.add(FastMath.exp(MathUtils.pow2((j + k) * step) / var2));
        }
        halfKernel.add(g[i]);
      }
    }
  }

  private void increaseKernel(int kradius) {
    if (halfKernel.size() < kradius) {
      final double step = 1.0 / currentScale;
      for (int i = halfKernel.size(); i < kradius; i++) {
        final double v = FastMath.exp(MathUtils.pow2(i * step) / var2);
        if (v == 0) {
          break;
        }
        halfKernel.add(v);
      }
    }
  }

  @SuppressWarnings("unused")
  private void increaseKernel(int scale, int kradius) {
    if (currentScale < scale) {
      currentScale = scale;
      halfKernel.resetQuick();
      halfKernel.add(1);
    }

    if (halfKernel.size() < kradius) {
      final double step = 1.0 / currentScale;
      for (int i = halfKernel.size(); i < kradius; i++) {
        final double v = FastMath.exp(MathUtils.pow2(i * step) / var2);
        if (v == 0) {
          break;
        }
        halfKernel.add(v);
      }
    }
  }

  private static double[] buildKernel(double[] kernel, int kradius, boolean edgeCorrection) {
    // Clip in the event that zeros occurred during computation
    if (kernel[kradius - 1] == 0) {
      while (kernel[--kradius] == 0) {
        /* decrement */
      }
      if (kradius == 1) {
        return new double[] {1};
      }
      kernel = Arrays.copyOf(kernel, 2 * kradius - 1);
    }

    // Edge correction
    if (edgeCorrection && kradius > 3) {
      double sqrtSlope = Double.MAX_VALUE;
      int radius = kradius;
      while (radius > kradius / 2) {
        radius--;
        final double a = Math.sqrt(kernel[radius]) / (kradius - radius);
        if (a < sqrtSlope) {
          sqrtSlope = a;
        } else {
          break;
        }
      }
      for (int r1 = radius + 2; r1 < kradius; r1++) {
        kernel[r1] = ((kradius - r1) * (kradius - r1) * sqrtSlope * sqrtSlope);
      }
    }

    // Normalise
    double sum = kernel[0];
    for (int i = 1; i < kradius; i++) {
      sum += 2 * kernel[i];
    }
    for (int i = 0; i < kradius; i++) {
      kernel[i] /= sum;
    }

    // Create symmetrical
    System.arraycopy(kernel, 0, kernel, kradius - 1, kradius);
    for (int i = kradius, j = i - 2; i < kernel.length; i++, j--) {
      kernel[j] = kernel[i];
    }
    return kernel;
  }

  /**
   * Create a 1-dimensional normalized Gaussian kernel with standard deviation sigma. To avoid a
   * step due to the cutoff at a finite value, the near-edge values are replaced by a 2nd-order
   * polynomial with its minimum=0 at the first out-of-kernel pixel. Thus, the kernel function has a
   * smooth 1st derivative in spite of finite length.
   *
   * @param sigma Standard deviation
   * @param range the range
   * @param edgeCorrection Set to true to perform the edge correction
   * @return The kernel, decaying towards zero, which would be reached at the first out of kernel
   *         index
   */
  public static double[] makeGaussianKernel(final double sigma, double range,
      boolean edgeCorrection) {
    // Limit range for the Gaussian
    if (range < 1) {
      range = 1;
    } else if (range > 38) {
      range = 38;
    }

    // Build half the kernel into the full kernel array. This is duplicated later.
    final int kradius = getGaussianHalfWidth(sigma, range) + 1;
    final double[] kernel = new double[2 * kradius - 1];

    kernel[0] = 1;
    final double s2 = sigma * sigma;
    for (int i = 1; i < kradius; i++) {
      // Gaussian function
      kernel[i] = FastMath.exp(-0.5 * i * i / s2);
      if (kernel[i] == 0) {
        break;
      }
    }

    return buildKernel(kernel, kradius, edgeCorrection);
  }

  /**
   * Create a 1-dimensional normalized Gaussian kernel with standard deviation sigma. The kernel is
   * constructed using the Error function (Erf) to compute the sum of the Gaussian from x-0.5 to
   * x+0.5 for each x sample point.
   *
   * @param sigma Standard deviation
   * @param range the range
   * @return The Erf kernel
   */
  public static double[] makeErfGaussianKernel(double sigma, double range) {
    // Limit range for the Gaussian
    if (range < 1) {
      range = 1;
    } else if (range > 38) {
      range = 38;
    }

    // Build half the kernel into the full kernel array. This is duplicated later.
    final int kradius = getGaussianHalfWidth(sigma, range) + 1;
    final double[] kernel = new double[2 * kradius - 1];

    if (kradius == 1) {
      kernel[0] = 1;
      return kernel;
    }

    // Use the error function to get the integral of the Gaussian.
    final double sqrt_var_by_2 = Math.sqrt(sigma * sigma * 2);

    double upper = org.apache.commons.math3.special.Erf.erf(-0.5 / sqrt_var_by_2);
    for (int i = 0; i < kradius; i++) {
      final double lower = upper;
      upper = org.apache.commons.math3.special.Erf.erf((i + 0.5) / sqrt_var_by_2);
      kernel[i] = (upper - lower) * 0.5;
      if (kernel[i] == 0) {
        break;
      }
    }

    return buildKernel(kernel, kradius, false);
  }
}
