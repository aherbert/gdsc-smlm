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

import java.util.Arrays;
import org.apache.commons.math3.util.FastMath;
import org.jtransforms.fft.DoubleFFT_1D;
import org.jtransforms.utils.CommonUtils;

/**
 * Simple class to perform convolution.
 */
public class Convolution {
  // Allow h as an input parameter name
  // CHECKSTYLE.OFF: ParameterName

  /** The maximum size supported for scaled convolution. */
  public static final int MAX = 1 << 30;

  /**
   * Interface to handle a convolution value.
   */
  @FunctionalInterface
  public interface ConvolutionValueProcedure {
    /**
     * Executes this procedure.
     *
     * @param value the value of the convolution
     * @return true, if further values should be computed
     */
    boolean execute(double value);
  }

  /**
   * Interface to handle two convolution valuess.
   */
  @FunctionalInterface
  public interface DoubleConvolutionValueProcedure {
    /**
     * Executes this procedure.
     *
     * @param value1 the value of the convolution of the first input
     * @param value2 the value of the convolution of the second input
     * @return true, if further values should be computed
     */
    boolean execute(double value1, double value2);
  }

  /**
   * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution"> convolution</a> between two
   * sequences.
   *
   * <p>The solution is obtained via straightforward computation of the convolution sum (and not via
   * FFT). Whenever the computation needs an element that would be located at an index outside the
   * input arrays, the value is assumed to be zero.
   *
   * <p>This has been taken from Apache Commons Math v3.3: org.apache.commons.math3.util.MathArrays
   *
   * @param x First sequence. Typically, this sequence will represent an input signal to a system.
   * @param h Second sequence. Typically, this sequence will represent the impulse response of the
   *        system.
   * @return the convolution of {@code x} and {@code h}. This array's length will be
   *         {@code x.length + h.length - 1}.
   * @throws IllegalArgumentException if either {@code x} or {@code h} is {@code null} or either
   *         {@code x} or {@code h} is empty.
   */
  public static double[] convolve(double[] x, double[] h) {
    checkInput(x, h);

    final int xLen = x.length;
    final int hLen = h.length;

    // Initialise the output array
    final int totalLength = xLen + hLen - 1;
    final double[] y = new double[totalLength];

    // Straightforward implementation of the convolution sum
    for (int n = 0; n < totalLength; n++) {
      double yn = 0;
      int hi = FastMath.max(0, n + 1 - xLen);
      int xi = n - hi;
      while (hi < hLen && xi >= 0) {
        yn += x[xi--] * h[hi++];
      }
      y[n] = yn;
    }

    return y;
  }

  /**
   * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution"> convolution</a> between two
   * sequences.
   *
   * <p>The solution is obtained via straightforward computation of the convolution sum (and not via
   * FFT). Whenever the computation needs an element that would be located at an index outside the
   * input arrays, the value is assumed to be zero.
   *
   * <p>The convolution is computed dynamically and can be stopped.
   *
   * @param x First sequence. Typically, this sequence will represent an input signal to a system.
   * @param h Second sequence. Typically, this sequence will represent the impulse response of the
   *        system.
   * @param result Output procedure for the convolution of {@code x} and {@code h}. This total
   *        number of times this is called will be {@code x.length + h.length - 1}.
   * @throws IllegalArgumentException if either {@code x} or {@code h} is {@code null} or either
   *         {@code x} or {@code h} is empty.
   */
  public static void convolve(double[] x, double[] h, ConvolutionValueProcedure result) {
    // As above but dynamically output the result

    checkInput(x, h);
    checkProcedure(result);

    final int xLen = x.length;
    final int hLen = h.length;

    // Initialise the output array
    final int totalLength = xLen + hLen - 1;
    if (totalLength <= 0) {
      throw new IllegalArgumentException("Unsupported size: " + ((long) xLen + hLen - 1));
    }

    // Straightforward implementation of the convolution sum
    for (int n = 0; n < totalLength; n++) {
      double yn = 0;
      int hi = FastMath.max(0, n + 1 - xLen);
      int xi = n - hi;
      while (hi < hLen && xi >= 0) {
        yn += x[xi--] * h[hi++];
      }
      if (!result.execute(yn)) {
        break;
      }
    }
  }

  /**
   * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution"> convolution</a> between one
   * sequence and two other sequences.
   *
   * <p>The solution is obtained via straightforward computation of the convolution sum (and not via
   * FFT). Whenever the computation needs an element that would be located at an index outside the
   * input arrays, the value is assumed to be zero.
   *
   * <p>This has been adapted from Apache Commons Math v3.3:
   * org.apache.commons.math3.util.MathArrays
   *
   * @param x First sequence.
   * @param h1 Second sequence 1.
   * @param h2 Second sequence 2.
   * @return the convolution of {@code x} and {@code h1} and {@code x} and {@code h2}. This array's
   *         length will be [2][{@code x.length + h1.length - 1}].
   * @throws IllegalArgumentException If any input is null or empty. If h1 and h2 are different
   *         lengths.
   */
  public static double[][] convolve(double[] x, double[] h1, double[] h2) {
    checkInput(x, h1, h2);

    final int xLen = x.length;
    final int hLen = h1.length;

    // Initialise the output array
    final int totalLength = xLen + hLen - 1;
    final double[][] y = new double[2][totalLength];

    // Straightforward implementation of the convolution sum
    for (int n = 0; n < totalLength; n++) {
      double yn1 = 0;
      double yn2 = 0;
      int hi = FastMath.max(0, n + 1 - xLen);
      int xi = n - hi;
      while (hi < hLen && xi >= 0) {
        yn1 += x[xi] * h1[hi];
        yn2 += x[xi] * h2[hi];
        xi--;
        hi++;
      }
      y[0][n] = yn1;
      y[1][n] = yn2;
    }

    return y;
  }

  /**
   * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution"> convolution</a> between one
   * sequence and two other sequences.
   *
   * <p>The solution is obtained via straightforward computation of the convolution sum (and not via
   * FFT). Whenever the computation needs an element that would be located at an index outside the
   * input arrays, the value is assumed to be zero.
   *
   * <p>This has been adapted from Apache Commons Math v3.3:
   * org.apache.commons.math3.util.MathArrays
   *
   * @param x First sequence.
   * @param h1 Second sequence 1.
   * @param h2 Second sequence 2.
   * @param result Output procedure for the convolution of {@code x} and {@code h1} or {@code h2}.
   *        This total number of times this is called will be {@code x.length + h1.length - 1}.
   * @throws IllegalArgumentException If any input is null or empty. If h1 and h2 are different
   *         lengths.
   */
  public static void convolve(double[] x, double[] h1, double[] h2,
      DoubleConvolutionValueProcedure result) {
    // As above but dynamically output the result

    checkInput(x, h1, h2);
    checkProcedure(result);

    final int xLen = x.length;
    final int hLen = h1.length;

    // Initialise the output array
    final int totalLength = xLen + hLen - 1;
    if (totalLength <= 0) {
      throw new IllegalArgumentException("Unsupported size: " + ((long) xLen + hLen - 1));
    }

    // Straightforward implementation of the convolution sum
    for (int n = 0; n < totalLength; n++) {
      double yn1 = 0;
      double yn2 = 0;
      int hi = FastMath.max(0, n + 1 - xLen);
      int xi = n - hi;
      while (hi < hLen && xi >= 0) {
        yn1 += x[xi] * h1[hi];
        yn2 += x[xi] * h2[hi];
        xi--;
        hi++;
      }
      if (!result.execute(yn1, yn2)) {
        break;
      }
    }
  }

  /**
   * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution"> convolution</a> between two
   * sequences.
   *
   * <p>The scale is used to increase the size of h dynamically to H with zero fill. The length of H
   * is thus ((h.length-1) * scale + 1);
   *
   * @param x First sequence.
   * @param h Second sequence.
   * @param scale the scale
   * @return the convolution of {@code x} and {@code h}. This array's length will be
   *         {@code x.length + H.length - 1}.
   * @throws IllegalArgumentException if either {@code x} or {@code h} is {@code null} or either
   *         {@code x} or {@code h} is empty.
   * @throws IllegalArgumentException if the scale is not strictly positive
   */
  public static double[] convolve(double[] x, double[] h, int scale) {
    checkInput(x, h);
    checkScale(scale);

    if (h.length == 1 || scale == 1) {
      // No scaling
      return convolve(x, h);
    }

    final int xLen = x.length;
    final int hLen = h.length;

    // Initialise the output array
    final int totalLength = checkLength(xLen, hLen, scale);
    final double[] y = new double[totalLength];

    // Convolution sum. x is reversed verses h.
    // h is scaled up with zeros.
    // This is equivalent to using x every interval of scale.
    for (int n = 0; n < totalLength; n++) {
      double yn = 0;
      // hi is the index in the scaled up distribution H
      final int hi = FastMath.max(0, n + 1 - xLen);
      // xi is the index in the input distribution x
      int xi = n - hi;

      // hi has to be scaled.
      // The modulus indicates how many values are zero
      // before the first non-zero value in H (in the descending direction).
      final int mod = hi % scale;
      // ihi is the index in input distribution h (in the descending direction).
      int ihi = hi / scale;
      // If there are non-zero value shift the indices
      if (mod != 0) {
        // Shift ihi by one for the next non-zero value (in the ascending direction)
        ihi++;
        // Shift xi by the number of zero values (in the descending direction)
        xi -= (scale - mod);
      }

      while (ihi < hLen && xi >= 0) {
        yn += x[xi] * h[ihi++];
        xi -= scale;
      }
      y[n] = yn;
    }

    return y;
  }

  /**
   * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution"> convolution</a> between two
   * sequences.
   *
   * <p>The scale is used to increase the size of h dynamically to H with zero fill. The length of H
   * is thus ((h.length-1) * scale + 1);
   *
   * <p>The convolution is computed dynamically and can be stopped.
   *
   * @param x First sequence.
   * @param h Second sequence.
   * @param scale the scale
   * @param result Output procedure for the convolution of {@code x} and {@code h}. This total
   *        number of times this is called will be {@code x.length + H.length - 1}.
   * @throws IllegalArgumentException if either {@code x} or {@code h} is {@code null} or either
   *         {@code x} or {@code h} is empty.
   * @throws IllegalArgumentException if the scale is not strictly positive
   */
  public static void convolve(double[] x, double[] h, int scale, ConvolutionValueProcedure result) {
    // As above but dynamically output the result

    checkInput(x, h);
    checkScale(scale);
    checkProcedure(result);

    final int xLen = x.length;
    final int hLen = h.length;

    // For consistency just support up to the max for integers.
    // This could be changed to use long for the index.
    final int totalLength = checkLength(xLen, hLen, scale);

    for (int n = 0; n < totalLength; n++) {
      double yn = 0;
      final int hi = FastMath.max(0, n + 1 - xLen);
      int xi = n - hi;
      final int mod = hi % scale;
      int ihi = hi / scale;
      if (mod != 0) {
        ihi++;
        xi -= (scale - mod);
      }
      while (ihi < hLen && xi >= 0) {
        yn += x[xi] * h[ihi++];
        xi -= scale;
      }
      if (!result.execute(yn)) {
        break;
      }
    }
  }

  /**
   * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution"> convolution</a> between two
   * sequences.
   *
   * <p>The scale is used to increase the size of h dynamically to H with zero fill. The length of H
   * is thus ((h.length-1) * scale + 1);
   *
   * @param x First sequence.
   * @param h1 Second sequence 1.
   * @param h2 Second sequence 2.
   * @param scale the scale
   * @return the convolution of {@code x} and {@code h1} and {@code x} and {@code h2}. This array's
   *         length will be [2][{@code x.length + h1.length - 1}].
   * @throws IllegalArgumentException If any input is null or empty. If h1 and h2 are different
   *         lengths.
   * @throws IllegalArgumentException if the scale is not strictly positive
   */
  public static double[][] convolve(double[] x, double[] h1, double[] h2, int scale) {
    checkInput(x, h1, h2);
    checkScale(scale);

    if (h1.length == 1 || scale == 1) {
      // No scaling
      return convolve(x, h1, h2);
    }

    final int xLen = x.length;
    final int hLen = h1.length;

    // Initialise the output array
    final int totalLength = checkLength(xLen, hLen, scale);
    final double[][] y = new double[2][totalLength];

    // Convolution sum. x is reversed verses h.
    // h is scaled up with zeros.
    // This is equivalent to using x every interval of scale.
    for (int n = 0; n < totalLength; n++) {
      double yn1 = 0;
      double yn2 = 0;
      // hi is the index in the scaled up distribution H
      final int hi = FastMath.max(0, n + 1 - xLen);
      // xi is the index in the input distribution x
      int xi = n - hi;

      // hi has to be scaled.
      // The modulus indicates how many values are zero
      // before the first non-zero value in H (in the descending direction).
      final int mod = hi % scale;
      // ihi is the index in input distribution h (in the descending direction).
      int ihi = hi / scale;
      // If there are non-zero value shift the indices
      if (mod != 0) {
        // Shift ihi by one for the next non-zero value (in the ascending direction)
        ihi++;
        // Shift xi by the number of zero values (in the descending direction)
        xi -= (scale - mod);
      }

      // int xi = n - ihi * scale
      while (ihi < hLen && xi >= 0) {
        yn1 += x[xi] * h1[ihi];
        yn2 += x[xi] * h2[ihi];
        ihi++;
        xi -= scale;
      }
      y[0][n] = yn1;
      y[1][n] = yn2;
    }

    return y;
  }

  /**
   * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution"> convolution</a> between two
   * sequences.
   *
   * <p>The scale is used to increase the size of h dynamically to H with zero fill. The length of H
   * is thus ((h1.length-1) * scale + 1);
   *
   * <p>The convolution is computed dynamically and can be stopped.
   *
   * @param x First sequence.
   * @param h1 Second sequence 1.
   * @param h2 Second sequence 2.
   * @param scale the scale
   * @param result Output procedure for the convolution of {@code x} and {@code h1} or {@code h2}.
   *        This total number of times this is called will be {@code x.length + h1.length - 1}.
   * @throws IllegalArgumentException If any input is null or empty. If h1 and h2 are different
   *         lengths.
   * @throws IllegalArgumentException if the scale is not strictly positive
   * @throws IllegalArgumentException if the output size is above the max size supported
   */
  public static void convolve(double[] x, double[] h1, double[] h2, int scale,
      DoubleConvolutionValueProcedure result) {
    // As above but dynamically output the result

    checkInput(x, h1, h2);
    checkScale(scale);
    checkProcedure(result);

    final int xLen = x.length;
    final int hLen = h1.length;

    // For consistency just support up to the max for integers.
    // This could be changed to use long for the index.
    final int totalLength = checkLength(xLen, hLen, scale);

    for (int n = 0; n < totalLength; n++) {
      double yn1 = 0;
      double yn2 = 0;
      final int hi = FastMath.max(0, n + 1 - xLen);
      int xi = n - hi;
      final int mod = hi % scale;
      int ihi = hi / scale;
      if (mod != 0) {
        ihi++;
        xi -= (scale - mod);
      }
      while (ihi < hLen && xi >= 0) {
        yn1 += x[xi] * h1[ihi];
        yn2 += x[xi] * h2[ihi];
        ihi++;
        xi -= scale;
      }
      if (!result.execute(yn1, yn2)) {
        break;
      }
    }
  }

  /**
   * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution"> convolution</a> between two
   * sequences.
   *
   * <p>The solution is obtained via multiplication in the frequency domain. To reduce edge wrap
   * artifacts the input signals should be windowed to zero at the ends.
   *
   * @param x First sequence. Typically, this sequence will represent an input signal to a system.
   * @param h Second sequence. Typically, this sequence will represent the impulse response of the
   *        system.
   * @return the convolution of {@code x} and {@code h}. This array's length will be
   *         {@code x.length + h.length - 1}.
   * @throws IllegalArgumentException if either {@code x} or {@code h} is {@code null} or either
   *         {@code x} or {@code h} is empty.
   */
  public static double[] convolveFft(double[] x, double[] h) {
    checkInput(x, h);

    final int xLen = x.length;
    final int hLen = h.length;
    final int totalLength = xLen + hLen - 1;

    // Get length to a power of 2
    final int newL = CommonUtils.nextPow2(totalLength);

    // Double the new length for complex values in DoubleFFT_1D
    x = Arrays.copyOf(x, 2 * newL);
    h = Arrays.copyOf(h, x.length);

    final DoubleFFT_1D fft = new DoubleFFT_1D(newL);

    // FFT
    fft.realForwardFull(x);
    fft.realForwardFull(h);

    // Complex multiply. Reuse data array
    for (int i = 0; i < x.length; i += 2) {
      final int j = i + 1;
      final double xi = x[i];
      final double xj = x[j];
      final double hi = h[i];
      final double hj = h[j];
      h[i] = hi * xi - hj * xj;
      h[j] = hi * xj + hj * xi;
    }

    // Inverse FFT
    fft.complexInverse(h, true);

    // Fill result with real part
    final double[] y = new double[totalLength];
    for (int i = 0; i < totalLength; i++) {
      y[i] = h[2 * i];
    }
    return y;
  }

  /**
   * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution"> convolution</a> between one
   * sequence and two other sequences.
   *
   * <p>The solution is obtained via multiplication in the frequency domain. To reduce edge wrap
   * artifacts the input signals should be windowed to zero at the ends.
   *
   * @param x First sequence.
   * @param h1 Second sequence 1.
   * @param h2 Second sequence 2.
   * @return the convolution of {@code x} and {@code h1} and {@code x} and {@code h2}. This array's
   *         length will be [2][{@code x.length + h1.length - 1}].
   * @throws IllegalArgumentException If any input is null or empty. If h1 and h2 are different
   *         lengths.
   */
  public static double[][] convolveFft(double[] x, double[] h1, double[] h2) {
    checkInput(x, h1, h2);

    final int xLen = x.length;
    final int hLen = h1.length;
    final int totalLength = xLen + hLen - 1;

    // Get length to a power of 2
    final int newL = CommonUtils.nextPow2(totalLength);

    // Double the new length for complex values in DoubleFFT_1D
    x = Arrays.copyOf(x, 2 * newL);
    h1 = Arrays.copyOf(h1, x.length);
    h2 = Arrays.copyOf(h2, x.length);

    final DoubleFFT_1D fft = new DoubleFFT_1D(newL);

    // FFT
    fft.realForwardFull(x);
    fft.realForwardFull(h1);
    fft.realForwardFull(h2);

    // Complex multiply. Reuse data array
    for (int i = 0; i < x.length; i += 2) {
      final int j = i + 1;
      final double xi = x[i];
      final double xj = x[j];
      double hi = h1[i];
      double hj = h1[j];
      h1[i] = hi * xi - hj * xj;
      h1[j] = hi * xj + hj * xi;
      hi = h2[i];
      hj = h2[j];
      h2[i] = hi * xi - hj * xj;
      h2[j] = hi * xj + hj * xi;
    }

    // Inverse FFT
    fft.complexInverse(h1, true);
    fft.complexInverse(h2, true);

    // Fill result with real part
    final double[][] y = new double[2][totalLength];
    for (int i = 0; i < totalLength; i++) {
      y[0][i] = h1[2 * i];
      y[1][i] = h2[2 * i];
    }
    return y;
  }

  /**
   * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution"> convolution</a> between two
   * sequences.
   *
   * <p>The solution is obtained using either the spatial or frequency domain depending on the size.
   * The switch is made when the min array length is above 127 and the product of the lengths is
   * above 40000. Speed tests have been performed for single threaded FFT computation. The FFT
   * library begins multi-threaded computation when the size of the array is above length 8192.
   *
   * @param x First sequence. Typically, this sequence will represent an input signal to a system.
   * @param h Second sequence. Typically, this sequence will represent the impulse response of the
   *        system.
   * @return the convolution of {@code x} and {@code h}. This array's length will be
   *         {@code x.length + h.length - 1}.
   * @throws IllegalArgumentException if either {@code x} or {@code h} is {@code null} or either
   *         {@code x} or {@code h} is empty.
   */
  public static double[] convolveFast(double[] x, double[] h) {
    checkInput(x, h);
    if (isFft(x.length, h.length)) {
      return convolveFft(x, h);
    }
    return convolve(x, h);
  }

  /**
   * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution"> convolution</a> between one
   * sequence and two other sequences.
   *
   * <p>The solution is obtained using either the spatial or frequency domain depending on the size.
   * The switch is made when the min array length is above 127 and the product of the lengths is
   * above 40000. Speed tests have been performed for single threaded FFT computation. The FFT
   * library begins multi-threaded computation when the size of the array is above length 8192.
   *
   * @param x First sequence.
   * @param h1 Second sequence 1.
   * @param h2 Second sequence 2.
   * @return the convolution of {@code x} and {@code h1} and {@code x} and {@code h2}. This array's
   *         length will be [2][{@code x.length + h1.length - 1}].
   * @throws IllegalArgumentException If any input is null or empty. If h1 and h2 are different
   *         lengths.
   */
  public static double[][] convolveFast(double[] x, double[] h1, double[] h2) {
    checkInput(x, h1, h2);
    if (isFft(x.length, h1.length)) {
      return convolveFft(x, h1, h2);
    }
    return convolve(x, h1, h2);
  }

  /**
   * Checks if convolution will use the FFT method.
   *
   * @param length1 the length 1
   * @param length2 the length 2
   * @return true, if using the FFT method
   */
  public static boolean isFft(int length1, int length2) {
    // See Junit class ConvolveTest to determine when to switch to the FFT method.
    // This is not perfect for all length combinations but the switch will happen
    // when the two methods are roughly the same speed.
    int min;
    int max;
    if (length1 < length2) {
      min = length1;
      max = length2;
    } else {
      min = length2;
      max = length1;
    }
    return (min >= 128 && (long) min * (long) max > 40000L);
  }

  private static void checkInput(double[] x, double[] h) {
    if (x == null) {
      throw new IllegalArgumentException("Input x is null");
    }
    if (h == null) {
      throw new IllegalArgumentException("Input h is null");
    }
    if (x.length == 0 || h.length == 0) {
      throw new IllegalArgumentException("Input x or h have no length");
    }
  }

  private static void checkInput(double[] x, double[] h1, double[] h2) {
    if (x == null) {
      throw new IllegalArgumentException("Input x is null");
    }
    if (h1 == null) {
      throw new IllegalArgumentException("Input h1 is null");
    }
    if (h2 == null) {
      throw new IllegalArgumentException("Input h2 is null");
    }
    if (x.length == 0 || h1.length == 0) {
      throw new IllegalArgumentException("Input x or h1 have no length");
    }
    if (h1.length != h2.length) {
      throw new IllegalArgumentException("Input h1 and h2 have different length");
    }
  }

  private static void checkProcedure(Object procedure) {
    if (procedure == null) {
      throw new IllegalArgumentException("Procedure is null");
    }
  }

  private static void checkScale(int scale) {
    if (scale < 1) {
      throw new IllegalArgumentException("Scale must be strictly positive");
    }
  }

  private static int checkLength(int xlength, int hlength, int scale) {
    final double scaledHlength = (double) (hlength - 1) * scale + 1;
    final double totalLength = xlength + scaledHlength - 1;
    if (totalLength > MAX) {
      throw new IllegalArgumentException("Scale creates unsupported size: " + totalLength);
    }
    return (int) totalLength;
  }
}
