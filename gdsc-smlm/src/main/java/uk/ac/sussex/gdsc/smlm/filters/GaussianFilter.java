/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.filters;

import java.awt.Rectangle;
import uk.ac.sussex.gdsc.core.utils.MathUtils;

/**
 * Computes a Gaussian convolution in the spatial domain for each point within the array.
 *
 * <p>Adapted from {@code ij.plugin.filter.GaussianBlur}.
 */
public class GaussianFilter extends BaseWeightedFilter {
  /** number of pixels to add for upscaling. */
  private static final int UPSCALE_K_RADIUS = 2;
  /** minimum standard deviation in the downscaled image. */
  private static final double MIN_DOWNSCALED_SIGMA = 4.;

  /**
   * The default accuracy. This is the highest acceptable value of the accuracy for maximum speed.
   * Lower values will increase the kernel size.
   */
  public static final double DEFAULT_ACCURACY = 0.02;

  private final double accuracy;

  private double lastSigma;
  private int lastMaxRadius;
  private float[][] kernel;
  private float[] downscaleKernel;
  private float[] upscaleKernel;
  private int lastUnitLength;

  private Normaliser normaliser;
  private double sx;
  private double sy;

  @Override
  protected void newWeights() {
    normaliser = null;
  }

  /**
   * Updates the weighted normaliser for the Gaussian.
   *
   * @param sigmaX The Gaussian standard deviation in X
   * @param sigmaY The Gaussian standard deviation in Y
   */
  private void updateWeightedNormaliser(final double sigmaX, final double sigmaY) {
    // Cache the normaliser
    if (normaliser == null || sx != sigmaX || sy != sigmaY) {
      final float[] normalisation = weights.clone();
      sx = sigmaX;
      sy = sigmaY;
      final GaussianFilter gf = new GaussianFilter(accuracy);
      gf.convolve(normalisation, weightWidth, weightHeight, sigmaX, sigmaY);
      normaliser = new PerPixelNormaliser(normalisation);
    }
  }

  /**
   * Instantiates a new gaussian filter.
   *
   * <p>Use the default accuracy of {@value #DEFAULT_ACCURACY}.
   */
  public GaussianFilter() {
    this(DEFAULT_ACCURACY);
  }

  /**
   * Instantiates a new gaussian filter.
   *
   * @param accuracy Accuracy of kernel, should not be above 0.02. Better (lower) accuracy needs
   *        slightly more computing time as the kernel size is increased.
   */
  public GaussianFilter(double accuracy) {
    this.accuracy = accuracy;
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected GaussianFilter(GaussianFilter source) {
    super(source);
    // These are thread safe
    accuracy = source.accuracy;

    lastSigma = source.lastSigma;
    lastMaxRadius = source.lastMaxRadius;
    kernel = source.kernel;
    downscaleKernel = source.downscaleKernel;
    upscaleKernel = source.upscaleKernel;
    lastUnitLength = source.lastUnitLength;

    normaliser = source.normaliser;
    sx = source.sx;
    sy = source.sy;
  }

  /**
   * Create a copy.
   *
   * @return the copy
   */
  public GaussianFilter copy() {
    return new GaussianFilter(this);
  }

  /**
   * Compute the Gaussian convolution. Pixels within border regions (defined by 3 sigma) are
   * unchanged.
   *
   * <p>Note: the input data is destructively modified.
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param sigma The Gaussian standard deviation
   */
  public void convolveInternal(float[] data, final int maxx, final int maxy, final double sigma) {
    final int border = getBorder(sigma);
    final Rectangle roi = new Rectangle(border, border, maxx - 2 * border, maxy - 2 * border);
    if (roi.width < 1 || roi.height < 1) {
      return;
    }
    // Q. Should the extra lines parameter be used here?
    convolve(data, roi, maxx, maxy, sigma, sigma, border);
  }

  /**
   * Get the border that will be ignored for the specified Gaussian standard deviation.
   *
   * @param sigma the Gaussian standard deviation
   * @return The border
   */
  public static int getBorder(double sigma) {
    return (int) (3 * sigma);
  }

  /**
   * Compute the Gaussian convolution.
   *
   * <p>Note: the input data is destructively modified.
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param sigma The Gaussian standard deviation
   */
  public void convolve(float[] data, final int maxx, final int maxy, final double sigma) {
    final Rectangle roi = new Rectangle(maxx, maxy);
    convolve(data, roi, maxx, maxy, sigma, sigma, 0);
  }

  /**
   * Compute the Gaussian convolution.
   *
   * <p>Note: the input data is destructively modified.
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param sigmaX The Gaussian standard deviation in X
   * @param sigmaY The Gaussian standard deviation in Y
   */
  public void convolve(float[] data, final int maxx, final int maxy, final double sigmaX,
      final double sigmaY) {
    final Rectangle roi = new Rectangle(maxx, maxy);
    convolve(data, roi, maxx, maxy, sigmaX, sigmaY, 0);
  }

  private void convolve(float[] data, Rectangle roi, final int maxx, final int maxy,
      final double sigmaX, final double sigmaY, final int extraLines) {
    // If not all data points will be over-written clone the input before blurring
    final float[] wdata = (extraLines != 0) ? data.clone() : data;

    if (hasWeights()) {
      final int size = data.length;
      if (weights.length != size || this.weightWidth != maxx || this.weightHeight != maxy) {
        throw new IllegalStateException("Weights are not the correct size");
      }

      updateWeightedNormaliser(sigmaX, sigmaY);

      // Apply weights
      for (int i = 0; i < size; i++) {
        wdata[i] *= weights[i];
      }

      blur1Direction(wdata, roi, maxx, maxy, sigmaX, true, extraLines);
      blur1Direction(wdata, roi, maxx, maxy, sigmaY, false, 0);

      // Normalise. This assumes the roi is only ever constructed as a border.
      if (extraLines != 0) {
        normaliser.normalise(wdata, data, maxx, maxy, extraLines);
      } else {
        // wdata == data
        normaliser.normalise(wdata, size);
      }
    } else {
      blur1Direction(wdata, roi, maxx, maxy, sigmaX, true, extraLines);
      blur1Direction(wdata, roi, maxx, maxy, sigmaY, false, 0);

      if (extraLines != 0) {
        // Copy back to the input array so that the extra blurred lines
        // are not present.
        NonNormaliser.INSTANCE.normalise(wdata, data, maxx, maxy, extraLines);
      }
    }
  }

  /**
   * Blur an image in one direction (x or y) by a Gaussian.
   *
   * @param pixels The input/output data (packed in YX order)
   * @param roi The region to blur
   * @param width The width of the data
   * @param height The height of the data
   * @param sigma Standard deviation of the Gaussian
   * @param xDirection True for blurring in x direction, false for y direction
   * @param extraLines Number of lines (parallel to the blurring direction) below and above the roi
   *        bounds that should be processed.
   */
  private void blur1Direction(final float[] pixels, Rectangle roi, final int width,
      final int height, final double sigma, final boolean xDirection, final int extraLines) {
    // number of points per line (line can be a row or column)
    final int length = xDirection ? width : height;
    // increment of the pixels array index to the next point in a line
    final int pointInc = xDirection ? 1 : width;
    // increment of the pixels array index to the next line
    final int lineInc = xDirection ? width : 1;
    // the first line to process
    final int lineFromA = (xDirection ? roi.y : roi.x) - extraLines;
    final int lineFrom;
    if (lineFromA < 0) {
      lineFrom = 0;
    } else {
      lineFrom = lineFromA;
    }
    // the last line+1 to process
    final int lineToA = (xDirection ? roi.y + roi.height : roi.x + roi.width) + extraLines;
    final int lineTo;
    if (lineToA > (xDirection ? height : width)) {
      lineTo = (xDirection ? height : width);
    } else {
      lineTo = lineToA;
    }
    // first point of a line that needs to be written
    final int writeFrom = xDirection ? roi.x : roi.y;
    final int writeTo = xDirection ? roi.x + roi.width : roi.y + roi.height;

    // large radius (sigma): scale down, then convolve, then scale up
    final boolean doDownscaling = sigma > 2 * MIN_DOWNSCALED_SIGMA + 0.5;
    final int reduceBy = doDownscaling
        // downscale by this factor
        ? Math.min((int) Math.floor(sigma / MIN_DOWNSCALED_SIGMA), length)
        : 1;
    // Downscaling and upscaling blur the image a bit - we have to correct the standard deviation
    // for this: Downscaling gives std deviation sigma = 1/sqrt(3); upscale gives sigma = 1/2 (in
    // downscaled pixels). All sigma^2 values add to full sigma^2, which should be the desired value
    final double sigmaGauss =
        doDownscaling ? Math.sqrt(sigma * sigma / (reduceBy * reduceBy) - 1. / 3. - 1. / 4.)
            : sigma;
    final int maxLength = doDownscaling
        // downscaled line can't be longer
        ? (length + reduceBy - 1) / reduceBy + 2 * (UPSCALE_K_RADIUS + 1)
        : length;
    final float[][] gaussKernel = makeGaussianKernel(sigmaGauss, maxLength);
    // Gaussian kernel radius after upscaling
    final int kRadius = gaussKernel[0].length * reduceBy;
    // not including broadening by downscale & upscale
    final int readFrom = (writeFrom - kRadius < 0) ? 0 : writeFrom - kRadius;
    final int readTo = (writeTo + kRadius > length) ? length : writeTo + kRadius;
    // line length for convolution
    final int newLength =
        doDownscaling ? (readTo - readFrom + reduceBy - 1) / reduceBy + 2 * (UPSCALE_K_RADIUS + 1)
            : length;
    // input point corresponding to cache index 0
    final int unscaled0 = readFrom - (UPSCALE_K_RADIUS + 1) * reduceBy;
    // the following is relevant for upscaling only
    if (doDownscaling) {
      createScalingKernels(reduceBy);
    }

    // holds data before convolution (after downscaling, if any)
    final float[] cache1 = new float[newLength];
    // holds data after convolution
    final float[] cache2 = doDownscaling ? new float[newLength] : null;

    int pixel0 = lineFrom * lineInc;
    for (int line = lineFrom; line < lineTo; line += 1, pixel0 += lineInc) {
      if (doDownscaling) {
        downscaleLine(pixels, cache1, downscaleKernel, reduceBy, pixel0, unscaled0, length,
            pointInc, newLength);
        convolveLine(cache1, cache2, gaussKernel, 1, newLength - 1, 0, 1);
        upscaleLine(cache2, pixels, upscaleKernel, reduceBy, pixel0, unscaled0, writeFrom, writeTo,
            pointInc);
      } else {
        int pi = pixel0 + readFrom * pointInc;
        for (int i = readFrom; i < readTo; i++, pi += pointInc) {
          cache1[i] = pixels[pi];
        }
        convolveLine(cache1, pixels, gaussKernel, writeFrom, writeTo, pixel0, pointInc);
      }
    }
  }

  private void createScalingKernels(int unitLength) {
    if (downscaleKernel == null || lastUnitLength != unitLength) {
      lastUnitLength = unitLength;
      downscaleKernel = makeDownscaleKernel(unitLength);
      upscaleKernel = makeUpscaleKernel(unitLength);
    }
  }

  /**
   * Scale a line (row or column or part thereof) down by a factor {@code reduceBy} and write
   * the result into {@code cache}. Input line pixel # {@code unscaled0} will correspond
   * to output line pixel # 0. {@code unscaled0} may be negative. Out-of-line pixels of the
   * input are replaced by the edge pixels.
   *
   * @param pixels input array
   * @param cache output array
   * @param kernel downscale kernel, runs form -1.5 to +1.5 in downscaled coordinates
   * @param reduceBy downscaling factor
   * @param pixel0 index in pixels array corresponding to start of line or column
   * @param unscaled0 index in input line corresponding to output line index 0, May be negative.
   * @param length length of full input line or column
   * @param pointInc spacing of values in input array (1 for lines, image width for columns)
   * @param newLength length of downscaled data
   */
  private static void downscaleLine(final float[] pixels, final float[] cache,
      final float[] kernel, final int reduceBy, final int pixel0, final int unscaled0,
      final int length, final int pointInc, final int newLength) {
    int pi = pixel0 + pointInc * (unscaled0 - reduceBy * 3 / 2); // pointer in pixels array
    final int pLast = pixel0 + pointInc * (length - 1);
    for (int xout = -1; xout <= newLength; xout++) {
      float sum0 = 0;
      float sum1 = 0;
      float sum2 = 0;
      for (int x = 0; x < reduceBy; x++, pi += pointInc) {
        final float v = pixels[MathUtils.clip(pixel0, pLast, pi)];
        sum0 += v * kernel[x + 2 * reduceBy];
        sum1 += v * kernel[x + reduceBy];
        sum2 += v * kernel[x];
      }
      if (xout > 0) {
        cache[xout - 1] += sum0;
      }
      if (xout >= 0 && xout < newLength) {
        cache[xout] += sum1;
      }
      if (xout + 1 < newLength) {
        cache[xout + 1] = sum2;
      }
    }
  }

  /*
   * Create a kernel for downscaling. The kernel function preserves norm and 1st moment (i.e.,
   * position) and has fixed 2nd moment, (in contrast to linear interpolation). In scaled space, the
   * length of the kernel runs from -1.5 to +1.5, and the standard deviation is 1/2. Array index
   * corresponding to the kernel center is unitLength*3/2
   */
  private static float[] makeDownscaleKernel(final int unitLength) {
    final int mid = unitLength * 3 / 2;
    final float[] kernel = new float[3 * unitLength];
    for (int i = 0; i <= unitLength / 2; i++) {
      final double x = i / (double) unitLength;
      final float v = (float) ((0.75 - x * x) / unitLength);
      kernel[mid - i] = v;
      kernel[mid + i] = v;
    }
    for (int i = unitLength / 2 + 1; i < (unitLength * 3 + 1) / 2; i++) {
      final double x = i / (double) unitLength;
      final float v = (float) ((0.125 + 0.5 * (x - 1) * (x - 2)) / unitLength);
      kernel[mid - i] = v;
      kernel[mid + i] = v;
    }
    return kernel;
  }

  /**
   * Scale a line up by factor {@code reduceBy} and write as a row or column (or part thereof)
   * to the pixels array of a FloatProcessor.
   */
  private static void upscaleLine(final float[] cache, final float[] pixels,
      final float[] kernel, final int reduceBy, final int pixel0, final int unscaled0,
      final int writeFrom, final int writeTo, final int pointInc) {
    int pi = pixel0 + pointInc * writeFrom;
    for (int xout = writeFrom; xout < writeTo; xout++, pi += pointInc) {
      // the corresponding point in the cache (if exact) or the one above
      final int xin = (xout - unscaled0 + reduceBy - 1) / reduceBy;
      final int x = reduceBy - 1 - (xout - unscaled0 + reduceBy - 1) % reduceBy;
      pixels[pi] = cache[xin - 2] * kernel[x] + cache[xin - 1] * kernel[x + reduceBy]
          + cache[xin] * kernel[x + 2 * reduceBy] + cache[xin + 1] * kernel[x + 3 * reduceBy];
    }
  }

  /**
   * Create a kernel for upscaling. The kernel function is a convolution of four unit squares, i.e.,
   * four uniform kernels with value +1 from -0.5 to +0.5 (in downscaled coordinates). The second
   * derivative of this kernel is smooth, the third is not. Its standard deviation is 1/sqrt(3) in
   * downscaled coordinates. The kernel runs from [-2 to +2[, corresponding to array index 0 ...
   * 4*unitLength (whereby the last point is not in the array any more).
   */
  private static float[] makeUpscaleKernel(final int unitLength) {
    final float[] kernel = new float[4 * unitLength];
    final int mid = 2 * unitLength;
    kernel[0] = 0;
    for (int i = 0; i < unitLength; i++) {
      final double x = i / (double) unitLength;
      final float v = (float) (2. / 3. - x * x * (1 - 0.5 * x));
      kernel[mid + i] = v;
      kernel[mid - i] = v;
    }
    for (int i = unitLength; i < 2 * unitLength; i++) {
      final double x = i / (double) unitLength;
      final float v = (float) ((2. - x) * (2. - x) * (2. - x) / 6.);
      kernel[mid + i] = v;
      kernel[mid - i] = v;
    }
    return kernel;
  }

  /**
   * Convolve a line with a symmetric kernel and write to a separate array, possibly the pixels
   * array of a FloatProcessor (as a row or column or part thereof).
   *
   * @param input Input array containing the line
   * @param pixels Float array for output, can be the pixels of a FloatProcessor
   * @param kernel "One-sided" kernel array, kernel[0][n] must contain the kernel itself,
   *        kernel[1][n] must contain the running sum over all kernel elements from kernel[0][n+1]
   *        to the periphery. The kernel must be normalized, i.e. sum(kernel[0][n]) = 1 where n runs
   *        from the kernel periphery (last element) to 0 and back. Normalization should include all
   *        kernel points, also these not calculated because they are not needed.
   * @param writeFrom Index of the first point in the line that should be written
   * @param writeTo Index+1 of the last point in the line that should be written
   * @param point0 Array index of first element of the 'line' in pixels (i.e., lineNumber * lineInc)
   * @param pointInc Increment of the pixels array index to the next point (for an ImageProcessor,
   *        it should be {@code 1} for a row, {@code width} for a column)
   */
  private static void convolveLine(final float[] input, final float[] pixels,
      final float[][] kernel, final int writeFrom, final int writeTo, final int point0,
      final int pointInc) {
    final int length = input.length;
    final float first = input[0]; // out-of-edge pixels are replaced by nearest edge pixels
    final float last = input[length - 1];
    final float[] kern = kernel[0]; // the kernel itself
    final float kern0 = kern[0];
    final float[] kernSum = kernel[1]; // the running sum over the kernel
    final int kRadius = kern.length;
    final int firstPart = kRadius < length ? kRadius : length;
    int pi = point0 + writeFrom * pointInc;
    int index = writeFrom;
    // while the sum would include pixels < 0
    for (; index < firstPart; index++, pi += pointInc) {
      float result = input[index] * kern0;
      result += kernSum[index] * first;
      if (index + kRadius > length) {
        result += kernSum[length - index - 1] * last;
      }
      for (int k = 1; k < kRadius; k++) {
        float value = 0;
        if (index - k >= 0) {
          value += input[index - k];
        }
        if (index + k < length) {
          value += input[index + k];
        }
        result += kern[k] * value;
      }
      pixels[pi] = result;
    }
    final int iEndInside = length - kRadius < writeTo ? length - kRadius : writeTo;
    // while only pixels within the line are be addressed (the easy case)
    for (; index < iEndInside; index++, pi += pointInc) {
      float result = input[index] * kern0;
      for (int k = 1; k < kRadius; k++) {
        result += kern[k] * (input[index - k] + input[index + k]);
      }
      pixels[pi] = result;
    }
    // while the sum would include pixels >= length
    for (; index < writeTo; index++, pi += pointInc) {
      float result = input[index] * kern0;
      if (index < kRadius) {
        result += kernSum[index] * first;
      }
      if (index + kRadius >= length) {
        result += kernSum[length - index - 1] * last;
      }
      for (int k = 1; k < kRadius; k++) {
        float value = 0;
        if (index - k >= 0) {
          value += input[index - k];
        }
        if (index + k < length) {
          value += input[index + k];
        }
        result += kern[k] * value;
      }
      pixels[pi] = result;
    }
  }

  /**
   * Create a 1-dimensional normalized Gaussian kernel with standard deviation sigma and the running
   * sum over the kernel.
   *
   * <p>Note: this is one side of the kernel only, not the full kernel as used by the Convolver
   * class of ImageJ. To avoid a step due to the cutoff at a finite value, the near-edge values are
   * replaced by a 2nd-order polynomial with its minimum=0 at the first out-of-kernel pixel. Thus,
   * the kernel function has a smooth 1st derivative in spite of finite length.
   *
   * @param sigma Standard deviation, i.e. radius of decay to 1/sqrt(e), in pixels.
   * @param maxRadius Maximum radius of the kernel: Limits kernel size in case of large sigma,
   *        should be set to image width or height. For small values of maxRadius, the kernel
   *        returned may have a larger radius, however.
   * @return A 2*n array. Array[0][n] is the kernel, decaying towards zero, which would be reached
   *         at kernel.length (unless kernel size is limited by maxRadius). Array[1][n] holds the
   *         sum over all kernel values > n, including non-calculated values in case the kernel size
   *         is limited by {@code maxRadius}.
   */
  private float[][] makeGaussianKernel(final double sigma, int maxRadius) {
    if (maxRadius < 50) {
      maxRadius = 50; // too small maxRadius would result in inaccurate sum.
    }

    // Use cached kernel
    if (kernel != null && sigma == lastSigma && maxRadius == lastMaxRadius) {
      return kernel;
    }

    int kradius = getHalfWidth(sigma) + 1;
    if (kradius > maxRadius) {
      kradius = maxRadius;
    }
    kernel = new float[2][kradius];
    lastSigma = sigma;
    lastMaxRadius = maxRadius;

    for (int i = 0; i < kradius; i++) {
      // Gaussian function
      kernel[0][i] = (float) (Math.exp(-0.5 * i * i / sigma / sigma));
    }
    if (kradius < maxRadius && kradius > 3) { // edge correction
      double sqrtSlope = Double.MAX_VALUE;
      int radius = kradius;
      while (radius > kradius / 2) {
        radius--;
        final double a = Math.sqrt(kernel[0][radius]) / (kradius - radius);
        if (a < sqrtSlope) {
          sqrtSlope = a;
        } else {
          break;
        }
      }
      for (int r1 = radius + 2; r1 < kradius; r1++) {
        kernel[0][r1] = (float) ((kradius - r1) * (kradius - r1) * sqrtSlope * sqrtSlope);
      }
    }
    double sum; // sum over all kernel elements for normalization
    if (kradius < maxRadius) {
      sum = kernel[0][0];
      for (int i = 1; i < kradius; i++) {
        sum += 2 * kernel[0][i];
      }
    } else {
      sum = sigma * Math.sqrt(2 * Math.PI);
    }

    double rsum = 0.5 + 0.5 * kernel[0][0] / sum;
    for (int i = 0; i < kradius; i++) {
      final double v = (kernel[0][i] / sum);
      kernel[0][i] = (float) v;
      rsum -= v;
      kernel[1][i] = (float) rsum;
    }
    return kernel;
  }

  /**
   * Get half the width of the region smoothed by the filter for the specified standard deviation.
   * The full region size is 2N + 1
   *
   * @param sigma the Gaussian standard deviation
   * @return The half width
   */
  public int getHalfWidth(double sigma) {
    return (int) Math.ceil(sigma * Math.sqrt(-2 * Math.log(accuracy)));
  }
}
