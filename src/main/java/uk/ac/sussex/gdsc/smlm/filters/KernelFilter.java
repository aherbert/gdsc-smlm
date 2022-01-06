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

package uk.ac.sussex.gdsc.smlm.filters;

import ij.process.FloatProcessor;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;

/**
 * Computes a convolution in the spatial domain for each point within the array. Pixels outside the
 * array are set to the value of the appropriate edge pixel.
 *
 * <p>Adapted from {@link ij.plugin.filter.Convolver}
 */
public class KernelFilter extends BaseWeightedFilter {
  /** The kernel. */
  protected float[] kernel;
  /** The kernel width. */
  protected final int kw;
  /** The kernel height. */
  protected final int kh;
  /** The kernel scale (1.0/[sum of the kernel]). */
  protected final double scale;

  private Normaliser normaliser;
  private boolean convolution;

  @Override
  protected void newWeights() {
    normaliser = null;
  }

  /**
   * Updates the weighted normaliser for the convolution.
   */
  private void updateWeightedNormaliser() {
    // Cache the normaliser
    if (normaliser == null) {
      final float[] normalisation = weights.clone();
      final KernelFilter kf = new KernelFilter(kernel, kw, kh, scale);
      kf.convolve(normalisation, weightWidth, weightHeight);
      normaliser = new PerPixelNormaliser(normalisation);
    }
  }

  /**
   * Instantiates a new kernel filter.
   *
   * @param kernel the kernel
   * @param kw the kernel width (must be odd)
   * @param kh the kernel height (must be odd)
   */
  public KernelFilter(float[] kernel, int kw, int kh) {
    if (kw < 1 || kh < 1) {
      throw new IllegalArgumentException("Kernel width & height must be positive");
    }
    if (kw * kh != kernel.length) {
      throw new IllegalArgumentException("Kernel width x height != kernel length");
    }
    if (isEven(kw) || isEven(kh)) {
      throw new IllegalArgumentException("Kernel width or height not odd (" + kw + "x" + kh + ")");
    }
    this.kernel = kernel.clone();
    this.kw = kw;
    this.kh = kh;
    this.scale = getScale(kernel);
  }

  /**
   * Instantiates a new kernel filter with no normaliser.
   *
   * @param kernel the kernel
   * @param kw the kernel width (must be odd)
   * @param kh the kernel height (must be odd)
   * @param scale the scale
   */
  private KernelFilter(float[] kernel, int kw, int kh, double scale) {
    this.kernel = kernel.clone();
    this.kw = kw;
    this.kh = kh;
    this.scale = scale;
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected KernelFilter(KernelFilter source) {
    super(source);
    // Copy the kernel as it can be rotated in place
    kernel = source.getKernal();
    // These are thread safe
    kw = source.kw;
    kh = source.kh;
    scale = source.scale;
    normaliser = source.normaliser;
    convolution = source.convolution;
  }

  /**
   * Create a copy.
   *
   * @return the copy
   */
  public KernelFilter copy() {
    return new KernelFilter(this);
  }

  /**
   * Gets the scale of the kernel (i.e. 1/sum). This is used to scale the convolution result.
   *
   * @param kernel the kernel
   * @return the scale
   */
  public static double getScale(float[] kernel) {
    double sum = 0.0;
    for (int i = 0; i < kernel.length; i++) {
      sum += kernel[i];
    }
    if (sum != 0.0) {
      return 1.0 / sum;
    }
    return 1.0;
  }

  /**
   * Compute the convolution.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   */
  public void convolve(float[] data, final int maxx, final int maxy) {
    convolveInternal(data, maxx, maxy, 0);
  }

  /**
   * Compute the convolution.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data the data
   * @param maxx the maxx
   * @param maxy the maxy
   * @param border the border
   */
  public void convolve(float[] data, final int maxx, final int maxy, int border) {
    if (border < 0) {
      border = 0;
    } else {
      final int size = 2 * border + 1;
      if (size > maxx || size > maxy) {
        return; // Nothing to do
      }
    }
    convolveInternal(data, maxx, maxy, border);
  }

  /**
   * Compute the convolution.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data the data
   * @param maxx the maxx
   * @param maxy the maxy
   * @param border the border
   */
  private void convolveInternal(float[] data, final int maxx, final int maxy, int border) {
    // Clone the input
    final float[] in = data.clone();

    if (hasWeights()) {
      final int size = data.length;
      if (weights.length != size || this.weightWidth != maxx || this.weightHeight != maxy) {
        throw new IllegalStateException("Weights are not the correct size");
      }

      updateWeightedNormaliser();

      // Apply weights
      for (int i = 0; i < size; i++) {
        in[i] *= weights[i];
      }

      convolveData(in, data, maxx, maxy, border);

      // Normalise. This assumes the roi is only ever constructed as a border.
      if (border != 0) {
        normaliser.normalise(data, data, maxx, maxy, border);
      } else {
        normaliser.normalise(data, size);
      }
    } else {
      convolveData(in, data, maxx, maxy, border);
    }
  }

  /**
   * Convolve the data.
   *
   * @param in the input data (not modified)
   * @param out the output data (modified)
   * @param width the width
   * @param height the height
   * @param border the border
   */
  protected void convolveData(float[] in, float[] out, final int width, final int height,
      int border) {
    final int x1 = border;
    final int y1 = border;
    final int x2 = width - border;
    final int y2 = height - border;
    final int uc = kw / 2;
    final int vc = kh / 2;
    final int xedge = width - uc;
    final int yedge = height - vc;
    for (int y = y1; y < y2; y++) {
      final boolean edgeY = y < vc || y >= yedge;
      for (int x = x1, c = x1 + y * width; x < x2; x++) {
        double sum = 0.0;
        int ki = 0;
        // Determine if at the edge
        if (edgeY || x < uc || x >= xedge) {
          for (int v = -vc; v <= vc; v++) {
            // Create a safe y-index
            int yIndex = y + v;
            if (yIndex < 0) {
              yIndex = 0;
            } else if (yIndex >= height) {
              yIndex = height - 1;
            }
            yIndex *= width;

            for (int u = -uc; u <= uc; u++) {
              // if (i >= kernel.length) // work around for JIT compiler bug on Linux
              // IJ.log("kernel index error: " + i);
              sum += getPixel(x + u, yIndex, in, width) * kernel[ki++];
            }
          }
        } else {
          // Internal
          for (int v = -vc; v <= vc; v++) {
            for (int u = -uc, offset = x - uc + (y + v) * width; u++ <= uc;) {
              sum += in[offset++] * kernel[ki++];
            }
          }
        }
        out[c++] = (float) (sum * scale);
      }
    }
  }

  /**
   * Gets the pixel respecting the image boundaries.
   *
   * @param x the x
   * @param yIndex the y index in the 2D array
   * @param pixels the pixels
   * @param width the width
   * @return the pixel
   */
  private static float getPixel(int x, int yIndex, float[] pixels, int width) {
    if (x < 0) {
      x = 0;
    } else if (x >= width) {
      x = width - 1;
    }
    return pixels[x + yIndex];
  }

  /**
   * Gets the kernal.
   *
   * @return the kernal
   */
  public float[] getKernal() {
    return kernel.clone();
  }

  /**
   * Gets the kernal width.
   *
   * @return the kernal width
   */
  public int getKernalWidth() {
    return kw;
  }

  /**
   * Gets the kernal height.
   *
   * @return the kernal height
   */
  public int getKernalHeight() {
    return kh;
  }

  /**
   * Gets the kernal scale.
   *
   * @return the kernal scale
   */
  public double getKernalScale() {
    return scale;
  }

  /**
   * Rotate 180 degrees around the middle.
   *
   * @param kernel the kernel
   */
  public static void rotate180(float[] kernel) {
    // A yx-packed 2D kernel can be just be reversed
    SimpleArrayUtils.reverse(kernel);
  }

  /**
   * Checks if the number is even. It must be positive.
   *
   * @param number the number
   * @return true if even
   */
  private static boolean isEven(int number) {
    return (number & 1) == 0;
  }

  /**
   * Pad an image processor to make the width and height an odd number.
   *
   * @param fp the fp
   * @return the float processor
   */
  public static FloatProcessor pad(FloatProcessor fp) {
    int kw = fp.getWidth();
    int kh = fp.getHeight();
    // Ensure odd size(to avoid exceptions)
    if ((kw & 1) != 1) {
      kw++;
    }
    if ((kh & 1) != 1) {
      kh++;
    }
    if (kw != fp.getWidth() || kh != fp.getHeight()) {
      final FloatProcessor fp2 = new FloatProcessor(kw, kh);
      fp2.insert(fp, 0, 0);
      fp = fp2;
    }
    return fp;
  }

  /**
   * Checks if is a convolution filter. The default is correlation.
   *
   * @return true, if is a convolution filter
   */
  public boolean isConvolution() {
    return convolution;
  }

  /**
   * Sets the convolution flag. The default is correlation. Changing the type just rotates the
   * kernel 180.
   *
   * @param convolution the new convolution flag
   */
  public void setConvolution(boolean convolution) {
    if (this.convolution != convolution) {
      // Reset the weighted normalisation
      newWeights();
      rotate180(kernel);
    }
    this.convolution = convolution;
  }
}
