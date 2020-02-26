/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.ij.utils;

import ij.ImageStack;
import ij.process.ImageProcessor;
import java.util.Arrays;
import java.util.logging.Logger;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunctionGradient;
import uk.ac.sussex.gdsc.core.math.interpolation.CubicSplinePosition;
import uk.ac.sussex.gdsc.core.math.interpolation.CustomTricubicFunction;
import uk.ac.sussex.gdsc.core.math.interpolation.CustomTricubicFunctionUtils;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.ImageWindow;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.function.cspline.CubicSplineCalculator;
import uk.ac.sussex.gdsc.smlm.math3.optim.PositionChecker;
import uk.ac.sussex.gdsc.smlm.math3.optim.nonlinear.scalar.gradient.BfgsOptimizer;

/**
 * Perform 3D image alignment using normalised cross-correlation.
 *
 * <p>Uses the following formula:
 *
 * <pre>
 *  ( Σ xiyi - nx̄ӯ ) / ( (Σ xi^2 - nx̄^2) (Σ yi^2 - nӯ^2) )^0.5
 * </pre>
 *
 * <p>The summation in the numerator is computed using a conjugate multiplication in the frequency
 * domain. The summation terms are computed using rolling sum tables. Images are converted to the
 * full range of an unsigned 16-bit integer before computation to avoid errors in the rolling sum
 * tables. This should have minimal impact on the correlation value since it is normalised.
 *
 * @see <a href=
 *      "https://en.wikipedia.org/wiki/Pearson_correlation_coefficient">https://en.wikipedia.org/wiki/Pearson_correlation_coefficient</a>
 * @see <a href="http://scribblethink.org/Work/nvisionInterface/nip.html">Fast Normalized
 *      Cross-Correlation by J.P. Lewis</a>
 */
public class Image3DAligner {
  /**
   * The limit for the range of the data as an integer.
   *
   * <p>When this it too high the sumXy from the DHT conjugate multiplication does not match the sum
   * from correlation in the spatial domain.
   *
   * <p>In theory the largest sumXy should be 2^bits * 2^bits * max integer (the size of the largest
   * array). 10-bit integer: 2^10 * 2^10 * 2^31 = 2^51. This is smaller than the mantissa of a
   * double (2^52) so should be represented correctly.
   */
  private static final double LIMIT = 1024;

  private static final int X = 0;
  private static final int XX = 1;
  private static final int Y = 0;
  private static final int YY = 1;

  // For optimisation

  /** Do not have a maximum evaluations as we will converge using the refinements parameter. */
  private static final MaxEval maxEvaluations = new MaxEval(Integer.MAX_VALUE);

  /** The bounds of the spline are always 0-1. */
  private static final SimpleBounds bounds =
      new SimpleBounds(new double[3], SimpleArrayUtils.newDoubleArray(3, 1));

  /** Set a maximum step length at 1 pixel scaled to the spline dimensions. */
  private static final BfgsOptimizer.StepLength stepLength =
      new BfgsOptimizer.StepLength(SimpleArrayUtils.newDoubleArray(3, 1.0 / 3));
  /**
   * This is the cut-off for the maximum gradient relative to the function value. When gradients are
   * too small then the optimisation will end.
   */
  private static final BfgsOptimizer.GradientTolerance gradientTolerance =
      new BfgsOptimizer.GradientTolerance(1e-6);

  /**
   * The search mode for sub-pixel refinement.
   */
  public enum SearchMode {
    /**
     * Perform a binary search by condensing the cube vertices around the highest value of a
     * tricubic interpolation.
     */
    BINARY,
    /** Use the local gradient of a tricubic interpolation to find the maximum. */
    GRADIENT
  }

  private double edgeWindow;
  private double relativeThreshold = 1e-6;
  private SearchMode searchMode = SearchMode.GRADIENT;
  private boolean checkCorrelation = true;
  private double minimumOverlap = 0.5;
  private double minimumDimensionOverlap = 0.75;
  private boolean fastMultiply = true;

  /** The number of slices (max z) of the discrete Hartley transform. */
  private int ns;
  /** The number of rows (max y) of the discrete Hartley transform. */
  private int nr;
  /** The number of columns (max x) of the discrete Hartley transform. */
  private int nc;
  /** The number of rows by columns of the discrete Hartley transform. */
  private int nrByNc;

  private DhtData reference;

  // Not thread safe as they are used for the target image
  private DhtData target;
  private double[] buffer;
  private double[] region;
  private double frequencyDomainCorrelationError;
  private int[] cropDimensions;

  // Allow cached window weights
  private double[] wx;
  private double[] wy;
  private double[] wz;

  private CubicSplineCalculator calc;

  private class DhtData {
    DoubleDht3D dht;
    double[] input;
    double[] sum;
    double[] sumSq;
    // Original dimensions and 3D size
    int width;
    int height;
    int depth;
    int size;
    // Insert position
    int ix;
    int iy;
    int iz;

    DhtData(DoubleDht3D dht, int width, int height, int depth) {
      setDht(dht, width, height, depth);
    }

    void setDht(DoubleDht3D dht, int width, int height, int depth) {
      this.dht = dht;
      sum = resize(sum);
      sumSq = resize(sumSq);
      this.width = width;
      this.height = height;
      this.depth = depth;
      size = width * height * depth;
      ix = getInsert(dht.nc, width);
      iy = getInsert(dht.nr, height);
      iz = getInsert(dht.ns, depth);
      // Make storage of the original data optional. It is just used for
      // the spatial domain correlation check
      if (isCheckCorrelation()) {
        input = resize(input);
      }
    }

    private double[] resize(double[] data) {
      return (data == null || data.length != dht.getDataLength()) ? new double[dht.getDataLength()]
          : data;
    }
  }

  /**
   * Instantiates a new image aligner with a default edge window of 0.25
   */
  public Image3DAligner() {
    this(0.25);
  }

  /**
   * Instantiates a new image aligner.
   *
   * @param edgeWindow the alpha value for the Tukey edge window
   */
  public Image3DAligner(double edgeWindow) {
    setEdgeWindow(edgeWindow);
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected Image3DAligner(Image3DAligner source) {
    edgeWindow = source.edgeWindow;
    relativeThreshold = source.relativeThreshold;
    searchMode = source.searchMode;
    checkCorrelation = source.checkCorrelation;
    minimumOverlap = source.minimumOverlap;
    minimumDimensionOverlap = source.minimumDimensionOverlap;
    fastMultiply = source.fastMultiply;
    ns = source.ns;
    nr = source.nr;
    nc = source.nc;
    nrByNc = source.nrByNc;
    reference = source.reference;
    wx = source.wx;
    wy = source.wy;
    wz = source.wz;
  }

  /**
   * Copy the aligner. This copies the initialised state for use in alignment on multiple threads
   * concurrently.
   *
   * @return the image aligner
   */
  public Image3DAligner copy() {
    return new Image3DAligner(this);
  }

  /**
   * Sets the reference image and assumes the target image will be the same size.
   *
   * <p>The dimension are converted to the next power of 2 for speed. The combined size must fit
   * within the maximum size of a single array.
   *
   * @param image the image (destructively modified)
   * @throws IllegalArgumentException If any dimension is less than 2, or if the combined target
   *         dimensions is too large for an array
   */
  public void setReference(ImageStack image) {
    setReference(image, image.getWidth(), image.getHeight(), image.getSize());
  }

  /**
   * Sets the reference image and assumes the target image will be the same size.
   *
   * <p>The dimension are converted to the next power of 2 for speed. The combined size must fit
   * within the maximum size of a single array.
   *
   * @param image the image (destructively modified)
   * @throws IllegalArgumentException If any dimension is less than 2, or if the combined target
   *         dimensions is too large for an array
   */
  public void setReference(Image3D image) {
    setReference(image, image.getWidth(), image.getHeight(), image.getSize());
  }

  /**
   * Sets the reference image and the size of the target image.
   *
   * <p>The dimension are converted to the next power of 2 for speed. The combined size must fit
   * within the maximum size of a single array.
   *
   * @param image the image (may be destructively modified)
   * @param width the width of the target image
   * @param height the height of the target image
   * @param depth the depth of the target image
   * @throws IllegalArgumentException If any dimension is less than 2, or if the combined target
   *         dimensions is too large for an array
   */
  public void setReference(ImageStack image, int width, int height, int depth) {
    check3D(image);
    if (width < 2 || height < 2 || depth < 2) {
      throw new IllegalArgumentException("Require a 3D target image");
    }
    nc = MathUtils.nextPow2(Math.max(width, image.getWidth()));
    nr = MathUtils.nextPow2(Math.max(height, image.getHeight()));
    ns = MathUtils.nextPow2(Math.max(depth, image.getSize()));
    // Check the image will fit in an Image3D
    Image3D.checkSize(nc, nr, ns, true);
    nrByNc = nr * nc;
    // Window and pad the reference
    setReference(createDht(image, reference));
  }


  /**
   * Sets the reference image and the size of the target image.
   *
   * <p>The dimension are converted to the next power of 2 for speed. The combined size must fit
   * within the maximum size of a single array.
   *
   * @param image the image (may be destructively modified)
   * @param width the width of the target image
   * @param height the height of the target image
   * @param depth the depth of the target image
   * @throws IllegalArgumentException If any dimension is less than 2, or if the combined target
   *         dimensions is too large for an array
   */
  public void setReference(Image3D image, int width, int height, int depth) {
    check3D(image);
    if (width < 2 || height < 2 || depth < 2) {
      throw new IllegalArgumentException("Require a 3D target image");
    }
    nc = MathUtils.nextPow2(Math.max(width, image.getWidth()));
    nr = MathUtils.nextPow2(Math.max(height, image.getHeight()));
    ns = MathUtils.nextPow2(Math.max(depth, image.getSize()));
    nrByNc = nr * nc;
    // Window and pad the reference
    setReference(createDht(image, reference));
  }

  private void setReference(DhtData dhtData) {
    reference = dhtData;
    if (fastMultiply) {
      reference.dht.initialiseFastMultiply();
    }
  }

  private static void check3D(ImageStack image) {
    if (image.getWidth() < 2 || image.getHeight() < 2 || image.getSize() < 2) {
      throw new IllegalArgumentException("Require a 3D image");
    }
    // Check for data
    final int size = image.getWidth() * image.getHeight();
    for (int s = 1; s <= image.getSize(); s++) {
      final ImageProcessor ip = image.getProcessor(s);
      for (int i = 0; i < size; i++) {
        if (ip.getf(i) != 0) {
          return;
        }
      }
    }
    throw new IllegalArgumentException("No data in 3D image");
  }

  private static void check3D(Image3D image) {
    if (image.getWidth() < 2 || image.getHeight() < 2 || image.getSize() < 2) {
      throw new IllegalArgumentException("Require a 3D image");
    }
    // Check for data
    for (int i = 0, size = image.getDataLength(); i < size; i++) {
      if (image.get(i) != 0) {
        return;
      }
    }
    throw new IllegalArgumentException("No data in 3D image");
  }

  private DhtData createDht(ImageStack image, DhtData dhtData) {
    if (image.getBitDepth() != 32) {
      return createDht(new FloatImage3D(image), dhtData);
    }

    // Shift mean to 0 with optional window
    final int width = image.getWidth();
    final int height = image.getHeight();
    final int depth = image.getSize();
    final double[] lwx = createXWindow(width);
    final double[] lwy = createYWindow(height);
    final double[] lwz = createZWindow(depth);

    // We need to compute the weighted centre
    final double[] sum = new double[2];

    for (int z = 0; z < depth; z++) {
      final float[] pixels = (float[]) image.getPixels(1 + z);
      if (lwz[z] == 0) {
        // Special case happens with Tukey window at the ends
      } else {
        calculateWeightedCentre(pixels, width, height, lwx, lwy, lwz[z], sum);
      }
    }

    final double shift = sum[0] / sum[1];

    for (int z = 0; z < depth; z++) {
      final float[] pixels = (float[]) image.getPixels(1 + z);
      if (lwz[z] == 0) {
        // Special case happens with Tukey window at the ends
        Arrays.fill(pixels, 0);
      } else {
        applyWindow(pixels, width, height, lwx, lwy, lwz[z], shift);
      }
    }

    DoubleDht3D dht;

    // Pad into the desired data size.
    // We always do this so the data is reused
    final int size = ns * nr * nc;
    double[] dest;
    if (dhtData == null || dhtData.dht.getDataLength() != size) {
      dest = new double[size];
    } else {
      // Re-use space
      dest = dhtData.dht.getData();
      Arrays.fill(dest, 0);
    }
    dht = new DoubleDht3D(nc, nr, ns, dest, false);
    final int ix = getInsert(nc, width);
    final int iy = getInsert(nr, height);
    final int iz = getInsert(ns, depth);
    dht.insert(ix, iy, iz, image);

    if (dhtData == null) {
      dhtData = new DhtData(dht, width, height, depth);
    } else {
      dhtData.setDht(dht, width, height, depth);
    }

    return prepareDht(dhtData);
  }

  private DhtData createDht(Image3D image, DhtData dhtData) {
    // Shift mean to 0 with optional window
    final int width = image.getWidth();
    final int height = image.getHeight();
    final int depth = image.getSize();
    final double[] lwx = createXWindow(width);
    final double[] lwy = createYWindow(height);
    final double[] lwz = createZWindow(depth);
    final int inc = image.nrByNc;

    // We need to compute the weighted centre
    final double[] sum = new double[2];

    for (int z = 0, i = 0; z < depth; z++) {
      if (lwz[z] == 0) {
        // Special case happens with Tukey window at the ends
      } else {
        calculateWeightedCentre(image, i, width, height, lwx, lwy, lwz[z], sum);
      }
      i += inc;
    }

    final double shift = sum[0] / sum[1];

    for (int z = 0, i = 0; z < depth; z++) {
      if (lwz[z] == 0) {
        // Special case happens with Tukey window at the ends
        for (int j = 0; j < inc; j++) {
          image.set(i++, 0);
        }
      } else {
        applyWindow(image, i, width, height, lwx, lwy, lwz[z], shift);
        i += inc;
      }
    }

    DoubleDht3D dht;

    // Pad into the desired data size.
    // We always do this to handle input of float/double Image3D data.
    final int size = ns * nr * nc;
    double[] dest;
    if (dhtData == null || dhtData.dht.getDataLength() != size) {
      dest = new double[size];
    } else {
      // Re-use space
      dest = dhtData.dht.getData();
      Arrays.fill(dest, 0);
    }
    dht = new DoubleDht3D(nc, nr, ns, dest, false);
    final int ix = getInsert(nc, width);
    final int iy = getInsert(nr, height);
    final int iz = getInsert(ns, depth);
    dht.insert(ix, iy, iz, image);

    if (dhtData == null) {
      dhtData = new DhtData(dht, width, height, depth);
    } else {
      dhtData.setDht(dht, width, height, depth);
    }

    return prepareDht(dhtData);
  }

  private double[] createXWindow(int n) {
    wx = createWindow(wx, n);
    return wx;
  }

  private double[] createYWindow(int n) {
    wy = createWindow(wy, n);
    return wy;
  }

  private double[] createZWindow(int n) {
    wz = createWindow(wz, n);
    return wz;
  }

  private double[] createWindow(double[] window, int n) {
    if (window == null || window.length != n) {
      return ImageWindow.tukey(n, edgeWindow);
    }
    return window;
  }

  private static void calculateWeightedCentre(float[] image, int maxx, int maxy, double[] wx,
      double[] wy, double wz, double[] sum) {
    calculateWeightedCentre(image, 0, maxx, maxy, wx, wy, wz, sum);
  }

  private static void calculateWeightedCentre(float[] image, int index, int maxx, int maxy,
      double[] wx, double[] wy, double wz, double[] sum) {
    for (int y = 0, i = index; y < maxy; y++) {
      final double wyz = wy[y] * wz;
      for (int x = 0; x < maxx; x++, i++) {
        final double weight = wx[x] * wyz;
        sum[0] += image[i] * weight;
        sum[1] += weight;
      }
    }
  }

  private static void calculateWeightedCentre(Image3D image, int index, int maxx, int maxy,
      double[] wx, double[] wy, double wz, double[] sum) {
    for (int y = 0, i = index; y < maxy; y++) {
      final double wyz = wy[y] * wz;
      for (int x = 0; x < maxx; x++, i++) {
        final double weight = wx[x] * wyz;
        sum[0] += image.get(i) * weight;
        sum[1] += weight;
      }
    }
  }

  private static void applyWindow(float[] image, int maxx, int maxy, double[] wx, double[] wy,
      double wz, double shift) {
    applyWindow(image, 0, maxx, maxy, wx, wy, wz, shift);
  }

  private static void applyWindow(float[] image, int index, int maxx, int maxy, double[] wx,
      double[] wy, double wz, double shift) {
    for (int y = 0, i = index; y < maxy; y++) {
      final double wyz = wy[y] * wz;
      for (int x = 0; x < maxx; x++, i++) {
        image[i] = (float) ((image[i] - shift) * wx[x] * wyz);
      }
    }
  }

  private static void applyWindow(Image3D image, int index, int maxx, int maxy, double[] wx,
      double[] wy, double wz, double shift) {
    for (int y = 0, i = index; y < maxy; y++) {
      final double wyz = wy[y] * wz;
      for (int x = 0; x < maxx; x++, i++) {
        image.set(i, (image.get(i) - shift) * wx[x] * wyz);
      }
    }
  }

  private static int getInsert(int maxN, int n) {
    // Note the FHT power spectrum centre is at n/2 of an even sized image.
    // So we must insert the centre at that point. To do this we check for odd/even
    // and offset if necessary.
    final int diff = maxN - n;
    return ((diff & 1) == 1) ? (diff + 1) / 2 : diff / 2;
  }

  /**
   * Prepare the DHT.
   *
   * <p>Converts the data to quantised data. Any zero value (from padding or weighting) remains
   * zero.
   *
   * <p>This may reduce the precision slightly but allows the computation of a rolling sum table
   * have minimal errors. The rolling sum and sum-of-squares table is computed and the DHT is
   * transformed to the frequency domain.
   *
   * @param dhtData the dht data
   * @return the DHT data
   */
  private static DhtData prepareDht(DhtData dhtData) {
    final DoubleDht3D dht = dhtData.dht;
    final double[] sum = dhtData.sum;
    final double[] sumSq = dhtData.sumSq;

    // Note previous versions converted to 10-bit integer data. However the 3D DHT creates very
    // large
    // output values and errors occurred when computing the conjugate multiple in the frequency
    // domain verses the spatial domain. A check has been added to compute the spatial domain
    // correlation for the corresponding max correlation in the frequency domain. This allow
    // the code to report when the correlation value is incorrect.

    final double[] data = dht.getData();
    final double[] limits = MathUtils.limits(data);
    final double min = limits[0];
    final double max = limits[1];

    // Note: The image has been shifted to a mean of 0 so that zero padding
    // for frequency domain transform does not add any information.
    // We need to maintain the sign information and ensure that zero is still
    // zero.

    final double scale = LIMIT / (max - min);

    // Compute the rolling sum tables
    final int nrByNc = dht.nrByNc;
    final int nc = dht.nc;
    final int nr = dht.nr;
    final int ns = dht.ns;

    // This has been adapted from Image3D to compute two rolling sum table at once

    // First build a table for each XY slice
    for (int s = 0; s < ns; s++) {
      double sum1 = 0;
      double sum2 = 0;
      int index = s * nrByNc;
      // Initialise first row sum
      // sum = rolling sum of (0 - column)
      for (int c = 0; c < nc; c++, index++) {
        final double v = transform(data[index], scale);
        data[index] = v;
        sum1 += v;
        sum2 += v * v;
        sum[index] = sum1;
        sumSq[index] = sum2;
      }
      // Remaining rows
      // sum = rolling sum of (0 - column) + sum of same position above
      for (int r = 1, ii = index - nc; r < nr; r++) {
        sum1 = 0;
        sum2 = 0;
        for (int c = 0; c < nc; c++, index++, ii++) {
          final double v = transform(data[index], scale);
          data[index] = v;
          sum1 += v;
          sum2 += v * v;
          // Add the sum from the previous row
          sum[index] = sum1 + sum[ii];
          sumSq[index] = sum2 + sumSq[ii];
        }
      }
    }

    // Now sum across slices
    // sum = rolling sum of (0,0 to column,row) + sum of same position above
    // => rolling sum of (0,0,0 to column,row,slice)
    for (int s = 1; s < ns; s++) {
      int i1 = s * nrByNc;
      int i2 = i1 - nrByNc;
      for (int j = 0; j < nrByNc; j++, i1++, i2++) {
        sum[i1] += sum[i2];
        sumSq[i1] += sumSq[i2];
      }
    }

    // Store after numerical transform
    if (dhtData.input != null && dhtData.input.length == sumSq.length) {
      System.arraycopy(dht.getData(), 0, dhtData.input, 0, sumSq.length);
    }

    // Transform the data
    dht.transform();
    return dhtData;
  }

  private static double transform(double value, double scale) {
    // Ensure zero is zero
    if (value == 0.0) {
      return 0.0;
    }

    // Maintain the sign information
    return Math.round(value * scale); // / scale
  }



  /**
   * Align the image with the reference. Compute the translation required to move the target image
   * onto the reference image for maximum correlation.
   *
   * @param image the image
   * @return [x,y,z,value]
   * @throws IllegalArgumentException If any dimension is less than 2, or if larger than the
   *         initialised reference
   */
  public double[] align(ImageStack image) {
    return align(image, 0, 0);
  }

  /**
   * Align the image with the reference with sub-pixel accuracy. Compute the translation required to
   * move the target image onto the reference image for maximum correlation.
   *
   * <p>Refinement uses a default sub-pixel accuracy of 1e-2;
   *
   * @param image the image
   * @param refinements the refinements for sub-pixel accuracy
   * @return [x,y,z,value]
   * @throws IllegalArgumentException If any dimension is less than 2, or if larger than the
   *         initialised reference
   */
  public double[] align(ImageStack image, int refinements) {
    check3D(image);
    final int width = image.getWidth();
    final int height = image.getHeight();
    final int depth = image.getSize();
    if (width > nc || height > nr || depth > ns) {
      throw new IllegalArgumentException("Image is larger than the initialised reference");
    }

    target = createDht(image, target);
    return align(target, refinements, 1e-2);
  }

  /**
   * Align the image with the reference with sub-pixel accuracy. Compute the translation required to
   * move the target image onto the reference image for maximum correlation.
   *
   * @param image the image
   * @param refinements the refinements for sub-pixel accuracy
   * @param error the error for sub-pixel accuracy (i.e. stop when improvements are less than this
   *        error)
   * @return [x,y,z,value]
   * @throws IllegalArgumentException If any dimension is less than 2, or if larger than the
   *         initialised reference
   */
  public double[] align(ImageStack image, int refinements, double error) {
    check3D(image);
    final int width = image.getWidth();
    final int height = image.getHeight();
    final int depth = image.getSize();
    if (width > nc || height > nr || depth > ns) {
      throw new IllegalArgumentException("Image is larger than the initialised reference");
    }

    target = createDht(image, target);
    return align(target, refinements, error);
  }

  /**
   * Align the image with the reference. Compute the translation required to move the target image
   * onto the reference image for maximum correlation.
   *
   * @param image the image
   * @return [x,y,z,value]
   * @throws IllegalArgumentException If any dimension is less than 2, or if larger than the
   *         initialised reference
   */
  public double[] align(Image3D image) {
    return align(image, 0, 0);
  }

  /**
   * Align the image with the reference with sub-pixel accuracy. Compute the translation required to
   * move the target image onto the reference image for maximum correlation.
   *
   * <p>Refinement uses a default sub-pixel accuracy of 1e-2;
   *
   * @param image the image
   * @param refinements the refinements for sub-pixel accuracy
   * @return [x,y,z,value]
   * @throws IllegalArgumentException If any dimension is less than 2, or if larger than the
   *         initialised reference
   */
  public double[] align(Image3D image, int refinements) {
    check3D(image);
    final int width = image.getWidth();
    final int height = image.getHeight();
    final int depth = image.getSize();
    if (width > nc || height > nr || depth > ns) {
      throw new IllegalArgumentException("Image is larger than the initialised reference");
    }

    target = createDht(image, target);
    return align(target, refinements, 1e-2);
  }

  /**
   * Align the image with the reference with sub-pixel accuracy. Compute the translation required to
   * move the target image onto the reference image for maximum correlation.
   *
   * @param image the image
   * @param refinements the maximum number of refinements for sub-pixel accuracy
   * @param error the error for sub-pixel accuracy (i.e. stop when improvements are less than this
   *        error)
   * @return [x,y,z,value]
   * @throws IllegalArgumentException If any dimension is less than 2, or if larger than the
   *         initialised reference
   */
  public double[] align(Image3D image, int refinements, double error) {
    check3D(image);
    final int width = image.getWidth();
    final int height = image.getHeight();
    final int depth = image.getSize();
    if (width > nc || height > nr || depth > ns) {
      throw new IllegalArgumentException("Image is larger than the initialised reference");
    }

    target = createDht(image, target);
    return align(target, refinements, error);
  }

  /**
   * Align the image with the reference with sub-pixel accuracy. Compute the translation required to
   * move the target image onto the reference image for maximum correlation.
   *
   * @param target the target
   * @param refinements the maximum number of refinements for sub-pixel accuracy
   * @param error the error for sub-pixel accuracy (i.e. stop when improvements are less than this
   *        error)
   * @return [x,y,z,value]
   * @throws IllegalArgumentException If any dimension is less than 2, or if larger than the
   *         initialised reference
   */
  private double[] align(DhtData target, int refinements, double error) {
    // Multiply by the reference. This allows the reference to be shared across threads.
    final DoubleDht3D correlation = target.dht.conjugateMultiply(reference.dht, buffer);
    buffer = correlation.getData(); // Store for reuse
    correlation.inverseTransform();
    correlation.swapOctants();

    // Normalise:
    // ( Σ xiyi - nx̄ӯ ) / ( (Σ xi^2 - nx̄^2) (Σ yi^2 - nӯ^2) )^0.5
    //
    // (sumXy - sumX*sumY/n) / sqrt( (sumXx - sumX^2 / n) * (sumYy - sumY^2 / n) )

    // Only do this over the range where at least half the original images overlap,
    // i.e. the insert point of one will be the middle of the other when shifted.
    int ix = Math.min(reference.ix, target.ix);
    int iy = Math.min(reference.iy, target.iy);
    int iz = Math.min(reference.iz, target.iz);
    int ixw = Math.max(reference.ix + reference.width, target.ix + target.width);
    int iyh = Math.max(reference.iy + reference.height, target.iy + target.height);
    int izd = Math.max(reference.iz + reference.depth, target.iz + target.depth);

    if (minimumDimensionOverlap > 0) {
      final double f = (1 - minimumDimensionOverlap) / 2;
      final int ux = (int) (Math.round(Math.min(reference.width, target.width) * f));
      final int uy = (int) (Math.round(Math.min(reference.height, target.height) * f));
      final int uz = (int) (Math.round(Math.min(reference.depth, target.depth) * f));
      ix += ux;
      ixw -= ux;
      iy += uy;
      iyh -= uy;
      iz += uz;
      izd -= uz;
    }

    cropDimensions = new int[] {ix, iy, iz, ixw - ix, iyh - iy, izd - iz};

    // The maximum correlation unnormalised. Since this is unnormalised
    // it will be biased towards the centre of the image. This is used
    // to restrict the bounds for finding the maximum of the normalised correlation
    // which should be close to this.
    int maxi = correlation.findMaxIndex(ix, iy, iz, cropDimensions[3], cropDimensions[4],
        cropDimensions[5]);

    // Check in the spatial domain
    checkCorrelation(target, correlation, maxi);

    // Compute sum from rolling sum using:
    // sum(x,y,z,w,h,d) =
    // + s(x+w-1,y+h-1,z+d-1)
    // - s(x-1,y+h-1,z+d-1)
    // - s(x+w-1,y-1,z+d-1)
    // + s(x-1,y-1,z+d-1)
    // /* Image above must be subtracted so reverse sign*/
    // - s(x+w-1,y+h-1,z-1)
    // + s(x-1,y+h-1,z-1)
    // + s(x+w-1,y-1,z-1)
    // - s(x-1,y-1,z-1)
    // Note:
    // s(i,j,k) = 0 when either i,j,k < 0
    // i = imax when i>imax
    // j = jmax when j>jmax
    // k = kmax when k>kmax

    // Note: The correlation is for the movement of the reference over the target
    final int nc_2 = nc / 2;
    final int nr_2 = nr / 2;
    final int ns_2 = ns / 2;
    final int[] centre = new int[] {nc_2, nr_2, ns_2};

    // Compute the shift from the centre
    final int dx = nc_2 - ix;
    final int dy = nr_2 - iy;
    final int dz = ns_2 - iz;

    // For the reference (moved -dx,-dy,-dz over the target)
    int rx = -dx;
    int ry = -dy;
    int rz = -dz;

    // For the target (moved dx,dy,dz over the reference)
    int tx = dx;
    int ty = dy;
    int tz = dz;

    // Precompute the x-1,x+w-1,y-1,y+h-1
    final int nx = cropDimensions[3];
    final int[] rx1 = new int[nx];
    final int[] rxw1 = new int[nx];
    final int[] tx1 = new int[nx];
    final int[] txw1 = new int[nx];
    final int[] width = new int[nx];
    for (int c = ix, i = 0; c < ixw; c++, i++) {
      rx1[i] = Math.max(-1, rx - 1);
      rxw1[i] = Math.min(nc, rx + nc) - 1;
      rx++;
      tx1[i] = Math.max(-1, tx - 1);
      txw1[i] = Math.min(nc, tx + nc) - 1;
      tx--;
      width[i] = rxw1[i] - rx1[i];
    }
    final int ny = cropDimensions[4];
    final int[] ry1 = new int[ny];
    final int[] ryh1 = new int[ny];
    final int[] ty1 = new int[ny];
    final int[] tyh1 = new int[ny];
    final int[] h = new int[ny];
    for (int r = iy, j = 0; r < iyh; r++, j++) {
      ry1[j] = Math.max(-1, ry - 1);
      ryh1[j] = Math.min(nr, ry + nr) - 1;
      ry++;
      ty1[j] = Math.max(-1, ty - 1);
      tyh1[j] = Math.min(nr, ty + nr) - 1;
      ty--;
      h[j] = ryh1[j] - ry1[j];
    }

    final double[] rs = reference.sum;
    final double[] rss = reference.sumSq;
    final double[] ts = target.sum;
    final double[] tss = target.sumSq;
    final double[] rsum = new double[2];
    final double[] tsum = new double[2];

    final int size = Math.min(reference.size, target.size);
    final int minimumN = (int) (Math.round(size * minimumOverlap));
    int maxj = -1;
    double max = 0;

    for (int s = iz; s < izd; s++) {
      // Compute the z-1,z+d-1
      final int rz_1 = Math.max(-1, rz - 1);
      final int rz_d_1 = Math.min(ns, rz + ns) - 1;
      rz++;
      final int tz_1 = Math.max(-1, tz - 1);
      final int tz_d_1 = Math.min(ns, tz + ns) - 1;
      tz--;
      final int d = rz_d_1 - rz_1;

      for (int r = iy, j = 0; r < iyh; r++, j++) {
        final int base = s * nrByNc + r * nc;
        final int hd = h[j] * d;
        for (int c = ix, i = 0; c < ixw; c++, i++) {
          final double sumXy = buffer[base + c];

          compute(rx1[i], ry1[j], rz_1, rxw1[i], ryh1[j], rz_d_1, width[i], h[j], d, rs, rss, rsum);
          compute(tx1[i], ty1[j], tz_1, txw1[i], tyh1[j], tz_d_1, width[i], h[j], d, ts, tss, tsum);

          // Compute the correlation
          // (sumXy - sumX*sumY/n) / sqrt( (sumXx - sumX^2 / n) * (sumYy - sumY^2 / n) )

          final int n = width[i] * hd;
          final double numerator = sumXy - (rsum[X] * tsum[Y] / n);
          final double denominator1 = rsum[XX] - (rsum[X] * rsum[X] / n);
          final double denominator2 = tsum[YY] - (tsum[Y] * tsum[Y] / n);

          double corr;
          if (denominator1 == 0 || denominator2 == 0) {
            // If there is data and all the variances are the same then correlation is perfect
            if (rsum[XX] == tsum[YY] && rsum[XX] == sumXy && rsum[XX] > 0) {
              corr = 1;
            } else {
              corr = 0;
            }
          } else {
            // Leave as raw for debugging, i.e. do not clip to range [-1:1]
            corr = numerator / Math.sqrt(denominator1 * denominator2);
          }

          buffer[base + c] = corr;

          if (n < minimumN) {
            continue;
          }

          // Check normalisation with some margin for error
          if (corr > 1.0001) {
            // Normalisation has failed.
            // This occurs when the correlation sum XY is incorrect.
            // The other terms are exact due to the quantisation to integer data.
            // It is likely to occur at the bounds.
            continue;
          }

          if (corr > max) {
            max = corr;
            maxj = base + c;
          } else if (corr == max) {
            // Get shift from centre
            final int[] xyz1 = correlation.getXyz(maxj);
            final int[] xyz2 = correlation.getXyz(base + c);
            int d1 = 0;
            int d2 = 0;
            for (int k = 0; k < 3; k++) {
              d1 += MathUtils.pow2(xyz1[k] - centre[k]);
              d2 += MathUtils.pow2(xyz2[k] - centre[k]);
            }
            if (d2 < d1) {
              max = corr;
              maxj = base + c;
            }
          }
        }
      }
    }

    // The maximum correlation with normalisation
    maxi = maxj; // correlation.findMaxIndex(ix, iy, iz, iw - ix, ih - iy, id - iz);
    final int[] xyz = correlation.getXyz(maxi);

    // Report the shift required to move from the centre of the target image to the reference
    // @formatter:off
    final double[] result = new double[] {
      nc_2 - xyz[0],
      nr_2 - xyz[1],
      ns_2 - xyz[2],
      buffer[maxi]
    };
    // @formatter:on

    if (refinements > 0) {
      // Perform sub-pixel alignment
      // Create a cubic spline using a small region of pixels around the maximum
      if (calc == null) {
        calc = new CubicSplineCalculator();
      }
      // Avoid out-of-bounds errors. Only use the range that was normalised
      final int x = MathUtils.clip(ix, ixw - 4, xyz[0] - 1);
      final int y = MathUtils.clip(iy, iyh - 4, xyz[1] - 1);
      final int z = MathUtils.clip(iz, izd - 4, xyz[2] - 1);
      final DoubleImage3D crop = correlation.crop(x, y, z, 4, 4, 4, region);
      region = crop.getData();
      final CustomTricubicFunction f = CustomTricubicFunctionUtils.create(calc.compute(region));

      // Find the maximum starting at the current origin
      final int ox = xyz[0] - x;
      final int oy = xyz[1] - y;
      final int oz = xyz[2] - z;

      // Scale to the cubic spline dimensions of 0-1
      final double[] origin = new double[] {ox / 3.0, oy / 3.0, oz / 3.0};

      // Simple condensing search
      if (searchMode == SearchMode.BINARY) {
        // Can this use the current origin as a start point?
        // Currently we evaluate 8-cube vertices. A better search
        // would evaluate 27 points around the optimum, pick the best then condense
        // the range.
        final double[] optimum = f.search(true, refinements, relativeThreshold, -1);
        final double value = optimum[3];
        if (value > result[3]) {
          result[3] = value;
          // Convert the maximum back with scaling
          for (int i = 0; i < 3; i++) {
            result[i] -= (optimum[i] - origin[i]) * 3.0;
          }
          return result;
        }
      } else {
        // Gradient search
        try {
          final SplineFunction sf = new SplineFunction(f, origin);

          final BfgsOptimizer optimiser = new BfgsOptimizer(
              // Use a simple check on the relative value change and
              // set the number of refinements
              new SimpleValueChecker(relativeThreshold, -1, refinements));

          final PointValuePair opt = optimiser.optimize(maxEvaluations, bounds, gradientTolerance,
              stepLength, new InitialGuess(origin),
              // Scale the error for the position check
              new PositionChecker(-1, error / 3.0), new ObjectiveFunction(sf::value),
              new ObjectiveFunctionGradient(point -> {
                // This must be new each time
                final double[] partialDerivative1 = new double[3];
                sf.value(point, partialDerivative1);
                return partialDerivative1;
              }));

          // Check it is higher. Invert since we did a minimisation.
          final double value = -opt.getValue();
          if (value > result[3]) {
            result[3] = value;
            // Convert the maximum back with scaling
            final double[] optimum = opt.getPointRef();
            for (int i = 0; i < 3; i++) {
              result[i] -= (optimum[i] - origin[i]) * 3.0;
            }
            return result;
          }
        } catch (final Exception ex) {
          // Ignore this
        }
      }
    }

    return result;
  }

  /**
   * Check the correlation in the spatial domain verses the maximum correlation in the frequency
   * domain.
   *
   * @param target the target
   * @param correlation the correlation
   * @param maxi the index of the maximum correlation
   */
  private void checkCorrelation(DhtData target, DoubleDht3D correlation, int maxi) {
    if (target.input == null || reference.input == null) {
      // No check possible
      return;
    }

    // The maximum correlation without normalisation
    final int[] xyz = correlation.getXyz(maxi);

    // Find the range for the target and reference
    final int nc_2 = nc / 2;
    final int nr_2 = nr / 2;
    final int ns_2 = ns / 2;
    final int tx = Math.max(0, xyz[0] - nc_2);
    final int ty = Math.max(0, xyz[1] - nr_2);
    final int tz = Math.max(0, xyz[2] - ns_2);
    final int width = Math.min(nc, xyz[0] + nc_2) - tx;
    final int height = Math.min(nr, xyz[1] + nr_2) - ty;
    final int depth = Math.min(ns, xyz[2] + ns_2) - tz;

    // For the reference we express as a shift relative to the centre
    // and subtract the half-width.
    // Formally: (nc_2 - xyz[0]) // shift
    // + nc_2 // centre
    // - nc_2 // Half width
    final int rx = Math.max(0, -xyz[0] + nc_2);
    final int ry = Math.max(0, -xyz[1] + nr_2);
    final int rz = Math.max(0, -xyz[2] + ns_2);

    final double[] tar = target.input;
    final double[] ref = reference.input;
    final double frequencyCorrelation = correlation.get(maxi);
    double spatialCorrelation = 0;
    for (int z = 0; z < depth; z++) {
      for (int y = 0; y < height; y++) {
        int ti = (tz + z) * nrByNc + (ty + y) * nc + tx;
        int ri = (rz + z) * nrByNc + (ry + y) * nc + rx;
        for (int x = 0; x < width; x++) {
          spatialCorrelation += tar[ti++] * ref[ri++];
        }
      }
    }

    frequencyDomainCorrelationError =
        DoubleEquality.relativeError(frequencyCorrelation, spatialCorrelation);
    if (frequencyDomainCorrelationError > 0.05) {
      final double finalSpatialCorrelation = spatialCorrelation;
      Logger.getLogger(getClass().getName())
          .warning(() -> String.format("3D Correlation Error = %s : Spatial = %s, Freq = %s",
              MathUtils.rounded(frequencyDomainCorrelationError), finalSpatialCorrelation,
              frequencyCorrelation));
    }
  }

  /**
   * Compute the sum from the rolling sum tables.
   *
   * @param x1 the x value -1
   * @param y1 the y value -1
   * @param z1 the z value -1
   * @param xw1 the x value +w -1
   * @param yh1 the y value +h -1
   * @param zd1 the z value +d -1
   * @param width the width
   * @param height the height
   * @param depth the depth
   * @param sum the sum table
   * @param sumSq the sum-of-squares table
   * @param result the sum (output = [sum, sum-of-squares])
   */
  private void compute(int x1, int y1, int z1, int xw1, int yh1, int zd1, int width, int height,
      int depth, double[] sum, double[] sumSq, double[] result) {
    // Compute sum from rolling sum using:
    // sum(x,y,z,w,h,d) =
    // + s(x+w-1,y+h-1,z+d-1)
    // - s(x-1,y+h-1,z+d-1)
    // - s(x+w-1,y-1,z+d-1)
    // + s(x-1,y-1,z+d-1)
    // /* Image above must be subtracted so reverse sign*/
    // - s(x+w-1,y+h-1,z-1)
    // + s(x-1,y+h-1,z-1)
    // + s(x+w-1,y-1,z-1)
    // - s(x-1,y-1,z-1)
    // Note:
    // s(i,j,k) = 0 when either i,j,k < 0
    // i = imax when i>imax
    // j = jmax when j>jmax
    // k = kmax when k>kmax

    // This has been adapted from Image3D to compute the twos sums together

    // int xw_yh_zd = reference.dht.getIndex(xw1, yh1, zd1);
    final int xw_yh_zd = zd1 * nrByNc + yh1 * nc + xw1;
    result[0] = 0;
    result[1] = 0;
    if (z1 >= 0) {
      final int xw_yh_z = xw_yh_zd - depth * nrByNc;
      if (y1 >= 0) {
        final int h_ = height * nc;
        if (x1 >= 0) {
          result[0] = sum[xw_yh_zd - width - h_] - sum[xw_yh_z - width - h_] - sum[xw_yh_zd - width]
              + sum[xw_yh_z - width];
          result[1] = sumSq[xw_yh_zd - width - h_] - sumSq[xw_yh_z - width - h_]
              - sumSq[xw_yh_zd - width] + sumSq[xw_yh_z - width];
        }
        result[0] = result[0] + sum[xw_yh_z - h_] - sum[xw_yh_zd - h_];
        result[1] = result[1] + sumSq[xw_yh_z - h_] - sumSq[xw_yh_zd - h_];
      } else if (x1 >= 0) {
        result[0] = sum[xw_yh_z - width] - sum[xw_yh_zd - width];
        result[1] = sumSq[xw_yh_z - width] - sumSq[xw_yh_zd - width];
      }
      result[0] = result[0] + sum[xw_yh_zd] - sum[xw_yh_z];
      result[1] = result[1] + sumSq[xw_yh_zd] - sumSq[xw_yh_z];
    } else {
      if (y1 >= 0) {
        final int h_ = height * nc;
        if (x1 >= 0) {
          result[0] = sum[xw_yh_zd - width - h_] - sum[xw_yh_zd - width];
          result[1] = sumSq[xw_yh_zd - width - h_] - sumSq[xw_yh_zd - width];
        }
        result[0] -= sum[xw_yh_zd - h_];
        result[1] -= sumSq[xw_yh_zd - h_];
      } else if (x1 >= 0) {
        result[0] = -sum[xw_yh_zd - width];
        result[1] = -sumSq[xw_yh_zd - width];
      }
      result[0] = result[0] + sum[xw_yh_zd];
      result[1] = result[1] + sumSq[xw_yh_zd];
    }
  }

  private static class SplineFunction {
    final CustomTricubicFunction function;
    CubicSplinePosition[] sp = new CubicSplinePosition[3];

    SplineFunction(CustomTricubicFunction function, double[] origin) {
      this.function = function;
      for (int i = 0; i < 3; i++) {
        sp[i] = new CubicSplinePosition(origin[i]);
      }
    }

    double value(double[] point) {
      initialise(point);
      // BFGS algorithm minimises so invert
      return -function.value(sp[0], sp[1], sp[2]);
    }

    void value(double[] point, double[] derivative) {
      initialise(point);
      function.value(sp[0], sp[1], sp[2], derivative);
      // BFGS algorithm minimises so invert
      for (int i = 0; i < 3; i++) {
        derivative[i] = -derivative[i];
      }
    }

    void initialise(double[] point) {
      // Allow caching the spline positions
      for (int i = 0; i < 3; i++) {
        if (sp[i].getX() != point[i]) {
          sp[i] = new CubicSplinePosition(point[i]);
        }
      }
    }
  }

  /**
   * Gets the correlation image from the last alignment.
   *
   * @return the correlation (or null)
   */
  public Image3D getCorrelation() {
    try {
      final DoubleImage3D image = new DoubleImage3D(nc, nr, ns, buffer);
      image.fillOutside(cropDimensions[0], cropDimensions[1], cropDimensions[2], cropDimensions[3],
          cropDimensions[5], cropDimensions[5], 0);
      return image;
    } catch (final IllegalArgumentException ex) {
      // Thrown when buffer is null or does not match the dimensions.
      return null;
    }
  }

  /**
   * Gets the frequency domain correlation error from the last correlation.
   *
   * @return the frequency domain correlation error
   */
  public double getFrequencyDomainCorrelationError() {
    return frequencyDomainCorrelationError;
  }

  /**
   * Gets the edge window.
   *
   * @return the edge window
   */
  public double getEdgeWindow() {
    return edgeWindow;
  }

  /**
   * Sets the edge window.
   *
   * @param edgeWindow the new edge window
   */
  public void setEdgeWindow(double edgeWindow) {
    this.edgeWindow = MathUtils.clip(0, 1, edgeWindow);
  }

  /**
   * Gets the relative threshold for change in the correlation value for halting refinement. If this
   * is negative it is disabled.
   *
   * @return the relative threshold
   */
  public double getRelativeThreshold() {
    return relativeThreshold;
  }

  /**
   * Sets the relative threshold for change in the correlation value for halting refinement. Set to
   * negative to disable. Refinement will then only be halted by the number of refinement steps or
   * the position error.
   *
   * @param relativeThreshold the new relative threshold
   */
  public void setRelativeThreshold(double relativeThreshold) {
    this.relativeThreshold = relativeThreshold;
  }

  /**
   * Gets the search mode.
   *
   * @return the search mode
   */
  public SearchMode getSearchMode() {
    return searchMode;
  }

  /**
   * Sets the search mode.
   *
   * @param searchMode the new search mode
   */
  public void setSearchMode(SearchMode searchMode) {
    this.searchMode = searchMode;
  }

  /**
   * Checks if the spatial domain correlation check is enabled.
   *
   * @return true, if the spatial domain correlation check is enabled
   */
  public boolean isCheckCorrelation() {
    return checkCorrelation;
  }

  /**
   * Sets the spatial domain correlation check flag. If true then the original untransformed data
   * will be stored in memory. The point of the highest correlation in the frequency domain will be
   * recomputed in the spatial domain. The error between the two can be returned using
   * {@link #getFrequencyDomainCorrelationError()}.
   *
   * @param checkCorrelation the new check correlation flag
   */
  public void setCheckCorrelation(boolean checkCorrelation) {
    this.checkCorrelation = checkCorrelation;
  }

  /**
   * Gets the minimum overlap between the smaller image and the other image.
   *
   * @return the minimum overlap
   */
  public double getMinimumOverlap() {
    return minimumOverlap;
  }

  /**
   * Sets the minimum overlap between the smaller image and the other image.
   *
   * @param minimumOverlap the new minimum overlap
   */
  public void setMinimumOverlap(double minimumOverlap) {
    this.minimumOverlap = MathUtils.clip(0, 1, minimumOverlap);
  }

  /**
   * Gets the minimum overlap between the smaller image and the other image in each dimension.
   *
   * @return the minimum dimension overlap
   */
  public double getMinimumDimensionOverlap() {
    return minimumDimensionOverlap;
  }

  /**
   * Sets the minimum overlap between the smaller image and the other image in each dimension.
   *
   * @param minimumDimensionOverlap the new minimum dimension overlap
   */
  public void setMinimumDimensionOverlap(double minimumDimensionOverlap) {
    this.minimumDimensionOverlap = MathUtils.clip(0, 1, minimumDimensionOverlap);
  }

  /**
   * Checks if is fast multiply.
   *
   * @return true, if is fast multiply
   */
  public boolean isFastMultiply() {
    return fastMultiply;
  }

  /**
   * Sets the fast multiply flag. This initialises the DHT for multiplication at the cost of extra
   * memory storage. The storage requirements are 2 double arrays and 1 integer array of the same
   * length at the FHT data.
   *
   * @param fastMultiply the new fast multiply flag
   */
  public void setFastMultiply(boolean fastMultiply) {
    this.fastMultiply = fastMultiply;
  }
}
