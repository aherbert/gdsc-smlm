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

package uk.ac.sussex.gdsc.smlm.model;

import uk.ac.sussex.gdsc.smlm.math3.distribution.CustomPoissonDistribution;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.Precision;

import java.util.Arrays;

/**
 * Contains methods for generating models of a Point Spread Function.
 */
public abstract class PSFModel {
  /** The rand.om data generator */
  protected RandomDataGenerator rand;
  private double[] psf;
  private int x0min;
  private int x1min;
  private int x0max;
  private int x1max;
  private int[] samplePositions;

  /**
   * Instantiates a new PSF model.
   */
  public PSFModel() {
    setRandomGenerator(new RandomDataGenerator());
  }

  /**
   * Instantiates a new PSF model.
   *
   * @param randomGenerator the random generator
   */
  public PSFModel(RandomGenerator randomGenerator) {
    setRandomGenerator(randomGenerator);
  }

  /**
   * Instantiates a new PSF model.
   *
   * @param randomDataGenerator the random data generator
   */
  public PSFModel(RandomDataGenerator randomDataGenerator) {
    setRandomGenerator(randomDataGenerator);
  }

  /**
   * Construct a PSF function on the provided data.
   *
   * @param data The data (can be null)
   * @param width The data width
   * @param height The data height
   * @param sum The integral
   * @param x0 The centre in dimension 0
   * @param x1 The centre in dimension 1
   * @param x2 The centre in dimension 2
   * @param poissonNoise Add Poisson noise
   * @return The total sum added to the image (useful when poissonNoise is added)
   */
  public abstract double create3D(float[] data, final int width, final int height, final double sum,
      double x0, double x1, double x2, boolean poissonNoise);

  /**
   * Construct a PSF function on the provided data.
   *
   * @param data The data (can be null)
   * @param width The data width
   * @param height The data height
   * @param sum The integral
   * @param x0 The centre in dimension 0
   * @param x1 The centre in dimension 1
   * @param x2 The centre in dimension 2
   * @param poissonNoise Add Poisson noise
   * @return The total sum added to the image (useful when poissonNoise is added)
   */
  public abstract double create3D(double[] data, final int width, final int height,
      final double sum, double x0, double x1, double x2, boolean poissonNoise);

  /**
   * Construct a PSF function on the provided data.
   *
   * @param data The data (can be null)
   * @param width The data width
   * @param height The data height
   * @param sum The integral
   * @param x0 The centre in dimension 0
   * @param x1 The centre in dimension 1
   * @param poissonNoise Add Poisson noise
   * @return The total sum added to the image (useful when poissonNoise is added)
   */
  public double create2D(float[] data, final int width, final int height, final double sum,
      double x0, double x1, boolean poissonNoise) {
    return create3D(data, width, height, sum, x0, x1, 0, poissonNoise);
  }

  /**
   * Construct a PSF function on the provided data.
   *
   * @param data The data (can be null)
   * @param width The data width
   * @param height The data height
   * @param sum The integral
   * @param x0 The centre in dimension 0
   * @param x1 The centre in dimension 1
   * @param poissonNoise Add Poisson noise
   * @return The total sum added to the image (useful when poissonNoise is added)
   */
  public double create2D(double[] data, final int width, final int height, final double sum,
      double x0, double x1, boolean poissonNoise) {
    return create3D(data, width, height, sum, x0, x1, 0, poissonNoise);
  }

  /**
   * Gets the last drawn PSF.
   *
   * @return The last drawn PSF
   */
  public double[] getPSF() {
    return psf;
  }

  /**
   * Gets the minimum position in dimension 0 for the last drawn/sampled PSF.
   *
   * @return The minimum position in dimension 0 for the last drawn/sampled PSF
   */
  public int getX0min() {
    return x0min;
  }

  /**
   * Gets the maximum position in dimension 0 for the last drawn/sampled PSF.
   *
   * @return The maximum position in dimension 0 for the last drawn/sampled PSF
   */
  public int getX0max() {
    return x0max;
  }

  /**
   * Gets the minimum position in dimension 1 for the last drawn/sampled PSF.
   *
   * @return The minimum position in dimension 1 for the last drawn/sampled PSF
   */
  public int getX1min() {
    return x1min;
  }

  /**
   * Gets the maximum position in dimension 1 for the last drawn/sampled PSF.
   *
   * @return The maximum position in dimension 1 for the last drawn/sampled PSF
   */
  public int getX1max() {
    return x1max;
  }

  /**
   * Insert the psf into the data.
   *
   * @param data The input data (width*height)
   * @param x0min The minimum position to insert in dimension 0
   * @param x1min The minimum position to insert in dimension 1
   * @param x0max The maximum position to insert in dimension 0
   * @param x1max The maximum position to insert in dimension 1
   * @param width The width of the input data
   * @param psf The PSF data
   * @param poissonNoise Set to true to add Poisson noise to the PSF
   * @return The sum of the PSF inserted
   */
  protected double insert(float[] data, int x0min, int x1min, int x0max, int x1max, int width,
      double[] psf, boolean poissonNoise) {
    final int x0range = x0max - x0min;
    final int x1range = x1max - x1min;

    if (x0range < 1 || x1range < 1) {
      this.psf = null;
      this.x0min = 0;
      this.x0max = 0;
      this.x1min = 0;
      this.x1max = 0;
      return 0;
    }

    this.psf = psf;
    this.x0min = x0min;
    this.x0max = x0max;
    this.x1min = x1min;
    this.x1max = x1max;

    if (poissonNoise) {
      final CustomPoissonDistribution pd =
          new CustomPoissonDistribution(rand.getRandomGenerator(), 1);
      for (int i = 0; i < psf.length; i++) {
        if (psf[i] > 0) {
          pd.setMeanUnsafe(psf[i]);
          psf[i] = pd.sample();
        }
      }
    }

    // Insert the function into the input data
    for (int y = 0; y < x1range; y++) {
      // Locate the insert location
      int indexTo = (y + x1min) * width + x0min;
      int indexFrom = y * x0range;
      for (int x = 0; x < x0range; x++) {
        data[indexTo++] += psf[indexFrom++];
      }
    }

    double total = 0;
    for (final double d : psf) {
      total += d;
    }
    return total;
  }

  /**
   * Insert the psf into the data.
   *
   * @param data The input data (width*height)
   * @param x0min The minimum position to insert in dimension 0
   * @param x1min The minimum position to insert in dimension 1
   * @param x0max The maximum position to insert in dimension 0
   * @param x1max The maximum position to insert in dimension 1
   * @param width The width of the input data
   * @param psf The PSF data
   * @param poissonNoise Set to true to add Poisson noise to the PSF
   * @return The sum of the PSF inserted
   */
  protected double insert(double[] data, int x0min, int x1min, int x0max, int x1max, int width,
      double[] psf, boolean poissonNoise) {
    final int x0range = x0max - x0min;
    final int x1range = x1max - x1min;

    if (x0range < 1 || x1range < 1) {
      this.psf = null;
      this.x0min = 0;
      this.x0max = 0;
      this.x1min = 0;
      this.x1max = 0;
      return 0;
    }

    this.psf = psf;
    this.x0min = x0min;
    this.x0max = x0max;
    this.x1min = x1min;
    this.x1max = x1max;

    if (poissonNoise) {
      final CustomPoissonDistribution pd =
          new CustomPoissonDistribution(rand.getRandomGenerator(), 1);
      for (int i = 0; i < psf.length; i++) {
        if (psf[i] > 0) {
          pd.setMeanUnsafe(psf[i]);
          psf[i] = pd.sample();
        }
      }
    }

    // Insert the function into the input data
    for (int y = 0; y < x1range; y++) {
      // Locate the insert location
      int indexTo = (y + x1min) * width + x0min;
      int indexFrom = y * x0range;
      for (int x = 0; x < x0range; x++) {
        data[indexTo++] += psf[indexFrom++];
      }
    }

    double total = 0;
    for (final double d : psf) {
      total += d;
    }
    return total;
  }

  /**
   * Remove the last added PSF from the data. This can be invoked after any call to draw a PSF into
   * an input data array.
   *
   * @param data the data
   * @param width the width
   * @param height the height
   */
  public void erase(float[] data, int width, int height) {
    erase(data, width, height, psf, x0min, x0max, x1min, x1max);
  }

  /**
   * Remove the last added PSF from the data. This can be invoked after any call to draw a PSF into
   * an input data array.
   *
   * @param data the data
   * @param width the width
   * @param height the height
   */
  public void erase(double[] data, int width, int height) {
    erase(data, width, height, psf, x0min, x0max, x1min, x1max);
  }

  /**
   * Remove the PSF from the data. Can be invoked using a saved copy of the PSF previously drawn by
   * the model obtained from the appropriate get() methods.
   *
   * @param data the data
   * @param width the width
   * @param height the height
   * @param psf the psf
   * @param x0min the minimum position in dimension 0 for the last drawn/sampled PSF
   * @param x0max the maximum position in dimension 0 for the last drawn/sampled PSF
   * @param x1min the minimum position in dimension 1 for the last drawn/sampled PSF
   * @param x1max the maximum position in dimension 1 for the last drawn/sampled PSF
   */
  public void erase(float[] data, int width, int height, double[] psf, int x0min, int x0max,
      int x1min, int x1max) {
    if (psf == null) {
      return;
    }

    final int x0range = x0max - x0min;
    final int x1range = x1max - x1min;

    // min should always be less than max
    if (x0range < 1) {
      throw new IllegalArgumentException("Dimension 0 range not within data bounds");
    }
    if (x1range < 1) {
      throw new IllegalArgumentException("Dimension 1 range not within data bounds");
    }

    // Remove from the input data
    for (int y = 0; y < x1range; y++) {
      // Locate the insert location
      int indexTo = (y + x1min) * width + x0min;
      int indexFrom = y * x0range;
      for (int x = 0; x < x0range; x++) {
        data[indexTo++] -= psf[indexFrom++];
      }
    }
  }

  /**
   * Remove the PSF from the data. Can be invoked using a saved copy of the PSF previously drawn by
   * the model obtained from the appropriate get() methods.
   *
   * @param data the data
   * @param width the width
   * @param height the height
   * @param psf the psf
   * @param x0min the minimum position in dimension 0 for the last drawn/sampled PSF
   * @param x0max the maximum position in dimension 0 for the last drawn/sampled PSF
   * @param x1min the minimum position in dimension 1 for the last drawn/sampled PSF
   * @param x1max the maximum position in dimension 1 for the last drawn/sampled PSF
   */
  public void erase(double[] data, int width, int height, double[] psf, int x0min, int x0max,
      int x1min, int x1max) {
    if (psf == null) {
      return;
    }

    final int x0range = x0max - x0min;
    final int x1range = x1max - x1min;

    // min should always be less than max
    if (x0range < 1) {
      throw new IllegalArgumentException("Dimension 0 range not within data bounds");
    }
    if (x1range < 1) {
      throw new IllegalArgumentException("Dimension 1 range not within data bounds");
    }

    // Remove from the input data
    for (int y = 0; y < x1range; y++) {
      // Locate the insert location
      int indexTo = (y + x1min) * width + x0min;
      int indexFrom = y * x0range;
      for (int x = 0; x < x0range; x++) {
        data[indexTo++] -= psf[indexFrom++];
      }
    }
  }

  /**
   * Produce a shallow copy of this object. This shares the pre-computed PSF data but will allow the
   * copy to store its own version of the most recently created PSF.
   *
   * @return A shallow copy of this object
   */
  public abstract PSFModel copy();

  /**
   * Sample a PSF function on the provided data.
   *
   * @param data The data (can be null)
   * @param width The data width
   * @param height The data height
   * @param n The number of samples
   * @param x0 The centre in dimension 0
   * @param x1 The centre in dimension 1
   * @param x2 The centre in dimension 2
   * @return The number of samples drawn on the image (useful to detect samples outside the image
   *         bounds)
   */
  public abstract int sample3D(float[] data, final int width, final int height, final int n,
      double x0, double x1, double x2);

  /**
   * Sample a PSF function on the provided data.
   *
   * @param data The data (can be null)
   * @param width The data width
   * @param height The data height
   * @param n The number of samples
   * @param x0 The centre in dimension 0
   * @param x1 The centre in dimension 1
   * @param x2 The centre in dimension 2
   * @return The number of samples drawn on the image (useful to detect samples outside the image
   *         bounds)
   */
  public abstract int sample3D(double[] data, final int width, final int height, final int n,
      double x0, double x1, double x2);

  /**
   * Sample a PSF function on the provided data.
   *
   * @param data The data (can be null)
   * @param width The data width
   * @param height The data height
   * @param n The number of samples
   * @param x0 The centre in dimension 0
   * @param x1 The centre in dimension 1
   * @return The number of samples drawn on the image (useful to detect samples outside the image
   *         bounds)
   */
  public double sample2D(float[] data, final int width, final int height, final int n, double x0,
      double x1) {
    return sample3D(data, width, height, n, x0, x1, 0);
  }

  /**
   * Sample a PSF function on the provided data.
   *
   * @param data The data (can be null)
   * @param width The data width
   * @param height The data height
   * @param n The number of samples
   * @param x0 The centre in dimension 0
   * @param x1 The centre in dimension 1
   * @return The number of samples drawn on the image (useful to detect samples outside the image
   *         bounds)
   */
  public double sample2D(double[] data, final int width, final int height, final int n, double x0,
      double x1) {
    return sample3D(data, width, height, n, x0, x1, 0);
  }

  /**
   * Insert a set of sampled XY positions into the data.
   *
   * @param data The data
   * @param width The data width
   * @param height The data height
   * @param x The x-positions
   * @param y The y-positions
   * @return The number of samples that were inside the data bounds
   */
  @SuppressWarnings("null")
  protected int insertSample(double[] data, final int width, final int height, double[] x,
      double[] y) {
    int n = 0;
    samplePositions = new int[(x == null) ? 0 : x.length];

    x0max = x1max = 0;
    x0min = x1min = width;
    for (int i = 0; i < samplePositions.length; i++) {
      if (x[i] < 0 || x[i] >= width || y[i] < 0 || y[i] >= height) {
        continue;
      }
      final int xp = (int) x[i];
      final int yp = (int) y[i];
      if (x0min > xp) {
        x0min = xp;
      }
      if (x0max < xp) {
        x0max = xp;
      }
      if (x1min > yp) {
        x1min = yp;
      }
      if (x1max < yp) {
        x1max = yp;
      }
      final int index = yp * width + xp;
      samplePositions[n++] = index;
      data[index] += 1;
    }

    if (n < samplePositions.length) {
      samplePositions = Arrays.copyOf(samplePositions, n);
    }
    return samplePositions.length;
  }

  /**
   * Insert a set of sampled XY positions into the data.
   *
   * @param data The data
   * @param width The data width
   * @param height The data height
   * @param x The x-positions
   * @param y The y-positions
   * @return The number of samples that were inside the data bounds
   */
  @SuppressWarnings("null")
  protected int insertSample(float[] data, final int width, final int height, double[] x,
      double[] y) {
    int n = 0;
    samplePositions = new int[(x == null) ? 0 : x.length];

    x0max = x1max = 0;
    x0min = x1min = width;
    for (int i = 0; i < samplePositions.length; i++) {
      if (x[i] < 0 || x[i] >= width || y[i] < 0 || y[i] >= height) {
        continue;
      }
      final int xp = (int) x[i];
      final int yp = (int) y[i];
      if (x0min > xp) {
        x0min = xp;
      }
      if (x0max < xp) {
        x0max = xp;
      }
      if (x1min > yp) {
        x1min = yp;
      }
      if (x1max < yp) {
        x1max = yp;
      }
      final int index = yp * width + xp;
      samplePositions[n++] = index;
      data[index] += 1;
    }

    if (n < samplePositions.length) {
      samplePositions = Arrays.copyOf(samplePositions, n);
    }
    return samplePositions.length;
  }

  /**
   * Return the positions where samples were added to the data. The size of the array should equal
   * the number of samples added by a sample(...) method.
   *
   * @return The positions in the data where samples where added
   */
  public int[] getSamplePositions() {
    return samplePositions;
  }

  /**
   * Remove the last added PSF from the data. This can be invoked after any call to sample a PSF
   * into an input data array.
   *
   * @param data the data
   * @param width the width
   * @param height the height
   */
  public void eraseSample(float[] data, int width, int height) {
    eraseSample(data, width, height, samplePositions);
  }

  /**
   * Remove the last added PSF from the data. This can be invoked after any call to sample a PSF
   * into an input data array.
   *
   * @param data the data
   * @param width the width
   * @param height the height
   */
  public void eraseSample(double[] data, int width, int height) {
    eraseSample(data, width, height, samplePositions);
  }

  /**
   * Remove the PSF from the data. Can be invoked using a saved copy of the PSF previously drawn by
   * the model obtained from the appropriate get() methods.
   *
   * @param data the data
   * @param width the width
   * @param height the height
   * @param samplePositions the sample positions
   */
  public void eraseSample(double[] data, int width, int height, int[] samplePositions) {
    if (samplePositions == null) {
      return;
    }

    // Remove from the input data
    for (final int i : samplePositions) {
      data[i] -= 1;
    }
  }

  /**
   * Remove the PSF from the data. Can be invoked using a saved copy of the PSF previously drawn by
   * the model obtained from the appropriate get() methods.
   *
   * @param data the data
   * @param width the width
   * @param height the height
   * @param samplePositions the sample positions
   */
  public void eraseSample(float[] data, int width, int height, int[] samplePositions) {
    if (samplePositions == null) {
      return;
    }

    // Remove from the input data
    for (final int i : samplePositions) {
      data[i] -= 1;
    }
  }

  /**
   * Set the random generator used for the random data generator to create data.
   *
   * @param randomGenerator the new random generator
   */
  public void setRandomGenerator(RandomGenerator randomGenerator) {
    if (randomGenerator == null) {
      throw new IllegalArgumentException("Random generator was null");
    }
    rand = new RandomDataGenerator(randomGenerator);
  }

  /**
   * Set the random data generator used to create data.
   *
   * @param randomDataGenerator the new random generator
   */
  public void setRandomGenerator(RandomDataGenerator randomDataGenerator) {
    if (randomDataGenerator == null) {
      throw new IllegalArgumentException("Random generator was null");
    }
    rand = randomDataGenerator;
  }

  /**
   * Get the value of the PSF function.
   *
   * @param width The data width
   * @param height The data height
   * @param x0 The centre in dimension 0
   * @param x1 The centre in dimension 1
   * @param x2 The centre in dimension 2
   * @param value the value
   * @return true, if successful; false, if no value was computed
   * @throws IllegalArgumentException If value or gradient are not length [width*height]
   */
  public boolean getValue(final int width, final int height, double x0, double x1, double x2,
      double[] value) {
    checkSize(width, height);
    final int size = width * height;
    if (value.length != size) {
      throw new IllegalArgumentException("Value is not the correct size");
    }
    Arrays.fill(value, 0);
    return computeValue(width, height, x0, x1, x2, value);
  }

  /**
   * Check size.
   *
   * @param width the width
   * @param height the height
   * @throws IllegalArgumentException If width * height is too large for an array
   */
  protected void checkSize(int width, int height) {
    if (width < 1) {
      throw new IllegalArgumentException("Width cannot be less than 1");
    }
    if (height < 1) {
      throw new IllegalArgumentException("Height cannot be less than 1");
    }
    if ((double) width * height > Integer.MAX_VALUE) {
      throw new IllegalArgumentException("width*height is too large");
    }
  }

  /**
   * Compute the value of the PSF function.
   *
   * <p>This can be over-ridden if the result is different from
   * {@link #create3D(double[], int, int, double, double, double, double, boolean)} using a sum of 1
   * and no Poisson noise.
   *
   * @param width The data width
   * @param height The data height
   * @param x0 The centre in dimension 0
   * @param x1 The centre in dimension 1
   * @param x2 The centre in dimension 2
   * @param value the value
   * @return true, if successful; false, if no value was computed
   */
  protected boolean computeValue(final int width, final int height, double x0, double x1, double x2,
      double[] value) {
    // Default implementation. Allow this to be overridden.
    return create3D(value, width, height, 1, x0, x1, x2, false) != 0;
  }

  /**
   * Get the value and partial gradient of the PSF function.
   *
   * <p>If gradient[i].length is not 3 then the correct array size will be created. Allows the
   * function to be called using <code>double[][] gradient = new double[width*height][]</code>.
   *
   * @param width The data width
   * @param height The data height
   * @param x0 The centre in dimension 0
   * @param x1 The centre in dimension 1
   * @param x2 The centre in dimension 2
   * @param value the value
   * @param gradient the partial gradient at each point for each dimension
   * @return the value and gradient true, if successful; false, if no value was computed
   * @throws IllegalArgumentException If value or gradient are not length [width*height]
   */
  public boolean getValueAndGradient(final int width, final int height, double x0, double x1,
      double x2, double[] value, double[][] gradient) {
    checkSize(width, height);
    final int size = width * height;
    if (value.length != size) {
      throw new IllegalArgumentException("Value is not the correct size");
    }
    if (gradient.length != size) {
      throw new IllegalArgumentException("Gradient is not the correct size");
    }
    for (int i = 0; i < gradient.length; i++) {
      if (gradient[i] == null || gradient[i].length != 3) {
        gradient[i] = new double[3];
      } else {
        gradient[i][0] = 0;
        gradient[i][1] = 0;
        gradient[i][2] = 0;
      }
    }
    Arrays.fill(value, 0);
    return computeValueAndGradient(width, height, x0, x1, x2, value, gradient);
  }

  /**
   * Compute the value and partial gradient of the PSF function.
   *
   * @param width The data width
   * @param height The data height
   * @param x0 The centre in dimension 0
   * @param x1 The centre in dimension 1
   * @param x2 The centre in dimension 2
   * @param value the value
   * @param jacobian the jacobian (partial gradient at each point for each dimension)
   * @return true, if successful; false, if no value was computed
   */
  protected abstract boolean computeValueAndGradient(final int width, final int height, double x0,
      double x1, double x2, double[] value, double[][] jacobian);

  /**
   * Compute the value and partial gradient of the PSF function using numerical gradients.
   *
   * <p>This is a helper function for sub-classes and arguments are unchecked.
   *
   * @param width The data width
   * @param height The data height
   * @param x0 The centre in dimension 0
   * @param x1 The centre in dimension 1
   * @param x2 The centre in dimension 2
   * @param value the value
   * @param jacobian the jacobian (partial gradient at each point for each dimension)
   * @param dx the delta for each dimension for numerical gradient
   * @return true, if successful; false, if no value was computed
   */
  protected boolean computeValueAndGradient(final int width, final int height, double x0, double x1,
      double x2, double[] value, double[][] jacobian, double[] dx) {
    final int size = width * height;
    final double[] v1 = new double[size];
    final double[] v2 = new double[size];
    // Compute the value
    Arrays.fill(value, 0);
    if (!computeValue(width, height, x0, x1, x2, value)) {
      return false;
    }
    final double[] x = {x0, x1, x2};
    for (int i = 0; i < 3; i++) {
      // Numerical gradient
      final double p = x[i];
      double delta = Precision.representableDelta(p, dx[i]);
      x[i] = p + delta;
      Arrays.fill(v1, 0);
      final boolean upper = computeValue(width, height, x[0], x[1], x[2], v1);
      x[i] = p - delta;
      Arrays.fill(v2, 0);
      final boolean lower = computeValue(width, height, x[0], x[1], x[2], v2);
      x[i] = p;
      final double[] u = (upper) ? v1 : value;
      final double[] l = (lower) ? v2 : value;
      if (u == l) {
        return false;
      }
      if (upper && lower) {
        delta *= 2;
      }
      for (int j = 0; j < size; j++) {
        jacobian[j][i] = (u[j] - l[j]) / delta;
      }
    }
    return true;
  }
}
