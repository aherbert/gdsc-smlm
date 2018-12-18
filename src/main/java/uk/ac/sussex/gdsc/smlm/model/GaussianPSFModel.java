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

package uk.ac.sussex.gdsc.smlm.model;

import uk.ac.sussex.gdsc.smlm.function.gaussian.AstigmatismZModel;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.NullAstigmatismZModel;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.SingleAstigmatismErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.utils.Pair;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.special.Erf;

/**
 * Contains methods for generating models of a Point Spread Function using a Gaussian approximation.
 */
public class GaussianPSFModel extends PSFModel {
  private double s0;
  private double s1;
  private AstigmatismZModel zModel;
  private double range = 5;

  /**
   * Instantiates a new gaussian PSF model.
   *
   * @param s0 The Gaussian standard deviation dimension 0
   * @param s1 The Gaussian standard deviation dimension 1
   */
  public GaussianPSFModel(double s0, double s1) {
    super();
    zModel = new NullAstigmatismZModel(s0, s1);
  }

  /**
   * Instantiates a new gaussian PSF model.
   *
   * @param zModel the z model
   */
  public GaussianPSFModel(AstigmatismZModel zModel) {
    super();
    if (zModel == null) {
      throw new IllegalArgumentException("Model must not be null");
    }
    this.zModel = zModel;
  }

  /**
   * Instantiates a new gaussian PSF model.
   *
   * @param randomGenerator the random generator
   * @param s0 The Gaussian standard deviation dimension 0
   * @param s1 The Gaussian standard deviation dimension 1
   */
  public GaussianPSFModel(RandomGenerator randomGenerator, double s0, double s1) {
    super(randomGenerator);
    zModel = new NullAstigmatismZModel(s0, s1);
  }

  /**
   * Instantiates a new gaussian PSF model.
   *
   * @param randomGenerator the random generator
   * @param zModel the z model
   */
  public GaussianPSFModel(RandomGenerator randomGenerator, AstigmatismZModel zModel) {
    super(randomGenerator);
    if (zModel == null) {
      throw new IllegalArgumentException("Model must not be null");
    }
    this.zModel = zModel;
  }

  /**
   * Instantiates a new gaussian PSF model.
   *
   * @param randomDataGenerator the random data generator
   * @param s0 The Gaussian standard deviation dimension 0
   * @param s1 The Gaussian standard deviation dimension 1
   */
  public GaussianPSFModel(RandomDataGenerator randomDataGenerator, double s0, double s1) {
    super(randomDataGenerator);
    zModel = new NullAstigmatismZModel(s0, s1);
  }

  /**
   * Instantiates a new gaussian PSF model.
   *
   * @param randomDataGenerator the random data generator
   * @param zModel the z model
   */
  public GaussianPSFModel(RandomDataGenerator randomDataGenerator, AstigmatismZModel zModel) {
    super(randomDataGenerator);
    if (zModel == null) {
      throw new IllegalArgumentException("Model must not be null");
    }
    this.zModel = zModel;
  }

  /**
   * Private constructor used in the {@link #copy()} method.
   */
  private GaussianPSFModel() {
    super();
  }

  /** {@inheritDoc} */
  @Override
  public double create3D(float[] data, final int width, final int height, final double sum,
      double x0, double x1, double x2, boolean poissonNoise) {
    if (sum == 0) {
      return 0;
    }
    try {
      final double d =
          gaussian2D(data, width, height, sum, x0, x1, getS0(x2), getS1(x2), poissonNoise);
      // if (d == 0)
      // {
      // System.out.printf("No data inserted: %f @ %f %f %f (%f x %f)\n", sum, x0, x1, x2, scale *
      // zeroS0,
      // scale * zeroS1);
      // }
      return d;
    } catch (final IllegalArgumentException ex) {
      // System.out.println(ex.getMessage());
      return 0;
    }
  }

  /** {@inheritDoc} */
  @Override
  public double create3D(double[] data, final int width, final int height, final double sum,
      double x0, double x1, double x2, boolean poissonNoise) {
    if (sum == 0) {
      return 0;
    }
    try {
      return gaussian2D(data, width, height, sum, x0, x1, getS0(x2), getS1(x2), poissonNoise);
    } catch (final IllegalArgumentException ex) {
      // System.out.println(ex.getMessage());
      return 0;
    }
  }

  /**
   * Gets the width in dimension 0 for the given z-depth.
   *
   * @param z the z
   * @return the s0
   */
  public double getS0(double z) {
    return zModel.getSx(z);
  }

  /**
   * Gets the width in dimension 1 for the given z-depth.
   *
   * @param z the z
   * @return the s1
   */
  public double getS1(double z) {
    return zModel.getSy(z);
  }

  /**
   * Construct a Gaussian 2D function on the provided data. Only evaluates the function within +/- 5
   * standard deviations in each direction from the centre (allows populating large images).
   *
   * <p>Builds the pixel approximation using the Gaussian error function as described in Smith et
   * al, (2010). Fast, single-molecule localisation that achieves theoretically minimum uncertainty.
   * Nature Methods 7, 373-375 (supplementary note).
   *
   * @param data The data (can be null)
   * @param width The data width
   * @param height The data height
   * @param sum The Gaussian integral
   * @param x0 The Gaussian centre in dimension 0
   * @param x1 The Gaussian centre in dimension 1
   * @param s0 The Gaussian standard deviation dimension 0
   * @param s1 The Gaussian standard deviation dimension 1
   * @param poissonNoise Add Poisson noise
   * @return The total sum added to the image (useful when poissonNoise is added)
   */
  public double gaussian2D(float[] data, final int width, final int height, final double sum,
      double x0, double x1, double s0, double s1, boolean poissonNoise) {
    if (sum == 0) {
      return 0;
    }
    // Parameter check
    checkSize(width, height);
    if (data == null) {
      data = new float[width * height];
    } else if (data.length < width * height) {
      throw new IllegalArgumentException("Data length cannot be smaller than width * height");
    }

    s0 = Math.abs(s0);
    s1 = Math.abs(s1);

    // Evaluate the Gaussian error function on a pixel grid covering +/- 5 SD
    final int x0min = clip((int) (x0 - range * s0), width);
    final int x1min = clip((int) (x1 - range * s1), height);
    final int x0max = clip((int) Math.ceil(x0 + range * s0), width);
    final int x1max = clip((int) Math.ceil(x1 + range * s1), height);

    final int x0range = x0max - x0min;
    final int x1range = x1max - x1min;

    // min should always be less than max
    if (x0range < 1) {
      throw new IllegalArgumentException("Gaussian dimension 0 range not within data bounds");
    }
    if (x1range < 1) {
      throw new IllegalArgumentException("Gaussian dimension 1 range not within data bounds");
    }

    // Shift centre to origin and compute gaussian
    final double[] gauss = gaussian2D(x0range, x1range, sum, x0 - x0min, x1 - x1min, s0, s1);

    return insert(data, x0min, x1min, x0max, x1max, width, gauss, poissonNoise);
  }

  /**
   * Construct a Gaussian 2D function on the provided data. Only evaluates the function within +/- 5
   * standard deviations in each direction from the centre (allows populating large images).
   *
   * <p>Builds the pixel approximation using the Gaussian error function as described in Smith et
   * al, (2010). Fast, single-molecule localisation that achieves theoretically minimum uncertainty.
   * Nature Methods 7, 373-375 (supplementary note).
   *
   * @param data The data (can be null)
   * @param width The data width
   * @param height The data height
   * @param sum The Gaussian integral
   * @param x0 The Gaussian centre in dimension 0
   * @param x1 The Gaussian centre in dimension 1
   * @param s0 The Gaussian standard deviation dimension 0
   * @param s1 The Gaussian standard deviation dimension 1
   * @param poissonNoise Add Poisson noise
   * @return The total sum added to the image (useful when poissonNoise is added)
   */
  public double gaussian2D(double[] data, final int width, final int height, final double sum,
      double x0, double x1, double s0, double s1, boolean poissonNoise) {
    if (sum == 0) {
      return 0;
    }
    // Parameter check
    checkSize(width, height);
    if (data == null) {
      data = new double[width * height];
    } else if (data.length < width * height) {
      throw new IllegalArgumentException("Data length cannot be smaller than width * height");
    }

    s0 = Math.abs(s0);
    s1 = Math.abs(s1);

    // Evaluate the Gaussian error function on a pixel grid covering +/- 5 SD
    final int x0min = clip((int) (x0 - range * s0), width);
    final int x1min = clip((int) (x1 - range * s1), height);
    final int x0max = clip((int) Math.ceil(x0 + range * s0), width);
    final int x1max = clip((int) Math.ceil(x1 + range * s1), height);

    final int x0range = x0max - x0min;
    final int x1range = x1max - x1min;

    // min should always be less than max
    if (x0range < 1) {
      throw new IllegalArgumentException("Gaussian dimension 0 range not within data bounds");
    }
    if (x1range < 1) {
      throw new IllegalArgumentException("Gaussian dimension 1 range not within data bounds");
    }

    // Shift centre to origin and compute gaussian
    final double[] gauss = gaussian2D(x0range, x1range, sum, x0 - x0min, x1 - x1min, s0, s1);

    return insert(data, x0min, x1min, x0max, x1max, width, gauss, poissonNoise);
  }

  private static final double ONE_OVER_ROOT2 = 1.0 / Math.sqrt(2);

  /**
   * Construct a Gaussian 2D function based at the origin using the specified range in each
   * dimension.
   *
   * <p>Builds the pixel approximation using the Gaussian error function as described in Smith et
   * al, (2010). Fast, single-molecule localisation that achieves theoretically minimum uncertainty.
   * Nature Methods 7, 373-375 (supplementary note).
   *
   * @param x0range The maximum range in dimension 0 (width)
   * @param x1range The maximum range in dimension 1 (height)
   * @param sum The Gaussian integral
   * @param x0 The Gaussian centre in dimension 0
   * @param x1 The Gaussian centre in dimension 1
   * @param s0 The Gaussian standard deviation dimension 0
   * @param s1 The Gaussian standard deviation dimension 1
   * @return The data (packed in yx order, length = x0range * x1range)
   */
  public double[] gaussian2D(int x0range, int x1range, double sum, double x0, double x1, double s0,
      double s1) {
    s0 = Math.abs(s0);
    s1 = Math.abs(s1);

    this.s0 = s0;
    this.s1 = s1;

    // Compute Gaussian error function grid up to and including the final grid position
    final double[] erf0 = new double[x0range + 1];
    final double[] erf1 = new double[x1range + 1];

    final double denom0 = ONE_OVER_ROOT2 / s0;
    final double denom1 = ONE_OVER_ROOT2 / s1;

    // Note: The 0.5 factors are moved to reduce computations
    for (int x = 0; x <= x0range; x++) {
      // erf0[x] = 0.5 * Erf.erf((x - x0) * denom0);
      erf0[x] = Erf.erf((x - x0) * denom0);
    }
    for (int y = 0; y <= x1range; y++) {
      // erf1[y] = 0.5 * Erf.erf((y - x1) * denom1);
      erf1[y] = Erf.erf((y - x1) * denom1);
    }
    sum *= 0.5; // Incorporate the 0.5 factor for Y into the sum

    // Pre-compute deltaE0
    final double[] deltaE0 = new double[x0range];
    for (int x = 0; x < x0range; x++) {
      // Incorporate the 0.5 factor for X into the delta
      deltaE0[x] = 0.5 * (erf0[x + 1] - erf0[x]);
    }

    // Compute Gaussian using the difference of the Gaussian error function
    final double[] data = new double[x0range * x1range];
    for (int y = 0, i = 0; y < x1range; y++) {
      // Include the sum into the first deltaE to get the Gaussian integral
      final double deltaE1 = sum * (erf1[y + 1] - erf1[y]);

      for (int x = 0; x < x0range; x++, i++) {
        data[i] = deltaE0[x] * deltaE1;
      }
    }

    return data;
  }

  private static int clip(int x, int max) {
    if (x < 0) {
      x = 0;
    }
    if (x > max) {
      x = max;
    }
    return x;
  }

  /**
   * Gets the standard deviation dimension 0 for the last drawn Gaussian.
   *
   * @return The standard deviation dimension 0 for the last drawn Gaussian.
   */
  public double getS0() {
    return s0;
  }

  /**
   * Gets the standard deviation dimension 1 for the last drawn Gaussian.
   *
   * @return The standard deviation dimension 1 for the last drawn Gaussian.
   */
  public double getS1() {
    return s1;
  }

  /** {@inheritDoc} */
  @Override
  public GaussianPSFModel copy() {
    final GaussianPSFModel model = new GaussianPSFModel();
    model.zModel = zModel;
    return model;
  }

  @Override
  public int sample3D(float[] data, int width, int height, int n, double x0, double x1, double x2) {
    if (n <= 0) {
      return insertSample(data, width, height, null, null);
    }
    final double[][] sample = sample(n, x0, x1, getS0(x2), getS1(x2));
    return insertSample(data, width, height, sample[0], sample[1]);
  }

  @Override
  public int sample3D(double[] data, int width, int height, int n, double x0, double x1,
      double x2) {
    if (n <= 0) {
      return insertSample(data, width, height, null, null);
    }
    final double[][] sample = sample(n, x0, x1, getS0(x2), getS1(x2));
    return insertSample(data, width, height, sample[0], sample[1]);
  }

  /**
   * Sample from a Gaussian distribution.
   *
   * @param n The number of samples
   * @param x0 The Gaussian centre in dimension 0
   * @param x1 The Gaussian centre in dimension 1
   * @param s0 The Gaussian standard deviation dimension 0
   * @param s1 The Gaussian standard deviation dimension 1
   * @return The sample x and y values
   */
  public double[][] sample(final int n, final double x0, final double x1, final double s0,
      final double s1) {
    this.s0 = s0;
    this.s1 = s1;
    final double[] x = sample(n, x0, s0);
    final double[] y = sample(n, x1, s1);
    return new double[][] {x, y};
  }

  private double[] sample(final int n, final double mu, final double sigma) {
    final double[] x = new double[n];
    final RandomGenerator random = rand.getRandomGenerator();
    for (int i = 0; i < n; i++) {
      x[i] = sigma * random.nextGaussian() + mu;
    }
    return x;
  }

  @Override
  protected boolean computeValue(int width, int height, double x0, double x1, double x2,
      double[] value) {
    final double s0 = getS0(x2);
    final double s1 = getS1(x2);

    // Evaluate the Gaussian error function on a pixel grid covering +/- 5 SD
    final int x0min = clip((int) (x0 - range * s0), width);
    final int x1min = clip((int) (x1 - range * s1), height);
    final int x0max = clip((int) Math.ceil(x0 + range * s0), width);
    final int x1max = clip((int) Math.ceil(x1 + range * s1), height);

    final int x0range = x0max - x0min;
    final int x1range = x1max - x1min;

    // min should always be less than max
    if (x0range < 1) {
      throw new IllegalArgumentException("Gaussian dimension 0 range not within data bounds");
    }
    if (x1range < 1) {
      throw new IllegalArgumentException("Gaussian dimension 1 range not within data bounds");
    }

    // Evaluate using a function.
    // This allows testing the function verses the default gaussian2D() method
    // as they should be nearly identical (the Erf function may be different).
    // Thus the gradient can be evaluated using the same function.
    final ErfGaussian2DFunction f = createGaussianFunction(x0range, x1range);
    final double[] p = new double[Gaussian2DFunction.PARAMETERS_PER_PEAK + 1];
    p[Gaussian2DFunction.SIGNAL] = 1;
    // The function computes the centre of the pixel as 0,0
    p[Gaussian2DFunction.X_POSITION] = x0 - x0min - 0.5;
    p[Gaussian2DFunction.Y_POSITION] = x1 - x1min - 0.5;
    p[Gaussian2DFunction.Z_POSITION] = x2;
    final double[] v = f.computeValues(p);
    for (int y = 0; y < x1range; y++) {
      // Locate the insert location
      int indexTo = (y + x1min) * width + x0min;
      int indexFrom = y * x0range;
      for (int x = 0; x < x0range; x++) {
        value[indexTo++] = v[indexFrom++];
      }
    }
    return true;
  }

  private ErfGaussian2DFunction createGaussianFunction(final int x0range, final int x1range) {
    final ErfGaussian2DFunction f =
        new SingleAstigmatismErfGaussian2DFunction(x0range, x1range, zModel);
    f.setErfFunction(ErfGaussian2DFunction.ErfFunction.COMMONS_MATH);
    return f;
  }

  @Override
  protected boolean computeValueAndGradient(int width, int height, double x0, double x1, double x2,
      double[] value, double[][] jacobian) {
    final double s0 = getS0(x2);
    final double s1 = getS1(x2);

    // Evaluate the Gaussian error function on a pixel grid covering +/- 5 SD
    final int x0min = clip((int) (x0 - range * s0), width);
    final int x1min = clip((int) (x1 - range * s1), height);
    final int x0max = clip((int) Math.ceil(x0 + range * s0), width);
    final int x1max = clip((int) Math.ceil(x1 + range * s1), height);

    final int x0range = x0max - x0min;
    final int x1range = x1max - x1min;

    // min should always be less than max
    if (x0range < 1) {
      throw new IllegalArgumentException("Gaussian dimension 0 range not within data bounds");
    }
    if (x1range < 1) {
      throw new IllegalArgumentException("Gaussian dimension 1 range not within data bounds");
    }

    // Evaluate using a function.
    // This allows testing the function verses the default gaussian2D() method
    // as they should be nearly identical (the Erf function may be different).
    // Thus the gradient can be evaluated using the same function.
    final ErfGaussian2DFunction f = createGaussianFunction(x0range, x1range);
    final double[] p = new double[Gaussian2DFunction.PARAMETERS_PER_PEAK + 1];
    p[Gaussian2DFunction.SIGNAL] = 1;
    // The function computes the centre of the pixel as 0,0.
    // The PSF sets the centre as 0.5,0.5.
    p[Gaussian2DFunction.X_POSITION] = x0 - x0min - 0.5;
    p[Gaussian2DFunction.Y_POSITION] = x1 - x1min - 0.5;
    p[Gaussian2DFunction.Z_POSITION] = x2;
    final int i0 = f.findGradientIndex(Gaussian2DFunction.X_POSITION);
    final int i1 = f.findGradientIndex(Gaussian2DFunction.Y_POSITION);
    final int i2 = f.findGradientIndex(Gaussian2DFunction.Z_POSITION);
    final Pair<double[], double[][]> pair = f.computeValuesAndJacobian(p);
    final double[] v = pair.a;
    final double[][] j = pair.b;
    for (int y = 0; y < x1range; y++) {
      // Locate the insert location
      int indexTo = (y + x1min) * width + x0min;
      int indexFrom = y * x0range;
      for (int x = 0; x < x0range; x++) {
        value[indexTo] = v[indexFrom];
        jacobian[indexTo][0] = j[indexFrom][i0];
        jacobian[indexTo][1] = j[indexFrom][i1];
        jacobian[indexTo][2] = j[indexFrom][i2];
        indexTo++;
        indexFrom++;
      }
    }
    return true;
  }

  /**
   * Gets the range over which the Gaussian is evaluated. This is in units of standard deviation.
   *
   * @return the range
   */
  public double getRange() {
    return range;
  }

  /**
   * Sets the range over which the Gaussian is evaluated. This is in units of standard deviation.
   *
   * @param range the new range
   * @throws IllegalArgumentException If the range is not strictly positive
   */
  public void setRange(double range) {
    if (!(range > 0)) {
      throw new IllegalArgumentException("Range must be strictly positive");
    }
    this.range = range;
  }
}
