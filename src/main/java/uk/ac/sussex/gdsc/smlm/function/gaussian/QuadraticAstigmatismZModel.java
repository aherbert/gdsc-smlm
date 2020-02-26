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

package uk.ac.sussex.gdsc.smlm.function.gaussian;

/**
 * Implements an astigmatism model of a 2D Gaussian function, where z-depth determines the x and y
 * width. This is a simple quadratic model where the z-depth for the width at 1.5 can be specified.
 */
public class QuadraticAstigmatismZModel implements AstigmatismZModel {
  /** The gamma parameter (half the distance between the focal planes). */
  public final double gamma;
  /** The z-depth where the width is 1.5. */
  public final double zDepth;
  /** Pre-computed factor for the second derivative of s given z. */
  private final double d2sOverDz2;

  /**
   * Constructor.
   *
   * @param gamma The gamma parameter (half the distance between the focal planes)
   * @param zDepth The z-depth where the width is 1.5
   */
  public QuadraticAstigmatismZModel(double gamma, double zDepth) {
    this.gamma = gamma;
    this.zDepth = zDepth;
    d2sOverDz2 = 1.0 / (zDepth * zDepth);
  }

  /**
   * Gets the standard deviation, first and second derivatives for the z-depth.
   *
   * @param z the z
   * @param zDepth The z-depth where the width is 1.5
   * @param dsdz the first and second derivative of s given z
   * @return the standard deviation
   */
  public static double getS2(double z, double zDepth, double[] dsdz) {
    z /= zDepth; // Scale so z=1 at the configured z-depth
    final double s = 1.0 + z * z * 0.5;
    dsdz[0] = z / zDepth;
    dsdz[1] = 1.0 / (zDepth * zDepth);
    return s;
  }

  /**
   * Gets the standard deviation, first and second derivatives for the z-depth.
   *
   * @param z the z
   * @param dsdz the first and second derivative of s given z
   * @return the standard deviation
   */
  private double getS2(double z, double[] dsdz) {
    final double s = getS1(z, zDepth, dsdz);
    dsdz[1] = d2sOverDz2; // Use cached value
    return s;
  }

  /**
   * Gets the standard deviation and first derivative for the z-depth.
   *
   * @param z the z
   * @param zDepth The z-depth where the width is 1.5
   * @param dsdz the first derivative of s given z
   * @return the standard deviation
   */
  public static double getS1(double z, double zDepth, double[] dsdz) {
    z /= zDepth; // Scale so z=1 at the configured z-depth
    final double s = 1.0 + z * z * 0.5;
    dsdz[0] = z / zDepth;
    return s;
  }

  /**
   * Gets the standard deviation for the z-depth.
   *
   * @param z the z
   * @param zDepth The z-depth where the width is 1.5
   * @return the standard deviation
   */
  public static double getS(double z, double zDepth) {
    z /= zDepth; // Scale so z=1 at the configured z-depth
    return 1.0 + z * z * 0.5;
  }

  @Override
  public double getSx(double z) {
    return getS(z - gamma, zDepth);
  }

  @Override
  public double getSx(double z, double[] dsdz) {
    return getS1(z - gamma, zDepth, dsdz);
  }

  @Override
  public double getSx2(double z, double[] dsdz) {
    return getS2(z - gamma, dsdz);
  }

  @Override
  public double getSy(double z) {
    return getS(z + gamma, zDepth);
  }

  @Override
  public double getSy(double z, double[] dsdz) {
    return getS1(z + gamma, zDepth, dsdz);
  }

  @Override
  public double getSy2(double z, double[] dsdz) {
    return getS2(z + gamma, dsdz);
  }
}
