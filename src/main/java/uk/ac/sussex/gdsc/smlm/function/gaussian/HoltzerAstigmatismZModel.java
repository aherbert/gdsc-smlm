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

package uk.ac.sussex.gdsc.smlm.function.gaussian;

/**
 * Implements a astigmatism model of a 2D Gaussian function, where z-depth determines the x and y
 * width.
 *
 * <p>Ref: Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically
 * minimum uncertainty. Nature Methods 7, 373-375 (supplementary note).
 *
 * <p>Ref: Holtzer, L., Meckel, T. &amp; Schmidt, T. Nanometric three-dimensional tracking of
 * individual quantum dots in cells. Applied Physics Letters 90, 1–3 (2007).
 */
public class HoltzerAstigmatismZModel implements AstigmatismZModel {
  /** The width in the x focal plane. */
  public final double s0x;
  /** The width in the y focal plane. */
  public final double s0y;
  /** the gamma parameter (half the distance between the focal planes). */
  public final double gamma;
  /** one over the depth of focus squared (1./d^2) */
  public final double oneOverD2;
  /** Empirical constant A for the x-astigmatism of the PSF. */
  public final double ax;
  /** Empirical constant B for the x-astigmatism of the PSF. */
  public final double bx;
  /** Empirical constant A for the y-astigmatism of the PSF. */
  public final double ay;
  /** Empirical constant B for the y-astigmatism of the PSF. */
  public final double by;

  /**
   * Static constructor.
   *
   * <p>Note that a positive gamma puts the focal plane for the X-dimension above the z-centre
   * (positive Z) and the focal plane for the Y-dimension below the z-centre (negative Z). If gamma
   * is negative then the orientation of the focal planes of X and Y are reversed.
   *
   * @param s0x The width in the x focal plane
   * @param s0y The width in the y focal plane
   * @param gamma the gamma parameter (half the distance between the focal planes)
   * @param d the depth of focus
   * @param ax Empirical constant A for the x-astigmatism of the PSF
   * @param bx Empirical constant B for the x-astigmatism of the PSF
   * @param ay Empirical constant A for the y-astigmatism of the PSF
   * @param by Empirical constant B for the y-astigmatism of the PSF
   * @return the holtzer astimatism Z model
   */
  public static HoltzerAstigmatismZModel create(double s0x, double s0y, double gamma, double d,
      double ax, double bx, double ay, double by) {
    final double d2 = d * d;
    return new HoltzerAstigmatismZModel(s0x, s0y, gamma, 1.0 / d2, ax, bx, ay, by);
  }

  /**
   * Constructor.
   *
   * <p>Note that a positive gamma puts the focal plane for the X-dimension above the z-centre
   * (positive Z) and the focal plane for the Y-dimension below the z-centre (negative Z). If gamma
   * is negative then the orientation of the focal planes of X and Y are reversed.
   *
   * @param s0x The width in the x focal plane
   * @param s0y The width in the y focal plane
   * @param gamma the gamma parameter (half the distance between the focal planes)
   * @param oneOverD2 one over the depth of focus squared (1/d^2)
   * @param ax Empirical constant A for the x-astigmatism of the PSF
   * @param bx Empirical constant B for the x-astigmatism of the PSF
   * @param ay Empirical constant A for the y-astigmatism of the PSF
   * @param by Empirical constant B for the y-astigmatism of the PSF
   */
  public HoltzerAstigmatismZModel(double s0x, double s0y, double gamma, double oneOverD2, double ax,
      double bx, double ay, double by) {
    this.s0x = s0x;
    this.s0y = s0y;
    this.gamma = gamma;
    this.oneOverD2 = oneOverD2;
    this.ax = ax;
    this.bx = bx;
    this.ay = ay;
    this.by = by;
  }

  /**
   * Gets the standard deviation, first and second derivatives for the z-depth.
   *
   * @param s0 the width in the focal plane
   * @param z the z
   * @param oneOverD2 one over the depth of focus squared (1/d^2)
   * @param a Empirical constant A for the astigmatism of the PSF
   * @param b Empirical constant B for the astigmatism of the PSF
   * @param dsdz the first and second derivative of s given z
   * @return the standard deviation
   */
  public static double getS2(double s0, double z, double oneOverD2, double a, double b,
      double[] dsdz) {
    final double z2 = z * z;
    final double z3 = z2 * z;
    final double z4 = z2 * z2;
    // Eq. 17a
    final double s = Math.sqrt(1 + oneOverD2 * (z2 + a * z3 + b * z4));
    // Eq. 19a
    dsdz[0] = s0 * (oneOverD2 * (2 * z + a * 3 * z2 + b * 4 * z3)) / (2 * s);
    // Eq. 19b
    dsdz[1] = s0 * ((oneOverD2 * (2 + a * 6 * z + b * 12 * z2)) / (2 * s)
        - pow2(oneOverD2 * (2 * z + a * 3 * z2 + b * 4 * z3)) / (4 * s * s * s));
    return s0 * s;
  }

  private static double pow2(final double d) {
    return d * d;
  }

  /**
   * Gets the standard deviation and first derivative for the z-depth.
   *
   * @param s0 the width in the focal plane
   * @param z the z
   * @param oneOverD2 one over the depth of focus squared (1/d^2)
   * @param a Empirical constant A for the astigmatism of the PSF
   * @param b Empirical constant B for the astigmatism of the PSF
   * @param dsdz the first derivative of s given z
   * @return the standard deviation
   */
  public static double getS1(double s0, double z, double oneOverD2, double a, double b,
      double[] dsdz) {
    final double z2 = z * z;
    final double z3 = z2 * z;
    final double z4 = z2 * z2;
    // Eq. 17a
    final double s = Math.sqrt(1 + oneOverD2 * (z2 + a * z3 + b * z4));
    // Eq. 19a
    dsdz[0] = s0 * (oneOverD2 * (2 * z + a * 3 * z2 + b * 4 * z3)) / (2 * s);
    return s0 * s;
  }

  /**
   * Gets the standard deviation for the z-depth.
   *
   * @param s0 the width in the focal plane
   * @param z the z
   * @param oneOverD2 one over the depth of focus squared (1/d^2)
   * @param a Empirical constant A for the astigmatism of the PSF
   * @param b Empirical constant B for the astigmatism of the PSF
   * @return the standard deviation
   */
  public static double getS(double s0, double z, double oneOverD2, double a, double b) {
    final double z2 = z * z;
    final double z3 = z2 * z;
    final double z4 = z2 * z2;
    // Eq. 17a
    return s0 * Math.sqrt(1 + oneOverD2 * (z2 + a * z3 + b * z4));
  }

  /** {@inheritDoc} */
  @Override
  public double getSx(double z) {
    return getS(s0x, z - gamma, oneOverD2, ax, bx);
  }

  /** {@inheritDoc} */
  @Override
  public double getSx(double z, double[] dsdz) {
    return getS1(s0x, z - gamma, oneOverD2, ax, bx, dsdz);
  }

  /** {@inheritDoc} */
  @Override
  public double getSx2(double z, double[] dsdz) {
    return getS2(s0x, z - gamma, oneOverD2, ax, bx, dsdz);
  }

  /** {@inheritDoc} */
  @Override
  public double getSy(double z) {
    return getS(s0y, z + gamma, oneOverD2, ay, by);
  }

  /** {@inheritDoc} */
  @Override
  public double getSy(double z, double[] dsdz) {
    return getS1(s0y, z + gamma, oneOverD2, ay, by, dsdz);
  }

  /** {@inheritDoc} */
  @Override
  public double getSy2(double z, double[] dsdz) {
    return getS2(s0y, z + gamma, oneOverD2, ay, by, dsdz);
  }

  /** {@inheritDoc} */
  @Override
  public String toString() {
    return String.format("s0x=%f s0y=%f gamma=%f 1/d^2=%f Ax=%f Bx=%f Ay=%f By=%f", s0x, s0y, gamma,
        oneOverD2, ax, bx, ay, by);
  }
}
