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
package uk.ac.sussex.gdsc.smlm.function.gaussian.erf;

// import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.util.FastMath;

import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.smlm.function.ExtendedGradient2Procedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Procedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient2Procedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;

/**
 * Evaluates a 2-dimensional Gaussian function for a single peak.
 */
public class SingleFreeCircularErfGaussian2DFunction extends SingleErfGaussian2DFunction {
  /** The gradient indices. */
  static final int[] gradientIndices;
  static {
    gradientIndices = createGradientIndices(1, new SingleFreeCircularErfGaussian2DFunction(1, 1));
  }

  /**
   * Constructor.
   *
   * @param maxx The maximum x value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   * @param maxy The maximum y value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   */
  public SingleFreeCircularErfGaussian2DFunction(int maxx, int maxy) {
    super(maxx, maxy);
  }

  @Override
  public ErfGaussian2DFunction copy() {
    return new SingleFreeCircularErfGaussian2DFunction(maxx, maxy);
  }

  /** {@inheritDoc} */
  @Override
  public void initialise0(double[] a) {
    tB = a[Gaussian2DFunction.BACKGROUND];
    tI = a[Gaussian2DFunction.SIGNAL];
    // Pre-compute the offset by 0.5
    final double tx = a[Gaussian2DFunction.X_POSITION] + 0.5;
    final double ty = a[Gaussian2DFunction.Y_POSITION] + 0.5;
    final double tsx = abs(a[Gaussian2DFunction.X_SD]);
    final double tsy = abs(a[Gaussian2DFunction.Y_SD]);

    createDeltaETable(ONE_OVER_ROOT2 / tsx, deltaEx, tx);
    createDeltaETable(ONE_OVER_ROOT2 / tsy, deltaEy, ty);
  }

  @Override
  public double integral(double[] a) {
    final double tB = a[Gaussian2DFunction.BACKGROUND];
    final double tI = a[Gaussian2DFunction.SIGNAL];
    // Pre-compute the offset by 0.5
    final double tx = a[Gaussian2DFunction.X_POSITION] + 0.5;
    final double ty = a[Gaussian2DFunction.Y_POSITION] + 0.5;
    final double tsx = abs(a[Gaussian2DFunction.X_SD]);
    final double tsy = abs(a[Gaussian2DFunction.Y_SD]);

    return tB * size() + tI * compute1DIntegral(ONE_OVER_ROOT2 / tsx, maxx, tx)
        * compute1DIntegral(ONE_OVER_ROOT2 / tsy, maxy, ty);
  }

  /** {@inheritDoc} */
  @Override
  public void initialise1(double[] a) {
    create1Arrays();
    tB = a[Gaussian2DFunction.BACKGROUND];
    tI = a[Gaussian2DFunction.SIGNAL];
    // Pre-compute the offset by 0.5
    final double tx = a[Gaussian2DFunction.X_POSITION] + 0.5;
    final double ty = a[Gaussian2DFunction.Y_POSITION] + 0.5;
    final double tsx = abs(a[Gaussian2DFunction.X_SD]);
    final double tsy = abs(a[Gaussian2DFunction.Y_SD]);

    // We can pre-compute part of the derivatives for position and sd in arrays
    // since the Gaussian is XY separable
    createFirstOrderTables(tI, deltaEx, du_dtx, du_dtsx, tx, tsx);
    createFirstOrderTables(tI, deltaEy, du_dty, du_dtsy, ty, tsy);
  }

  /** {@inheritDoc} */
  @Override
  public void initialise2(double[] a) {
    create2Arrays();
    tB = a[Gaussian2DFunction.BACKGROUND];
    tI = a[Gaussian2DFunction.SIGNAL];
    // Pre-compute the offset by 0.5
    final double tx = a[Gaussian2DFunction.X_POSITION] + 0.5;
    final double ty = a[Gaussian2DFunction.Y_POSITION] + 0.5;
    final double tsx = abs(a[Gaussian2DFunction.X_SD]);
    final double tsy = abs(a[Gaussian2DFunction.Y_SD]);

    // We can pre-compute part of the derivatives for position and sd in arrays
    // since the Gaussian is XY separable
    createSecondOrderTables(tI, deltaEx, du_dtx, du_dtsx, d2u_dtx2, d2u_dtsx2, tx, tsx);
    createSecondOrderTables(tI, deltaEy, du_dty, du_dtsy, d2u_dty2, d2u_dtsy2, ty, tsy);
  }

  /** {@inheritDoc} */
  @Override
  public void initialiseExtended2(double[] a) {
    createEx2Arrays();
    tB = a[Gaussian2DFunction.BACKGROUND];
    tI = a[Gaussian2DFunction.SIGNAL];
    // Pre-compute the offset by 0.5
    final double tx = a[Gaussian2DFunction.X_POSITION] + 0.5;
    final double ty = a[Gaussian2DFunction.Y_POSITION] + 0.5;
    final double tsx = abs(a[Gaussian2DFunction.X_SD]);
    final double tsy = abs(a[Gaussian2DFunction.Y_SD]);

    // We can pre-compute part of the derivatives for position and sd in arrays
    // since the Gaussian is XY separable
    createExSecondOrderTables(tI, deltaEx, du_dtx, du_dtsx, d2u_dtx2, d2u_dtsx2, d2deltaEx_dtsxdx,
        tx, tsx);
    createExSecondOrderTables(tI, deltaEy, du_dty, du_dtsy, d2u_dty2, d2u_dtsy2, d2deltaEy_dtsydy,
        ty, tsy);
  }

  /**
   * Creates the delta E array. This is the sum of the Gaussian function using the error function
   * for each of the pixels from 0 to n.
   *
   * @param one_sSqrt2 one over (s times sqrt(2))
   * @param deltaE the delta E for dimension 0 (difference between the error function at the start
   *        and end of each pixel)
   * @param u the mean of the Gaussian for dimension 0
   */
  protected void createDeltaETable(double one_sSqrt2, double[] deltaE, double u) {
    // For documentation see SingleFreeCircularErfGaussian2DFunction.createSecondOrderTables(...)

    double x_u_p12 = -u;
    double erf_x_minus = 0.5 * erf(x_u_p12 * one_sSqrt2);
    for (int i = 0, n = deltaE.length; i < n; i++) {
      x_u_p12 += 1.0;
      final double erf_x_plus = 0.5 * erf(x_u_p12 * one_sSqrt2);
      deltaE[i] = erf_x_plus - erf_x_minus;
      erf_x_minus = erf_x_plus;
    }
  }

  /**
   * Creates the first order derivatives.
   *
   * @param tI the target intensity
   * @param deltaE the delta E for dimension 0 (difference between the error function at the start
   *        and end of each pixel)
   * @param du_dx the first order x derivative for dimension 0
   * @param du_ds the first order s derivative for dimension 0
   * @param u the mean of the Gaussian for dimension 0
   * @param s the standard deviation of the Gaussian for dimension 0
   */
  protected void createFirstOrderTables(double tI, double[] deltaE, double[] du_dx, double[] du_ds,
      double u, double s) {
    createFirstOrderTables(ONE_OVER_ROOT2 / s, 0.5 / (s * s), tI * ONE_OVER_ROOT2PI / s,
        tI * ONE_OVER_ROOT2PI / (s * s), deltaE, du_dx, du_ds, u);
  }

  /**
   * Creates the first order derivatives.
   *
   * @param one_sSqrt2 one over (s times sqrt(2))
   * @param one_2ss one over (2 * s^2)
   * @param I_sSqrt2pi the intensity over (s * sqrt(2*pi))
   * @param I_ssSqrt2pi the intensity over (s^2 * sqrt(2*pi))
   * @param deltaE the delta E for dimension 0 (difference between the error function at the start
   *        and end of each pixel)
   * @param du_dx the first order x derivative for dimension 0
   * @param du_ds the first order s derivative for dimension 0
   * @param u the mean of the Gaussian for dimension 0
   */
  protected void createFirstOrderTables(double one_sSqrt2, double one_2ss, double I_sSqrt2pi,
      double I_ssSqrt2pi, double[] deltaE, double[] du_dx, double[] du_ds, double u) {
    // For documentation see SingleFreeCircularErfGaussian2DFunction.createSecondOrderTables(...)

    double x_u_p12 = -u;
    double erf_x_minus = 0.5 * erf(x_u_p12 * one_sSqrt2);
    double exp_x_minus = FastMath.exp(-(x_u_p12 * x_u_p12 * one_2ss));
    for (int i = 0, n = deltaE.length; i < n; i++) {
      final double x_u_m12 = x_u_p12;
      x_u_p12 += 1.0;
      final double erf_x_plus = 0.5 * erf(x_u_p12 * one_sSqrt2);
      deltaE[i] = erf_x_plus - erf_x_minus;
      erf_x_minus = erf_x_plus;

      final double exp_x_plus = FastMath.exp(-(x_u_p12 * x_u_p12 * one_2ss));
      du_dx[i] = I_sSqrt2pi * (exp_x_minus - exp_x_plus);
      // Compute: I0 * G21(xk)
      du_ds[i] = I_ssSqrt2pi * (x_u_m12 * exp_x_minus - x_u_p12 * exp_x_plus);

      exp_x_minus = exp_x_plus;
    }
  }

  /**
   * Creates the first and second order derivatives.
   *
   * @param tI the target intensity
   * @param deltaE the delta E for dimension 0 (difference between the error function at the start
   *        and end of each pixel)
   * @param du_dx the first order x derivative for dimension 0
   * @param du_ds the first order s derivative for dimension 0
   * @param d2u_dx2 the second order x derivative for dimension 0
   * @param d2u_ds2 the second order s derivative for dimension 0
   * @param u the mean of the Gaussian for dimension 0
   * @param s the standard deviation of the Gaussian for dimension 0
   */
  protected void createSecondOrderTables(double tI, double[] deltaE, double[] du_dx, double[] du_ds,
      double[] d2u_dx2, double[] d2u_ds2, double u, double s) {
    final double ss = s * s;
    final double one_sSqrt2pi = ONE_OVER_ROOT2PI / s;
    final double one_ssSqrt2pi = ONE_OVER_ROOT2PI / ss;
    final double one_sssSqrt2pi = one_sSqrt2pi / ss;
    final double one_sssssSqrt2pi = one_sssSqrt2pi / ss;
    createSecondOrderTables(tI, ONE_OVER_ROOT2 / s, 0.5 / ss, tI * one_sSqrt2pi, tI * one_ssSqrt2pi,
        tI * one_sssSqrt2pi, ss, one_sssSqrt2pi, one_sssssSqrt2pi, deltaE, du_dx, du_ds, d2u_dx2,
        d2u_ds2, u);
  }

  /**
   * Creates the first and second order derivatives.
   *
   * @param tI the target intensity
   * @param one_sSqrt2 one over (s times sqrt(2))
   * @param one_2ss one over (2 * s^2)
   * @param I_sSqrt2pi the intensity over (s * sqrt(2*pi))
   * @param I_ssSqrt2pi the intensity over (s^2 * sqrt(2*pi))
   * @param I_sssSqrt2pi the intensity over (s^3 * sqrt(2*pi))
   * @param ss the standard deviation squared
   * @param one_sssSqrt2pi one over (s^3 * sqrt(2*pi))
   * @param one_sssssSqrt2pi one over (s^5 * sqrt(2*pi))
   * @param deltaE the delta E for dimension 0 (difference between the error function at the start
   *        and end of each pixel)
   * @param du_dx the first order x derivative for dimension 0
   * @param du_ds the first order s derivative for dimension 0
   * @param d2u_dx2 the second order x derivative for dimension 0
   * @param d2u_ds2 the second order s derivative for dimension 0
   * @param u the mean of the Gaussian for dimension 0
   */
  protected void createSecondOrderTables(double tI, double one_sSqrt2, double one_2ss,
      double I_sSqrt2pi, double I_ssSqrt2pi, double I_sssSqrt2pi, double ss, double one_sssSqrt2pi,
      double one_sssssSqrt2pi, double[] deltaE, double[] du_dx, double[] du_ds, double[] d2u_dx2,
      double[] d2u_ds2, double u) {
    // Note: The paper by Smith, et al computes the integral for the kth pixel centred at (x,y)
    // If x=u then the Erf will be evaluated at x-u+0.5 - x-u-0.5 => integral from -0.5 to 0.5.
    // This code sets the first pixel at (0,0).

    // All computations for pixel k (=(x,y)) that require the exponential can use x,y indices for
    // the
    // lower boundary value and x+1,y+1 indices for the upper value.

    // Working example of this in GraspJ source code:
    // https://github.com/isman7/graspj/blob/master/graspj/src/main/java/eu/brede/graspj/opencl/src/functions/
    // I have used the same notation for clarity

    // The first position:
    // Offset x by the position and get the pixel lower bound.
    // (x - u - 0.5) with x=0 and u offset by +0.5
    double x_u_p12 = -u;
    double erf_x_minus = 0.5 * erf(x_u_p12 * one_sSqrt2);
    double exp_x_minus = FastMath.exp(-(x_u_p12 * x_u_p12 * one_2ss));
    for (int i = 0, n = deltaE.length; i < n; i++) {
      final double x_u_m12 = x_u_p12;
      x_u_p12 += 1.0;
      final double erf_x_plus = 0.5 * erf(x_u_p12 * one_sSqrt2);
      deltaE[i] = erf_x_plus - erf_x_minus;
      erf_x_minus = erf_x_plus;

      final double exp_x_plus = FastMath.exp(-(x_u_p12 * x_u_p12 * one_2ss));
      du_dx[i] = I_sSqrt2pi * (exp_x_minus - exp_x_plus);
      // Compute: I0 * G21(xk)
      final double pre2 = (x_u_m12 * exp_x_minus - x_u_p12 * exp_x_plus);
      du_ds[i] = I_ssSqrt2pi * pre2;

      // Second derivatives
      d2u_dx2[i] = I_sssSqrt2pi * pre2;

      // Compute G31(xk)
      final double G31 = one_sssSqrt2pi * pre2;

      // Compute G53(xk)
      final double G53 = one_sssssSqrt2pi
          * (MathUtils.pow3(x_u_m12) * exp_x_minus - MathUtils.pow3(x_u_p12) * exp_x_plus);
      d2u_ds2[i] = tI * (G53 - 2 * G31);

      exp_x_minus = exp_x_plus;
    }
  }

  /**
   * Creates the first and extended second order derivatives.
   *
   * @param tI the target intensity
   * @param deltaE the delta E for dimension 0 (difference between the error function at the start
   *        and end of each pixel)
   * @param du_dx the first order x derivative for dimension 0
   * @param du_ds the first order s derivative for dimension 0
   * @param d2u_dx2 the second order x derivative for dimension 0
   * @param d2u_ds2 the second order s derivative for dimension 0
   * @param d2deltaE_dsdx the second order deltaE s,x derivative for dimension 0
   * @param u the mean of the Gaussian for dimension 0
   * @param s the standard deviation of the Gaussian for dimension 0
   */
  protected void createExSecondOrderTables(double tI, double[] deltaE, double[] du_dx,
      double[] du_ds, double[] d2u_dx2, double[] d2u_ds2, double[] d2deltaE_dsdx, double u,
      double s) {
    final double ss = s * s;
    final double one_sSqrt2pi = ONE_OVER_ROOT2PI / s;
    final double one_ssSqrt2pi = ONE_OVER_ROOT2PI / ss;
    final double one_sssSqrt2pi = one_sSqrt2pi / ss;
    final double one_sssssSqrt2pi = one_sssSqrt2pi / ss;
    createExSecondOrderTables(tI, ONE_OVER_ROOT2 / s, 0.5 / ss, tI * one_sSqrt2pi,
        tI * one_ssSqrt2pi, tI * one_sssSqrt2pi, ss, one_sssSqrt2pi, one_sssssSqrt2pi, deltaE,
        du_dx, du_ds, d2u_dx2, d2u_ds2, d2deltaE_dsdx, u);
  }

  /**
   * Creates the first and second order derivatives.
   *
   * @param tI the target intensity
   * @param one_sSqrt2 one over (s times sqrt(2))
   * @param one_2ss one over (2 * s^2)
   * @param I_sSqrt2pi the intensity over (s * sqrt(2*pi))
   * @param I_ssSqrt2pi the intensity over (s^2 * sqrt(2*pi))
   * @param I_sssSqrt2pi the intensity over (s^3 * sqrt(2*pi))
   * @param ss the standard deviation squared
   * @param one_sssSqrt2pi one over (s^3 * sqrt(2*pi))
   * @param one_sssssSqrt2pi one over (s^5 * sqrt(2*pi))
   * @param deltaE the delta E for dimension 0 (difference between the error function at the start
   *        and end of each pixel)
   * @param du_dx the first order x derivative for dimension 0
   * @param du_ds the first order s derivative for dimension 0
   * @param d2u_dx2 the second order x derivative for dimension 0
   * @param d2u_ds2 the second order s derivative for dimension 0
   * @param d2deltaE_dsdx the second order deltaE s,x derivative for dimension 0
   * @param u the mean of the Gaussian for dimension 0
   */
  protected void createExSecondOrderTables(double tI, double one_sSqrt2, double one_2ss,
      double I_sSqrt2pi, double I_ssSqrt2pi, double I_sssSqrt2pi, double ss, double one_sssSqrt2pi,
      double one_sssssSqrt2pi, double[] deltaE, double[] du_dx, double[] du_ds, double[] d2u_dx2,
      double[] d2u_ds2, double[] d2deltaE_dsdx, double u) {
    // Note: The paper by Smith, et al computes the integral for the kth pixel centred at (x,y)
    // If x=u then the Erf will be evaluated at x-u+0.5 - x-u-0.5 => integral from -0.5 to 0.5.
    // This code sets the first pixel at (0,0).

    // All computations for pixel k (=(x,y)) that require the exponential can use x,y indices for
    // the
    // lower boundary value and x+1,y+1 indices for the upper value.

    // Working example of this in GraspJ source code:
    // https://github.com/isman7/graspj/blob/master/graspj/src/main/java/eu/brede/graspj/opencl/src/functions/
    // I have used the same notation for clarity

    // The first position:
    // Offset x by the position and get the pixel lower bound.
    // (x - u - 0.5) with x=0 and u offset by +0.5
    double x_u_p12 = -u;
    double erf_x_minus = 0.5 * erf(x_u_p12 * one_sSqrt2);
    double exp_x_minus = FastMath.exp(-(x_u_p12 * x_u_p12 * one_2ss));
    for (int i = 0, n = deltaE.length; i < n; i++) {
      final double x_u_m12 = x_u_p12;
      x_u_p12 += 1.0;
      final double erf_x_plus = 0.5 * erf(x_u_p12 * one_sSqrt2);
      deltaE[i] = erf_x_plus - erf_x_minus;
      erf_x_minus = erf_x_plus;

      final double exp_x_plus = FastMath.exp(-(x_u_p12 * x_u_p12 * one_2ss));
      du_dx[i] = I_sSqrt2pi * (exp_x_minus - exp_x_plus);
      // Compute: I0 * G21(xk)
      final double pre2 = (x_u_m12 * exp_x_minus - x_u_p12 * exp_x_plus);
      du_ds[i] = I_ssSqrt2pi * pre2;

      // Second derivatives
      d2u_dx2[i] = I_sssSqrt2pi * pre2;

      // Compute G31(xk)
      final double G31 = one_sssSqrt2pi * pre2;

      d2deltaE_dsdx[i] = I_ssSqrt2pi * (x_u_m12 * x_u_m12 * exp_x_minus / ss - exp_x_minus
          + exp_x_plus - x_u_p12 * x_u_p12 * exp_x_plus / ss);

      // Compute G53(xk)
      final double G53 = one_sssssSqrt2pi
          * (MathUtils.pow3(x_u_m12) * exp_x_minus - MathUtils.pow3(x_u_p12) * exp_x_plus);
      d2u_ds2[i] = tI * (G53 - 2 * G31);

      exp_x_minus = exp_x_plus;
    }
  }

  /** {@inheritDoc} */
  @Override
  public double eval(final int i, final double[] duda) {
    // Unpack the predictor into the dimensions
    final int y = i / maxx;
    final int x = i % maxx;

    // Return in order of Gaussian2DFunction.createGradientIndices().
    // Use pre-computed gradients
    duda[0] = 1.0;
    duda[1] = deltaEx[x] * deltaEy[y];
    duda[2] = du_dtx[x] * deltaEy[y];
    duda[3] = du_dty[y] * deltaEx[x];
    duda[4] = du_dtsx[x] * deltaEy[y];
    duda[5] = du_dtsy[y] * deltaEx[x];

    return tB + tI * duda[1];
  }

  /** {@inheritDoc} */
  @Override
  public double eval(final int i, final double[] duda, final double[] d2uda2) {
    // Unpack the predictor into the dimensions
    final int y = i / maxx;
    final int x = i % maxx;

    // Return in order of Gaussian2DFunction.createGradientIndices().
    // Use pre-computed gradients
    duda[0] = 1.0;
    duda[1] = deltaEx[x] * deltaEy[y];
    duda[2] = du_dtx[x] * deltaEy[y];
    duda[3] = du_dty[y] * deltaEx[x];
    duda[4] = du_dtsx[x] * deltaEy[y];
    duda[5] = du_dtsy[y] * deltaEx[x];
    d2uda2[0] = 0;
    d2uda2[1] = 0;
    d2uda2[2] = d2u_dtx2[x] * deltaEy[y];
    d2uda2[3] = d2u_dty2[y] * deltaEx[x];
    d2uda2[4] = d2u_dtsx2[x] * deltaEy[y];
    d2uda2[5] = d2u_dtsy2[y] * deltaEx[x];

    return tB + tI * duda[1];
  }

  @Override
  public boolean evaluatesBackground() {
    return true;
  }

  @Override
  public boolean evaluatesSignal() {
    return true;
  }

  @Override
  public boolean evaluatesAngle() {
    return false;
  }

  @Override
  public boolean evaluatesPosition() {
    return true;
  }

  @Override
  public boolean evaluatesSD0() {
    return true;
  }

  @Override
  public boolean evaluatesSD1() {
    return true;
  }

  @Override
  public int getGradientParametersPerPeak() {
    return 5;
  }

  /** {@inheritDoc} */
  @Override
  public int[] gradientIndices() {
    return gradientIndices;
  }

  /** {@inheritDoc} */
  @Override
  public int getNumberOfGradients() {
    return 6;
  }

  /** {@inheritDoc} */
  @Override
  public void forEach(Gradient1Procedure procedure) {
    final double[] duda = new double[getNumberOfGradients()];
    duda[0] = 1.0;
    for (int y = 0; y < maxy; y++) {
      final double du_dty = this.du_dty[y];
      final double deltaEy = this.deltaEy[y];
      final double du_dtsy = this.du_dtsy[y];
      for (int x = 0; x < maxx; x++) {
        duda[1] = deltaEx[x] * deltaEy;
        duda[2] = du_dtx[x] * deltaEy;
        duda[3] = du_dty * deltaEx[x];
        duda[4] = du_dtsx[x] * deltaEy;
        duda[5] = du_dtsy * deltaEx[x];
        procedure.execute(tB + tI * duda[1], duda);
      }
    }
  }

  /** {@inheritDoc} */
  @Override
  public void forEach(Gradient2Procedure procedure) {
    final double[] duda = new double[getNumberOfGradients()];
    final double[] d2uda2 = new double[getNumberOfGradients()];
    duda[0] = 1.0;
    for (int y = 0; y < maxy; y++) {
      final double du_dty = this.du_dty[y];
      final double deltaEy = this.deltaEy[y];
      final double du_dtsy = this.du_dtsy[y];
      final double d2u_dty2 = this.d2u_dty2[y];
      final double d2u_dtsy2 = this.d2u_dtsy2[y];
      for (int x = 0; x < maxx; x++) {
        duda[1] = deltaEx[x] * deltaEy;
        duda[2] = du_dtx[x] * deltaEy;
        duda[3] = du_dty * deltaEx[x];
        duda[4] = du_dtsx[x] * deltaEy;
        duda[5] = du_dtsy * deltaEx[x];
        d2uda2[2] = d2u_dtx2[x] * deltaEy;
        d2uda2[3] = d2u_dty2 * deltaEx[x];
        d2uda2[4] = d2u_dtsx2[x] * deltaEy;
        d2uda2[5] = d2u_dtsy2 * deltaEx[x];
        procedure.execute(tB + tI * duda[1], duda, d2uda2);
      }
    }
  }

  /** {@inheritDoc} */
  @Override
  public void forEach(ExtendedGradient2Procedure procedure) {
    final int n = getNumberOfGradients();
    final double[] duda = new double[n];
    final double[] d2udadb = new double[n * n];
    duda[0] = 1.0;
    final double[] du_dtsx_tI = new double[maxx];
    for (int x = 0; x < maxx; x++) {
      du_dtsx_tI[x] = du_dtsx[x] / tI;
    }
    for (int y = 0; y < maxy; y++) {
      final double du_dty = this.du_dty[y];
      final double du_dty_tI = du_dty / tI;
      final double deltaEy = this.deltaEy[y];
      final double du_dtsy = this.du_dtsy[y];
      final double du_dtsy_tI = du_dtsy / tI;
      final double d2u_dty2 = this.d2u_dty2[y];
      final double d2u_dtsy2 = this.d2u_dtsy2[y];
      final double d2deltaEy_dtsydy = this.d2deltaEy_dtsydy[y];
      for (int x = 0; x < maxx; x++) {
        duda[1] = deltaEx[x] * deltaEy;
        duda[2] = du_dtx[x] * deltaEy;
        duda[3] = du_dty * deltaEx[x];
        duda[4] = du_dtsx[x] * deltaEy;
        duda[5] = du_dtsy * deltaEx[x];

        // Compute all the partial second order derivatives

        // Background are all 0

        // Signal,X
        d2udadb[8] = duda[2] / tI;
        // Signal,Y
        d2udadb[9] = duda[3] / tI;
        // Signal,X SD
        d2udadb[10] = duda[4] / tI;
        // Signal,Y SD
        d2udadb[11] = duda[5] / tI;

        // X,Signal
        d2udadb[13] = d2udadb[8];
        // X,X
        d2udadb[14] = d2u_dtx2[x] * deltaEy;
        // X,Y
        d2udadb[15] = du_dtx[x] * du_dty_tI;
        // X,X SD
        d2udadb[16] = deltaEy * d2deltaEx_dtsxdx[x];
        // X,Y SD
        d2udadb[17] = du_dtx[x] * du_dtsy_tI;

        // Y,Signal
        d2udadb[19] = d2udadb[9];
        // Y,X
        d2udadb[20] = d2udadb[15];
        // Y,Y
        d2udadb[21] = d2u_dty2 * deltaEx[x];
        // Y,X SD
        d2udadb[22] = du_dty * du_dtsx_tI[x];
        // Y,Y SD
        d2udadb[23] = deltaEx[x] * d2deltaEy_dtsydy;

        // X SD,Signal
        d2udadb[25] = d2udadb[10];
        // X SD,X
        d2udadb[26] = d2udadb[16];
        // X SD,Y
        d2udadb[27] = d2udadb[22];
        // X SD,X SD
        d2udadb[28] = d2u_dtsx2[x] * deltaEy;
        // X SD,Y SD
        d2udadb[29] = du_dtsy * du_dtsx_tI[x];

        // Y SD,Signal
        d2udadb[31] = d2udadb[11];
        // Y SD,X
        d2udadb[32] = d2udadb[17];
        // Y SD,Y
        d2udadb[33] = d2udadb[23];
        // Y SD,X SD
        d2udadb[34] = d2udadb[29];
        // Y SD,Y SD
        d2udadb[35] = d2u_dtsy2 * deltaEx[x];

        procedure.executeExtended(tB + tI * duda[1], duda, d2udadb);
      }
    }
  }
}
