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

package uk.ac.sussex.gdsc.smlm.function.gaussian.erf;

import org.apache.commons.math3.util.FastMath;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.smlm.function.ExtendedGradient2Procedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Procedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient2Procedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;

/**
 * Evaluates a 2-dimensional Gaussian function for a single peak.
 */
public class MultiFreeCircularErfGaussian2DFunction extends MultiErfGaussian2DFunction {
  // Allow underscores in the variables used during computation
  // CHECKSTYLE.OFF: ParameterName|LocalVariableName

  /**
   * Constructor.
   *
   * @param numberOfPeaks The number of peaks
   * @param maxx The maximum x value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   * @param maxy The maximum y value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   */
  public MultiFreeCircularErfGaussian2DFunction(int numberOfPeaks, int maxx, int maxy) {
    super(numberOfPeaks, maxx, maxy);
  }

  @Override
  protected int[] createGradientIndices() {
    return replicateGradientIndices(SingleFreeCircularErfGaussian2DFunction.gradientIndices);
  }

  @Override
  public ErfGaussian2DFunction copy() {
    return new MultiFreeCircularErfGaussian2DFunction(numberOfPeaks, maxx, maxy);
  }

  @Override
  public void initialise0(double[] a) {
    tb = a[Gaussian2DFunction.BACKGROUND];

    for (int n = 0, i = 0; n < numberOfPeaks; n++, i += PARAMETERS_PER_PEAK) {
      tI[n] = a[i + Gaussian2DFunction.SIGNAL];
      // Pre-compute the offset by 0.5
      final double tx = a[i + Gaussian2DFunction.X_POSITION] + 0.5;
      final double ty = a[i + Gaussian2DFunction.Y_POSITION] + 0.5;
      final double tsx = Math.abs(a[i + Gaussian2DFunction.X_SD]);
      final double tsy = Math.abs(a[i + Gaussian2DFunction.Y_SD]);

      createDeltaETable(n, maxx, ONE_OVER_ROOT2 / tsx, deltaEx, tx);
      createDeltaETable(n, maxy, ONE_OVER_ROOT2 / tsy, deltaEy, ty);
    }
  }

  @Override
  public double integral(double[] a) {
    double sum = a[Gaussian2DFunction.BACKGROUND] * size();
    for (int n = 0, i = 0; n < numberOfPeaks; n++, i += PARAMETERS_PER_PEAK) {
      final double tI = a[i + Gaussian2DFunction.SIGNAL];
      // Pre-compute the offset by 0.5
      final double tx = a[i + Gaussian2DFunction.X_POSITION] + 0.5;
      final double ty = a[i + Gaussian2DFunction.Y_POSITION] + 0.5;
      final double tsx = Math.abs(a[i + Gaussian2DFunction.X_SD]);
      final double tsy = Math.abs(a[i + Gaussian2DFunction.Y_SD]);
      sum += tI * compute1DIntegral(ONE_OVER_ROOT2 / tsx, maxx, tx)
          * compute1DIntegral(ONE_OVER_ROOT2 / tsy, maxy, ty);
    }
    return sum;
  }

  @Override
  public void initialise1(double[] a) {
    create1Arrays();
    tb = a[Gaussian2DFunction.BACKGROUND];
    for (int n = 0, i = 0; n < numberOfPeaks; n++, i += PARAMETERS_PER_PEAK) {
      tI[n] = a[i + Gaussian2DFunction.SIGNAL];
      // Pre-compute the offset by 0.5
      final double tx = a[i + Gaussian2DFunction.X_POSITION] + 0.5;
      final double ty = a[i + Gaussian2DFunction.Y_POSITION] + 0.5;
      final double tsx = Math.abs(a[i + Gaussian2DFunction.X_SD]);
      final double tsy = Math.abs(a[i + Gaussian2DFunction.Y_SD]);

      // We can pre-compute part of the derivatives for position and sd in arrays
      // since the Gaussian is XY separable
      createFirstOrderTables(n, maxx, tI[n], deltaEx, duDtx, duDtsx, tx, tsx);
      createFirstOrderTables(n, maxy, tI[n], deltaEy, duDty, duDtsy, ty, tsy);
    }
  }

  @Override
  public void initialise2(double[] a) {
    create2Arrays();
    tb = a[Gaussian2DFunction.BACKGROUND];
    for (int n = 0, i = 0; n < numberOfPeaks; n++, i += PARAMETERS_PER_PEAK) {
      tI[n] = a[i + Gaussian2DFunction.SIGNAL];
      // Pre-compute the offset by 0.5
      final double tx = a[i + Gaussian2DFunction.X_POSITION] + 0.5;
      final double ty = a[i + Gaussian2DFunction.Y_POSITION] + 0.5;
      final double tsx = Math.abs(a[i + Gaussian2DFunction.X_SD]);
      final double tsy = Math.abs(a[i + Gaussian2DFunction.Y_SD]);

      // We can pre-compute part of the derivatives for position and sd in arrays
      // since the Gaussian is XY separable
      createSecondOrderTables(n, maxx, tI[n], deltaEx, duDtx, duDtsx, d2uDtx2, d2uDtsx2, tx, tsx);
      createSecondOrderTables(n, maxy, tI[n], deltaEy, duDty, duDtsy, d2uDty2, d2uDtsy2, ty, tsy);
    }
  }

  @Override
  public void initialiseExtended2(double[] a) {
    createEx2Arrays();
    tb = a[Gaussian2DFunction.BACKGROUND];
    for (int n = 0, i = 0; n < numberOfPeaks; n++, i += PARAMETERS_PER_PEAK) {
      tI[n] = a[i + Gaussian2DFunction.SIGNAL];
      // Pre-compute the offset by 0.5
      final double tx = a[i + Gaussian2DFunction.X_POSITION] + 0.5;
      final double ty = a[i + Gaussian2DFunction.Y_POSITION] + 0.5;
      final double tsx = Math.abs(a[i + Gaussian2DFunction.X_SD]);
      final double tsy = Math.abs(a[i + Gaussian2DFunction.Y_SD]);

      // We can pre-compute part of the derivatives for position and sd in arrays
      // since the Gaussian is XY separable
      createExSecondOrderTables(n, maxx, tI[n], deltaEx, duDtx, duDtsx, d2uDtx2, d2uDtsx2,
          d2deltaExDtsxDx, tx, tsx);
      createExSecondOrderTables(n, maxy, tI[n], deltaEy, duDty, duDtsy, d2uDty2, d2uDtsy2,
          d2deltaEyDtsyDy, ty, tsy);
    }
  }

  /**
   * Creates the delta E array. This is the sum of the Gaussian function using the error function
   * for each of the pixels from 0 to n.
   *
   * @param n the peak number
   * @param max the maximum for the dimension
   * @param one_sSqrt2 one over (s times sqrt(2))
   * @param deltaE the delta E for dimension 0 (difference between the error function at the start
   *        and end of each pixel)
   * @param u the mean of the Gaussian for dimension 0
   */
  protected void createDeltaETable(int n, int max, double one_sSqrt2, double[] deltaE, double u) {
    // For documentation see SingleFreeCircularErfGaussian2DFunction.createSecondOrderTables(...)

    double x_u_p12 = -u;
    double erf_x_minus = 0.5 * erf(x_u_p12 * one_sSqrt2);
    for (int i = 0, j = n * max; i < max; i++, j++) {
      x_u_p12 += 1.0;
      final double erf_x_plus = 0.5 * erf(x_u_p12 * one_sSqrt2);
      deltaE[j] = erf_x_plus - erf_x_minus;
      erf_x_minus = erf_x_plus;
    }
  }

  /**
   * Creates the first order derivatives.
   *
   * @param n the peak number
   * @param max the maximum for the dimension
   * @param tI the target intensity
   * @param deltaE the delta E for dimension 0 (difference between the error function at the start
   *        and end of each pixel)
   * @param duDx the first order x derivative for dimension 0
   * @param duDs the first order s derivative for dimension 0
   * @param u the mean of the Gaussian for dimension 0
   * @param s the standard deviation of the Gaussian for dimension 0
   */
  protected void createFirstOrderTables(int n, int max, double tI, double[] deltaE, double[] duDx,
      double[] duDs, double u, double s) {
    createFirstOrderTables(n, max, ONE_OVER_ROOT2 / s, 0.5 / (s * s), tI * ONE_OVER_ROOT2PI / s,
        tI * ONE_OVER_ROOT2PI / (s * s), deltaE, duDx, duDs, u);
  }

  /**
   * Creates the first order derivatives.
   *
   * @param n the peak number
   * @param max the maximum for the dimension
   * @param one_sSqrt2 one over (s times sqrt(2))
   * @param one_2ss one over (2 * s^2)
   * @param I_sSqrt2pi the intensity over (s * sqrt(2*pi))
   * @param I_ssSqrt2pi the intensity over (s^2 * sqrt(2*pi))
   * @param deltaE the delta E for dimension 0 (difference between the error function at the start
   *        and end of each pixel)
   * @param duDx the first order x derivative for dimension 0
   * @param duDs the first order s derivative for dimension 0
   * @param u the mean of the Gaussian for dimension 0
   */
  protected void createFirstOrderTables(int n, int max, double one_sSqrt2, double one_2ss,
      double I_sSqrt2pi, double I_ssSqrt2pi, double[] deltaE, double[] duDx, double[] duDs,
      double u) {
    // For documentation see SingleFreeCircularErfGaussian2DFunction.createSecondOrderTables(...)

    double x_u_p12 = -u;
    double erf_x_minus = 0.5 * erf(x_u_p12 * one_sSqrt2);
    double exp_x_minus = FastMath.exp(-(x_u_p12 * x_u_p12 * one_2ss));
    for (int i = 0, j = n * max; i < max; i++, j++) {
      final double x_u_m12 = x_u_p12;
      x_u_p12 += 1.0;
      final double erf_x_plus = 0.5 * erf(x_u_p12 * one_sSqrt2);
      deltaE[j] = erf_x_plus - erf_x_minus;
      erf_x_minus = erf_x_plus;

      final double exp_x_plus = FastMath.exp(-(x_u_p12 * x_u_p12 * one_2ss));
      duDx[j] = I_sSqrt2pi * (exp_x_minus - exp_x_plus);
      // Compute: I0 * G21(xk)
      duDs[j] = I_ssSqrt2pi * (x_u_m12 * exp_x_minus - x_u_p12 * exp_x_plus);

      exp_x_minus = exp_x_plus;
    }
  }

  /**
   * Creates the first and second order derivatives.
   *
   * @param n the peak number
   * @param max the maximum for the dimension
   * @param tI the target intensity
   * @param deltaE the delta E for dimension 0 (difference between the error function at the start
   *        and end of each pixel)
   * @param duDx the first order x derivative for dimension 0
   * @param duDs the first order s derivative for dimension 0
   * @param d2uDx2 the second order x derivative for dimension 0
   * @param d2uDs2 the second order s derivative for dimension 0
   * @param u the mean of the Gaussian for dimension 0
   * @param s the standard deviation of the Gaussian for dimension 0
   */
  protected void createSecondOrderTables(int n, int max, double tI, double[] deltaE, double[] duDx,
      double[] duDs, double[] d2uDx2, double[] d2uDs2, double u, double s) {
    final double ss = s * s;
    final double one_sSqrt2pi = ONE_OVER_ROOT2PI / s;
    final double one_ssSqrt2pi = ONE_OVER_ROOT2PI / ss;
    final double one_sssSqrt2pi = one_sSqrt2pi / ss;
    final double one_sssssSqrt2pi = one_sssSqrt2pi / ss;
    createSecondOrderTables(n, max, tI, ONE_OVER_ROOT2 / s, 0.5 / ss, tI * one_sSqrt2pi,
        tI * one_ssSqrt2pi, tI * one_sssSqrt2pi, ss, one_sssSqrt2pi, one_sssssSqrt2pi, deltaE, duDx,
        duDs, d2uDx2, d2uDs2, u);
  }

  /**
   * Creates the first and second order derivatives.
   *
   * @param n the peak number
   * @param max the maximum for the dimension
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
   * @param duDx the first order x derivative for dimension 0
   * @param duDs the first order s derivative for dimension 0
   * @param d2uDx2 the second order x derivative for dimension 0
   * @param d2uDs2 the second order s derivative for dimension 0
   * @param u the mean of the Gaussian for dimension 0
   */
  protected void createSecondOrderTables(int n, int max, double tI, double one_sSqrt2,
      double one_2ss, double I_sSqrt2pi, double I_ssSqrt2pi, double I_sssSqrt2pi, double ss,
      double one_sssSqrt2pi, double one_sssssSqrt2pi, double[] deltaE, double[] duDx, double[] duDs,
      double[] d2uDx2, double[] d2uDs2, double u) {
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
    for (int i = 0, j = n * max; i < max; i++, j++) {
      final double x_u_m12 = x_u_p12;
      x_u_p12 += 1.0;
      final double erf_x_plus = 0.5 * erf(x_u_p12 * one_sSqrt2);
      deltaE[j] = erf_x_plus - erf_x_minus;
      erf_x_minus = erf_x_plus;

      final double exp_x_plus = FastMath.exp(-(x_u_p12 * x_u_p12 * one_2ss));
      duDx[j] = I_sSqrt2pi * (exp_x_minus - exp_x_plus);
      // Compute: I0 * G21(xk)
      final double pre2 = (x_u_m12 * exp_x_minus - x_u_p12 * exp_x_plus);
      duDs[j] = I_ssSqrt2pi * pre2;

      // Second derivatives
      d2uDx2[j] = I_sssSqrt2pi * pre2;

      // Compute G31(xk)
      final double G31 = one_sssSqrt2pi * pre2;

      // Compute G53(xk)
      final double G53 = one_sssssSqrt2pi
          * (MathUtils.pow3(x_u_m12) * exp_x_minus - MathUtils.pow3(x_u_p12) * exp_x_plus);
      d2uDs2[j] = tI * (G53 - 2 * G31);

      exp_x_minus = exp_x_plus;
    }
  }

  /**
   * Creates the first and second order derivatives.
   *
   * @param n the peak number
   * @param max the maximum for the dimension
   * @param tI the target intensity
   * @param deltaE the delta E for dimension 0 (difference between the error function at the start
   *        and end of each pixel)
   * @param duDx the first order x derivative for dimension 0
   * @param duDs the first order s derivative for dimension 0
   * @param d2uDx2 the second order x derivative for dimension 0
   * @param d2uDs2 the second order s derivative for dimension 0
   * @param d2deltaE_dsdx the second order deltaE s,x derivative for dimension 0
   * @param u the mean of the Gaussian for dimension 0
   * @param s the standard deviation of the Gaussian for dimension 0
   */
  protected void createExSecondOrderTables(int n, int max, double tI, double[] deltaE,
      double[] duDx, double[] duDs, double[] d2uDx2, double[] d2uDs2, double[] d2deltaE_dsdx,
      double u, double s) {
    final double ss = s * s;
    final double one_sSqrt2pi = ONE_OVER_ROOT2PI / s;
    final double one_ssSqrt2pi = ONE_OVER_ROOT2PI / ss;
    final double one_sssSqrt2pi = one_sSqrt2pi / ss;
    final double one_sssssSqrt2pi = one_sssSqrt2pi / ss;
    createExSecondOrderTables(n, max, tI, ONE_OVER_ROOT2 / s, 0.5 / ss, tI * one_sSqrt2pi,
        tI * one_ssSqrt2pi, tI * one_sssSqrt2pi, ss, one_sssSqrt2pi, one_sssssSqrt2pi, deltaE, duDx,
        duDs, d2uDx2, d2uDs2, d2deltaE_dsdx, u);
  }

  /**
   * Creates the first and second order derivatives.
   *
   * @param n the peak number
   * @param max the maximum for the dimension
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
   * @param duDx the first order x derivative for dimension 0
   * @param duDs the first order s derivative for dimension 0
   * @param d2uDx2 the second order x derivative for dimension 0
   * @param d2uDs2 the second order s derivative for dimension 0
   * @param d2deltaE_dsdx the second order deltaE s,x derivative for dimension 0
   * @param u the mean of the Gaussian for dimension 0
   */
  protected void createExSecondOrderTables(int n, int max, double tI, double one_sSqrt2,
      double one_2ss, double I_sSqrt2pi, double I_ssSqrt2pi, double I_sssSqrt2pi, double ss,
      double one_sssSqrt2pi, double one_sssssSqrt2pi, double[] deltaE, double[] duDx, double[] duDs,
      double[] d2uDx2, double[] d2uDs2, double[] d2deltaE_dsdx, double u) {
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
    for (int i = 0, j = n * max; i < max; i++, j++) {
      final double x_u_m12 = x_u_p12;
      x_u_p12 += 1.0;
      final double erf_x_plus = 0.5 * erf(x_u_p12 * one_sSqrt2);
      deltaE[j] = erf_x_plus - erf_x_minus;
      erf_x_minus = erf_x_plus;

      final double exp_x_plus = FastMath.exp(-(x_u_p12 * x_u_p12 * one_2ss));
      duDx[j] = I_sSqrt2pi * (exp_x_minus - exp_x_plus);
      // Compute: I0 * G21(xk)
      final double pre2 = (x_u_m12 * exp_x_minus - x_u_p12 * exp_x_plus);
      duDs[j] = I_ssSqrt2pi * pre2;

      // Second derivatives
      d2uDx2[j] = I_sssSqrt2pi * pre2;

      // Compute G31(xk)
      final double G31 = one_sssSqrt2pi * pre2;

      d2deltaE_dsdx[j] = I_ssSqrt2pi * (x_u_m12 * x_u_m12 * exp_x_minus / ss - exp_x_minus
          + exp_x_plus - x_u_p12 * x_u_p12 * exp_x_plus / ss);

      // Compute G53(xk)
      final double G53 = one_sssssSqrt2pi
          * (MathUtils.pow3(x_u_m12) * exp_x_minus - MathUtils.pow3(x_u_p12) * exp_x_plus);
      d2uDs2[j] = tI * (G53 - 2 * G31);

      exp_x_minus = exp_x_plus;
    }
  }

  @Override
  public double eval(final int i, final double[] duda) {
    // Unpack the predictor into the dimensions
    int yy = i / maxx;
    int xx = i % maxx;

    // Return in order of Gaussian2DFunction.createGradientIndices().
    // Use pre-computed gradients
    duda[0] = 1.0;
    double I = tb;
    for (int n = 0, a = 1; n < numberOfPeaks; n++, xx += maxx, yy += maxy) {
      duda[a] = deltaEx[xx] * deltaEy[yy];
      I += tI[n] * duda[a++];
      duda[a++] = duDtx[xx] * deltaEy[yy];
      duda[a++] = duDty[yy] * deltaEx[xx];
      duda[a++] = duDtsx[xx] * deltaEy[yy];
      duda[a++] = duDtsy[yy] * deltaEx[xx];
    }
    return I;
  }

  @Override
  public double eval2(final int i, final double[] duda, final double[] d2uda2) {
    // Unpack the predictor into the dimensions
    int yy = i / maxx;
    int xx = i % maxx;

    // Return in order of Gaussian2DFunction.createGradientIndices().
    // Use pre-computed gradients
    duda[0] = 1.0;
    d2uda2[0] = 0;
    double I = tb;
    for (int n = 0, a = 1; n < numberOfPeaks; n++, xx += maxx, yy += maxy) {
      duda[a] = deltaEx[xx] * deltaEy[yy];
      I += tI[n] * duda[a];
      d2uda2[a++] = 0;
      duda[a] = duDtx[xx] * deltaEy[yy];
      d2uda2[a++] = d2uDtx2[xx] * deltaEy[yy];
      duda[a] = duDty[yy] * deltaEx[xx];
      d2uda2[a++] = d2uDty2[yy] * deltaEx[xx];
      duda[a] = duDtsx[xx] * deltaEy[yy];
      d2uda2[a++] = d2uDtsx2[xx] * deltaEy[yy];
      duda[a] = duDtsy[yy] * deltaEx[xx];
      d2uda2[a++] = d2uDtsy2[yy] * deltaEx[xx];
    }
    return I;
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

  @Override
  public void forEach(Gradient1Procedure procedure) {
    final double[] duda = new double[getNumberOfGradients()];
    duda[0] = 1.0;
    // Note: This unrolling does not perform better in JUnit speed test
    // // Unroll for the number of peaks
    // if (numberOfPeaks == 2)
    // {
    // for (int y = 0; y < maxy; y++)
    // {
    // for (int x = 0, xx = maxx, yy = maxy; x < maxx; x++, xx++, yy++)
    // {
    // duda[1] = deltaEx[x] * deltaEy[y];
    // duda[2] = du_dtx[x] * deltaEy[y];
    // duda[3] = du_dty[y] * deltaEx[x];
    // duda[4] = du_dtsx[x] * deltaEy[y];
    // duda[5] = du_dtsy[y] * deltaEx[x];
    // duda[6] = deltaEx[xx] * deltaEy[yy];
    // duda[7] = du_dtx[xx] * deltaEy[yy];
    // duda[8] = du_dty[yy] * deltaEx[xx];
    // duda[9] = du_dtsx[xx] * deltaEy[yy];
    // duda[10] = du_dtsy[yy] * deltaEx[xx];
    // procedure.execute(tb + tI[0] * duda[1] + tI[1] * duda[6], duda);
    // }
    // }
    // }
    // else
    // {
    for (int y = 0; y < maxy; y++) {
      for (int x = 0; x < maxx; x++) {
        double I = tb;
        for (int n = 0, xx = x, yy = y, a = 1; n < numberOfPeaks; n++, xx += maxx, yy += maxy) {
          duda[a] = deltaEx[xx] * deltaEy[yy];
          I += tI[n] * duda[a++];
          duda[a++] = duDtx[xx] * deltaEy[yy];
          duda[a++] = duDty[yy] * deltaEx[xx];
          duda[a++] = duDtsx[xx] * deltaEy[yy];
          duda[a++] = duDtsy[yy] * deltaEx[xx];
        }
        procedure.execute(I, duda);
      }
    }
  }

  @Override
  public void forEach(Gradient2Procedure procedure) {
    final double[] duda = new double[getNumberOfGradients()];
    final double[] d2uda2 = new double[getNumberOfGradients()];
    duda[0] = 1.0;
    for (int y = 0; y < maxy; y++) {
      for (int x = 0; x < maxx; x++) {
        double I = tb;
        for (int n = 0, xx = x, yy = y, a = 1; n < numberOfPeaks; n++, xx += maxx, yy += maxy) {
          duda[a] = deltaEx[xx] * deltaEy[yy];
          I += tI[n] * duda[a++];
          duda[a] = duDtx[xx] * deltaEy[yy];
          d2uda2[a++] = d2uDtx2[xx] * deltaEy[yy];
          duda[a] = duDty[yy] * deltaEx[xx];
          d2uda2[a++] = d2uDty2[yy] * deltaEx[xx];
          duda[a] = duDtsx[xx] * deltaEy[yy];
          d2uda2[a++] = d2uDtsx2[xx] * deltaEy[yy];
          duda[a] = duDtsy[yy] * deltaEx[xx];
          d2uda2[a++] = d2uDtsy2[yy] * deltaEx[xx];
        }
        procedure.execute(I, duda, d2uda2);
      }
    }
  }

  @Override
  public void forEach(ExtendedGradient2Procedure procedure) {
    final int ng = getNumberOfGradients();
    final double[] duda = new double[ng];
    final double[] d2udadb = new double[ng * ng];
    duda[0] = 1.0;
    final double[] du_dtsx_tI = new double[duDtsx.length];
    for (int x = 0; x < maxx; x++) {
      for (int n = 0, xx = x; n < numberOfPeaks; n++, xx += maxx) {
        du_dtsx_tI[xx] = duDtsx[xx] / tI[n];
      }
    }
    final double[] du_dty_tI = new double[numberOfPeaks];
    final double[] du_dtsy_tI = new double[numberOfPeaks];
    for (int y = 0; y < maxy; y++) {
      for (int n = 0, yy = y; n < numberOfPeaks; n++, yy += maxy) {
        du_dty_tI[n] = duDty[yy] / tI[n];
        du_dtsy_tI[n] = duDtsy[yy] / tI[n];
      }
      for (int x = 0; x < maxx; x++) {
        double I = tb;
        for (int n = 0, xx = x, yy = y, a = 1; n < numberOfPeaks; n++, xx += maxx, yy += maxy) {
          duda[a] = deltaEx[xx] * deltaEy[yy];
          I += tI[n] * duda[a];
          duda[a + 1] = duDtx[xx] * deltaEy[yy];
          duda[a + 2] = duDty[yy] * deltaEx[xx];
          duda[a + 3] = duDtsx[xx] * deltaEy[yy];
          duda[a + 4] = duDtsy[yy] * deltaEx[xx];

          // Compute all the partial second order derivatives
          final double tI = this.tI[n];

          // Background are all 0

          final int k = a * ng + a;
          // Signal,X
          d2udadb[k + 1] = duda[a + 1] / tI;
          // Signal,Y
          d2udadb[k + 2] = duda[a + 2] / tI;
          // Signal,X SD
          d2udadb[k + 3] = duda[a + 3] / tI;
          // Signal,Y SD
          d2udadb[k + 4] = duda[a + 4] / tI;

          a += 5;

          final int kk = k + ng;
          // X,Signal
          d2udadb[kk] = d2udadb[k + 1];
          // X,X
          d2udadb[kk + 1] = d2uDtx2[xx] * deltaEy[yy];
          // X,Y
          d2udadb[kk + 2] = duDtx[xx] * du_dty_tI[n];
          // X,X SD
          d2udadb[kk + 3] = deltaEy[yy] * d2deltaExDtsxDx[xx];
          // X,Y SD
          d2udadb[kk + 4] = duDtx[xx] * du_dtsy_tI[n];

          final int kkk = kk + ng;
          // Y,Signal
          d2udadb[kkk] = d2udadb[k + 2];
          // Y,X
          d2udadb[kkk + 1] = d2udadb[kk + 2];
          // Y,Y
          d2udadb[kkk + 2] = d2uDty2[yy] * deltaEx[xx];
          // Y,X SD
          d2udadb[kkk + 3] = duDty[yy] * du_dtsx_tI[xx];
          // Y,Y SD
          d2udadb[kkk + 4] = deltaEx[xx] * d2deltaEyDtsyDy[yy];

          final int kkkk = kkk + ng;
          // X SD,Signal
          d2udadb[kkkk] = d2udadb[k + 3];
          // X SD,X
          d2udadb[kkkk + 1] = d2udadb[kk + 3];
          // X SD,Y
          d2udadb[kkkk + 2] = d2udadb[kkk + 3];
          // X SD,X SD
          d2udadb[kkkk + 3] = d2uDtsx2[xx] * deltaEy[yy];
          // X SD,Y SD
          d2udadb[kkkk + 4] = duDtsy[yy] * du_dtsx_tI[xx];

          final int kkkkk = kkkk + ng;
          // Y SD,Signal
          d2udadb[kkkkk] = d2udadb[k + 4];
          // Y SD,X
          d2udadb[kkkkk + 1] = d2udadb[kk + 4];
          // Y SD,Y
          d2udadb[kkkkk + 2] = d2udadb[kkk + 4];
          // Y SD,X SD
          d2udadb[kkkkk + 3] = d2udadb[kkkk + 4];
          // Y SD,Y SD
          d2udadb[kkkkk + 4] = d2uDtsy2[yy] * deltaEx[xx];

        }
        procedure.executeExtended(I, duda, d2udadb);
      }
    }
  }
}
