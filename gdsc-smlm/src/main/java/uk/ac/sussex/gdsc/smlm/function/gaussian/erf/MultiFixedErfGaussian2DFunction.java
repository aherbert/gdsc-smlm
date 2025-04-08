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

package uk.ac.sussex.gdsc.smlm.function.gaussian.erf;

import uk.ac.sussex.gdsc.smlm.function.ExtendedGradient2Procedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Procedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient2Procedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.utils.StdMath;

/**
 * Evaluates a 2-dimensional Gaussian function for a single peak.
 */
public class MultiFixedErfGaussian2DFunction extends MultiCircularErfGaussian2DFunction {
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
  public MultiFixedErfGaussian2DFunction(int numberOfPeaks, int maxx, int maxy) {
    super(numberOfPeaks, maxx, maxy);
  }

  @Override
  protected void create1Arrays() {
    if (duDtx != null) {
      return;
    }
    duDtx = new double[deltaEx.length];
    duDty = new double[deltaEy.length];
  }

  @Override
  protected void create2Arrays() {
    if (d2uDtx2 != null) {
      return;
    }
    d2uDtx2 = new double[deltaEx.length];
    d2uDty2 = new double[deltaEy.length];
    create1Arrays();
  }

  @Override
  protected int[] createGradientIndices() {
    return replicateGradientIndices(SingleFixedErfGaussian2DFunction.gradientIndices);
  }

  @Override
  public ErfGaussian2DFunction copy() {
    return new MultiFixedErfGaussian2DFunction(numberOfPeaks, maxx, maxy);
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
      final double s = Math.abs(a[i + Gaussian2DFunction.X_SD]);

      // We can pre-compute part of the derivatives for position and sd in arrays
      // since the Gaussian is XY separable

      final double one_sSqrt2 = ONE_OVER_ROOT2 / s;
      final double one_2ss = 0.5 / (s * s);
      final double I_sSqrt2pi = tI[n] * ONE_OVER_ROOT2PI / s;
      createFirstOrderTables(n, maxx, one_sSqrt2, one_2ss, I_sSqrt2pi, deltaEx, duDtx, tx);
      createFirstOrderTables(n, maxy, one_sSqrt2, one_2ss, I_sSqrt2pi, deltaEy, duDty, ty);
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
      final double s = Math.abs(a[i + Gaussian2DFunction.X_SD]);

      // We can pre-compute part of the derivatives for position and sd in arrays
      // since the Gaussian is XY separable
      final double ss = s * s;
      final double one_sSqrt2 = ONE_OVER_ROOT2 / s;
      final double one_2ss = 0.5 / ss;
      final double I_sSqrt2pi = tI[n] * ONE_OVER_ROOT2PI / s;
      final double I_sssSqrt2pi = I_sSqrt2pi / ss;
      createSecondOrderTables(n, maxx, one_sSqrt2, one_2ss, I_sSqrt2pi, I_sssSqrt2pi, deltaEx,
          duDtx, d2uDtx2, tx);
      createSecondOrderTables(n, maxy, one_sSqrt2, one_2ss, I_sSqrt2pi, I_sssSqrt2pi, deltaEy,
          duDty, d2uDty2, ty);
    }
  }

  @Override
  public void initialiseExtended2(double[] a) {
    initialise2(a);
  }

  /**
   * Creates the first order derivatives.
   *
   * @param n the peak number
   * @param max the maximum for the dimension
   * @param one_sSqrt2 one over (s times sqrt(2))
   * @param one_2ss one over (2 * s^2)
   * @param I_sSqrt2pi the intensity over (s * sqrt(2*pi))
   * @param deltaE the delta E for dimension 0 (difference between the error function at the start
   *        and end of each pixel)
   * @param duDx the first order x derivative for dimension 0
   * @param u the mean of the Gaussian for dimension 0
   */
  protected void createFirstOrderTables(int n, int max, double one_sSqrt2, double one_2ss,
      double I_sSqrt2pi, double[] deltaE, double[] duDx, double u) {
    // For documentation see SingleFreeCircularErfGaussian2DFunction.createSecondOrderTables(...)

    double x_u_p12 = -u;
    double erf_x_minus = 0.5 * erf(x_u_p12 * one_sSqrt2);
    double exp_x_minus = StdMath.exp(-(x_u_p12 * x_u_p12 * one_2ss));
    for (int i = 0, j = n * max; i < max; i++, j++) {
      x_u_p12 += 1.0;
      final double erf_x_plus = 0.5 * erf(x_u_p12 * one_sSqrt2);
      deltaE[j] = erf_x_plus - erf_x_minus;
      erf_x_minus = erf_x_plus;

      final double exp_x_plus = StdMath.exp(-(x_u_p12 * x_u_p12 * one_2ss));
      duDx[j] = I_sSqrt2pi * (exp_x_minus - exp_x_plus);

      exp_x_minus = exp_x_plus;
    }
  }

  /**
   * Creates the first and second order derivatives.
   *
   * @param n the peak number
   * @param max the maximum for the dimension
   * @param one_sSqrt2 one over (s times sqrt(2))
   * @param one_2ss one over (2 * s^2)
   * @param I_sSqrt2pi the intensity over (s * sqrt(2*pi))
   * @param I_sssSqrt2pi the intensity over (s^3 * sqrt(2*pi))
   * @param deltaE the delta E for dimension 0 (difference between the error function at the start
   *        and end of each pixel)
   * @param duDx the first order x derivative for dimension 0
   * @param d2uDx2 the second order x derivative for dimension 0
   * @param u the mean of the Gaussian for dimension 0
   */
  protected void createSecondOrderTables(int n, int max, double one_sSqrt2, double one_2ss,
      double I_sSqrt2pi, double I_sssSqrt2pi, double[] deltaE, double[] duDx, double[] d2uDx2,
      double u) {
    // For documentation see SingleFreeCircularErfGaussian2DFunction.createSecondOrderTables(...)

    double x_u_p12 = -u;
    double erf_x_minus = 0.5 * erf(x_u_p12 * one_sSqrt2);
    double exp_x_minus = StdMath.exp(-(x_u_p12 * x_u_p12 * one_2ss));
    for (int i = 0, j = n * max; i < max; i++, j++) {
      final double x_u_m12 = x_u_p12;
      x_u_p12 += 1.0;
      final double erf_x_plus = 0.5 * erf(x_u_p12 * one_sSqrt2);
      deltaE[j] = erf_x_plus - erf_x_minus;
      erf_x_minus = erf_x_plus;

      final double exp_x_plus = StdMath.exp(-(x_u_p12 * x_u_p12 * one_2ss));
      duDx[j] = I_sSqrt2pi * (exp_x_minus - exp_x_plus);
      d2uDx2[j] = I_sssSqrt2pi * (x_u_m12 * exp_x_minus - x_u_p12 * exp_x_plus);

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
    return false;
  }

  @Override
  public boolean evaluatesSD1() {
    return false;
  }

  @Override
  public int getGradientParametersPerPeak() {
    return 3;
  }

  @Override
  public void forEach(Gradient1Procedure procedure) {
    final double[] duda = new double[getNumberOfGradients()];
    duda[0] = 1.0;
    for (int y = 0; y < maxy; y++) {
      for (int x = 0; x < maxx; x++) {
        double I = tb;
        for (int n = 0, xx = x, yy = y, a = 1; n < numberOfPeaks; n++, xx += maxx, yy += maxy) {
          duda[a] = deltaEx[xx] * deltaEy[yy];
          I += tI[n] * duda[a++];
          duda[a++] = duDtx[xx] * deltaEy[yy];
          duda[a++] = duDty[yy] * deltaEx[xx];
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
    final double[] du_dty_tI = new double[numberOfPeaks];
    for (int y = 0; y < maxy; y++) {
      for (int n = 0, yy = y; n < numberOfPeaks; n++, yy += maxy) {
        du_dty_tI[n] = duDty[yy] / tI[n];
      }
      for (int x = 0; x < maxx; x++) {
        double I = tb;
        for (int n = 0, xx = x, yy = y, a = 1; n < numberOfPeaks; n++, xx += maxx, yy += maxy) {
          duda[a] = deltaEx[xx] * deltaEy[yy];
          I += tI[n] * duda[a];
          duda[a + 1] = duDtx[xx] * deltaEy[yy];
          duda[a + 2] = duDty[yy] * deltaEx[xx];

          // Compute all the partial second order derivatives
          final double tI = this.tI[n];

          // Background are all 0

          final int k = a * ng + a;
          // Signal,X
          d2udadb[k + 1] = duda[a + 1] / tI;
          // Signal,Y
          d2udadb[k + 2] = duda[a + 2] / tI;

          a += 3;

          final int kk = k + ng;
          // X,Signal
          d2udadb[kk] = d2udadb[k + 1];
          // X,X
          d2udadb[kk + 1] = d2uDtx2[xx] * deltaEy[yy];
          // X,Y
          d2udadb[kk + 2] = duDtx[xx] * du_dty_tI[n];

          final int kkk = kk + ng;
          // Y,Signal
          d2udadb[kkk] = d2udadb[k + 2];
          // Y,X
          d2udadb[kkk + 1] = d2udadb[kk + 2];
          // Y,Y
          d2udadb[kkk + 2] = d2uDty2[yy] * deltaEx[xx];
        }
        procedure.executeExtended(I, duda, d2udadb);
      }
    }
  }
}
