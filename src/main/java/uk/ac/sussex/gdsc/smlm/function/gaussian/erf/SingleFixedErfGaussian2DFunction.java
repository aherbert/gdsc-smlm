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

import uk.ac.sussex.gdsc.smlm.function.ExtendedGradient2Procedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Procedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient2Procedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.utils.StdMath;

/**
 * Evaluates a 2-dimensional Gaussian function for a single peak.
 */
public class SingleFixedErfGaussian2DFunction extends SingleCircularErfGaussian2DFunction {
  // Allow underscores in the variables used during computation
  // CHECKSTYLE.OFF: ParameterName|LocalVariableName

  /** The gradient indices. */
  static final int[] gradientIndices;

  static {
    gradientIndices = createGradientIndices(1, new SingleFixedErfGaussian2DFunction(1, 1));
  }

  /**
   * Constructor.
   *
   * @param maxx The maximum x value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   * @param maxy The maximum y value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   */
  public SingleFixedErfGaussian2DFunction(int maxx, int maxy) {
    super(maxx, maxy);
  }

  @Override
  protected void create1Arrays() {
    if (duDtx != null) {
      return;
    }
    duDtx = new double[this.maxx];
    duDty = new double[this.maxy];
  }

  @Override
  protected void create2Arrays() {
    if (d2uDtx2 != null) {
      return;
    }
    d2uDtx2 = new double[this.maxx];
    d2uDty2 = new double[this.maxy];
    create1Arrays();
  }

  @Override
  public ErfGaussian2DFunction copy() {
    return new SingleFixedErfGaussian2DFunction(maxx, maxy);
  }

  @Override
  public void initialise1(double[] a) {
    create1Arrays();
    tb = a[Gaussian2DFunction.BACKGROUND];
    tI = a[Gaussian2DFunction.SIGNAL];
    // Pre-compute the offset by 0.5
    final double tx = a[Gaussian2DFunction.X_POSITION] + 0.5;
    final double ty = a[Gaussian2DFunction.Y_POSITION] + 0.5;
    final double s = Math.abs(a[Gaussian2DFunction.X_SD]);

    // We can pre-compute part of the derivatives for position and sd in arrays
    // since the Gaussian is XY separable

    final double one_sSqrt2 = ONE_OVER_ROOT2 / s;
    final double one_2ss = 0.5 / (s * s);
    final double I_sSqrt2pi = tI * ONE_OVER_ROOT2PI / s;
    createFirstOrderTables(one_sSqrt2, one_2ss, I_sSqrt2pi, deltaEx, duDtx, tx);
    createFirstOrderTables(one_sSqrt2, one_2ss, I_sSqrt2pi, deltaEy, duDty, ty);
  }

  @Override
  public void initialise2(double[] a) {
    create2Arrays();
    tb = a[Gaussian2DFunction.BACKGROUND];
    tI = a[Gaussian2DFunction.SIGNAL];
    // Pre-compute the offset by 0.5
    final double tx = a[Gaussian2DFunction.X_POSITION] + 0.5;
    final double ty = a[Gaussian2DFunction.Y_POSITION] + 0.5;
    final double s = Math.abs(a[Gaussian2DFunction.X_SD]);

    // We can pre-compute part of the derivatives for position and sd in arrays
    // since the Gaussian is XY separable
    final double ss = s * s;
    final double one_sSqrt2 = ONE_OVER_ROOT2 / s;
    final double one_2ss = 0.5 / ss;
    final double I_sSqrt2pi = tI * ONE_OVER_ROOT2PI / s;
    final double I_sssSqrt2pi = I_sSqrt2pi / ss;
    createSecondOrderTables(one_sSqrt2, one_2ss, I_sSqrt2pi, I_sssSqrt2pi, deltaEx, duDtx, d2uDtx2,
        tx);
    createSecondOrderTables(one_sSqrt2, one_2ss, I_sSqrt2pi, I_sssSqrt2pi, deltaEy, duDty, d2uDty2,
        ty);
  }

  @Override
  public void initialiseExtended2(double[] a) {
    initialise2(a);
  }

  /**
   * Creates the first order derivatives.
   *
   * @param one_sSqrt2 one over (s times sqrt(2))
   * @param one_2ss one over (2 * s^2)
   * @param I_sSqrt2pi the intensity over (s * sqrt(2*pi))
   * @param deltaE the delta E for dimension 0 (difference between the error function at the start
   *        and end of each pixel)
   * @param duDx the first order x derivative for dimension 0
   * @param u the mean of the Gaussian for dimension 0
   */
  protected void createFirstOrderTables(double one_sSqrt2, double one_2ss, double I_sSqrt2pi,
      double[] deltaE, double[] duDx, double u) {
    // For documentation see SingleFreeCircularErfGaussian2DFunction.createSecondOrderTables(...)

    double x_u_p12 = -u;
    double erf_x_minus = 0.5 * erf(x_u_p12 * one_sSqrt2);
    double exp_x_minus = StdMath.exp(-(x_u_p12 * x_u_p12 * one_2ss));
    for (int i = 0, n = deltaE.length; i < n; i++) {
      x_u_p12 += 1.0;
      final double erf_x_plus = 0.5 * erf(x_u_p12 * one_sSqrt2);
      deltaE[i] = erf_x_plus - erf_x_minus;
      erf_x_minus = erf_x_plus;

      final double exp_x_plus = StdMath.exp(-(x_u_p12 * x_u_p12 * one_2ss));
      duDx[i] = I_sSqrt2pi * (exp_x_minus - exp_x_plus);

      exp_x_minus = exp_x_plus;
    }
  }

  /**
   * Creates the first and second order derivatives.
   *
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
  protected void createSecondOrderTables(double one_sSqrt2, double one_2ss, double I_sSqrt2pi,
      double I_sssSqrt2pi, double[] deltaE, double[] duDx, double[] d2uDx2, double u) {
    // For documentation see SingleFreeCircularErfGaussian2DFunction.createSecondOrderTables(...)

    double x_u_p12 = -u;
    double erf_x_minus = 0.5 * erf(x_u_p12 * one_sSqrt2);
    double exp_x_minus = StdMath.exp(-(x_u_p12 * x_u_p12 * one_2ss));
    for (int i = 0, n = deltaE.length; i < n; i++) {
      final double x_u_m12 = x_u_p12;
      x_u_p12 += 1.0;
      final double erf_x_plus = 0.5 * erf(x_u_p12 * one_sSqrt2);
      deltaE[i] = erf_x_plus - erf_x_minus;
      erf_x_minus = erf_x_plus;

      final double exp_x_plus = StdMath.exp(-(x_u_p12 * x_u_p12 * one_2ss));
      duDx[i] = I_sSqrt2pi * (exp_x_minus - exp_x_plus);
      d2uDx2[i] = I_sssSqrt2pi * (x_u_m12 * exp_x_minus - x_u_p12 * exp_x_plus);

      exp_x_minus = exp_x_plus;
    }
  }

  /**
   * Evaluates a 2-dimensional Gaussian function for a single peak.
   *
   * @param i Input predictor
   * @param duda Partial gradient of function with respect to each coefficient
   * @return The predicted value
   */
  @Override
  public double eval(final int i, final double[] duda) {
    // Unpack the predictor into the dimensions
    final int y = i / maxx;
    final int x = i % maxx;

    // Return in order of Gaussian2DFunction.createGradientIndices().
    // Use pre-computed gradients
    duda[0] = 1.0;
    duda[1] = deltaEx[x] * deltaEy[y];
    duda[2] = duDtx[x] * deltaEy[y];
    duda[3] = duDty[y] * deltaEx[x];

    return tb + tI * duda[1];
  }

  /**
   * Evaluates a 2-dimensional Gaussian function for a single peak.
   *
   * @param i Input predictor
   * @param duda Partial first gradient of function with respect to each coefficient
   * @param d2uda2 Partial second gradient of function with respect to each coefficient
   * @return The predicted value
   */
  @Override
  public double eval2(final int i, final double[] duda, final double[] d2uda2) {
    // Unpack the predictor into the dimensions
    final int y = i / maxx;
    final int x = i % maxx;

    // Return in order of Gaussian2DFunction.createGradientIndices().
    // Use pre-computed gradients
    duda[0] = 1.0;
    duda[1] = deltaEx[x] * deltaEy[y];
    duda[2] = duDtx[x] * deltaEy[y];
    duda[3] = duDty[y] * deltaEx[x];
    d2uda2[0] = 0;
    d2uda2[1] = 0;
    d2uda2[2] = d2uDtx2[x] * deltaEy[y];
    d2uda2[3] = d2uDty2[y] * deltaEx[x];

    return tb + tI * duda[1];
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
  public int[] gradientIndices() {
    return gradientIndices;
  }

  @Override
  public int getNumberOfGradients() {
    return 4;
  }

  @Override
  public void forEach(Gradient1Procedure procedure) {
    final double[] duda = new double[getNumberOfGradients()];
    duda[0] = 1.0;
    for (int y = 0; y < maxy; y++) {
      final double deltaEy = this.deltaEy[y];
      final double du_dty = this.duDty[y];
      for (int x = 0; x < maxx; x++) {
        duda[1] = deltaEx[x] * deltaEy;
        duda[2] = duDtx[x] * deltaEy;
        duda[3] = du_dty * deltaEx[x];
        procedure.execute(tb + tI * duda[1], duda);
      }
    }
  }

  @Override
  public void forEach(Gradient2Procedure procedure) {
    final double[] duda = new double[getNumberOfGradients()];
    final double[] d2uda2 = new double[getNumberOfGradients()];
    duda[0] = 1.0;
    for (int y = 0; y < maxy; y++) {
      final double deltaEy = this.deltaEy[y];
      final double du_dty = this.duDty[y];
      final double d2u_dty2 = this.d2uDty2[y];
      for (int x = 0; x < maxx; x++) {
        duda[1] = deltaEx[x] * deltaEy;
        duda[2] = duDtx[x] * deltaEy;
        duda[3] = du_dty * deltaEx[x];
        d2uda2[2] = d2uDtx2[x] * deltaEy;
        d2uda2[3] = d2u_dty2 * deltaEx[x];
        procedure.execute(tb + tI * duda[1], duda, d2uda2);
      }
    }
  }

  @Override
  public void forEach(ExtendedGradient2Procedure procedure) {
    final int n = getNumberOfGradients();
    final double[] duda = new double[n];
    final double[] d2udadb = new double[n * n];
    duda[0] = 1.0;
    for (int y = 0; y < maxy; y++) {
      final double deltaEy = this.deltaEy[y];
      final double du_dty = this.duDty[y];
      final double du_dty_tI = du_dty / tI;
      final double d2u_dty2 = this.d2uDty2[y];
      for (int x = 0; x < maxx; x++) {
        duda[1] = deltaEx[x] * deltaEy;
        duda[2] = duDtx[x] * deltaEy;
        duda[3] = du_dty * deltaEx[x];

        // Compute all the partial second order derivatives

        // Background are all 0

        // Signal,X
        d2udadb[6] = duda[2] / tI;
        // Signal,Y
        d2udadb[7] = duda[3] / tI;

        // X,Signal
        d2udadb[9] = d2udadb[6];
        // X,X
        d2udadb[10] = d2uDtx2[x] * deltaEy;
        // X,Y
        d2udadb[11] = duDtx[x] * du_dty_tI;

        // Y,Signal
        d2udadb[13] = d2udadb[7];
        // Y,X
        d2udadb[14] = d2udadb[11];
        // Y,Y
        d2udadb[15] = d2u_dty2 * deltaEx[x];

        procedure.executeExtended(tb + tI * duda[1], duda, d2udadb);
      }
    }
  }
}
