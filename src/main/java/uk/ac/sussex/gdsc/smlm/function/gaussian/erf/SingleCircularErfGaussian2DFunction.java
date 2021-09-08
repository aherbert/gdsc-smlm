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

package uk.ac.sussex.gdsc.smlm.function.gaussian.erf;

import uk.ac.sussex.gdsc.smlm.function.ExtendedGradient2Procedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Procedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient2Procedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;

/**
 * Evaluates a 2-dimensional Gaussian function for a single peak.
 */
public class SingleCircularErfGaussian2DFunction extends SingleFreeCircularErfGaussian2DFunction {
  /** The gradient indices. */
  static final int[] gradientIndices;

  static {
    gradientIndices = createGradientIndices(1, new SingleCircularErfGaussian2DFunction(1, 1));
  }

  /**
   * Constructor.
   *
   * @param maxx The maximum x value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   * @param maxy The maximum y value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   */
  public SingleCircularErfGaussian2DFunction(int maxx, int maxy) {
    super(maxx, maxy);
  }

  @Override
  public ErfGaussian2DFunction copy() {
    return new SingleCircularErfGaussian2DFunction(maxx, maxy);
  }

  @Override
  public void initialise0(double[] a) {
    tb = a[Gaussian2DFunction.BACKGROUND];
    tI = a[Gaussian2DFunction.SIGNAL];
    // Pre-compute the offset by 0.5
    final double tx = a[Gaussian2DFunction.X_POSITION] + 0.5;
    final double ty = a[Gaussian2DFunction.Y_POSITION] + 0.5;
    final double s = Math.abs(a[Gaussian2DFunction.X_SD]);

    final double one_sSqrt2 = ONE_OVER_ROOT2 / s;
    createDeltaETable(one_sSqrt2, deltaEx, tx);
    createDeltaETable(one_sSqrt2, deltaEy, ty);
  }

  @Override
  public double integral(double[] a) {
    final double tb = a[Gaussian2DFunction.BACKGROUND];
    final double tI = a[Gaussian2DFunction.SIGNAL];
    // Pre-compute the offset by 0.5
    final double tx = a[Gaussian2DFunction.X_POSITION] + 0.5;
    final double ty = a[Gaussian2DFunction.Y_POSITION] + 0.5;
    final double s = Math.abs(a[Gaussian2DFunction.X_SD]);

    final double one_sSqrt2 = ONE_OVER_ROOT2 / s;
    return tb * size()
        + tI * compute1DIntegral(one_sSqrt2, maxx, tx) * compute1DIntegral(one_sSqrt2, maxy, ty);
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
    final double I_ssSqrt2pi = tI * ONE_OVER_ROOT2PI / (s * s);
    createFirstOrderTables(one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi, deltaEx, duDtx, duDtsx,
        tx);
    createFirstOrderTables(one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi, deltaEy, duDty, duDtsy,
        ty);
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
    final double one_sSqrt2pi = ONE_OVER_ROOT2PI / s;
    final double ss = s * s;
    final double one_sSqrt2 = ONE_OVER_ROOT2 / s;
    final double one_2ss = 0.5 / ss;
    final double I_sSqrt2pi = tI * ONE_OVER_ROOT2PI / s;
    final double I_ssSqrt2pi = tI * ONE_OVER_ROOT2PI / ss;
    final double I_sssSqrt2pi = I_sSqrt2pi / ss;
    final double one_sssSqrt2pi = one_sSqrt2pi / ss;
    final double one_sssssSqrt2pi = one_sssSqrt2pi / ss;
    createSecondOrderTables(tI, one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi, I_sssSqrt2pi, ss,
        one_sssSqrt2pi, one_sssssSqrt2pi, deltaEx, duDtx, duDtsx, d2uDtx2, d2uDtsx2, tx);
    createSecondOrderTables(tI, one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi, I_sssSqrt2pi, ss,
        one_sssSqrt2pi, one_sssssSqrt2pi, deltaEy, duDty, duDtsy, d2uDty2, d2uDtsy2, ty);
  }

  @Override
  public void initialiseExtended2(double[] a) {
    createEx2Arrays();
    tb = a[Gaussian2DFunction.BACKGROUND];
    tI = a[Gaussian2DFunction.SIGNAL];
    // Pre-compute the offset by 0.5
    final double tx = a[Gaussian2DFunction.X_POSITION] + 0.5;
    final double ty = a[Gaussian2DFunction.Y_POSITION] + 0.5;
    final double s = Math.abs(a[Gaussian2DFunction.X_SD]);

    // We can pre-compute part of the derivatives for position and sd in arrays
    // since the Gaussian is XY separable
    final double one_sSqrt2pi = ONE_OVER_ROOT2PI / s;
    final double ss = s * s;
    final double one_sSqrt2 = ONE_OVER_ROOT2 / s;
    final double one_2ss = 0.5 / ss;
    final double I_sSqrt2pi = tI * ONE_OVER_ROOT2PI / s;
    final double I_ssSqrt2pi = tI * ONE_OVER_ROOT2PI / ss;
    final double I_sssSqrt2pi = I_sSqrt2pi / ss;
    final double one_sssSqrt2pi = one_sSqrt2pi / ss;
    final double one_sssssSqrt2pi = one_sssSqrt2pi / ss;
    createExSecondOrderTables(tI, one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi, I_sssSqrt2pi, ss,
        one_sssSqrt2pi, one_sssssSqrt2pi, deltaEx, duDtx, duDtsx, d2uDtx2, d2uDtsx2,
        d2deltaExDtsxDx, tx);
    createExSecondOrderTables(tI, one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi, I_sssSqrt2pi, ss,
        one_sssSqrt2pi, one_sssssSqrt2pi, deltaEy, duDty, duDtsy, d2uDty2, d2uDtsy2,
        d2deltaEyDtsyDy, ty);
  }

  @Override
  public double eval(final int x, final double[] duda) {
    // Unpack the predictor into the dimensions
    final int yy = x / maxx;
    final int xx = x % maxx;

    // Return in order of Gaussian2DFunction.createGradientIndices().
    // Use pre-computed gradients
    duda[0] = 1.0;
    duda[1] = deltaEx[xx] * deltaEy[yy];
    duda[2] = duDtx[xx] * deltaEy[yy];
    duda[3] = duDty[yy] * deltaEx[xx];
    duda[4] = duDtsx[xx] * deltaEy[yy] + duDtsy[yy] * deltaEx[xx];

    return tb + tI * duda[1];
  }

  @Override
  public double eval2(final int x, final double[] duda, final double[] d2uda2) {
    // Unpack the predictor into the dimensions
    final int yy = x / maxx;
    final int xx = x % maxx;

    // Return in order of Gaussian2DFunction.createGradientIndices().
    // Use pre-computed gradients
    duda[0] = 1.0;
    duda[1] = deltaEx[xx] * deltaEy[yy];
    duda[2] = duDtx[xx] * deltaEy[yy];
    duda[3] = duDty[yy] * deltaEx[xx];
    duda[4] = duDtsx[xx] * deltaEy[yy] + duDtsy[yy] * deltaEx[xx];
    d2uda2[0] = 0;
    d2uda2[1] = 0;
    d2uda2[2] = d2uDtx2[xx] * deltaEy[yy];
    d2uda2[3] = d2uDty2[yy] * deltaEx[xx];
    // Working example of this in GraspJ source code:
    // https://github.com/isman7/graspj/blob/master/graspj/src/main/java/eu/brede/graspj/opencl/src/functions/psfmodel_derivatives_sigma.cl
    //@formatter:off
    d2uda2[4] =
        d2uDtsx2[xx] * deltaEy[yy] +
        d2uDtsy2[yy] * deltaEx[xx] +
        2 * duDtsx[xx] * duDtsy[yy] / tI;
    //@formatter:on

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
    return true;
  }

  @Override
  public boolean evaluatesSD1() {
    return false;
  }

  @Override
  public int getGradientParametersPerPeak() {
    return 4;
  }

  @Override
  public int[] gradientIndices() {
    return gradientIndices;
  }

  @Override
  public int getNumberOfGradients() {
    return 5;
  }

  @Override
  public void forEach(Gradient1Procedure procedure) {
    final double[] duda = new double[getNumberOfGradients()];
    duda[0] = 1.0;
    for (int y = 0; y < maxy; y++) {
      final double du_dty = this.duDty[y];
      final double deltaEy = this.deltaEy[y];
      final double du_dtsy = this.duDtsy[y];
      for (int x = 0; x < maxx; x++) {
        duda[1] = deltaEx[x] * deltaEy;
        duda[2] = duDtx[x] * deltaEy;
        duda[3] = du_dty * deltaEx[x];
        duda[4] = duDtsx[x] * deltaEy + du_dtsy * deltaEx[x];
        // invalidGradients(duda);
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
      final double du_dty = this.duDty[y];
      final double deltaEy = this.deltaEy[y];
      final double du_dtsy = this.duDtsy[y];
      final double two_du_dtsy_tI = 2 * this.duDtsy[y] / tI;
      final double d2u_dty2 = this.d2uDty2[y];
      final double d2u_dtsy2 = this.d2uDtsy2[y];
      for (int x = 0; x < maxx; x++) {
        duda[1] = deltaEx[x] * deltaEy;
        duda[2] = duDtx[x] * deltaEy;
        duda[3] = du_dty * deltaEx[x];
        duda[4] = duDtsx[x] * deltaEy + du_dtsy * deltaEx[x];
        d2uda2[2] = d2uDtx2[x] * deltaEy;
        d2uda2[3] = d2u_dty2 * deltaEx[x];
        //@formatter:off
        d2uda2[4] =
            d2uDtsx2[x] * deltaEy +
            d2u_dtsy2 * deltaEx[x] +
            duDtsx[x] * two_du_dtsy_tI;
        //@formatter:on
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
    final double[] du_dtsx_tI = new double[maxx];
    for (int x = 0; x < maxx; x++) {
      du_dtsx_tI[x] = duDtsx[x] / tI;
    }
    for (int y = 0; y < maxy; y++) {
      final double du_dty = this.duDty[y];
      final double du_dty_tI = du_dty / tI;
      final double deltaEy = this.deltaEy[y];
      final double du_dtsy = this.duDtsy[y];
      final double du_dtsy_tI = du_dtsy / tI;
      final double two_du_dtsy_tI = 2 * this.duDtsy[y] / tI;
      final double d2u_dty2 = this.d2uDty2[y];
      final double d2u_dtsy2 = this.d2uDtsy2[y];
      final double d2deltaEy_dtsydy = this.d2deltaEyDtsyDy[y];
      for (int x = 0; x < maxx; x++) {
        duda[1] = deltaEx[x] * deltaEy;
        duda[2] = duDtx[x] * deltaEy;
        duda[3] = du_dty * deltaEx[x];
        duda[4] = duDtsx[x] * deltaEy + du_dtsy * deltaEx[x];

        // Compute all the partial second order derivatives

        // Background are all 0

        // Signal,X
        d2udadb[7] = duda[2] / tI;
        // Signal,Y
        d2udadb[8] = duda[3] / tI;
        // Signal,X SD
        d2udadb[9] = duda[4] / tI;

        // X,Signal
        d2udadb[11] = d2udadb[7];
        // X,X
        d2udadb[12] = d2uDtx2[x] * deltaEy;
        // X,Y
        d2udadb[13] = duDtx[x] * du_dty_tI;
        // X,X SD
        d2udadb[14] = deltaEy * d2deltaExDtsxDx[x] + duDtx[x] * du_dtsy_tI;

        // Y,Signal
        d2udadb[16] = d2udadb[8];
        // Y,X
        d2udadb[17] = d2udadb[13];
        // Y,Y
        d2udadb[18] = d2u_dty2 * deltaEx[x];
        // Y,X SD
        d2udadb[19] = du_dty * du_dtsx_tI[x] + deltaEx[x] * d2deltaEy_dtsydy;

        // X SD,Signal
        d2udadb[21] = d2udadb[9];
        // X SD,X
        d2udadb[22] = d2udadb[14];
        // X SD,Y
        d2udadb[23] = d2udadb[19];
        // X SD,X SD
        //@formatter:off
        d2udadb[24] =
            d2uDtsx2[x] * deltaEy +
            d2u_dtsy2 * deltaEx[x] +
            duDtsx[x] * two_du_dtsy_tI;
        //@formatter:on

        procedure.executeExtended(tb + tI * duda[1], duda, d2udadb);
      }
    }
  }
}
