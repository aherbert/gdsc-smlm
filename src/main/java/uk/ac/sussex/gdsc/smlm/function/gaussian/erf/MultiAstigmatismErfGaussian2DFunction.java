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
import uk.ac.sussex.gdsc.smlm.function.gaussian.AstigmatismZModel;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;

/**
 * Evaluates a 2-dimensional Gaussian function for a single peak.
 */
public class MultiAstigmatismErfGaussian2DFunction extends MultiFreeCircularErfGaussian2DFunction {
  /** The z model. */
  protected final AstigmatismZModel zModel;

  // Required for the z-depth gradients

  /** The x|z pre-factors for first-order partial derivatives. */
  protected double[] dtsxDtz;
  /** The x|z pre-factors for first-order partial derivatives. */
  protected double[] d2tsxDtz2;
  /** The y|z pre-factors for second-order partial derivatives. */
  protected double[] dtsyDtz;
  /** The y|z pre-factors for second-order partial derivatives. */
  protected double[] d2tsyDtz2;

  /**
   * Constructor.
   *
   * @param numberOfPeaks The number of peaks
   * @param maxx The maximum x value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   * @param maxy The maximum y value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   * @param zModel the z model
   */
  public MultiAstigmatismErfGaussian2DFunction(int numberOfPeaks, int maxx, int maxy,
      AstigmatismZModel zModel) {
    super(numberOfPeaks, maxx, maxy);
    this.zModel = zModel;
  }

  @Override
  protected void create1Arrays() {
    if (duDtx != null) {
      return;
    }
    duDtx = new double[deltaEx.length];
    duDty = new double[deltaEy.length];
    duDtsx = new double[deltaEx.length];
    duDtsy = new double[deltaEy.length];
    dtsxDtz = new double[deltaEx.length];
    dtsyDtz = new double[deltaEy.length];
  }

  @Override
  protected void create2Arrays() {
    if (d2uDtx2 != null) {
      return;
    }
    d2uDtx2 = new double[deltaEx.length];
    d2uDty2 = new double[deltaEy.length];
    d2uDtsx2 = new double[deltaEx.length];
    d2uDtsy2 = new double[deltaEy.length];
    d2tsxDtz2 = new double[deltaEx.length];
    d2tsyDtz2 = new double[deltaEy.length];
    create1Arrays();
  }

  @Override
  protected int[] createGradientIndices() {
    return replicateGradientIndices(SingleAstigmatismErfGaussian2DFunction.gradientIndices);
  }

  @Override
  public ErfGaussian2DFunction copy() {
    return new MultiAstigmatismErfGaussian2DFunction(numberOfPeaks, maxx, maxy, zModel);
  }

  @Override
  public void initialise0(double[] a) {
    tb = a[Gaussian2DFunction.BACKGROUND];
    for (int n = 0, i = 0; n < numberOfPeaks; n++, i += PARAMETERS_PER_PEAK) {
      tI[n] = a[i + Gaussian2DFunction.SIGNAL];
      // Pre-compute the offset by 0.5
      final double tx = a[i + Gaussian2DFunction.X_POSITION] + 0.5;
      final double ty = a[i + Gaussian2DFunction.Y_POSITION] + 0.5;
      final double tz = a[i + Gaussian2DFunction.Z_POSITION];

      final double sx = zModel.getSx(tz);
      final double sy = zModel.getSy(tz);
      createDeltaETable(n, maxx, ONE_OVER_ROOT2 / sx, deltaEx, tx);
      createDeltaETable(n, maxy, ONE_OVER_ROOT2 / sy, deltaEy, ty);
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
      final double tz = a[i + Gaussian2DFunction.Z_POSITION];

      final double sx = zModel.getSx(tz);
      final double sy = zModel.getSy(tz);
      sum += tI * compute1DIntegral(ONE_OVER_ROOT2 / sx, maxx, tx)
          * compute1DIntegral(ONE_OVER_ROOT2 / sy, maxy, ty);
    }
    return sum;
  }

  @Override
  public void initialise1(double[] a) {
    create1Arrays();
    final double[] dsDz = new double[1];
    tb = a[Gaussian2DFunction.BACKGROUND];
    for (int n = 0, i = 0; n < numberOfPeaks; n++, i += PARAMETERS_PER_PEAK) {
      tI[n] = a[i + Gaussian2DFunction.SIGNAL];
      // Pre-compute the offset by 0.5
      final double tx = a[i + Gaussian2DFunction.X_POSITION] + 0.5;
      final double ty = a[i + Gaussian2DFunction.Y_POSITION] + 0.5;
      final double tz = a[i + Gaussian2DFunction.Z_POSITION];

      // We can pre-compute part of the derivatives for position and sd in arrays
      // since the Gaussian is XY separable
      final double sx = zModel.getSx(tz, dsDz);
      dtsxDtz[n] = dsDz[0];
      final double sy = zModel.getSy(tz, dsDz);
      dtsyDtz[n] = dsDz[0];
      createFirstOrderTables(n, maxx, tI[n], deltaEx, duDtx, duDtsx, tx, sx);
      createFirstOrderTables(n, maxy, tI[n], deltaEy, duDty, duDtsy, ty, sy);
    }
  }

  @Override
  public void initialise2(double[] a) {
    create2Arrays();
    final double[] dsDz = new double[2];
    tb = a[Gaussian2DFunction.BACKGROUND];
    for (int n = 0, i = 0; n < numberOfPeaks; n++, i += PARAMETERS_PER_PEAK) {
      tI[n] = a[i + Gaussian2DFunction.SIGNAL];
      // Pre-compute the offset by 0.5
      final double tx = a[i + Gaussian2DFunction.X_POSITION] + 0.5;
      final double ty = a[i + Gaussian2DFunction.Y_POSITION] + 0.5;
      final double tz = a[i + Gaussian2DFunction.Z_POSITION];

      // We can pre-compute part of the derivatives for position and sd in arrays
      // since the Gaussian is XY separable
      final double sx = zModel.getSx2(tz, dsDz);
      dtsxDtz[n] = dsDz[0];
      d2tsxDtz2[n] = dsDz[1];
      final double sy = zModel.getSy2(tz, dsDz);
      dtsyDtz[n] = dsDz[0];
      d2tsyDtz2[n] = dsDz[1];
      createSecondOrderTables(n, maxx, tI[n], deltaEx, duDtx, duDtsx, d2uDtx2, d2uDtsx2, tx, sx);
      createSecondOrderTables(n, maxy, tI[n], deltaEy, duDty, duDtsy, d2uDty2, d2uDtsy2, ty, sy);
    }
  }

  @Override
  public void initialiseExtended2(double[] a) {
    createEx2Arrays();
    final double[] dsDz = new double[2];
    tb = a[Gaussian2DFunction.BACKGROUND];
    for (int n = 0, i = 0; n < numberOfPeaks; n++, i += PARAMETERS_PER_PEAK) {
      tI[n] = a[i + Gaussian2DFunction.SIGNAL];
      // Pre-compute the offset by 0.5
      final double tx = a[i + Gaussian2DFunction.X_POSITION] + 0.5;
      final double ty = a[i + Gaussian2DFunction.Y_POSITION] + 0.5;
      final double tz = a[i + Gaussian2DFunction.Z_POSITION];

      // We can pre-compute part of the derivatives for position and sd in arrays
      // since the Gaussian is XY separable
      final double sx = zModel.getSx2(tz, dsDz);
      dtsxDtz[n] = dsDz[0];
      d2tsxDtz2[n] = dsDz[1];
      final double sy = zModel.getSy2(tz, dsDz);
      dtsyDtz[n] = dsDz[0];
      d2tsyDtz2[n] = dsDz[1];
      createExSecondOrderTables(n, maxx, tI[n], deltaEx, duDtx, duDtsx, d2uDtx2, d2uDtsx2,
          d2deltaExDtsxDx, tx, sx);
      createExSecondOrderTables(n, maxy, tI[n], deltaEy, duDty, duDtsy, d2uDty2, d2uDtsy2,
          d2deltaEyDtsyDy, ty, sy);
    }
    // Pre-apply the gradient mapping from width to z
    for (int x = 0; x < maxx; x++) {
      for (int n = 0, xx = x; n < numberOfPeaks; n++, xx += maxx) {
        d2deltaExDtsxDx[xx] *= dtsxDtz[n];
      }
    }
    for (int y = 0; y < maxy; y++) {
      for (int n = 0, yy = y; n < numberOfPeaks; n++, yy += maxy) {
        d2deltaEyDtsyDy[yy] *= dtsyDtz[n];
      }
    }
  }

  @Override
  public double eval(final int x, final double[] duda) {
    // Unpack the predictor into the dimensions
    int yy = x / maxx;
    int xx = x % maxx;

    // Return in order of Gaussian2DFunction.createGradientIndices().
    // Use pre-computed gradients
    duda[0] = 1.0;
    double value = tb;
    for (int n = 0, a = 1; n < numberOfPeaks; n++, xx += maxx, yy += maxy) {
      duda[a] = deltaEx[xx] * deltaEy[yy];
      value += tI[n] * duda[a++];
      duda[a++] = duDtx[xx] * deltaEy[yy];
      duda[a++] = duDty[yy] * deltaEx[xx];
      duda[a++] = duDtsx[xx] * deltaEy[yy] * dtsxDtz[n] + duDtsy[yy] * deltaEx[xx] * dtsyDtz[n];
    }
    return value;
  }

  @Override
  public double eval2(final int x, final double[] duda, final double[] d2uda2) {
    // Unpack the predictor into the dimensions
    int yy = x / maxx;
    int xx = x % maxx;

    // Return in order of Gaussian2DFunction.createGradientIndices().
    // Use pre-computed gradients
    duda[0] = 1.0;
    d2uda2[0] = 0;
    double value = tb;
    for (int n = 0, a = 1; n < numberOfPeaks; n++, xx += maxx, yy += maxy) {
      final double du_dsx = duDtsx[xx] * deltaEy[yy];
      final double du_dsy = duDtsy[yy] * deltaEx[xx];

      duda[a] = deltaEx[xx] * deltaEy[yy];
      value += tI[n] * duda[a];
      d2uda2[a++] = 0;
      duda[a] = duDtx[xx] * deltaEy[yy];
      d2uda2[a++] = d2uDtx2[xx] * deltaEy[yy];
      duda[a] = duDty[yy] * deltaEx[xx];
      d2uda2[a++] = d2uDty2[yy] * deltaEx[xx];
      duda[a] = du_dsx * dtsxDtz[n] + du_dsy * dtsyDtz[n];
      //@formatter:off
      d2uda2[a++] =
          d2uDtsx2[xx] * deltaEy[yy] * dtsxDtz[n] * dtsxDtz[n] +
          du_dsx * d2tsxDtz2[n] +
          d2uDtsy2[yy] * deltaEx[xx] * dtsyDtz[n] * dtsyDtz[n] +
          du_dsy * d2tsyDtz2[n] +
          // Add the equivalent term we add in the circular version.
          // Note: this is not in the Smith, et al (2010) paper but is
          // in the GraspJ source code and it works in JUnit tests.
          2 * duDtsx[xx] * dtsxDtz[n] * duDtsy[yy] * dtsyDtz[n] / tI[n];
      //@formatter:on
    }
    return value;
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
  public boolean evaluatesPosition() {
    return true;
  }

  @Override
  public boolean evaluatesZ() {
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
    return 4;
  }

  @Override
  public void forEach(Gradient1Procedure procedure) {
    final double[] duda = new double[getNumberOfGradients()];
    duda[0] = 1.0;
    final double[] deltaEy_by_dtsx_dtz = new double[numberOfPeaks];
    final double[] du_dtsy_by_dtsy_dtz = new double[numberOfPeaks];
    for (int y = 0; y < maxy; y++) {
      for (int n = 0, yy = y; n < numberOfPeaks; n++, yy += maxy) {
        deltaEy_by_dtsx_dtz[n] = deltaEy[yy] * dtsxDtz[n];
        du_dtsy_by_dtsy_dtz[n] = duDtsy[yy] * dtsyDtz[n];
      }

      for (int x = 0; x < maxx; x++) {
        double value = tb;
        for (int n = 0, xx = x, yy = y, a = 1; n < numberOfPeaks; n++, xx += maxx, yy += maxy) {
          duda[a] = deltaEx[xx] * deltaEy[yy];
          value += tI[n] * duda[a++];
          duda[a++] = duDtx[xx] * deltaEy[yy];
          duda[a++] = duDty[yy] * deltaEx[xx];
          duda[a++] = duDtsx[xx] * deltaEy_by_dtsx_dtz[n] + du_dtsy_by_dtsy_dtz[n] * deltaEx[xx];
        }
        procedure.execute(value, duda);
      }
    }
  }

  @Override
  public void forEach(Gradient2Procedure procedure) {
    final double[] duda = new double[getNumberOfGradients()];
    final double[] d2uda2 = new double[getNumberOfGradients()];
    duda[0] = 1.0;
    final double[] dtsx_dtz_2 = new double[numberOfPeaks];
    final double[] dtsy_dtz_2 = new double[numberOfPeaks];
    final double[] two_dtsx_dtz_by_dtsy_dtz_tI = new double[numberOfPeaks];
    for (int n = 0; n < numberOfPeaks; n++) {
      dtsx_dtz_2[n] = dtsxDtz[n] * dtsxDtz[n];
      dtsy_dtz_2[n] = dtsyDtz[n] * dtsyDtz[n];
      two_dtsx_dtz_by_dtsy_dtz_tI[n] = 2 * dtsxDtz[n] * dtsyDtz[n] / tI[n];
    }
    final double[] deltaEy_by_dtsx_dtz_2 = new double[numberOfPeaks];
    final double[] d2u_dtsy2_by_dtsy_dtz_2 = new double[numberOfPeaks];
    final double[] two_dtsx_dtz_by_du_dtsy_by_dtsy_dtz_tI = new double[numberOfPeaks];
    for (int y = 0; y < maxy; y++) {
      for (int n = 0, yy = y; n < numberOfPeaks; n++, yy += maxy) {
        deltaEy_by_dtsx_dtz_2[n] = deltaEy[yy] * dtsx_dtz_2[n];
        d2u_dtsy2_by_dtsy_dtz_2[n] = d2uDtsy2[yy] * dtsy_dtz_2[n];
        two_dtsx_dtz_by_du_dtsy_by_dtsy_dtz_tI[n] = two_dtsx_dtz_by_dtsy_dtz_tI[n] * duDtsy[yy];
      }

      for (int x = 0; x < maxx; x++) {
        double value = tb;
        for (int n = 0, xx = x, yy = y, a = 1; n < numberOfPeaks; n++, xx += maxx, yy += maxy) {
          final double du_dsx = duDtsx[xx] * deltaEy[yy];
          final double du_dsy = duDtsy[yy] * deltaEx[xx];

          duda[a] = deltaEx[xx] * deltaEy[yy];
          value += tI[n] * duda[a++];
          duda[a] = duDtx[xx] * deltaEy[yy];
          d2uda2[a++] = d2uDtx2[xx] * deltaEy[yy];
          duda[a] = duDty[yy] * deltaEx[xx];
          d2uda2[a++] = d2uDty2[yy] * deltaEx[xx];
          duda[a] = du_dsx * dtsxDtz[n] + du_dsy * dtsyDtz[n];
          //@formatter:off
          d2uda2[a++] =
              d2uDtsx2[xx] * deltaEy_by_dtsx_dtz_2[n] +
              du_dsx * d2tsxDtz2[n] +
              d2u_dtsy2_by_dtsy_dtz_2[n] * deltaEx[xx] +
              du_dsy * d2tsyDtz2[n] +
              // Add the equivalent term we add in the circular version.
              // Note: this is not in the Smith, et al (2010) paper but is
              // in the GraspJ source code and it works in JUnit tests.
              // 2 * du_dtsx[x] * dtsx_dtz * du_dtsy * dtsy_dtz / tI
              two_dtsx_dtz_by_du_dtsy_by_dtsy_dtz_tI[n] * duDtsx[xx];
          //@formatter:on
        }
        procedure.execute(value, duda, d2uda2);
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
    final double[] du_dtsy_by_dtsy_dtz_tI = new double[numberOfPeaks];
    final double[] du_dty_by_dtsx_dtz_tI = new double[numberOfPeaks];
    final double[] deltaEy_by_dtsx_dtz_2 = new double[numberOfPeaks];
    final double[] d2u_dtsy2_by_dtsy_dtz_2 = new double[numberOfPeaks];
    final double[] two_dtsx_dtz_by_du_dtsy_by_dtsy_dtz_tI = new double[numberOfPeaks];

    final double[] dtsx_dtz_2 = new double[numberOfPeaks];
    final double[] dtsy_dtz_2 = new double[numberOfPeaks];
    final double[] two_dtsx_dtz_by_dtsy_dtz_tI = new double[numberOfPeaks];
    final double[] dtsx_dtz_tI = new double[numberOfPeaks];
    final double[] dtsy_dtz_tI = new double[numberOfPeaks];
    for (int n = 0; n < numberOfPeaks; n++) {
      dtsx_dtz_2[n] = dtsxDtz[n] * dtsxDtz[n];
      dtsy_dtz_2[n] = dtsyDtz[n] * dtsyDtz[n];
      two_dtsx_dtz_by_dtsy_dtz_tI[n] = 2 * dtsxDtz[n] * dtsyDtz[n] / tI[n];
      dtsx_dtz_tI[n] = dtsxDtz[n] / tI[n];
      dtsy_dtz_tI[n] = dtsyDtz[n] / tI[n];
    }

    for (int y = 0; y < maxy; y++) {
      for (int n = 0, yy = y; n < numberOfPeaks; n++, yy += maxy) {
        du_dty_tI[n] = duDty[yy] / tI[n];
        du_dtsy_by_dtsy_dtz_tI[n] = duDtsy[yy] * dtsy_dtz_tI[n];
        du_dty_by_dtsx_dtz_tI[n] = duDty[yy] * dtsx_dtz_tI[n];
        deltaEy_by_dtsx_dtz_2[n] = deltaEy[yy] * dtsx_dtz_2[n];
        d2u_dtsy2_by_dtsy_dtz_2[n] = d2uDtsy2[yy] * dtsy_dtz_2[n];
        two_dtsx_dtz_by_du_dtsy_by_dtsy_dtz_tI[n] = two_dtsx_dtz_by_dtsy_dtz_tI[n] * duDtsy[yy];
      }
      for (int x = 0; x < maxx; x++) {
        double value = tb;
        for (int n = 0, xx = x, yy = y, a = 1; n < numberOfPeaks; n++, xx += maxx, yy += maxy) {
          final double du_dsx = duDtsx[xx] * deltaEy[yy];
          final double du_dsy = duDtsy[yy] * deltaEx[xx];

          duda[a] = deltaEx[xx] * deltaEy[yy];
          value += tI[n] * duda[a];
          duda[a + 1] = duDtx[xx] * deltaEy[yy];
          duda[a + 2] = duDty[yy] * deltaEx[xx];
          duda[a + 3] = du_dsx * dtsxDtz[n] + du_dsy * dtsyDtz[n];

          // Compute all the partial second order derivatives
          final double tI = this.tI[n];

          // Background are all 0

          final int k = a * ng + a;
          // Signal,X
          d2udadb[k + 1] = duda[a + 1] / tI;
          // Signal,Y
          d2udadb[k + 2] = duda[a + 2] / tI;
          // Signal,Z
          d2udadb[k + 3] = duda[a + 3] / tI;

          a += 4;

          final int kk = k + ng;
          // X,Signal
          d2udadb[kk] = d2udadb[k + 1];
          // X,X
          d2udadb[kk + 1] = d2uDtx2[xx] * deltaEy[yy];
          // X,Y
          d2udadb[kk + 2] = duDtx[xx] * du_dty_tI[n];
          // X,Z
          d2udadb[kk + 3] =
              deltaEy[yy] * d2deltaExDtsxDx[xx] + duDtx[xx] * du_dtsy_by_dtsy_dtz_tI[n];

          final int kkk = kk + ng;
          // Y,Signal
          d2udadb[kkk] = d2udadb[k + 2];
          // Y,X
          d2udadb[kkk + 1] = d2udadb[kk + 2];
          // Y,Y
          d2udadb[kkk + 2] = d2uDty2[yy] * deltaEx[xx];
          // X,Z
          d2udadb[kkk + 3] =
              duDtsx[xx] * du_dty_by_dtsx_dtz_tI[n] + deltaEx[xx] * d2deltaEyDtsyDy[yy];

          final int kkkk = kkk + ng;
          // Z,Signal
          d2udadb[kkkk] = d2udadb[k + 3];
          // Z,X
          d2udadb[kkkk + 1] = d2udadb[kk + 3];
          // Z,Y
          d2udadb[kkkk + 2] = d2udadb[kkk + 3];
          // Z,Z
          //@formatter:off
          d2udadb[kkkk + 3] =
              d2uDtsx2[xx] * deltaEy_by_dtsx_dtz_2[n] +
              du_dsx * d2tsxDtz2[n] +
              d2u_dtsy2_by_dtsy_dtz_2[n] * deltaEx[xx] +
              du_dsy * d2tsyDtz2[n] +
              // Add the equivalent term we add in the circular version.
              // Note: this is not in the Smith, et al (2010) paper but is
              // in the GraspJ source code and it works in JUnit tests.
              // 2 * du_dtsx[x] * dtsx_dtz * du_dtsy * dtsy_dtz / tI
              two_dtsx_dtz_by_du_dtsy_by_dtsy_dtz_tI[n] * duDtsx[xx];
          //@formatter:on
        }
        procedure.executeExtended(value, duda, d2udadb);
      }
    }
  }
}
