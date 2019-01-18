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

package uk.ac.sussex.gdsc.smlm.function.gaussian.erf;

import uk.ac.sussex.gdsc.smlm.function.ExtendedGradient2Procedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Procedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient2Procedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;

/**
 * Evaluates a 2-dimensional Gaussian function for a single peak.
 */
public class MultiCircularErfGaussian2DFunction extends MultiFreeCircularErfGaussian2DFunction {
  /**
   * Constructor.
   *
   * @param numberOfPeaks The number of peaks
   * @param maxx The maximum x value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   * @param maxy The maximum y value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   */
  public MultiCircularErfGaussian2DFunction(int numberOfPeaks, int maxx, int maxy) {
    super(numberOfPeaks, maxx, maxy);
  }

  @Override
  protected int[] createGradientIndices() {
    return replicateGradientIndices(SingleCircularErfGaussian2DFunction.gradientIndices);
  }

  @Override
  public ErfGaussian2DFunction copy() {
    return new MultiCircularErfGaussian2DFunction(numberOfPeaks, maxx, maxy);
  }

  @Override
  public void initialise0(double[] a) {
    tb = a[Gaussian2DFunction.BACKGROUND];
    for (int n = 0, i = 0; n < numberOfPeaks; n++, i += PARAMETERS_PER_PEAK) {
      tI[n] = a[i + Gaussian2DFunction.SIGNAL];
      // Pre-compute the offset by 0.5
      final double tx = a[i + Gaussian2DFunction.X_POSITION] + 0.5;
      final double ty = a[i + Gaussian2DFunction.Y_POSITION] + 0.5;
      final double s = Math.abs(a[i + Gaussian2DFunction.X_SD]);
      final double one_sSqrt2 = ONE_OVER_ROOT2 / s;

      createDeltaETable(n, maxx, one_sSqrt2, deltaEx, tx);
      createDeltaETable(n, maxy, one_sSqrt2, deltaEy, ty);
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
      final double s = Math.abs(a[i + Gaussian2DFunction.X_SD]);
      final double one_sSqrt2 = ONE_OVER_ROOT2 / s;
      sum += tI * compute1DIntegral(one_sSqrt2, maxx, tx) * compute1DIntegral(one_sSqrt2, maxy, ty);
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
      final double s = Math.abs(a[i + Gaussian2DFunction.X_SD]);

      // We can pre-compute part of the derivatives for position and sd in arrays
      // since the Gaussian is XY separable
      final double one_sSqrt2 = ONE_OVER_ROOT2 / s;
      final double one_2ss = 0.5 / (s * s);
      final double I_sSqrt2pi = tI[n] * ONE_OVER_ROOT2PI / s;
      final double I_ssSqrt2pi = tI[n] * ONE_OVER_ROOT2PI / (s * s);

      // We can pre-compute part of the derivatives for position and sd in arrays
      // since the Gaussian is XY separable
      createFirstOrderTables(n, maxx, one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi, deltaEx, duDtx,
          duDtsx, tx);
      createFirstOrderTables(n, maxy, one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi, deltaEy, duDty,
          duDtsy, ty);
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
      final double one_sSqrt2pi = ONE_OVER_ROOT2PI / s;
      final double ss = s * s;
      final double one_sSqrt2 = ONE_OVER_ROOT2 / s;
      final double one_2ss = 0.5 / ss;
      final double I_sSqrt2pi = tI[n] * ONE_OVER_ROOT2PI / s;
      final double I_ssSqrt2pi = tI[n] * ONE_OVER_ROOT2PI / ss;
      final double I_sssSqrt2pi = I_sSqrt2pi / ss;
      final double one_sssSqrt2pi = one_sSqrt2pi / ss;
      final double one_sssssSqrt2pi = one_sssSqrt2pi / ss;

      // We can pre-compute part of the derivatives for position and sd in arrays
      // since the Gaussian is XY separable
      createSecondOrderTables(n, maxx, tI[n], one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi,
          I_sssSqrt2pi, ss, one_sssSqrt2pi, one_sssssSqrt2pi, deltaEx, duDtx, duDtsx, d2uDtx2,
          d2uDtsx2, tx);
      createSecondOrderTables(n, maxy, tI[n], one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi,
          I_sssSqrt2pi, ss, one_sssSqrt2pi, one_sssssSqrt2pi, deltaEy, duDty, duDtsy, d2uDty2,
          d2uDtsy2, ty);
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
      final double s = Math.abs(a[i + Gaussian2DFunction.X_SD]);

      // We can pre-compute part of the derivatives for position and sd in arrays
      // since the Gaussian is XY separable
      final double one_sSqrt2pi = ONE_OVER_ROOT2PI / s;
      final double ss = s * s;
      final double one_sSqrt2 = ONE_OVER_ROOT2 / s;
      final double one_2ss = 0.5 / ss;
      final double I_sSqrt2pi = tI[n] * ONE_OVER_ROOT2PI / s;
      final double I_ssSqrt2pi = tI[n] * ONE_OVER_ROOT2PI / ss;
      final double I_sssSqrt2pi = I_sSqrt2pi / ss;
      final double one_sssSqrt2pi = one_sSqrt2pi / ss;
      final double one_sssssSqrt2pi = one_sssSqrt2pi / ss;

      // We can pre-compute part of the derivatives for position and sd in arrays
      // since the Gaussian is XY separable
      createExSecondOrderTables(n, maxx, tI[n], one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi,
          I_sssSqrt2pi, ss, one_sssSqrt2pi, one_sssssSqrt2pi, deltaEx, duDtx, duDtsx, d2uDtx2,
          d2uDtsx2, d2deltaExDtsxDx, tx);
      createExSecondOrderTables(n, maxy, tI[n], one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi,
          I_sssSqrt2pi, ss, one_sssSqrt2pi, one_sssssSqrt2pi, deltaEy, duDty, duDtsy, d2uDty2,
          d2uDtsy2, d2deltaEyDtsyDy, ty);
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
      duda[a++] = duDtsx[xx] * deltaEy[yy] + duDtsy[yy] * deltaEx[xx];
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
      duda[a] = deltaEx[xx] * deltaEy[yy];
      value += tI[n] * duda[a];
      d2uda2[a++] = 0;
      duda[a] = duDtx[xx] * deltaEy[yy];
      d2uda2[a++] = d2uDtx2[xx] * deltaEy[yy];
      duda[a] = duDty[yy] * deltaEx[xx];
      d2uda2[a++] = d2uDty2[yy] * deltaEx[xx];
      duda[a] = duDtsx[xx] * deltaEy[yy] + duDtsy[yy] * deltaEx[xx];
      //@formatter:off
      d2uda2[a++] = d2uDtsx2[xx] * deltaEy[yy] +
                d2uDtsy2[yy] * deltaEx[xx] +
                2 * duDtsx[xx] * duDtsy[yy] / tI[n];
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
  public void forEach(Gradient1Procedure procedure) {
    final double[] duda = new double[getNumberOfGradients()];
    duda[0] = 1.0;
    for (int y = 0; y < maxy; y++) {
      for (int x = 0; x < maxx; x++) {
        double value = tb;
        for (int n = 0, xx = x, yy = y, a = 1; n < numberOfPeaks; n++, xx += maxx, yy += maxy) {
          duda[a] = deltaEx[xx] * deltaEy[yy];
          value += tI[n] * duda[a++];
          duda[a++] = duDtx[xx] * deltaEy[yy];
          duda[a++] = duDty[yy] * deltaEx[xx];
          duda[a++] = duDtsx[xx] * deltaEy[yy] + duDtsy[yy] * deltaEx[xx];
        }
        // invalidGradients(duda);
        procedure.execute(value, duda);
      }
    }
  }

  @Override
  public void forEach(Gradient2Procedure procedure) {
    final double[] duda = new double[getNumberOfGradients()];
    final double[] d2uda2 = new double[getNumberOfGradients()];
    final double[] two_du_dtsy_tI = new double[numberOfPeaks];
    duda[0] = 1.0;
    for (int y = 0; y < maxy; y++) {
      for (int n = 0, yy = y; n < numberOfPeaks; n++, yy += maxy) {
        two_du_dtsy_tI[n] = 2 * this.duDtsy[yy] / tI[n];
      }
      for (int x = 0; x < maxx; x++) {
        double value = tb;
        for (int n = 0, xx = x, yy = y, a = 1; n < numberOfPeaks; n++, xx += maxx, yy += maxy) {
          duda[a] = deltaEx[xx] * deltaEy[yy];
          value += tI[n] * duda[a++];
          duda[a] = duDtx[xx] * deltaEy[yy];
          d2uda2[a++] = d2uDtx2[xx] * deltaEy[yy];
          duda[a] = duDty[yy] * deltaEx[xx];
          d2uda2[a++] = d2uDty2[yy] * deltaEx[xx];
          duda[a] = duDtsx[xx] * deltaEy[yy] + duDtsy[yy] * deltaEx[xx];
          //@formatter:off
          d2uda2[a++] = d2uDtsx2[xx] * deltaEy[yy] +
                  d2uDtsy2[yy] * deltaEx[xx] +
                    duDtsx[xx] * two_du_dtsy_tI[n];
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
    final double[] du_dtsx_tI = new double[duDtsx.length];
    for (int x = 0; x < maxx; x++) {
      for (int n = 0, xx = x; n < numberOfPeaks; n++, xx += maxx) {
        du_dtsx_tI[xx] = duDtsx[xx] / tI[n];
      }
    }
    final double[] du_dty_tI = new double[numberOfPeaks];
    final double[] du_dtsy_tI = new double[numberOfPeaks];
    final double[] two_du_dtsy_tI = new double[numberOfPeaks];
    for (int y = 0; y < maxy; y++) {
      for (int n = 0, yy = y; n < numberOfPeaks; n++, yy += maxy) {
        du_dty_tI[n] = duDty[yy] / tI[n];
        du_dtsy_tI[n] = duDtsy[yy] / tI[n];
        two_du_dtsy_tI[n] = 2 * duDtsy[yy] / tI[n];
      }
      for (int x = 0; x < maxx; x++) {
        double value = tb;
        for (int n = 0, xx = x, yy = y, a = 1; n < numberOfPeaks; n++, xx += maxx, yy += maxy) {
          duda[a] = deltaEx[xx] * deltaEy[yy];
          value += tI[n] * duda[a];
          duda[a + 1] = duDtx[xx] * deltaEy[yy];
          duda[a + 2] = duDty[yy] * deltaEx[xx];
          duda[a + 3] = duDtsx[xx] * deltaEy[yy] + duDtsy[yy] * deltaEx[xx];

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

          a += 4;

          final int kk = k + ng;
          // X,Signal
          d2udadb[kk] = d2udadb[k + 1];
          // X,X
          d2udadb[kk + 1] = d2uDtx2[xx] * deltaEy[yy];
          // X,Y
          d2udadb[kk + 2] = duDtx[xx] * du_dty_tI[n];
          // X,X SD
          d2udadb[kk + 3] = deltaEy[yy] * d2deltaExDtsxDx[xx] + duDtx[xx] * du_dtsy_tI[n];

          final int kkk = kk + ng;
          // Y,Signal
          d2udadb[kkk] = d2udadb[k + 2];
          // Y,X
          d2udadb[kkk + 1] = d2udadb[kk + 2];
          // Y,Y
          d2udadb[kkk + 2] = d2uDty2[yy] * deltaEx[xx];
          // Y,X SD
          d2udadb[kkk + 3] = duDty[yy] * du_dtsx_tI[xx] + deltaEx[xx] * d2deltaEyDtsyDy[yy];

          final int kkkk = kkk + ng;
          // X SD,Signal
          d2udadb[kkkk] = d2udadb[k + 3];
          // X SD,X
          d2udadb[kkkk + 1] = d2udadb[kk + 3];
          // X SD,Y
          d2udadb[kkkk + 2] = d2udadb[kkk + 3];
          // X SD,X SD
          //@formatter:off
          d2udadb[kkkk + 3] = d2uDtsx2[xx] * deltaEy[yy] +
                            d2uDtsy2[yy] * deltaEx[xx] +
                            duDtsx[xx] * two_du_dtsy_tI[n];
            //@formatter:on
        }
        procedure.executeExtended(value, duda, d2udadb);
      }
    }
  }
}
