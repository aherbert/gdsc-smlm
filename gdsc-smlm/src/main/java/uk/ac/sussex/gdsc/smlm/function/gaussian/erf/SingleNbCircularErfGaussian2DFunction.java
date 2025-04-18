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
import uk.ac.sussex.gdsc.smlm.function.ValueProcedure;

/**
 * Evaluates a 2-dimensional Gaussian function for a single peak.
 */
public class SingleNbCircularErfGaussian2DFunction extends SingleCircularErfGaussian2DFunction {
  /** The gradient indices. */
  static final int[] gradientIndices;

  static {
    gradientIndices = createGradientIndices(1, new SingleNbCircularErfGaussian2DFunction(1, 1));
  }

  /**
   * Constructor.
   *
   * @param maxx The maximum x value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   * @param maxy The maximum y value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   */
  public SingleNbCircularErfGaussian2DFunction(int maxx, int maxy) {
    super(maxx, maxy);
  }

  @Override
  public ErfGaussian2DFunction copy() {
    return new SingleNbCircularErfGaussian2DFunction(maxx, maxy);
  }

  @Override
  public double eval(final int x, final double[] duda) {
    // Unpack the predictor into the dimensions
    final int yy = x / maxx;
    final int xx = x % maxx;

    // Return in order of Gaussian2DFunction.createGradientIndices().
    // Use pre-computed gradients
    duda[0] = deltaEx[xx] * deltaEy[yy];
    duda[1] = duDtx[xx] * deltaEy[yy];
    duda[2] = duDty[yy] * deltaEx[xx];
    duda[3] = duDtsx[xx] * deltaEy[yy] + duDtsy[yy] * deltaEx[xx];

    return tb + tI * duda[0];
  }

  @Override
  public double eval2(final int x, final double[] duda, final double[] d2uda2) {
    // Unpack the predictor into the dimensions
    final int yy = x / maxx;
    final int xx = x % maxx;

    // Return in order of Gaussian2DFunction.createGradientIndices().
    // Use pre-computed gradients
    duda[0] = deltaEx[xx] * deltaEy[yy];
    duda[1] = duDtx[xx] * deltaEy[yy];
    duda[2] = duDty[yy] * deltaEx[xx];
    duda[3] = duDtsx[xx] * deltaEy[yy] + duDtsy[yy] * deltaEx[xx];
    d2uda2[0] = 0;
    d2uda2[1] = d2uDtx2[xx] * deltaEy[yy];
    d2uda2[2] = d2uDty2[yy] * deltaEx[xx];
    // Working example of this in GraspJ source code:
    // https://github.com/isman7/graspj/blob/master/graspj/src/main/java/eu/brede/graspj/opencl/src/functions/psfmodel_derivatives_sigma.cl
    //@formatter:off
    d2uda2[3] =
        d2uDtsx2[xx] * deltaEy[yy] +
        d2uDtsy2[yy] * deltaEx[xx] +
        2 * duDtsx[xx] * duDtsy[yy] / tI;
    //@formatter:on

    return tb + tI * duda[0];
  }

  @Override
  public boolean evaluatesBackground() {
    return false;
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
    return 4;
  }

  @Override
  public void forEach(ValueProcedure procedure) {
    if (tb == 0) {
      // Specialised implementation without a background.
      // (This function is likely to be used to compute the Gaussian integral
      // without a background.)
      for (int y = 0; y < maxy; y++) {
        final double tI_deltaEy = tI * deltaEy[y];
        for (int x = 0; x < maxx; x++) {
          procedure.execute(tI_deltaEy * deltaEx[x]);
        }
      }
    } else {
      super.forEach(procedure);
    }
  }

  @Override
  public void forEach(Gradient1Procedure procedure) {
    final double[] duda = new double[getNumberOfGradients()];
    for (int y = 0; y < maxy; y++) {
      final double du_dty = this.duDty[y];
      final double deltaEy = this.deltaEy[y];
      final double du_dtsy = this.duDtsy[y];
      for (int x = 0; x < maxx; x++) {
        duda[0] = deltaEx[x] * deltaEy;
        duda[1] = duDtx[x] * deltaEy;
        duda[2] = du_dty * deltaEx[x];
        duda[3] = duDtsx[x] * deltaEy + du_dtsy * deltaEx[x];
        // invalidGradients(duda);
        procedure.execute(tb + tI * duda[0], duda);
      }
    }
  }

  @Override
  public void forEach(Gradient2Procedure procedure) {
    final double[] duda = new double[getNumberOfGradients()];
    final double[] d2uda2 = new double[getNumberOfGradients()];
    for (int y = 0; y < maxy; y++) {
      final double du_dty = this.duDty[y];
      final double deltaEy = this.deltaEy[y];
      final double du_dtsy = this.duDtsy[y];
      final double two_du_dtsy_tI = 2 * this.duDtsy[y] / tI;
      final double d2u_dty2 = this.d2uDty2[y];
      final double d2u_dtsy2 = this.d2uDtsy2[y];
      for (int x = 0; x < maxx; x++) {
        duda[0] = deltaEx[x] * deltaEy;
        duda[1] = duDtx[x] * deltaEy;
        duda[2] = du_dty * deltaEx[x];
        duda[3] = duDtsx[x] * deltaEy + du_dtsy * deltaEx[x];
        d2uda2[1] = d2uDtx2[x] * deltaEy;
        d2uda2[2] = d2u_dty2 * deltaEx[x];
        //@formatter:off
        d2uda2[3] =
            d2uDtsx2[x] * deltaEy +
            d2u_dtsy2 * deltaEx[x] +
            duDtsx[x] * two_du_dtsy_tI;
        //@formatter:on
        procedure.execute(tb + tI * duda[0], duda, d2uda2);
      }
    }
  }

  @Override
  public void forEach(ExtendedGradient2Procedure procedure) {
    final int n = getNumberOfGradients();
    final double[] duda = new double[n];
    final double[] d2udadb = new double[n * n];
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
        duda[0] = deltaEx[x] * deltaEy;
        duda[1] = duDtx[x] * deltaEy;
        duda[2] = du_dty * deltaEx[x];
        duda[3] = duDtsx[x] * deltaEy + du_dtsy * deltaEx[x];

        // Compute all the partial second order derivatives

        // Signal,X
        d2udadb[1] = duda[1] / tI;
        // Signal,Y
        d2udadb[2] = duda[2] / tI;
        // Signal,X SD
        d2udadb[3] = duda[3] / tI;

        // X,Signal
        d2udadb[4] = d2udadb[1];
        // X,X
        d2udadb[5] = d2uDtx2[x] * deltaEy;
        // X,Y
        d2udadb[6] = duDtx[x] * du_dty_tI;
        // X,X SD
        d2udadb[7] = deltaEy * d2deltaExDtsxDx[x] + duDtx[x] * du_dtsy_tI;

        // Y,Signal
        d2udadb[8] = d2udadb[2];
        // Y,X
        d2udadb[9] = d2udadb[6];
        // Y,Y
        d2udadb[10] = d2u_dty2 * deltaEx[x];
        // Y,X SD
        d2udadb[11] = du_dty * du_dtsx_tI[x] + deltaEx[x] * d2deltaEy_dtsydy;

        // X SD,Signal
        d2udadb[12] = d2udadb[3];
        // X SD,X
        d2udadb[13] = d2udadb[7];
        // X SD,Y
        d2udadb[14] = d2udadb[11];
        // X SD,X SD
        //@formatter:off
        d2udadb[15] =
            d2uDtsx2[x] * deltaEy +
            d2u_dtsy2 * deltaEx[x] +
            duDtsx[x] * two_du_dtsy_tI;
        //@formatter:on

        procedure.executeExtended(tb + tI * duda[0], duda, d2udadb);
      }
    }
  }
}
