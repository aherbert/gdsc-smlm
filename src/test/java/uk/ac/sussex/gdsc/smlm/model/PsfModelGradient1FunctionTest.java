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

package uk.ac.sussex.gdsc.smlm.model;

import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Procedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.AstigmatismZModel;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.HoltzerAstigmatismZModel;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction.ErfFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.SingleAstigmatismErfGaussian2DFunction;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;

@SuppressWarnings({"javadoc"})
class PsfModelGradient1FunctionTest {
  @Test
  void canComputeValueAndGradient() {
    // Use a reasonable z-depth function from the Smith, et al (2010) paper (page 377)
    final double sx = 1.08;
    final double sy = 1.01;
    final double gamma = 0.389;
    final double d = 0.531;
    final double Ax = -0.0708;
    final double Bx = -0.073;
    final double Ay = 0.164;
    final double By = 0.0417;
    final AstigmatismZModel zModel =
        HoltzerAstigmatismZModel.create(sx, sy, gamma, d, Ax, Bx, Ay, By);

    // Small size ensure the PSF model covers the entire image
    final int maxx = 11;
    final int maxy = 11;
    final double[] ve = new double[maxx * maxy];
    final double[] vo = new double[maxx * maxy];
    final double[][] ge = new double[maxx * maxy][];
    final double[][] go = new double[maxx * maxy][];

    final PsfModelGradient1Function psf =
        new PsfModelGradient1Function(new GaussianPsfModel(zModel), maxx, maxy);
    final ErfGaussian2DFunction f = new SingleAstigmatismErfGaussian2DFunction(maxx, maxy, zModel);
    f.setErfFunction(ErfFunction.COMMONS_MATH);
    final double[] a2 = new double[Gaussian2DFunction.PARAMETERS_PER_PEAK + 1];
    final DoubleDoubleBiPredicate equality = TestHelper.doublesAreClose(1e-8, 0);

    final double c = maxx * 0.5;
    for (int i = -1; i <= 1; i++) {
      final double x0 = c + i * 0.33;
      for (int j = -1; j <= 1; j++) {
        final double x1 = c + j * 0.33;
        for (int k = -1; k <= 1; k++) {
          final double x2 = k * 0.33;
          for (final double in : new double[] {23.2, 405.67}) {
            // Background is constant for gradients so just use 1 value
            final double[] a = new double[] {2.2, in, x0, x1, x2};
            psf.initialise1(a);
            psf.forEach(new Gradient1Procedure() {
              int index = 0;

              @Override
              public void execute(double value, double[] dyDa) {
                vo[index] = value;
                go[index] = dyDa.clone();
                index++;
              }
            });
            a2[Gaussian2DFunction.BACKGROUND] = a[0];
            a2[Gaussian2DFunction.SIGNAL] = a[1];
            a2[Gaussian2DFunction.X_POSITION] = a[2] - 0.5;
            a2[Gaussian2DFunction.Y_POSITION] = a[3] - 0.5;
            a2[Gaussian2DFunction.Z_POSITION] = a[4];
            f.initialise1(a2);
            f.forEach(new Gradient1Procedure() {
              int index = 0;

              @Override
              public void execute(double value, double[] dyDa) {
                ve[index] = value;
                ge[index] = dyDa.clone();
                index++;
              }
            });

            for (int ii = 0; ii < ve.length; ii++) {
              TestAssertions.assertTest(ve[ii], vo[ii], equality);
              TestAssertions.assertArrayTest(ge[ii], go[ii], equality);
            }
          }
        }
      }
    }
  }
}
