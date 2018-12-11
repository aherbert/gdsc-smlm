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

import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.function.gaussian.HoltzerAstigmatismZModel;

@SuppressWarnings({"javadoc"})
public class MultiAstigmatismErfGaussian2DFunctionTest extends ErfGaussian2DFunctionTest {
  @Override
  protected void init() {
    flags = GaussianFunctionFactory.FIT_ERF_ASTIGMATISM;
    // Use a reasonable z-depth function from the Smith, et al (2010) paper (page 377)
    final double sx = 1.08;
    final double sy = 1.01;
    final double gamma = 0.389;
    final double d = 0.531;
    final double Ax = -0.0708;
    final double Bx = -0.073;
    final double Ay = 0.164;
    final double By = 0.0417;
    zModel = HoltzerAstigmatismZModel.create(sx, sy, gamma, d, Ax, Bx, Ay, By);
    f1 = new MultiAstigmatismErfGaussian2DFunction(1, maxx, maxy, zModel);
    f2 = new MultiAstigmatismErfGaussian2DFunction(2, maxx, maxy, zModel);
  }

  @Override
  protected void postInit() {
    // Even though the function does not evaluate the widths it can use them
    // to construct independent widths.
    // Test with different X and Y SD
    testw1 = new double[][] {{1.1, 1.1}, {1.1, 1.2}, {1.1, 1.4}};
    testw2 = new double[][] {{1.2, 1.2}, {1.2, 1.3}, {1.2, 1.5}};
  }
}
