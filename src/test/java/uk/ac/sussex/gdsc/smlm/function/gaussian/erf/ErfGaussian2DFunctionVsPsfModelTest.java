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

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.ContinuousUniformSampler;
import org.junit.jupiter.api.Assertions;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.smlm.function.StandardValueProcedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.model.GaussianPsfModel;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;

@SuppressWarnings({"javadoc"})
public class ErfGaussian2DFunctionVsPsfModelTest {
  private final int width = 10;
  private final int height = 9;

  @SeededTest
  public void computesSameAsPsfModel(RandomSeed seed) {
    final UniformRandomProvider rng = RngUtils.create(seed.getSeed());
    for (int i = 0; i < 10; i++) {
      //@formatter:off
      computesSameAsPsfModel(
          nextUniform(rng, 50, 100),
          nextUniform(rng, (width-1)/2.0, (width+1)/2.0),
          nextUniform(rng, (height-1)/2.0, (height+1)/2.0),
          nextUniform(rng, 0.5, 2),
          nextUniform(rng, 0.5, 2));
      //@formatter:on
    }
  }

  private void computesSameAsPsfModel(double sum, double x0, double x1, double s0, double s1) {
    final Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, width, height,
        GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
    final double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    a[Gaussian2DFunction.SIGNAL] = sum;
    a[Gaussian2DFunction.X_POSITION] = x0;
    a[Gaussian2DFunction.Y_POSITION] = x1;
    a[Gaussian2DFunction.X_SD] = s0;
    a[Gaussian2DFunction.Y_SD] = s1;
    final double[] o = new StandardValueProcedure().getValues(f, a);

    final GaussianPsfModel m = new GaussianPsfModel(s0, s1);
    final double[] e = new double[o.length];
    // Note that the Gaussian2DFunction has 0,0 at the centre of a pixel.
    // The model has 0.5,0.5 at the centre so add an offset.
    m.create2D(e, width, height, sum, x0 + 0.5, x1 + 0.5, null);

    // Since the model only computes within +/- 5 sd only check for equality
    // when the model is not zero (and there is a reasonable amount of signal)

    for (int i = 0; i < e.length; i++) {
      // Only check where there is a reasonable amount of signal
      if (e[i] > 1e-2) {
        final double error = DoubleEquality.relativeError(e[i], o[i]);
        // We expect a small error since the ErfGaussian2DFunction uses a
        // fast approximation of the Erf(..) (the error function). The PSFModel
        // uses the Apache commons implementation.
        if (error > 5e-4) {
          Assertions.fail(String.format("[%d] %s != %s  error = %f", i, Double.toString(e[i]),
              Double.toString(o[i]), error));
        }
      }
    }
  }

  private static double nextUniform(UniformRandomProvider rng, double min, double max) {
    return new ContinuousUniformSampler(rng, min, max).sample();
  }
}
