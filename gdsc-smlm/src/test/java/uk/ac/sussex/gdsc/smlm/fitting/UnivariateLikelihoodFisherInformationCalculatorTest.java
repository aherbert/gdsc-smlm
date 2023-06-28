/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.fitting;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.PoissonGradientProcedure;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.PoissonGradientProcedureUtils;
import uk.ac.sussex.gdsc.smlm.function.FisherInformation;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;
import uk.ac.sussex.gdsc.smlm.function.HalfPoissonFisherInformation;
import uk.ac.sussex.gdsc.smlm.function.OffsetFunctionFactory;
import uk.ac.sussex.gdsc.smlm.function.PoissonFisherInformation;
import uk.ac.sussex.gdsc.smlm.function.PoissonGaussianApproximationFisherInformation;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.test.api.Predicates;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngFactory;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;

@SuppressWarnings({"javadoc"})
class UnivariateLikelihoodFisherInformationCalculatorTest {
  enum Model {
    POISSON, HALF_POISSON, POISSON_GAUSSIAN
  }

  @SeededTest
  void canComputePoissonFisherInformation(RandomSeed seed) {
    final UniformRandomProvider r = RngFactory.create(seed.get());
    for (int n = 1; n < 10; n++) {
      computePoissonFisherInformation(r, Model.POISSON);
    }
  }

  @SeededTest
  void canComputeHalfPoissonFisherInformation(RandomSeed seed) {
    final UniformRandomProvider r = RngFactory.create(seed.get());
    for (int n = 1; n < 10; n++) {
      computePoissonFisherInformation(r, Model.HALF_POISSON);
    }
  }

  @SeededTest
  void canComputePoissonGaussianApproximationFisherInformation(RandomSeed seed) {
    final UniformRandomProvider r = RngFactory.create(seed.get());
    for (int n = 1; n < 10; n++) {
      computePoissonFisherInformation(r, Model.POISSON_GAUSSIAN);
    }
  }

  private static void computePoissonFisherInformation(UniformRandomProvider rng, Model model) {
    // Create function
    final Gaussian2DFunction func =
        GaussianFunctionFactory.create2D(1, 10, 10, GaussianFunctionFactory.FIT_ERF_CIRCLE, null);
    final double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    params[Gaussian2DFunction.BACKGROUND] = nextUniform(rng, 0.1, 0.3);
    params[Gaussian2DFunction.SIGNAL] = nextUniform(rng, 100, 300);
    params[Gaussian2DFunction.X_POSITION] = nextUniform(rng, 4, 6);
    params[Gaussian2DFunction.Y_POSITION] = nextUniform(rng, 4, 6);
    params[Gaussian2DFunction.X_SD] = nextUniform(rng, 1, 1.3);

    Gradient1Function f1 = func;
    FisherInformation fi;

    switch (model) {
      // Get a variance
      case POISSON_GAUSSIAN:
        final double var = 0.9 + 0.2 * rng.nextDouble();
        fi = new PoissonGaussianApproximationFisherInformation(Math.sqrt(var));
        f1 = (Gradient1Function) OffsetFunctionFactory.wrapFunction(func,
            SimpleArrayUtils.newDoubleArray(func.size(), var));
        break;
      case POISSON:
        fi = new PoissonFisherInformation();
        break;
      case HALF_POISSON:
        fi = new HalfPoissonFisherInformation();
        break;
      default:
        throw new IllegalStateException();
    }

    // This introduces a dependency on a different package, and relies on that
    // computing the correct answer. However that code predates this and so the
    // test ensures that the FisherInformationCalculator functions correctly.
    final PoissonGradientProcedure p1 = PoissonGradientProcedureUtils.create(f1);
    p1.computeFisherInformation(params);
    final double[] e = p1.getLinear();

    final FisherInformationCalculator calc =
        new UnivariateLikelihoodFisherInformationCalculator(func, fi);
    final FisherInformationMatrix I = calc.compute(params);
    final double[] o = I.getMatrix().data;

    final boolean emCcd = model == Model.HALF_POISSON;

    if (emCcd) {
      // Assumes half the poisson fisher information
      SimpleArrayUtils.multiply(e, 0.5);
    }

    Assertions.assertArrayEquals(e, o, 1e-6);
    final DoubleDoubleBiPredicate predicate = Predicates.doublesAreClose(5e-2, 0);

    if (model == Model.POISSON || model == Model.HALF_POISSON) {
      // Get the Mortensen approximation for fitting Poisson data with a Gaussian.
      // Set a to 100 for the square pixel adjustment.
      final double a = 100;
      final double s = params[Gaussian2DFunction.X_SD] * a;
      final double N = params[Gaussian2DFunction.SIGNAL];
      final double b2 = params[Gaussian2DFunction.BACKGROUND];
      double var = Gaussian2DPeakResultHelper.getMLVarianceX(a, s, N, b2, emCcd);

      // Convert expected variance to pixels
      var /= (a * a);

      // Get the limits by inverting the Fisher information
      final double[] crlb = I.crlb();

      TestAssertions.assertTest(var, crlb[2], predicate);
      TestAssertions.assertTest(var, crlb[3], predicate);
    }
  }

  @SeededTest
  void canComputePerPixelPoissonGaussianApproximationFisherInformation(RandomSeed seed) {
    final UniformRandomProvider r = RngFactory.create(seed.get());
    for (int n = 1; n < 10; n++) {
      canComputePerPixelPoissonGaussianApproximationFisherInformation(r);
    }
  }

  private static void
      canComputePerPixelPoissonGaussianApproximationFisherInformation(UniformRandomProvider rng) {
    // Create function
    final Gaussian2DFunction func =
        GaussianFunctionFactory.create2D(1, 10, 10, GaussianFunctionFactory.FIT_ERF_CIRCLE, null);
    final double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    params[Gaussian2DFunction.BACKGROUND] = nextUniform(rng, 0.1, 0.3);
    params[Gaussian2DFunction.SIGNAL] = nextUniform(rng, 100, 300);
    params[Gaussian2DFunction.X_POSITION] = nextUniform(rng, 4, 6);
    params[Gaussian2DFunction.Y_POSITION] = nextUniform(rng, 4, 6);
    params[Gaussian2DFunction.X_SD] = nextUniform(rng, 1, 1.3);

    Gradient1Function f1 = func;
    FisherInformation[] fi;

    // Get a per-pixel variance
    final double[] var = new double[func.size()];

    fi = new FisherInformation[var.length];
    for (int i = var.length; i-- > 0;) {
      var[i] = 0.9 + 0.2 * rng.nextDouble();
      fi[i] = new PoissonGaussianApproximationFisherInformation(Math.sqrt(var[i]));
    }

    f1 = (Gradient1Function) OffsetFunctionFactory.wrapFunction(func, var);

    // This introduces a dependency on a different package, and relies on that
    // computing the correct answer. However that code predates this and so the
    // test ensures that the FisherInformationCalculator functions correctly.
    final PoissonGradientProcedure p1 = PoissonGradientProcedureUtils.create(f1);
    p1.computeFisherInformation(params);
    final double[] e = p1.getLinear();

    final FisherInformationCalculator calc =
        new UnivariateLikelihoodFisherInformationCalculator(func, fi);
    final FisherInformationMatrix I = calc.compute(params);
    final double[] o = I.getMatrix().data;

    TestAssertions.assertArrayTest(e, o, Predicates.doublesAreClose(1e-6, 0));
  }

  private static double nextUniform(UniformRandomProvider rng, double min, double max) {
    return min + rng.nextDouble() * (max - min);
  }
}
