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

package uk.ac.sussex.gdsc.smlm.math3.analysis.integration;

import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;

@SuppressWarnings({"javadoc"})
public class CustomSimpsonIntegratorTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(CustomSimpsonIntegratorTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  private interface TestUnivariateFunction extends UnivariateFunction {
    /**
     * Sum the function between the lower and upper limit.
     *
     * @param ax the lower limit on x
     * @param bx the upper limit on x
     * @return the sum
     */
    public double sum(double ax, double bx);
  }

  private class LinearTestUnivariateFunction implements TestUnivariateFunction {
    @Override
    public double value(double x) {
      return x;
    }

    @Override
    public double sum(double ax, double bx) {
      // y=x => x^2/2
      return (bx * bx - ax * ax) / 2;
    }
  }

  private class QuadraticTestUnivariateFunction implements TestUnivariateFunction {
    @Override
    public double value(double x) {
      return x * x - x;
    }

    @Override
    public double sum(double ax, double bx) {
      // y=x^2 - x => x^3/3 - x^2 / 2
      final double u = MathUtils.pow3(bx) / 3 - MathUtils.pow2(bx) / 2;
      final double l = MathUtils.pow3(ax) / 3 - MathUtils.pow2(ax) / 2;
      return u - l;
    }
  }

  private class CubicTestUnivariateFunction implements TestUnivariateFunction {
    @Override
    public double value(double x) {
      return x * x * x - x * x;
    }

    @Override
    public double sum(double ax, double bx) {
      // y=x^3 - x^2 => x^4/4 - x^3 / 3
      final double u = MathUtils.pow4(bx) / 4 - MathUtils.pow3(bx) / 3;
      final double l = MathUtils.pow4(ax) / 4 - MathUtils.pow3(ax) / 3;
      return u - l;
    }
  }

  private class ReciprocalTestUnivariateFunction implements TestUnivariateFunction {
    @Override
    public double value(double x) {
      return 1.0 / x;
    }

    @Override
    public double sum(double ax, double bx) {
      // y=1/x => log(x)
      return Math.log(bx) - Math.log(ax);
    }
  }

  @Test
  public void canIntegrateFunction() {
    TestUnivariateFunction func;

    func = new LinearTestUnivariateFunction();
    canIntegrateFunction(func, 0.5, 2, 1);
    canIntegrateFunction(func, 0.5, 2, 2);
    canIntegrateFunction(func, 0.5, 2, 3);

    func = new QuadraticTestUnivariateFunction();
    canIntegrateFunction(func, 0.5, 2, 1);
    canIntegrateFunction(func, 0.5, 2, 2);
    canIntegrateFunction(func, 0.5, 2, 3);

    func = new CubicTestUnivariateFunction();
    canIntegrateFunction(func, 0.5, 2, 1);
    canIntegrateFunction(func, 0.5, 2, 2);
    canIntegrateFunction(func, 0.5, 2, 3);

    // Harder so use more iterations
    func = new ReciprocalTestUnivariateFunction();
    canIntegrateFunction(func, 0.5, 2, 5);
  }

  private static void canIntegrateFunction(TestUnivariateFunction func, double ax, double bx,
      int iter) {
    final double relativeAccuracy = 1e-4;
    // So it stops at the min iterations
    final double absoluteAccuracy = Double.POSITIVE_INFINITY;

    final CustomSimpsonIntegrator in = new CustomSimpsonIntegrator(
        // new SimpsonIntegrator(
        relativeAccuracy, absoluteAccuracy, iter,
        CustomSimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);

    final double e = func.sum(ax, bx);
    final double ee = simpson(func, ax, bx, iter);
    final double o = in.integrate(Integer.MAX_VALUE, func, ax, bx);

    logger.log(TestLogUtils.getRecord(Level.INFO, "%s iter=%d  %g-%g  e=%g  ee=%g  o=%g",
        func.getClass().getSimpleName(), iter, ax, bx, e, ee, o));

    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-6, 0);
    TestAssertions.assertTest(e, ee, predicate);
    TestAssertions.assertTest(e, o, predicate);

    // These should be the same within numeric tolerance
    TestAssertions.assertTest(ee, o, TestHelper.doublesAreClose(1e-12, 0));
  }

  private static double simpson(UnivariateFunction func, double ax, double bx, int iter) {
    // Simple Simpson integration:
    // https://en.wikipedia.org/wiki/Simpson%27s_rule

    // Number of sub intervals
    final int n = 1 << iter + 1;
    final double h = (bx - ax) / n; // sub-interval width
    double sum2 = 0;
    double sum4 = 0;
    for (int j = 1; j <= n / 2 - 1; j++) {
      sum2 += func.value(ax + (2 * j) * h);
    }
    for (int j = 1; j <= n / 2; j++) {
      sum4 += func.value(ax + (2 * j - 1) * h);
    }
    final double sum = (h / 3) * (func.value(ax) + 2 * sum2 + 4 * sum4 + func.value(bx));
    return sum;
  }
}
