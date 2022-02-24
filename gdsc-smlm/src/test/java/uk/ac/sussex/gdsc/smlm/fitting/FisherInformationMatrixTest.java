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

package uk.ac.sussex.gdsc.smlm.fitting;

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.rng.UniformRandomProvider;
import org.ejml.data.DenseMatrix64F;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.opentest4j.AssertionFailedError;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.rng.RandomUtils;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.GradientCalculatorUtils;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils.TestLevel;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;

@SuppressWarnings({"javadoc"})
class FisherInformationMatrixTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(FisherInformationMatrixTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  /** The level for logging output. */
  private final Level level = TestLevel.TEST_DEBUG;

  @SeededTest
  void canComputeCrlb(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.get());
    for (int n = 1; n < 10; n++) {
      testComputeCrlb(rg, n, 0, true);
    }
  }

  @SeededTest
  void canComputeCrlbWithZeros(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.get());
    for (int n = 2; n < 10; n++) {
      testComputeCrlb(rg, n, 1, true);
      testComputeCrlb(rg, n, n / 2, true);
    }
  }

  @SeededTest
  void canComputeCrlbWithReciprocal(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.get());
    for (int n = 1; n < 10; n++) {
      testComputeCrlb(rg, n, 0, false);
    }
  }

  @SeededTest
  void canComputeCrlbWithReciprocalWithZeros(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.get());
    for (int n = 2; n < 10; n++) {
      testComputeCrlb(rg, n, 1, false);
      testComputeCrlb(rg, n, n / 2, false);
    }
  }

  @SeededTest
  void inversionDoesNotMatchReciprocal(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.get());
    for (int n = 1; n < 10; n++) {
      final FisherInformationMatrix m = createFisherInformationMatrix(rg, n, 0);
      final double[] crlb = m.crlb();
      final double[] crlb2 = m.crlbReciprocal();
      // These increasingly do not match with increasing number of parameters.
      if (logger.isLoggable(level)) {
        logger.log(level,
            FunctionUtils.getSupplier("%s =? %s", Arrays.toString(crlb), Arrays.toString(crlb2)));
      }
      if (n > 1) {
        // Just do a sum so we have a test
        Assertions.assertThrows(AssertionFailedError.class, () -> {
          Assertions.assertEquals(MathUtils.sum(crlb), MathUtils.sum(crlb2));
        });
      }
    }
  }

  private double[] testComputeCrlb(UniformRandomProvider rg, int columns, int zeroColumns,
      boolean invert) {
    final FisherInformationMatrix m = createFisherInformationMatrix(rg, columns, zeroColumns);

    // Invert for Crlb
    final double[] crlb = (invert) ? m.crlb() : m.crlbReciprocal();
    if (logger.isLoggable(level)) {
      logger.log(level, FunctionUtils.getSupplier("columns=%d, zeroColumns=%d : %s", columns,
          zeroColumns, Arrays.toString(crlb)));
    }
    Assertions.assertNotNull(crlb,
        () -> String.format("Crlb failed: columns=%d, zeroColumns=%d", columns, zeroColumns));
    return crlb;
  }

  private static FisherInformationMatrix createFisherInformationMatrix(UniformRandomProvider rg,
      int columns, int zeroColumns) {
    final int maxx = 10;
    final int size = maxx * maxx;

    // Use a real Gaussian function here to compute the Fisher information.
    // The matrix may be sensitive to the type of equation used.
    int npeazeroColumnss = 1;
    Gaussian2DFunction fun = createFunction(maxx, npeazeroColumnss);
    while (fun.getNumberOfGradients() < columns) {
      npeazeroColumnss++;
      fun = createFunction(maxx, npeazeroColumnss);
    }

    final double[] a = new double[1 + npeazeroColumnss * Gaussian2DFunction.PARAMETERS_PER_PEAK];
    a[Gaussian2DFunction.BACKGROUND] = nextUniform(rg, 1, 5);
    for (int i = 0, j = 0; i < npeazeroColumnss; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK) {
      a[j + Gaussian2DFunction.SIGNAL] = nextUniform(rg, 100, 300);
      // Non-overlapping peazeroColumnss otherwise the Crlb are poor
      a[j + Gaussian2DFunction.X_POSITION] = nextUniform(rg, 2 + i * 2, 4 + i * 2);
      a[j + Gaussian2DFunction.Y_POSITION] = nextUniform(rg, 2 + i * 2, 4 + i * 2);
      a[j + Gaussian2DFunction.X_SD] = nextUniform(rg, 1.5, 2);
      a[j + Gaussian2DFunction.Y_SD] = nextUniform(rg, 1.5, 2);
    }
    fun.initialise(a);

    final GradientCalculator calc =
        GradientCalculatorUtils.newCalculator(fun.getNumberOfGradients());
    double[][] matrixI = calc.fisherInformationMatrix(size, a, fun);

    // Reduce to the desired size
    matrixI = Arrays.copyOf(matrixI, columns);
    for (int i = 0; i < columns; i++) {
      matrixI[i] = Arrays.copyOf(matrixI[i], columns);
    }

    // Zero selected columns
    if (zeroColumns > 0) {
      final int[] zero = RandomUtils.sample(zeroColumns, columns, rg);
      for (final int i : zero) {
        for (int j = 0; j < columns; j++) {
          matrixI[i][j] = matrixI[j][i] = 0;
        }
      }
    }

    // TestLog.fine(logger,"columns=%d, zeroColumns=%d", columns, zeroColumns);
    // for (int i = 0; i < columns; i++)
    // TestLog.fine(logger,Arrays.toString(I[i]));

    // Create matrix
    return new FisherInformationMatrix(matrixI, 1e-3);
  }

  private static double nextUniform(UniformRandomProvider rng, double min, double max) {
    return min + rng.nextDouble() * (max - min);
  }

  private static Gaussian2DFunction createFunction(int maxx, int npeaks) {
    final Gaussian2DFunction f = GaussianFunctionFactory.create2D(npeaks, maxx, maxx,
        GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
    return f;
  }

  @SeededTest
  void canProduceSubset(RandomSeed seed) {
    final int k = 5;
    final int n = 10;

    final UniformRandomProvider UniformRandomProvider = RngUtils.create(seed.get());
    final FisherInformationMatrix m = createRandomMatrix(UniformRandomProvider, n);
    final DenseMatrix64F e = m.getMatrix();
    if (logger.isLoggable(level)) {
      logger.log(level, String.valueOf(e));
    }

    for (int run = 1; run < 10; run++) {
      final int[] indices = RandomUtils.sample(k, n, UniformRandomProvider);
      Arrays.sort(indices);
      final DenseMatrix64F o = m.subset(indices).getMatrix();
      if (logger.isLoggable(level)) {
        logger.log(level, FunctionUtils.getSupplier(Arrays.toString(indices)));
        logger.log(level, String.valueOf(o));
      }
      for (int i = 0; i < indices.length; i++) {
        for (int j = 0; j < indices.length; j++) {
          Assertions.assertEquals(e.get(indices[i], indices[j]), o.get(i, j));
        }
      }
    }
  }

  private static FisherInformationMatrix createRandomMatrix(UniformRandomProvider rng, int size) {
    final double[] data = new double[size * size];
    for (int i = 0; i < data.length; i++) {
      data[i] = rng.nextDouble();
    }
    return new FisherInformationMatrix(data, size);
  }

  @SeededTest
  void computeWithSubsetReducesTheCrlb(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.get());
    final Gaussian2DFunction f = createFunction(10, 1);
    final int perPeak = f.getGradientParametersPerPeak();
    // Create a matrix with 2 peaks + background
    final FisherInformationMatrix m = createFisherInformationMatrix(rg, 1 + 2 * perPeak, 0);
    // Subset each peak
    final int[] indices = SimpleArrayUtils.natural(1 + perPeak);
    final FisherInformationMatrix m1 = m.subset(indices);
    for (int i = 1; i < indices.length; i++) {
      indices[i] += perPeak;
    }
    final FisherInformationMatrix m2 = m.subset(indices);

    // TestLog.fine(logger,m.getMatrix());
    // TestLog.fine(logger,m1.getMatrix());
    // TestLog.fine(logger,m2.getMatrix());

    final double[] crlb = m.crlb();
    final double[] crlb1 = m1.crlb();
    final double[] crlb2 = m2.crlb();
    final double[] crlbB = Arrays.copyOf(crlb1, crlb.length);
    System.arraycopy(crlb2, 1, crlbB, crlb1.length, perPeak);

    // TestLog.fine(logger,Arrays.toString(crlb));
    // TestLog.fine(logger,Arrays.toString(crlb1));
    // TestLog.fine(logger,Arrays.toString(crlb2));
    // TestLog.fine(logger,Arrays.toString(crlbB));

    // Removing the interaction between fit parameters lowers the bounds
    for (int i = 0; i < crlb.length; i++) {
      Assertions.assertTrue(crlbB[i] < crlb[i]);
    }
  }
}
