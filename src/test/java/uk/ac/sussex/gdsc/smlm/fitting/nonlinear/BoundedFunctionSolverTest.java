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

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear;

import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.stop.ErrorStoppingCriteria;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestSettings;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;

/**
 * Test that a bounded fitter can return the same results with and without bounds.
 */
@SuppressWarnings({"javadoc"})
public class BoundedFunctionSolverTest extends BaseFunctionSolverTest {
  // The following tests ensure that the Lvm can fit data without
  // requiring a bias (i.e. an offset to the background).
  // In a previous version the Lvm fitter was stable only if a bias existed.
  // The exact source of this instability is unknown as it could be due to
  // how the data was processed before or after fitting, or within the Lvm
  // fitter itself. However the process should be the same without a bias
  // and these tests ensure that is true.

  @SeededTest
  public void fitSingleGaussianLvmWithoutBias(RandomSeed seed) {
    fitSingleGaussianWithoutBias(seed, false, 0);
  }

  @SeededTest
  public void fitSingleGaussianCLvmWithoutBias(RandomSeed seed) {
    fitSingleGaussianWithoutBias(seed, false, 1);
  }

  @SeededTest
  public void fitSingleGaussianDcLvmWithoutBias(RandomSeed seed) {
    fitSingleGaussianWithoutBias(seed, false, 2);
  }

  @SeededTest
  public void fitSingleGaussianBLvmWithoutBias(RandomSeed seed) {
    fitSingleGaussianWithoutBias(seed, true, 0);
  }

  @SeededTest
  public void fitSingleGaussianBcLvmWithoutBias(RandomSeed seed) {
    fitSingleGaussianWithoutBias(seed, true, 1);
  }

  @SeededTest
  public void fitSingleGaussianBdcLvmWithoutBias(RandomSeed seed) {
    fitSingleGaussianWithoutBias(seed, true, 2);
  }

  private void fitSingleGaussianWithoutBias(RandomSeed seed, boolean applyBounds, int clamping) {
    final double bias = 100;

    final NonLinearFit solver = getLvm((applyBounds) ? 2 : 1, clamping, false);
    final NonLinearFit solver2 = getLvm((applyBounds) ? 2 : 1, clamping, false);

    final String name = getLvmName(applyBounds, clamping, false);

    final int loops = 5;
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final StoredDataStatistics[] stats = new StoredDataStatistics[6];

    for (final double s : signal) {
      final double[] expected = createParams(1, s, 0, 0, 1);
      double[] lower = null;
      double[] upper = null;
      if (applyBounds) {
        lower = createParams(0, s * 0.5, -0.2, -0.2, 0.8);
        upper = createParams(3, s * 2, 0.2, 0.2, 1.2);
        solver.setBounds(lower, upper);
      }

      final double[] expected2 = addBiasToParams(expected, bias);
      if (applyBounds) {
        final double[] lower2 = addBiasToParams(lower, bias);
        final double[] upper2 = addBiasToParams(upper, bias);
        solver2.setBounds(lower2, upper2);
      }

      for (int loop = loops; loop-- > 0;) {
        final double[] data = drawGaussian(expected, rg);
        final double[] data2 = data.clone();
        for (int i = 0; i < data.length; i++) {
          data2[i] += bias;
        }

        for (int i = 0; i < stats.length; i++) {
          stats[i] = new StoredDataStatistics();
        }

        for (final double db : base) {
          for (final double dx : shift) {
            for (final double dy : shift) {
              for (final double dsx : factor) {
                final double[] p = createParams(db, s, dx, dy, dsx);
                final double[] p2 = addBiasToParams(p, bias);

                final double[] fp = fitGaussian(solver, data, p, expected);
                final double[] fp2 = fitGaussian(solver2, data2, p2, expected2);

                // The result should be the same without a bias
                Assertions.assertEquals(solver.getEvaluations(), solver2.getEvaluations(),
                    () -> name + " Iterations");
                fp2[0] -= bias;
                Assertions.assertArrayEquals(fp, fp2, 1e-6, () -> name + " Solution");
              }
            }
          }
        }
      }
    }
  }

  // Standard Lvm
  @SeededTest
  public void canFitSingleGaussianLvm(RandomSeed seed) {
    fitSingleGaussianLvm(seed, 0, 0, false);
  }

  // Bounded/Clamped Lvm

  @SeededTest
  public void canFitSingleGaussianBLvmNoBounds(RandomSeed seed) {
    fitSingleGaussianLvm(seed, 1, 0, false);
  }

  @SeededTest
  public void canFitSingleGaussianBLvm(RandomSeed seed) {
    fitSingleGaussianLvm(seed, 2, 0, false);
  }

  @SeededTest
  public void canFitSingleGaussianCLvm(RandomSeed seed) {
    fitSingleGaussianLvm(seed, 0, 1, false);
  }

  @SeededTest
  public void canFitSingleGaussianDcLvm(RandomSeed seed) {
    fitSingleGaussianLvm(seed, 0, 2, false);
  }

  @SeededTest
  public void canFitSingleGaussianBcLvm(RandomSeed seed) {
    fitSingleGaussianLvm(seed, 2, 1, false);
  }

  @SeededTest
  public void canFitSingleGaussianBdcLvm(RandomSeed seed) {
    fitSingleGaussianLvm(seed, 2, 2, false);
  }

  // Mle Lvm

  @SeededTest
  public void canFitSingleGaussianLvmMle(RandomSeed seed) {
    fitSingleGaussianLvm(seed, 0, 0, true);
  }

  @SeededTest
  public void canFitSingleGaussianBLvmMleNoBounds(RandomSeed seed) {
    fitSingleGaussianLvm(seed, 1, 0, true);
  }

  @SeededTest
  public void canFitSingleGaussianBLvmMle(RandomSeed seed) {
    fitSingleGaussianLvm(seed, 2, 0, true);
  }

  private void fitSingleGaussianLvm(RandomSeed seed, int bounded, int clamping, boolean mle) {
    canFitSingleGaussian(seed, getLvm(bounded, clamping, mle), bounded == 2);
  }

  // Is Bounded/Clamped Lvm better?

  @SeededTest
  public void fitSingleGaussianBLvmBetterThanLvm(RandomSeed seed) {
    fitSingleGaussianBetterLvm(seed, true, 0, false, false, 0, false);
  }

  @SeededTest
  public void fitSingleGaussianCLvmBetterThanLvm(RandomSeed seed) {
    fitSingleGaussianBetterLvm(seed, false, 1, false, false, 0, false);
  }

  @SeededTest
  public void fitSingleGaussianBcLvmBetterThanLvm(RandomSeed seed) {
    fitSingleGaussianBetterLvm(seed, true, 1, false, false, 0, false);
  }

  @SeededTest
  public void fitSingleGaussianDcLvmBetterThanLvm(RandomSeed seed) {
    fitSingleGaussianBetterLvm(seed, false, 2, false, false, 0, false);
  }

  @SeededTest
  public void fitSingleGaussianBdcLvmBetterThanLvm(RandomSeed seed) {
    fitSingleGaussianBetterLvm(seed, true, 2, false, false, 0, false);
  }

  @SeededTest
  public void fitSingleGaussianLvmMleBetterThanLvm(RandomSeed seed) {
    fitSingleGaussianBetterLvm(seed, false, 0, true, false, 0, false);
  }

  @SeededTest
  public void fitSingleGaussianBLvmMleBetterThanLvm(RandomSeed seed) {
    fitSingleGaussianBetterLvm(seed, true, 0, true, false, 0, false);
  }

  @SeededTest
  public void fitSingleGaussianCLvmMleBetterThanLvm(RandomSeed seed) {
    fitSingleGaussianBetterLvm(seed, false, 1, true, false, 0, false);
  }

  @SeededTest
  public void fitSingleGaussianBcLvmMleBetterThanLvm(RandomSeed seed) {
    fitSingleGaussianBetterLvm(seed, true, 1, true, false, 0, false);
  }

  @SeededTest
  public void fitSingleGaussianDcLvmMleBetterThanLvm(RandomSeed seed) {
    fitSingleGaussianBetterLvm(seed, false, 2, true, false, 0, false);
  }

  @SeededTest
  public void fitSingleGaussianBdcLvmMleBetterThanLvm(RandomSeed seed) {
    fitSingleGaussianBetterLvm(seed, true, 2, true, false, 0, false);
  }

  @SeededTest
  public void fitSingleGaussianBLvmMleBetterThanLvmMle(RandomSeed seed) {
    fitSingleGaussianBetterLvm(seed, true, 0, true, false, 0, true);
  }

  @SeededTest
  public void fitSingleGaussianCLvmMleBetterThanLvmMle(RandomSeed seed) {
    fitSingleGaussianBetterLvm(seed, false, 1, true, false, 0, true);
  }

  @SeededTest
  public void fitSingleGaussianDcLvmMleBetterThanLvmMle(RandomSeed seed) {
    fitSingleGaussianBetterLvm(seed, false, 2, true, false, 0, true);
  }

  @SeededTest
  public void fitSingleGaussianBdcLvmMleBetterThanLvmMle(RandomSeed seed) {
    fitSingleGaussianBetterLvm(seed, true, 2, true, false, 0, true);
  }

  @SeededTest
  public void fitSingleGaussianBLvmMleBetterThanBLvm(RandomSeed seed) {
    fitSingleGaussianBetterLvm(seed, true, 0, true, true, 0, false);
  }

  @SeededTest
  public void fitSingleGaussianBcLvmMleBetterThanBcLvm(RandomSeed seed) {
    fitSingleGaussianBetterLvm(seed, true, 1, true, true, 1, false);
  }

  @SeededTest
  public void fitSingleGaussianBdcLvmMleBetterThanBdcLvm(RandomSeed seed) {
    fitSingleGaussianBetterLvm(seed, true, 2, true, true, 2, false);
  }

  private void fitSingleGaussianBetterLvm(RandomSeed seed, boolean bounded2, int clamping2,
      boolean mle2, boolean bounded, int clamping, boolean mle) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    final NonLinearFit solver = getLvm((bounded) ? 2 : 1, clamping, mle);
    final NonLinearFit solver2 = getLvm((bounded2) ? 2 : 1, clamping2, mle2);
    canFitSingleGaussianBetter(seed, solver, bounded, solver2, bounded2,
        getLvmName(bounded, clamping, mle), getLvmName(bounded2, clamping2, mle2));
  }

  private static NonLinearFit getLvm(int bounded, int clamping, boolean mle) {
    final Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, size, flags, null);
    final StoppingCriteria sc = new ErrorStoppingCriteria(5);
    sc.setMaximumIterations(100);
    final NonLinearFit solver =
        (bounded != 0 || clamping != 0) ? new BoundedNonLinearFit(f, sc, null)
            : new NonLinearFit(f, sc);
    if (clamping != 0) {
      final BoundedNonLinearFit bsolver = (BoundedNonLinearFit) solver;
      final ParameterBounds bounds = new ParameterBounds(f);
      bounds.setClampValues(defaultClampValues);
      bounds.setDynamicClamp(clamping == 2);
      bsolver.setBounds(bounds);
    }
    solver.setMle(mle);
    solver.setInitialLambda(1);
    return solver;
  }

  private static String getLvmName(boolean bounded, int clamping, boolean mle) {
    return ((bounded) ? "B" : "") + ((clamping == 0) ? "" : ((clamping == 1) ? "C" : "DC")) + "LVM"
        + ((mle) ? " MLE" : "");
  }

  private class CheatingStoppingCriteria extends StoppingCriteria {
    @Override
    public void evaluate(double oldError, double newError, double[] a) {
      // Do nothing
    }

    @Override
    public boolean areAchieved() {
      return true;
    }

    @Override
    public boolean areNotSatisfied() {
      return false;
    }
  }

  @SeededTest
  public void canFitAndComputeDeviationsLvm(RandomSeed seed) {
    canFitAndComputeDeviations(seed, false);
  }

  @SeededTest
  public void canFitAndComputeDeviationsLvmMle(RandomSeed seed) {
    canFitAndComputeDeviations(seed, true);
  }

  private void canFitAndComputeDeviations(RandomSeed seed, boolean mle) {
    final NonLinearFit solver1 = getLvm(0, 0, mle);
    final NonLinearFit solver2 = getLvm(0, 0, mle);
    solver1.setStoppingCriteria(new CheatingStoppingCriteria());
    fitAndComputeDeviationsMatch(seed, solver1, solver2, NoiseModel.EMCCD, false);
  }
}
