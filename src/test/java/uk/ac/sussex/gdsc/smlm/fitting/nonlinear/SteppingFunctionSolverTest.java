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

import org.junit.jupiter.api.Assumptions;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestSettings;

/**
 * Test that a stepping solver can fit a function.
 */
@SuppressWarnings({"javadoc"})
public class SteppingFunctionSolverTest extends BaseSteppingFunctionSolverTest {
  @SeededTest
  public void canFitSingleGaussianEmCcd_x_x__LseLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_x_C__LseLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, CLAMP, LSELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_x_DC_LseLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, LSELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_B_x__LseLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_B_C__LseLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, CLAMP, LSELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_B_DC_LseLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, LSELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_x_x__WlseLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, WLSELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_x_C__WlseLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, CLAMP, WLSELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_x_DC_WlseLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, WLSELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_B_x__WlseLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, NO_CLAMP, WLSELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_B_C__WlseLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, CLAMP, WLSELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_B_DC_WlseLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, WLSELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_x_x__MleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, MLELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_x_C__MleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, CLAMP, MLELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_x_DC_MleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, MLELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_B_x__MleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, NO_CLAMP, MLELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_B_C__MleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, CLAMP, MLELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_B_DC_MleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, MLELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_x_x__FastLogMleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_x_C__FastLogMleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_x_DC_FastLogMleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_B_x__FastLogMleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, NO_CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_B_C__FastLogMleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_B_DC_FastLogMleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_x_x__FastMle(RandomSeed seed) {
    // The FastMle method can generate very big steps that make the method unstable.
    // This test may fail depending on the random number generator.
    try {
      fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, FastMLE, NoiseModel.EMCCD);
    } catch (final AssertionError ex) {
      logger.log(TestLogUtils.getFailRecord(ex));
    }
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_x_C__FastMle(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, CLAMP, FastMLE, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_x_DC_FastMle(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, FastMLE, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_B_x__FastMle(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, NO_CLAMP, FastMLE, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_B_C__FastMle(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, CLAMP, FastMLE, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_B_DC_FastMle(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, FastMLE, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_x_x__BtFastMle(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, BtFastMLE, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_x_C__BtFastMle(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, CLAMP, BtFastMLE, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_x_DC_BtFastMle(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, BtFastMLE, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_B_x__BtFastMle(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, NO_CLAMP, BtFastMLE, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_B_C__BtFastMle(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, CLAMP, BtFastMLE, NoiseModel.EMCCD);
  }

  @SeededTest
  public void canFitSingleGaussianEmCcd_B_DC_BtFastMle(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, BtFastMLE, NoiseModel.EMCCD);
  }

  // Weighted solvers for sCMOS

  @SeededTest
  public void canFitSingleGaussianScmos_x_x__WlseLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, WLSELVM, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_x_C__WlseLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, CLAMP, WLSELVM, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_x_DC_WlseLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, WLSELVM, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_B_x__WlseLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, NO_CLAMP, WLSELVM, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_B_C__WlseLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, CLAMP, WLSELVM, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_B_DC_WlseLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, WLSELVM, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_x_x__MleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, MLELVM, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_x_C__MleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, CLAMP, MLELVM, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_x_DC_MleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, MLELVM, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_B_x__MleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, NO_CLAMP, MLELVM, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_B_C__MleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, CLAMP, MLELVM, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_B_DC_MleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, MLELVM, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_x_x__FastLogMleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, FastLogMLELVM, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_x_C__FastLogMleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, CLAMP, FastLogMLELVM, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_x_DC_FastLogMleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, FastLogMLELVM, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_B_x__FastLogMleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, NO_CLAMP, FastLogMLELVM, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_B_C__FastLogMleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, CLAMP, FastLogMLELVM, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_B_DC_FastLogMleLvm(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, FastLogMLELVM, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_x_x__FastMle(RandomSeed seed) {
    // The FastMle method can generate very big steps that make the method unstable
    // This test may fail depending on the random number generator.
    try {
      fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, FastMLE, NoiseModel.SCMOS);
    } catch (final AssertionError ex) {
      logger.log(TestLogUtils.getFailRecord(ex));
    }
  }

  @SeededTest
  public void canFitSingleGaussianScmos_x_C__FastMle(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, CLAMP, FastMLE, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_x_DC_FastMle(RandomSeed seed) {
    fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, FastMLE, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_B_x__FastMle(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, NO_CLAMP, FastMLE, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_B_C__FastMle(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, CLAMP, FastMLE, NoiseModel.SCMOS);
  }

  @SeededTest
  public void canFitSingleGaussianScmos_B_DC_FastMle(RandomSeed seed) {
    fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, FastMLE, NoiseModel.SCMOS);
  }

  private void fitSingleGaussian(RandomSeed seed, boolean bounded,
      SteppingFunctionSolverClamp clamp, SteppingFunctionSolverType type, NoiseModel noiseModel) {
    // org.junit.Assumptions.assumeTrue(false);
    final SteppingFunctionSolver solver = getSolver(clamp, type);
    canFitSingleGaussian(seed, solver, bounded, noiseModel);
  }

  // Is Bounded/Clamped better?

  @SeededTest
  public void fitSingleGaussianEmCcd_B_x__LseLvmBetterThanLseLvm(RandomSeed seed) {
    fitSingleGaussianBetter(seed, BOUNDED, NO_CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM,
        NoiseModel.EMCCD);
  }

  @SeededTest
  public void fitSingleGaussianEmCcd_x_C__LseLvmBetterThanLseLvm(RandomSeed seed) {
    fitSingleGaussianBetter(seed, NO_BOUND, CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM,
        NoiseModel.EMCCD);
  }

  @SeededTest
  public void fitSingleGaussianEmCcd_B_C__LseLvmBetterThanLseLvm(RandomSeed seed) {
    fitSingleGaussianBetter(seed, BOUNDED, CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM,
        NoiseModel.EMCCD);
  }

  @SeededTest
  public void fitSingleGaussianEmCcd_x_DC_LseLvmBetterThanLseLvm(RandomSeed seed) {
    fitSingleGaussianBetter(seed, NO_BOUND, DYNAMIC_CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM,
        NoiseModel.EMCCD);
  }

  @SeededTest
  public void fitSingleGaussianEmCcd_B_DC_LseLvmBetterThanLseLvm(RandomSeed seed) {
    fitSingleGaussianBetter(seed, BOUNDED, DYNAMIC_CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM,
        NoiseModel.EMCCD);
  }

  @SeededTest
  public void fitSingleGaussianEmCcd_x_x__MleLvmBetterThanLseLvm(RandomSeed seed) {
    fitSingleGaussianBetter(seed, NO_BOUND, NO_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM,
        NoiseModel.EMCCD);
  }

  @SeededTest
  public void fitSingleGaussianEmCcd_B_x__MleLvmBetterThanLseLvm(RandomSeed seed) {
    fitSingleGaussianBetter(seed, BOUNDED, NO_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM,
        NoiseModel.EMCCD);
  }

  @SeededTest
  public void fitSingleGaussianEmCcd_x_C__MleLvmBetterThanLseLvm(RandomSeed seed) {
    fitSingleGaussianBetter(seed, NO_BOUND, CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM,
        NoiseModel.EMCCD);
  }

  @SeededTest
  public void fitSingleGaussianEmCcd_B_C__MleLvmBetterThanLseLvm(RandomSeed seed) {
    fitSingleGaussianBetter(seed, BOUNDED, CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM,
        NoiseModel.EMCCD);
  }

  @SeededTest
  public void fitSingleGaussianEmCcd_x_DC_MleLvmBetterThanLseLvm(RandomSeed seed) {
    fitSingleGaussianBetter(seed, NO_BOUND, DYNAMIC_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM,
        NoiseModel.EMCCD);
  }

  @SeededTest
  public void fitSingleGaussianEmCcd_B_DC_MleLvmBetterThanLseLvm(RandomSeed seed) {
    fitSingleGaussianBetter(seed, BOUNDED, DYNAMIC_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM,
        NoiseModel.EMCCD);
  }

  @SeededTest
  public void fitSingleGaussianEmCcd_B_x__MleLvmBetterThanMleLvm(RandomSeed seed) {
    fitSingleGaussianBetter(seed, BOUNDED, NO_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, MLELVM,
        NoiseModel.EMCCD);
  }

  @SeededTest
  public void fitSingleGaussianEmCcd_x_C__MleLvmBetterThanMleLvm(RandomSeed seed) {
    fitSingleGaussianBetter(seed, NO_BOUND, CLAMP, MLELVM, NO_BOUND, NO_CLAMP, MLELVM,
        NoiseModel.EMCCD);
  }

  @SeededTest
  public void fitSingleGaussianEmCcd_x_DC_MleLvmBetterThanMleLvm(RandomSeed seed) {
    fitSingleGaussianBetter(seed, NO_BOUND, DYNAMIC_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, MLELVM,
        NoiseModel.EMCCD);
  }

  @SeededTest
  public void fitSingleGaussianEmCcd_B_DC_MleLvmBetterThanMleLvm(RandomSeed seed) {
    fitSingleGaussianBetter(seed, BOUNDED, DYNAMIC_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, MLELVM,
        NoiseModel.EMCCD);
  }

  @SeededTest
  public void fitSingleGaussianEmCcd_B_x__MleLvmBetterThanBLseLvm(RandomSeed seed) {
    fitSingleGaussianBetter(seed, BOUNDED, NO_CLAMP, MLELVM, BOUNDED, NO_CLAMP, LSELVM,
        NoiseModel.EMCCD);
  }

  @SeededTest
  public void fitSingleGaussianEmCcd_B_C__MleLvmBetterThanBcLseLvm(RandomSeed seed) {
    fitSingleGaussianBetter(seed, BOUNDED, CLAMP, MLELVM, BOUNDED, CLAMP, LSELVM, NoiseModel.EMCCD);
  }

  @SeededTest
  public void fitSingleGaussianEmCcd_B_DC_MleLvmBetterThanBdcLseLvm(RandomSeed seed) {
    fitSingleGaussianBetter(seed, BOUNDED, DYNAMIC_CLAMP, MLELVM, BOUNDED, DYNAMIC_CLAMP, LSELVM,
        NoiseModel.EMCCD);
  }

  // Note: The FastLogMLELVM converges too fast when there is still some
  // optimisation to do. The following tests have far fewer iterations with the fastLog version.

  @SeededTest
  public void fitSingleGaussianEmCcd_x_x__LseLvmBetterThanFastLogMleLvm(RandomSeed seed) {
    // This is actually a tie with the current fixed random number generator (50/50 for each)
    fitSingleGaussianBetter(seed, NO_BOUND, NO_CLAMP, LSELVM, NO_BOUND, NO_CLAMP, FastLogMLELVM,
        NoiseModel.EMCCD);
  }

  @SeededTest
  public void fitSingleGaussianEmCcd_x_x__MleLvmBetterThanFastLogMleLvm(RandomSeed seed) {
    fitSingleGaussianBetter(seed, NO_BOUND, NO_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, FastLogMLELVM,
        NoiseModel.EMCCD);
  }

  private void fitSingleGaussianBetter(RandomSeed seed, boolean bounded2,
      SteppingFunctionSolverClamp clamp2, SteppingFunctionSolverType type2, boolean bounded,
      SteppingFunctionSolverClamp clamp, SteppingFunctionSolverType type, NoiseModel noiseModel) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    final SteppingFunctionSolver solver = getSolver(clamp, type);
    final SteppingFunctionSolver solver2 = getSolver(clamp2, type2);
    canFitSingleGaussianBetter(seed, solver, bounded, solver2, bounded2,
        getName(bounded, clamp, type), getName(bounded2, clamp2, type2), noiseModel);
  }

  @SeededTest
  public void canFitAndComputeDeviationsEmCcd_LseLvm(RandomSeed seed) {
    canFitAndComputeDeviations(seed, SteppingFunctionSolverType.LSELVM, NoiseModel.EMCCD, false);
  }

  @SeededTest
  public void canFitAndComputeDeviationsScmos_LseLvm(RandomSeed seed) {
    canFitAndComputeDeviations(seed, SteppingFunctionSolverType.LSELVM, NoiseModel.SCMOS, false);
  }

  @SeededTest
  public void canFitAndComputeDeviationsEmCcd_WlseLvm(RandomSeed seed) {
    canFitAndComputeDeviations(seed, SteppingFunctionSolverType.WLSELVM, NoiseModel.EMCCD, false);
  }

  @SeededTest
  public void canFitAndComputeDeviationsScmos_WlseLvm(RandomSeed seed) {
    canFitAndComputeDeviations(seed, SteppingFunctionSolverType.WLSELVM, NoiseModel.SCMOS, false);
  }

  @SeededTest
  public void canFitAndComputeDeviationsScmos_WlseLvm_Weighted(RandomSeed seed) {
    canFitAndComputeDeviations(seed, SteppingFunctionSolverType.WLSELVM, NoiseModel.SCMOS, true);
  }

  @SeededTest
  public void canFitAndComputeDeviationsEmCcd_MleLvm(RandomSeed seed) {
    canFitAndComputeDeviations(seed, SteppingFunctionSolverType.MLELVM, NoiseModel.EMCCD, false);
  }

  @SeededTest
  public void canFitAndComputeDeviationsScmos_MleLvm(RandomSeed seed) {
    canFitAndComputeDeviations(seed, SteppingFunctionSolverType.MLELVM, NoiseModel.SCMOS, false);
  }

  @SeededTest
  public void canFitAndComputeDeviationsScmos_MleLvm_Weighted(RandomSeed seed) {
    canFitAndComputeDeviations(seed, SteppingFunctionSolverType.MLELVM, NoiseModel.SCMOS, true);
  }

  @SeededTest
  public void canFitAndComputeDeviationsEmCcd_FastMle(RandomSeed seed) {
    canFitAndComputeDeviations(seed, SteppingFunctionSolverType.FastMLE, NoiseModel.EMCCD, false);
  }

  @SeededTest
  public void canFitAndComputeDeviationsScmos_FastMle(RandomSeed seed) {
    canFitAndComputeDeviations(seed, SteppingFunctionSolverType.FastMLE, NoiseModel.SCMOS, false);
  }

  @SeededTest
  public void canFitAndComputeDeviationsScmos_FastMle_Weighted(RandomSeed seed) {
    canFitAndComputeDeviations(seed, SteppingFunctionSolverType.FastMLE, NoiseModel.SCMOS, true);
  }

  @SeededTest
  public void canFitAndComputeDeviationsEmCcd_BtFastMle(RandomSeed seed) {
    canFitAndComputeDeviations(seed, SteppingFunctionSolverType.BtFastMLE, NoiseModel.EMCCD, false);
  }

  @SeededTest
  public void canFitAndComputeDeviationsScmos_BtFastMle(RandomSeed seed) {
    canFitAndComputeDeviations(seed, SteppingFunctionSolverType.BtFastMLE, NoiseModel.SCMOS, false);
  }

  @SeededTest
  public void canFitAndComputeDeviationsScmos_BtFastMle_Weighted(RandomSeed seed) {
    canFitAndComputeDeviations(seed, SteppingFunctionSolverType.BtFastMLE, NoiseModel.SCMOS, true);
  }

  private void canFitAndComputeDeviations(RandomSeed seed, SteppingFunctionSolverType type,
      NoiseModel noiseModel, boolean useWeights) {
    final SteppingFunctionSolver solver1 =
        getSolver(SteppingFunctionSolverClamp.NO_CLAMP, type, noToleranceChecker);
    final SteppingFunctionSolver solver2 =
        getSolver(SteppingFunctionSolverClamp.NO_CLAMP, type, noToleranceChecker);
    fitAndComputeDeviationsMatch(seed, solver1, solver2, noiseModel, useWeights);
  }

  @SeededTest
  public void canFitAndComputeValueEmCcd_LseLvm(RandomSeed seed) {
    canFitAndComputeValue(seed, SteppingFunctionSolverType.LSELVM, NoiseModel.EMCCD, false);
  }

  @SeededTest
  public void canFitAndComputeValueScmos_LseLvm(RandomSeed seed) {
    canFitAndComputeValue(seed, SteppingFunctionSolverType.LSELVM, NoiseModel.SCMOS, false);
  }

  @SeededTest
  public void canFitAndComputeValueEmCcd_WlseLvm(RandomSeed seed) {
    canFitAndComputeValue(seed, SteppingFunctionSolverType.WLSELVM, NoiseModel.EMCCD, false);
  }

  @SeededTest
  public void canFitAndComputeValueScmos_WlseLvm(RandomSeed seed) {
    canFitAndComputeValue(seed, SteppingFunctionSolverType.WLSELVM, NoiseModel.SCMOS, false);
  }

  @SeededTest
  public void canFitAndComputeValueScmos_WlseLvm_Weighted(RandomSeed seed) {
    canFitAndComputeValue(seed, SteppingFunctionSolverType.WLSELVM, NoiseModel.SCMOS, true);
  }

  @SeededTest
  public void canFitAndComputeValueEmCcd_MleLvm(RandomSeed seed) {
    canFitAndComputeValue(seed, SteppingFunctionSolverType.MLELVM, NoiseModel.EMCCD, false);
  }

  @SeededTest
  public void canFitAndComputeValueScmos_MleLvm(RandomSeed seed) {
    canFitAndComputeValue(seed, SteppingFunctionSolverType.MLELVM, NoiseModel.SCMOS, false);
  }

  @SeededTest
  public void canFitAndComputeValueScmos_MleLvm_Weighted(RandomSeed seed) {
    canFitAndComputeValue(seed, SteppingFunctionSolverType.MLELVM, NoiseModel.SCMOS, true);
  }

  @SeededTest
  public void canFitAndComputeValueEmCcd_FastMle(RandomSeed seed) {
    canFitAndComputeValue(seed, SteppingFunctionSolverType.FastMLE, NoiseModel.EMCCD, false);
  }

  @SeededTest
  public void canFitAndComputeValueScmos_FastMle(RandomSeed seed) {
    canFitAndComputeValue(seed, SteppingFunctionSolverType.FastMLE, NoiseModel.SCMOS, false);
  }

  @SeededTest
  public void canFitAndComputeValueScmos_FastMle_Weighted(RandomSeed seed) {
    canFitAndComputeValue(seed, SteppingFunctionSolverType.FastMLE, NoiseModel.SCMOS, true);
  }

  @SeededTest
  public void canFitAndComputeValueEmCcd_BtFastMle(RandomSeed seed) {
    canFitAndComputeValue(seed, SteppingFunctionSolverType.BtFastMLE, NoiseModel.EMCCD, false);
  }

  @SeededTest
  public void canFitAndComputeValueScmos_BtFastMle(RandomSeed seed) {
    canFitAndComputeValue(seed, SteppingFunctionSolverType.BtFastMLE, NoiseModel.SCMOS, false);
  }

  @SeededTest
  public void canFitAndComputeValueScmos_BtFastMle_Weighted(RandomSeed seed) {
    canFitAndComputeValue(seed, SteppingFunctionSolverType.BtFastMLE, NoiseModel.SCMOS, true);
  }

  private void canFitAndComputeValue(RandomSeed seed, SteppingFunctionSolverType type,
      NoiseModel noiseModel, boolean useWeights) {
    final SteppingFunctionSolver solver1 =
        getSolver(SteppingFunctionSolverClamp.NO_CLAMP, type, noToleranceChecker);
    final SteppingFunctionSolver solver2 =
        getSolver(SteppingFunctionSolverClamp.NO_CLAMP, type, noToleranceChecker);
    fitAndComputeValueMatch(seed, solver1, solver2, noiseModel, useWeights);
  }
}
