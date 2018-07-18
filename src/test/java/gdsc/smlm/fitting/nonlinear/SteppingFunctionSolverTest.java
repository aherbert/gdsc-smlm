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
package gdsc.smlm.fitting.nonlinear;

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Test;

import gdsc.test.TestSettings;

/**
 * Test that a stepping solver can fit a function.
 */
@SuppressWarnings({ "javadoc" })
public class SteppingFunctionSolverTest extends BaseSteppingFunctionSolverTest
{
	@Test
	public void canFitSingleGaussianEMCCD_x_x__LSELVM()
	{
		fitSingleGaussian(NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_x_C__LSELVM()
	{
		fitSingleGaussian(NO_BOUND, CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_x_DC_LSELVM()
	{
		fitSingleGaussian(NO_BOUND, DYNAMIC_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_B_x__LSELVM()
	{
		fitSingleGaussian(BOUNDED, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_B_C__LSELVM()
	{
		fitSingleGaussian(BOUNDED, CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_B_DC_LSELVM()
	{
		fitSingleGaussian(BOUNDED, DYNAMIC_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_x_x__WLSELVM()
	{
		fitSingleGaussian(NO_BOUND, NO_CLAMP, WLSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_x_C__WLSELVM()
	{
		fitSingleGaussian(NO_BOUND, CLAMP, WLSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_x_DC_WLSELVM()
	{
		fitSingleGaussian(NO_BOUND, DYNAMIC_CLAMP, WLSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_B_x__WLSELVM()
	{
		fitSingleGaussian(BOUNDED, NO_CLAMP, WLSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_B_C__WLSELVM()
	{
		fitSingleGaussian(BOUNDED, CLAMP, WLSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_B_DC_WLSELVM()
	{
		fitSingleGaussian(BOUNDED, DYNAMIC_CLAMP, WLSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_x_x__MLELVM()
	{
		fitSingleGaussian(NO_BOUND, NO_CLAMP, MLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_x_C__MLELVM()
	{
		fitSingleGaussian(NO_BOUND, CLAMP, MLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_x_DC_MLELVM()
	{
		fitSingleGaussian(NO_BOUND, DYNAMIC_CLAMP, MLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_B_x__MLELVM()
	{
		fitSingleGaussian(BOUNDED, NO_CLAMP, MLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_B_C__MLELVM()
	{
		fitSingleGaussian(BOUNDED, CLAMP, MLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_B_DC_MLELVM()
	{
		fitSingleGaussian(BOUNDED, DYNAMIC_CLAMP, MLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_x_x__FastLogMLELVM()
	{
		fitSingleGaussian(NO_BOUND, NO_CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_x_C__FastLogMLELVM()
	{
		fitSingleGaussian(NO_BOUND, CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_x_DC_FastLogMLELVM()
	{
		fitSingleGaussian(NO_BOUND, DYNAMIC_CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_B_x__FastLogMLELVM()
	{
		fitSingleGaussian(BOUNDED, NO_CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_B_C__FastLogMLELVM()
	{
		fitSingleGaussian(BOUNDED, CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_B_DC_FastLogMLELVM()
	{
		fitSingleGaussian(BOUNDED, DYNAMIC_CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_x_x__FastMLE()
	{
		// The FastMLE method can generate very big steps that make the method unstable.
		// This test may fail depending on the random number generator.
		try
		{
			fitSingleGaussian(NO_BOUND, NO_CLAMP, FastMLE, NoiseModel.EMCCD);
		}
		catch (final AssertionError e)
		{
			TestSettings.logFailure(e);
		}
	}

	@Test
	public void canFitSingleGaussianEMCCD_x_C__FastMLE()
	{
		fitSingleGaussian(NO_BOUND, CLAMP, FastMLE, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_x_DC_FastMLE()
	{
		fitSingleGaussian(NO_BOUND, DYNAMIC_CLAMP, FastMLE, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_B_x__FastMLE()
	{
		fitSingleGaussian(BOUNDED, NO_CLAMP, FastMLE, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_B_C__FastMLE()
	{
		fitSingleGaussian(BOUNDED, CLAMP, FastMLE, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_B_DC_FastMLE()
	{
		fitSingleGaussian(BOUNDED, DYNAMIC_CLAMP, FastMLE, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_x_x__BTFastMLE()
	{
		fitSingleGaussian(NO_BOUND, NO_CLAMP, BTFastMLE, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_x_C__BTFastMLE()
	{
		fitSingleGaussian(NO_BOUND, CLAMP, BTFastMLE, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_x_DC_BTFastMLE()
	{
		fitSingleGaussian(NO_BOUND, DYNAMIC_CLAMP, BTFastMLE, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_B_x__BTFastMLE()
	{
		fitSingleGaussian(BOUNDED, NO_CLAMP, BTFastMLE, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_B_C__BTFastMLE()
	{
		fitSingleGaussian(BOUNDED, CLAMP, BTFastMLE, NoiseModel.EMCCD);
	}

	@Test
	public void canFitSingleGaussianEMCCD_B_DC_BTFastMLE()
	{
		fitSingleGaussian(BOUNDED, DYNAMIC_CLAMP, BTFastMLE, NoiseModel.EMCCD);
	}

	@Test(expected = AssertionError.class)
	public void cannotFitSingleGaussianEMCCD_x_x__JFastMLE()
	{
		// The JFastMLE method was built using a misinterpretation of the Newton
		// method in Numerical Recipes, 2nd Ed. This test is just here to prove that.
		TestSettings.assumeMaximumComplexity();

		// The JFastMLE method does not work
		fitSingleGaussian(NO_BOUND, NO_CLAMP, JFastMLE, NoiseModel.EMCCD);
	}

	//	@Test
	//	public void canFitSingleGaussianEMCCD_x_C__JFastMLE()
	//	{
	//		fitSingleGaussian(NO_BOUND, CLAMP, JFastMLE, NoiseModel.EMCCD);
	//	}
	//
	//	@Test
	//	public void canFitSingleGaussianEMCCD_x_DC_JFastMLE()
	//	{
	//		fitSingleGaussian(NO_BOUND, DYNAMIC_CLAMP, JFastMLE, NoiseModel.EMCCD);
	//	}
	//
	//	@Test
	//	public void canFitSingleGaussianEMCCD_B_x__JFastMLE()
	//	{
	//		fitSingleGaussian(BOUNDED, NO_CLAMP, JFastMLE, NoiseModel.EMCCD);
	//	}
	//
	//	@Test
	//	public void canFitSingleGaussianEMCCD_B_C__JFastMLE()
	//	{
	//		fitSingleGaussian(BOUNDED, CLAMP, JFastMLE, NoiseModel.EMCCD);
	//	}
	//
	//	@Test
	//	public void canFitSingleGaussianEMCCD_B_DC_JFastMLE()
	//	{
	//		fitSingleGaussian(BOUNDED, DYNAMIC_CLAMP, JFastMLE, NoiseModel.EMCCD);
	//	}

	// Weighted solvers for sCMOS

	@Test
	public void canFitSingleGaussianSCMOS_x_x__WLSELVM()
	{
		fitSingleGaussian(NO_BOUND, NO_CLAMP, WLSELVM, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_x_C__WLSELVM()
	{
		fitSingleGaussian(NO_BOUND, CLAMP, WLSELVM, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_x_DC_WLSELVM()
	{
		fitSingleGaussian(NO_BOUND, DYNAMIC_CLAMP, WLSELVM, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_B_x__WLSELVM()
	{
		fitSingleGaussian(BOUNDED, NO_CLAMP, WLSELVM, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_B_C__WLSELVM()
	{
		fitSingleGaussian(BOUNDED, CLAMP, WLSELVM, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_B_DC_WLSELVM()
	{
		fitSingleGaussian(BOUNDED, DYNAMIC_CLAMP, WLSELVM, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_x_x__MLELVM()
	{
		fitSingleGaussian(NO_BOUND, NO_CLAMP, MLELVM, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_x_C__MLELVM()
	{
		fitSingleGaussian(NO_BOUND, CLAMP, MLELVM, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_x_DC_MLELVM()
	{
		fitSingleGaussian(NO_BOUND, DYNAMIC_CLAMP, MLELVM, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_B_x__MLELVM()
	{
		fitSingleGaussian(BOUNDED, NO_CLAMP, MLELVM, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_B_C__MLELVM()
	{
		fitSingleGaussian(BOUNDED, CLAMP, MLELVM, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_B_DC_MLELVM()
	{
		fitSingleGaussian(BOUNDED, DYNAMIC_CLAMP, MLELVM, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_x_x__FastLogMLELVM()
	{
		fitSingleGaussian(NO_BOUND, NO_CLAMP, FastLogMLELVM, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_x_C__FastLogMLELVM()
	{
		fitSingleGaussian(NO_BOUND, CLAMP, FastLogMLELVM, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_x_DC_FastLogMLELVM()
	{
		fitSingleGaussian(NO_BOUND, DYNAMIC_CLAMP, FastLogMLELVM, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_B_x__FastLogMLELVM()
	{
		fitSingleGaussian(BOUNDED, NO_CLAMP, FastLogMLELVM, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_B_C__FastLogMLELVM()
	{
		fitSingleGaussian(BOUNDED, CLAMP, FastLogMLELVM, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_B_DC_FastLogMLELVM()
	{
		fitSingleGaussian(BOUNDED, DYNAMIC_CLAMP, FastLogMLELVM, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_x_x__FastMLE()
	{
		// The FastMLE method can generate very big steps that make the method unstable
		// This test may fail depending on the random number generator.
		try
		{
			fitSingleGaussian(NO_BOUND, NO_CLAMP, FastMLE, NoiseModel.SCMOS);
		}
		catch (final AssertionError e)
		{
			TestSettings.logFailure(e);
		}
	}

	@Test
	public void canFitSingleGaussianSCMOS_x_C__FastMLE()
	{
		fitSingleGaussian(NO_BOUND, CLAMP, FastMLE, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_x_DC_FastMLE()
	{
		fitSingleGaussian(NO_BOUND, DYNAMIC_CLAMP, FastMLE, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_B_x__FastMLE()
	{
		fitSingleGaussian(BOUNDED, NO_CLAMP, FastMLE, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_B_C__FastMLE()
	{
		fitSingleGaussian(BOUNDED, CLAMP, FastMLE, NoiseModel.SCMOS);
	}

	@Test
	public void canFitSingleGaussianSCMOS_B_DC_FastMLE()
	{
		fitSingleGaussian(BOUNDED, DYNAMIC_CLAMP, FastMLE, NoiseModel.SCMOS);
	}

	private void fitSingleGaussian(boolean bounded, SteppingFunctionSolverClamp clamp, SteppingFunctionSolverType type,
			NoiseModel noiseModel)
	{
		//org.junit.Assume.assumeTrue(false);
		final SteppingFunctionSolver solver = getSolver(clamp, type);
		canFitSingleGaussian(solver, bounded, noiseModel);
	}

	// Is Bounded/Clamped better?

	@Test
	public void fitSingleGaussianEMCCD_B_x__LSELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, NO_CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_x_C__LSELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_B_C__LSELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_x_DC_LSELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, DYNAMIC_CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_B_DC_LSELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, DYNAMIC_CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_x_x__MLELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, NO_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_B_x__MLELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, NO_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_x_C__MLELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_B_C__MLELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_x_DC_MLELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, DYNAMIC_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_B_DC_MLELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, DYNAMIC_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_B_x__MLELVMBetterThanMLELVM()
	{
		fitSingleGaussianBetter(BOUNDED, NO_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, MLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_x_C__MLELVMBetterThanMLELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, CLAMP, MLELVM, NO_BOUND, NO_CLAMP, MLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_x_DC_MLELVMBetterThanMLELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, DYNAMIC_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, MLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_B_DC_MLELVMBetterThanMLELVM()
	{
		fitSingleGaussianBetter(BOUNDED, DYNAMIC_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, MLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_B_x__MLELVMBetterThanBLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, NO_CLAMP, MLELVM, BOUNDED, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_B_C__MLELVMBetterThanBCLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, CLAMP, MLELVM, BOUNDED, CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_B_DC_MLELVMBetterThanBDCLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, DYNAMIC_CLAMP, MLELVM, BOUNDED, DYNAMIC_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	// Note: The FastLogMLELVM converges too fast when there is still some
	// optimisation to do. The following tests have far fewer iterations with the fastLog version.

	@Test
	public void fitSingleGaussianEMCCD_x_x__LSELVMBetterThanFastLogMLELVM()
	{
		// This is actually a tie with the current fixed random number generator (50/50 for each)
		fitSingleGaussianBetter(NO_BOUND, NO_CLAMP, LSELVM, NO_BOUND, NO_CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_x_x__MLELVMBetterThanFastLogMLELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, NO_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
	}

	private void fitSingleGaussianBetter(boolean bounded2, SteppingFunctionSolverClamp clamp2,
			SteppingFunctionSolverType type2, boolean bounded, SteppingFunctionSolverClamp clamp,
			SteppingFunctionSolverType type, NoiseModel noiseModel)
	{
		TestSettings.assumeMediumComplexity();
		final SteppingFunctionSolver solver = getSolver(clamp, type);
		final SteppingFunctionSolver solver2 = getSolver(clamp2, type2);
		canFitSingleGaussianBetter(solver, bounded, solver2, bounded2, getName(bounded, clamp, type),
				getName(bounded2, clamp2, type2), noiseModel);
	}

	@Test
	public void canFitAndComputeDeviationsEMCCD_LSELVM()
	{
		canFitAndComputeDeviations(SteppingFunctionSolverType.LSELVM, NoiseModel.EMCCD, false);
	}

	@Test
	public void canFitAndComputeDeviationsSCMOS_LSELVM()
	{
		canFitAndComputeDeviations(SteppingFunctionSolverType.LSELVM, NoiseModel.SCMOS, false);
	}

	@Test
	public void canFitAndComputeDeviationsEMCCD_WLSELVM()
	{
		canFitAndComputeDeviations(SteppingFunctionSolverType.WLSELVM, NoiseModel.EMCCD, false);
	}

	@Test
	public void canFitAndComputeDeviationsSCMOS_WLSELVM()
	{
		canFitAndComputeDeviations(SteppingFunctionSolverType.WLSELVM, NoiseModel.SCMOS, false);
	}

	@Test
	public void canFitAndComputeDeviationsSCMOS_WLSELVM_Weighted()
	{
		canFitAndComputeDeviations(SteppingFunctionSolverType.WLSELVM, NoiseModel.SCMOS, true);
	}

	@Test
	public void canFitAndComputeDeviationsEMCCD_MLELVM()
	{
		canFitAndComputeDeviations(SteppingFunctionSolverType.MLELVM, NoiseModel.EMCCD, false);
	}

	@Test
	public void canFitAndComputeDeviationsSCMOS_MLELVM()
	{
		canFitAndComputeDeviations(SteppingFunctionSolverType.MLELVM, NoiseModel.SCMOS, false);
	}

	@Test
	public void canFitAndComputeDeviationsSCMOS_MLELVM_Weighted()
	{
		canFitAndComputeDeviations(SteppingFunctionSolverType.MLELVM, NoiseModel.SCMOS, true);
	}

	@Test
	public void canFitAndComputeDeviationsEMCCD_FastMLE()
	{
		canFitAndComputeDeviations(SteppingFunctionSolverType.FastMLE, NoiseModel.EMCCD, false);
	}

	@Test
	public void canFitAndComputeDeviationsSCMOS_FastMLE()
	{
		canFitAndComputeDeviations(SteppingFunctionSolverType.FastMLE, NoiseModel.SCMOS, false);
	}

	@Test
	public void canFitAndComputeDeviationsSCMOS_FastMLE_Weighted()
	{
		canFitAndComputeDeviations(SteppingFunctionSolverType.FastMLE, NoiseModel.SCMOS, true);
	}

	@Test
	public void canFitAndComputeDeviationsEMCCD_BTFastMLE()
	{
		canFitAndComputeDeviations(SteppingFunctionSolverType.BTFastMLE, NoiseModel.EMCCD, false);
	}

	@Test
	public void canFitAndComputeDeviationsSCMOS_BTFastMLE()
	{
		canFitAndComputeDeviations(SteppingFunctionSolverType.BTFastMLE, NoiseModel.SCMOS, false);
	}

	@Test
	public void canFitAndComputeDeviationsSCMOS_BTFastMLE_Weighted()
	{
		canFitAndComputeDeviations(SteppingFunctionSolverType.BTFastMLE, NoiseModel.SCMOS, true);
	}

	private void canFitAndComputeDeviations(SteppingFunctionSolverType type, NoiseModel noiseModel, boolean useWeights)
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SteppingFunctionSolver solver1 = getSolver(SteppingFunctionSolverClamp.NO_CLAMP, type, noToleranceChecker);
		final SteppingFunctionSolver solver2 = getSolver(SteppingFunctionSolverClamp.NO_CLAMP, type, noToleranceChecker);
		fitAndComputeDeviationsMatch(rg, solver1, solver2, noiseModel, useWeights);
	}

	@Test
	public void canFitAndComputeValueEMCCD_LSELVM()
	{
		canFitAndComputeValue(SteppingFunctionSolverType.LSELVM, NoiseModel.EMCCD, false);
	}

	@Test
	public void canFitAndComputeValueSCMOS_LSELVM()
	{
		canFitAndComputeValue(SteppingFunctionSolverType.LSELVM, NoiseModel.SCMOS, false);
	}

	@Test
	public void canFitAndComputeValueEMCCD_WLSELVM()
	{
		canFitAndComputeValue(SteppingFunctionSolverType.WLSELVM, NoiseModel.EMCCD, false);
	}

	@Test
	public void canFitAndComputeValueSCMOS_WLSELVM()
	{
		canFitAndComputeValue(SteppingFunctionSolverType.WLSELVM, NoiseModel.SCMOS, false);
	}

	@Test
	public void canFitAndComputeValueSCMOS_WLSELVM_Weighted()
	{
		canFitAndComputeValue(SteppingFunctionSolverType.WLSELVM, NoiseModel.SCMOS, true);
	}

	@Test
	public void canFitAndComputeValueEMCCD_MLELVM()
	{
		canFitAndComputeValue(SteppingFunctionSolverType.MLELVM, NoiseModel.EMCCD, false);
	}

	@Test
	public void canFitAndComputeValueSCMOS_MLELVM()
	{
		canFitAndComputeValue(SteppingFunctionSolverType.MLELVM, NoiseModel.SCMOS, false);
	}

	@Test
	public void canFitAndComputeValueSCMOS_MLELVM_Weighted()
	{
		canFitAndComputeValue(SteppingFunctionSolverType.MLELVM, NoiseModel.SCMOS, true);
	}

	@Test
	public void canFitAndComputeValueEMCCD_FastMLE()
	{
		canFitAndComputeValue(SteppingFunctionSolverType.FastMLE, NoiseModel.EMCCD, false);
	}

	@Test
	public void canFitAndComputeValueSCMOS_FastMLE()
	{
		canFitAndComputeValue(SteppingFunctionSolverType.FastMLE, NoiseModel.SCMOS, false);
	}

	@Test
	public void canFitAndComputeValueSCMOS_FastMLE_Weighted()
	{
		canFitAndComputeValue(SteppingFunctionSolverType.FastMLE, NoiseModel.SCMOS, true);
	}

	@Test
	public void canFitAndComputeValueEMCCD_BTFastMLE()
	{
		canFitAndComputeValue(SteppingFunctionSolverType.BTFastMLE, NoiseModel.EMCCD, false);
	}

	@Test
	public void canFitAndComputeValueSCMOS_BTFastMLE()
	{
		canFitAndComputeValue(SteppingFunctionSolverType.BTFastMLE, NoiseModel.SCMOS, false);
	}

	@Test
	public void canFitAndComputeValueSCMOS_BTFastMLE_Weighted()
	{
		canFitAndComputeValue(SteppingFunctionSolverType.BTFastMLE, NoiseModel.SCMOS, true);
	}

	private void canFitAndComputeValue(SteppingFunctionSolverType type, NoiseModel noiseModel, boolean useWeights)
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SteppingFunctionSolver solver1 = getSolver(SteppingFunctionSolverClamp.NO_CLAMP, type, noToleranceChecker);
		final SteppingFunctionSolver solver2 = getSolver(SteppingFunctionSolverClamp.NO_CLAMP, type, noToleranceChecker);
		fitAndComputeValueMatch(rg, solver1, solver2, noiseModel, useWeights);
	}
}
