package gdsc.smlm.fitting.nonlinear;

import org.junit.Test;

/**
 * Test that a stepping solver can fit a function.
 */
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

	@Test(expected = AssertionError.class)
	public void cannotFitSingleGaussianEMCCD_x_x__FastMLE()
	{
		// The FastMLE method can generate very big steps that make the method unstable.
		// Currently we expect this test to fail.
		fitSingleGaussian(NO_BOUND, NO_CLAMP, FastMLE, NoiseModel.EMCCD);
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
		// The BTFastMLE method can generate very big steps that make the method unstable
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

	@Test(expected = AssertionError.class)
	public void cannotFitSingleGaussianSCMOS_x_x__FastMLE()
	{
		// The FastMLE method can generate very big steps that make the method unstable
		fitSingleGaussian(NO_BOUND, NO_CLAMP, FastMLE, NoiseModel.SCMOS);
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
		SteppingFunctionSolver solver = getSolver(clamp, type);
		canFitSingleGaussian(solver, bounded, noiseModel);
	}

	// Is Bounded/Clamped better?

	@Test
	public void fitSingleGaussianEMCCD_B_LSELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, NO_CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_C__LSELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_B_CLSELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_DC_LSELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, DYNAMIC_CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_B_DC_LSELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, DYNAMIC_CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_MLELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, NO_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_B_MLELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, NO_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_C__MLELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_B_CMLELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_DC_MLELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, DYNAMIC_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_B_DC_MLELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, DYNAMIC_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_B_MLELVMBetterThanMLELVM()
	{
		fitSingleGaussianBetter(BOUNDED, NO_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, MLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_C__MLELVMBetterThanMLELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, CLAMP, MLELVM, NO_BOUND, NO_CLAMP, MLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_DC_MLELVMBetterThanMLELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, DYNAMIC_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, MLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_B_DC_MLELVMBetterThanMLELVM()
	{
		fitSingleGaussianBetter(BOUNDED, DYNAMIC_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, MLELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_B_MLELVMBetterThanBLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, NO_CLAMP, MLELVM, BOUNDED, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_B_CMLELVMBetterThanBCLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, CLAMP, MLELVM, BOUNDED, CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	@Test
	public void fitSingleGaussianEMCCD_B_DC_MLELVMBetterThanBDCLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, DYNAMIC_CLAMP, MLELVM, BOUNDED, DYNAMIC_CLAMP, LSELVM, NoiseModel.EMCCD);
	}

	private void fitSingleGaussianBetter(boolean bounded2, SteppingFunctionSolverClamp clamp2,
			SteppingFunctionSolverType type2, boolean bounded, SteppingFunctionSolverClamp clamp,
			SteppingFunctionSolverType type, NoiseModel noiseModel)
	{
		//org.junit.Assume.assumeTrue(false);

		SteppingFunctionSolver solver = getSolver(clamp, type);
		SteppingFunctionSolver solver2 = getSolver(clamp2, type2);
		canFitSingleGaussianBetter(solver, bounded, solver2, bounded2, getName(bounded, clamp, type),
				getName(bounded2, clamp2, type2), noiseModel);
	}
}
