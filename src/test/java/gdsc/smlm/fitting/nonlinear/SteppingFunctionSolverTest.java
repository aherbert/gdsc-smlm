package gdsc.smlm.fitting.nonlinear;

import org.junit.Test;

/**
 * Test that a stepping solver can fit a function.
 */
public class SteppingFunctionSolverTest extends BaseSteppingFunctionSolverTest
{
	@Test
	public void canFitSingleGaussianLSELSELVM()
	{
		fitSingleGaussian(NO_BOUND, NO_CLAMP, LSELVM);
	}

	@Test
	public void canFitSingleGaussianCLSELSELVM()
	{
		fitSingleGaussian(NO_BOUND, CLAMP, LSELVM);
	}

	@Test
	public void canFitSingleGaussianDCLSELSELVM()
	{
		fitSingleGaussian(NO_BOUND, DYNAMIC_CLAMP, LSELVM);
	}

	@Test
	public void canFitSingleGaussianBLSELSELVM()
	{
		fitSingleGaussian(BOUNDED, NO_CLAMP, LSELVM);
	}

	@Test
	public void canFitSingleGaussianBCLSELSELVM()
	{
		fitSingleGaussian(BOUNDED, CLAMP, LSELVM);
	}

	@Test
	public void canFitSingleGaussianBDCLSELSELVM()
	{
		fitSingleGaussian(BOUNDED, DYNAMIC_CLAMP, LSELVM);
	}

	@Test
	public void canFitSingleGaussianLSEWLSELVM()
	{
		fitSingleGaussian(NO_BOUND, NO_CLAMP, WLSELVM);
	}

	@Test
	public void canFitSingleGaussianCLSEWLSELVM()
	{
		fitSingleGaussian(NO_BOUND, CLAMP, WLSELVM);
	}

	@Test
	public void canFitSingleGaussianDCLSEWLSELVM()
	{
		fitSingleGaussian(NO_BOUND, DYNAMIC_CLAMP, WLSELVM);
	}

	@Test
	public void canFitSingleGaussianBLSEWLSELVM()
	{
		fitSingleGaussian(BOUNDED, NO_CLAMP, WLSELVM);
	}

	@Test
	public void canFitSingleGaussianBCLSEWLSELVM()
	{
		fitSingleGaussian(BOUNDED, CLAMP, WLSELVM);
	}

	@Test
	public void canFitSingleGaussianBDCLSEWLSELVM()
	{
		fitSingleGaussian(BOUNDED, DYNAMIC_CLAMP, WLSELVM);
	}

	@Test
	public void canFitSingleGaussianLSEMLELVM()
	{
		fitSingleGaussian(NO_BOUND, NO_CLAMP, MLELVM);
	}

	@Test
	public void canFitSingleGaussianCLSEMLELVM()
	{
		fitSingleGaussian(NO_BOUND, CLAMP, MLELVM);
	}

	@Test
	public void canFitSingleGaussianDCLSEMLELVM()
	{
		fitSingleGaussian(NO_BOUND, DYNAMIC_CLAMP, MLELVM);
	}

	@Test
	public void canFitSingleGaussianBLSEMLELVM()
	{
		fitSingleGaussian(BOUNDED, NO_CLAMP, MLELVM);
	}

	@Test
	public void canFitSingleGaussianBCLSEMLELVM()
	{
		fitSingleGaussian(BOUNDED, CLAMP, MLELVM);
	}

	@Test
	public void canFitSingleGaussianBDCLSEMLELVM()
	{
		fitSingleGaussian(BOUNDED, DYNAMIC_CLAMP, MLELVM);
	}

	@Test
	public void canFitSingleGaussianMLENR()
	{
		fitSingleGaussian(NO_BOUND, NO_CLAMP, FastMLE);
	}

	@Test
	public void canFitSingleGaussianCMLENR()
	{
		fitSingleGaussian(NO_BOUND, CLAMP, FastMLE);
	}

	@Test
	public void canFitSingleGaussianDCMLENR()
	{
		fitSingleGaussian(NO_BOUND, DYNAMIC_CLAMP, FastMLE);
	}

	@Test
	public void canFitSingleGaussianBMLENR()
	{
		fitSingleGaussian(BOUNDED, NO_CLAMP, FastMLE);
	}

	@Test
	public void canFitSingleGaussianBCMLENR()
	{
		fitSingleGaussian(BOUNDED, CLAMP, FastMLE);
	}

	@Test
	public void canFitSingleGaussianBDCMLENR()
	{
		fitSingleGaussian(BOUNDED, DYNAMIC_CLAMP, FastMLE);
	}

	private void fitSingleGaussian(boolean bounded, SteppingFunctionSolverClamp clamp, SteppingFunctionSolverType type)
	{
		canFitSingleGaussian(getSolver(clamp, type), bounded);
	}

	// Is Bounded/Clamped better?

	@Test
	public void fitSingleGaussianBLSELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, NO_CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM);
	}

	@Test
	public void fitSingleGaussianCLSELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM);
	}

	@Test
	public void fitSingleGaussianBCLSELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM);
	}

	@Test
	public void fitSingleGaussianDCLSELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, DYNAMIC_CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM);
	}

	@Test
	public void fitSingleGaussianBDCLSELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, DYNAMIC_CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM);
	}

	@Test
	public void fitSingleGaussianMLELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, NO_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM);
	}

	@Test
	public void fitSingleGaussianBMLELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, NO_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM);
	}

	@Test
	public void fitSingleGaussianCMLELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM);
	}

	@Test
	public void fitSingleGaussianBCMLELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM);
	}

	@Test
	public void fitSingleGaussianDCMLELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, DYNAMIC_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM);
	}

	@Test
	public void fitSingleGaussianBDCMLELVMBetterThanLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, DYNAMIC_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM);
	}

	@Test
	public void fitSingleGaussianBMLELVMBetterThanMLELVM()
	{
		fitSingleGaussianBetter(BOUNDED, NO_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, MLELVM);
	}

	@Test
	public void fitSingleGaussianCMLELVMBetterThanMLELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, CLAMP, MLELVM, NO_BOUND, NO_CLAMP, MLELVM);
	}

	@Test
	public void fitSingleGaussianDCMLELVMBetterThanMLELVM()
	{
		fitSingleGaussianBetter(NO_BOUND, DYNAMIC_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, MLELVM);
	}

	@Test
	public void fitSingleGaussianBDCMLELVMBetterThanMLELVM()
	{
		fitSingleGaussianBetter(BOUNDED, DYNAMIC_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, MLELVM);
	}

	@Test
	public void fitSingleGaussianBMLELVMBetterThanBLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, NO_CLAMP, MLELVM, BOUNDED, NO_CLAMP, LSELVM);
	}

	@Test
	public void fitSingleGaussianBCMLELVMBetterThanBCLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, CLAMP, MLELVM, BOUNDED, CLAMP, LSELVM);
	}

	@Test
	public void fitSingleGaussianBDCMLELVMBetterThanBDCLSELVM()
	{
		fitSingleGaussianBetter(BOUNDED, DYNAMIC_CLAMP, MLELVM, BOUNDED, DYNAMIC_CLAMP, LSELVM);
	}

	private void fitSingleGaussianBetter(boolean bounded2, SteppingFunctionSolverClamp clamp2,
			SteppingFunctionSolverType type2, boolean bounded, SteppingFunctionSolverClamp clamp,
			SteppingFunctionSolverType type)
	{
		org.junit.Assume.assumeTrue(false);
		
		SteppingFunctionSolver solver = getSolver(clamp, type);
		SteppingFunctionSolver solver2 = getSolver(clamp2, type2);
		canFitSingleGaussianBetter(solver, bounded, solver2, bounded2, getName(bounded, clamp, type),
				getName(bounded2, clamp2, type2));
	}
}
