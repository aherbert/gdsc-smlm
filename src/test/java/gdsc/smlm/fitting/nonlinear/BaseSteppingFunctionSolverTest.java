package gdsc.smlm.fitting.nonlinear;

import gdsc.core.utils.NotImplementedException;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;

/**
 * Test that a stepping solver can return the same results with and without bounds.
 */
public abstract class BaseSteppingFunctionSolverTest extends BaseFunctionSolverTest
{
	enum SteppingFunctionSolverType
	{
		MLELVM, LSELVM, WLSELVM, MLENR
	}

	// For convenience declare vairable of the enum type
	SteppingFunctionSolverType MLELVM = SteppingFunctionSolverType.MLELVM;
	SteppingFunctionSolverType LSELVM = SteppingFunctionSolverType.LSELVM;
	SteppingFunctionSolverType WlSELVM = SteppingFunctionSolverType.WLSELVM;
	SteppingFunctionSolverType MLENR = SteppingFunctionSolverType.MLENR;

	SteppingFunctionSolver getSolver(int clamping, SteppingFunctionSolverType type)
	{
		ErfGaussian2DFunction f = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, size, size,
				GaussianFunctionFactory.FIT_ERF_CIRCLE, null);
		ToleranceChecker tc = new ToleranceChecker(1e-5, 1e-5, 0, 0, 100);
		ParameterBounds bounds = new ParameterBounds(f);
		if (clamping != 0)
		{
			bounds.setClampValues(defaultClampValues);
			bounds.setDynamicClamp(clamping == 2);
		}
		SteppingFunctionSolver solver;
		switch (type)
		{
			case LSELVM:
				solver = new LSELVMSteppingFunctionSolver(f, tc, bounds);
				break;
			case MLELVM:
				solver = new MLELVMSteppingFunctionSolver(f, tc, bounds);
				break;
			case WLSELVM:
				solver = new WLSELVMSteppingFunctionSolver(f, tc, bounds);
				break;
			case MLENR:
				solver = new NewtonRaphsonSteppingFunctionSolver(f, tc, bounds);
				break;
			default:
				throw new NotImplementedException();
		}
		if (solver instanceof LVMSteppingFunctionSolver)
			((LVMSteppingFunctionSolver) solver).setInitialLambda(1);
		return solver;
	}

	String getName(boolean bounded, int clamping, SteppingFunctionSolverType type)
	{
		return ((bounded) ? "B" : "") + ((clamping == 0) ? "" : ((clamping == 1) ? "C" : "DC")) + type;
	}
}
