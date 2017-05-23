package gdsc.smlm.fitting.nonlinear;

import gdsc.core.utils.NotImplementedException;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;

/**
 * Test that a stepping solver can return the same results with and without bounds.
 */
public abstract class BaseSteppingFunctionSolverTest extends BaseFunctionSolverTest
{
	enum SteppingFunctionSolverClamp
	{
		NO_CLAMP("x"), CLAMP("C"), DYNAMIC_CLAMP("DC");
		final String name;

		SteppingFunctionSolverClamp(String name)
		{
			this.name = name;
		}

		@Override
		public String toString()
		{
			return name;
		}
	}

	enum SteppingFunctionSolverType
	{
		MLELVM, LSELVM, WLSELVM, FastMLE
	}

	// For convenience declare variables of the enum type
	static final SteppingFunctionSolverClamp NO_CLAMP = SteppingFunctionSolverClamp.NO_CLAMP;
	static final SteppingFunctionSolverClamp CLAMP = SteppingFunctionSolverClamp.CLAMP;
	static final SteppingFunctionSolverClamp DYNAMIC_CLAMP = SteppingFunctionSolverClamp.DYNAMIC_CLAMP;
	static final SteppingFunctionSolverType MLELVM = SteppingFunctionSolverType.MLELVM;
	static final SteppingFunctionSolverType LSELVM = SteppingFunctionSolverType.LSELVM;
	static final SteppingFunctionSolverType WLSELVM = SteppingFunctionSolverType.WLSELVM;
	static final SteppingFunctionSolverType FastMLE = SteppingFunctionSolverType.FastMLE;
	static final boolean BOUNDED = true;
	static final boolean NO_BOUND = false;

	SteppingFunctionSolver getSolver(SteppingFunctionSolverClamp clamp, SteppingFunctionSolverType type)
	{
		ErfGaussian2DFunction f = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, size, size,
				GaussianFunctionFactory.FIT_ERF_CIRCLE, null);
		ToleranceChecker tc = new ToleranceChecker(1e-5, 1e-5, 0, 0, 100);
		ParameterBounds bounds = new ParameterBounds(f);
		switch (clamp)
		{
			case DYNAMIC_CLAMP:
				bounds.setDynamicClamp(true);
			case CLAMP:
				bounds.setClampValues(defaultClampValues);
			case NO_CLAMP:
			default:
				break;
		}
		SteppingFunctionSolver solver;
		switch (type)
		{
			case LSELVM:
				solver = new LSELVMSteppingFunctionSolver(f, tc, bounds);
				break;
			case MLELVM:
				solver = new MLELVMSteppingFunctionSolver(f, tc, bounds);
				// MLE requires a positive function value so use a lower bound
				solver.setBounds(new double[7], null);
				break;
			case WLSELVM:
				solver = new WLSELVMSteppingFunctionSolver(f, tc, bounds);
				break;
			case FastMLE:
				solver = new FastMLESteppingFunctionSolver(f, tc, bounds);
				// MLE requires a positive function value so use a lower bound
				solver.setBounds(new double[7], null);
				break;
			default:
				throw new NotImplementedException();
		}
		if (solver instanceof LVMSteppingFunctionSolver)
			((LVMSteppingFunctionSolver) solver).setInitialLambda(1);
		return solver;
	}

	String getName(boolean bounded, SteppingFunctionSolverClamp clamp, SteppingFunctionSolverType type)
	{
		return String.format("%s %-2s %s", (bounded) ? "B" : "x", clamp, type);
	}
}
