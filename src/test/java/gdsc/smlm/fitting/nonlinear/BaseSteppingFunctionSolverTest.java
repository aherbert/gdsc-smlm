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

import gdsc.core.utils.NotImplementedException;
import gdsc.smlm.function.FastLogFactory;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
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
		// Enum names should al be uppercase but this is just for a test so ignore that convention
		MLELVM, FastLogMLELVM, LSELVM, WLSELVM, FastMLE, JFastMLE, BTFastMLE
	}

	// For convenience declare variables of the enum type
	static final SteppingFunctionSolverClamp NO_CLAMP = SteppingFunctionSolverClamp.NO_CLAMP;
	static final SteppingFunctionSolverClamp CLAMP = SteppingFunctionSolverClamp.CLAMP;
	static final SteppingFunctionSolverClamp DYNAMIC_CLAMP = SteppingFunctionSolverClamp.DYNAMIC_CLAMP;
	static final SteppingFunctionSolverType MLELVM = SteppingFunctionSolverType.MLELVM;
	static final SteppingFunctionSolverType FastLogMLELVM = SteppingFunctionSolverType.FastLogMLELVM;
	static final SteppingFunctionSolverType LSELVM = SteppingFunctionSolverType.LSELVM;
	static final SteppingFunctionSolverType WLSELVM = SteppingFunctionSolverType.WLSELVM;
	static final SteppingFunctionSolverType FastMLE = SteppingFunctionSolverType.FastMLE;
	static final SteppingFunctionSolverType JFastMLE = SteppingFunctionSolverType.JFastMLE;
	static final SteppingFunctionSolverType BTFastMLE = SteppingFunctionSolverType.BTFastMLE;
	static final boolean BOUNDED = true;
	static final boolean NO_BOUND = false;

	static class NoToleranceChecker extends ToleranceChecker
	{
		public NoToleranceChecker()
		{
			super(0, 0);
		}

		@Override
		public int converged(double previousValue, double[] previousParameters, double currentValue,
				double[] currentParameters)
		{
			return STATUS_VALUE;
		}
	}

	static NoToleranceChecker noToleranceChecker = new NoToleranceChecker();

	SteppingFunctionSolver getSolver(SteppingFunctionSolverClamp clamp, SteppingFunctionSolverType type)
	{
		return getSolver(clamp, type, new ToleranceChecker(1e-5, 1e-5, 0, 0, 100));
	}

	@SuppressWarnings("deprecation")
	SteppingFunctionSolver getSolver(SteppingFunctionSolverClamp clamp, SteppingFunctionSolverType type,
			ToleranceChecker tc)
	{
		ErfGaussian2DFunction f = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, size, size, flags, null);
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
			case FastLogMLELVM:
				MLELVMSteppingFunctionSolver mleSolver = new MLELVMSteppingFunctionSolver(f, tc, bounds);
				solver = mleSolver;
				// MLE requires a positive function value so use a lower bound
				solver.setBounds(getLB(), null);
				// For testing the fast log version
				if (type == FastLogMLELVM)
					mleSolver.setFastLog(FastLogFactory.getFastLog());
				break;
			case WLSELVM:
				solver = new WLSELVMSteppingFunctionSolver(f, tc, bounds);
				break;
			case FastMLE:
				solver = new FastMLESteppingFunctionSolver(f, tc, bounds);
				// MLE requires a positive function value so use a lower bound
				solver.setBounds(getLB(), null);
				break;
			case BTFastMLE:
				solver = new BacktrackingFastMLESteppingFunctionSolver(f, tc, bounds);
				// MLE requires a positive function value so use a lower bound
				solver.setBounds(getLB(), null);
				break;
			case JFastMLE:
				ExtendedFastMLESteppingFunctionSolver efmSolver = new ExtendedFastMLESteppingFunctionSolver(f, tc,
						bounds);
				efmSolver.enableJacobianSolution(true);
				// MLE requires a positive function value so use a lower bound
				solver = efmSolver;
				solver.setBounds(getLB(), null);
				break;
			default:
				throw new NotImplementedException();
		}
		if (solver instanceof LVMSteppingFunctionSolver)
			((LVMSteppingFunctionSolver) solver).setInitialLambda(1);
		return solver;
	}

	double[] getLB()
	{
		double[] lb = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		lb[Gaussian2DFunction.Z_POSITION] = Double.NEGATIVE_INFINITY;
		return lb;
	}

	String getName(boolean bounded, SteppingFunctionSolverClamp clamp, SteppingFunctionSolverType type)
	{
		return String.format("%s %-2s %s", (bounded) ? "B" : "x", clamp, type);
	}
}
