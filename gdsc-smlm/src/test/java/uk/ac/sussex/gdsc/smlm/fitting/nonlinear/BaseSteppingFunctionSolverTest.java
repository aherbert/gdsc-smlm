/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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

import uk.ac.sussex.gdsc.core.data.NotImplementedException;
import uk.ac.sussex.gdsc.smlm.function.FastLogFactory;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;

/**
 * Test that a stepping solver can return the same results with and without bounds.
 */
@SuppressWarnings("javadoc")
public abstract class BaseSteppingFunctionSolverTest extends BaseFunctionSolverTest {
  enum SteppingFunctionSolverClamp {
    NO_CLAMP("x"), CLAMP("C"), DYNAMIC_CLAMP("DC");

    final String name;

    SteppingFunctionSolverClamp(String name) {
      this.name = name;
    }

    @Override
    public String toString() {
      return name;
    }
  }

  enum SteppingFunctionSolverType {
    // Enum names should all be uppercase but this is just for a test so ignore that convention
    MLELVM, FastLogMLELVM, LSELVM, WLSELVM, FastMLE
  }

  // For convenience declare variables of the enum type
  static final SteppingFunctionSolverClamp NO_CLAMP = SteppingFunctionSolverClamp.NO_CLAMP;
  static final SteppingFunctionSolverClamp CLAMP = SteppingFunctionSolverClamp.CLAMP;
  static final SteppingFunctionSolverClamp DYNAMIC_CLAMP =
      SteppingFunctionSolverClamp.DYNAMIC_CLAMP;
  static final SteppingFunctionSolverType MLELVM = SteppingFunctionSolverType.MLELVM;
  static final SteppingFunctionSolverType FastLogMLELVM = SteppingFunctionSolverType.FastLogMLELVM;
  static final SteppingFunctionSolverType LSELVM = SteppingFunctionSolverType.LSELVM;
  static final SteppingFunctionSolverType WLSELVM = SteppingFunctionSolverType.WLSELVM;
  static final SteppingFunctionSolverType FastMLE = SteppingFunctionSolverType.FastMLE;
  static final boolean BOUNDED = true;
  static final boolean NO_BOUND = false;

  static class NoToleranceChecker extends ToleranceChecker {
    public NoToleranceChecker() {
      super(0, 0);
    }

    @Override
    public int converged(double previousValue, double[] previousParameters, double currentValue,
        double[] currentParameters) {
      return STATUS_VALUE;
    }
  }

  static NoToleranceChecker noToleranceChecker = new NoToleranceChecker();

  SteppingFunctionSolver getSolver(SteppingFunctionSolverClamp clamp,
      SteppingFunctionSolverType type) {
    return getSolver(clamp, type, new ToleranceChecker(1e-5, 1e-5, 0, 0, 100));
  }

  SteppingFunctionSolver getSolver(SteppingFunctionSolverClamp clamp,
      SteppingFunctionSolverType type, ToleranceChecker tc) {
    final ErfGaussian2DFunction f =
        (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, size, size, flags, null);
    final ParameterBounds bounds = ParameterBounds.create(f);
    switch (clamp) {
      case DYNAMIC_CLAMP:
        bounds.setDynamicClamp(true);
        bounds.setClampValues(defaultClampValues);
        break;
      case CLAMP:
        bounds.setClampValues(defaultClampValues);
        break;
      case NO_CLAMP:
      default:
        break;
    }
    SteppingFunctionSolver solver;
    switch (type) {
      case LSELVM:
        solver = new LseLvmSteppingFunctionSolver(f, tc, bounds);
        break;
      case MLELVM:
      case FastLogMLELVM:
        final MleLvmSteppingFunctionSolver mleSolver =
            new MleLvmSteppingFunctionSolver(f, tc, bounds);
        solver = mleSolver;
        // MLE requires a positive function value so use a lower bound
        solver.setBounds(getLb(), null);
        // For testing the fast log version
        if (type == FastLogMLELVM) {
          mleSolver.setFastLog(FastLogFactory.getFastLog());
        }
        break;
      case WLSELVM:
        solver = new WLseLvmSteppingFunctionSolver(f, tc, bounds);
        break;
      case FastMLE:
        solver = new FastMleSteppingFunctionSolver(f, tc, bounds);
        // MLE requires a positive function value so use a lower bound
        solver.setBounds(getLb(), null);
        break;
      default:
        throw new NotImplementedException();
    }
    if (solver instanceof LvmSteppingFunctionSolver) {
      ((LvmSteppingFunctionSolver) solver).setInitialLambda(1);
    }
    return solver;
  }

  double[] getLb() {
    final double[] lb = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    lb[Gaussian2DFunction.Z_POSITION] = Double.NEGATIVE_INFINITY;
    return lb;
  }

  String getName(boolean bounded, SteppingFunctionSolverClamp clamp,
      SteppingFunctionSolverType type) {
    return String.format("%s %-2s %s", (bounded) ? "B" : "x", clamp, type);
  }
}
