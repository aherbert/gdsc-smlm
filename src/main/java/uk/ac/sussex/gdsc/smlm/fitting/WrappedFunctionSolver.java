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
package uk.ac.sussex.gdsc.smlm.fitting;

/**
 * Wrap a function solver. The default implementation passes all calls to the FunctioSolver
 * interface to the inner function solver.
 */
public class WrappedFunctionSolver implements FunctionSolver {
  /** The solver. */
  protected final FunctionSolver solver;

  /**
   * Instantiates a new wrapped function solver.
   *
   * @param solver the solver
   */
  public WrappedFunctionSolver(FunctionSolver solver) {
    if (solver == null) {
      throw new NullPointerException("FunctionSolver is null");
    }
    this.solver = solver;
  }

  // Pass through all interface calls to the inner function solver.

  @Override
  public FunctionSolverType getType() {
    return solver.getType();
  }

  @Override
  public FitStatus fit(double[] y, double[] f, double[] a, double[] aDev) {
    return solver.fit(y, f, a, aDev);
  }

  @Override
  public int getNumberOfFittedParameters() {
    return solver.getNumberOfFittedParameters();
  }

  @Override
  public int getNumberOfFittedPoints() {
    return solver.getNumberOfFittedPoints();
  }

  @Override
  public int getIterations() {
    return solver.getIterations();
  }

  @Override
  public int getEvaluations() {
    return solver.getEvaluations();
  }

  @Override
  public boolean isBounded() {
    return solver.isBounded();
  }

  @Override
  public boolean isConstrained() {
    return solver.isConstrained();
  }

  @Override
  public boolean isWeighted() {
    return solver.isWeighted();
  }

  @Override
  public boolean isStrictlyPositiveFunction() {
    return solver.isStrictlyPositiveFunction();
  }

  @Override
  public void setBounds(double[] lower, double[] upper) {
    solver.setBounds(lower, upper);
  }

  @Override
  public void setConstraints(double[] lower, double[] upper) {
    solver.setConstraints(lower, upper);
  }

  @Override
  public void setWeights(double[] weights) {
    solver.setWeights(weights);
  }

  @Override
  public double getValue() {
    return solver.getValue();
  }

  @Override
  public boolean evaluate(double[] y, double[] f, double[] a) {
    return solver.evaluate(y, f, a);
  }

  @Override
  public boolean computeDeviations(double[] y, double[] a, double[] aDev) {
    return solver.computeDeviations(y, a, aDev);
  }

  @Override
  public String getName(int i) {
    return solver.getName(i);
  }
}
