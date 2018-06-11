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
package gdsc.smlm.fitting;

/**
 * Wrap a function solver to scale the function value. Parameters that are scaled must be provided in the constructor.
 * It is assumed that a linear scale can be applied to all these parameters with the effect that the output function
 * value is reduced by the scale factor. This will be the case for offset parameters (if the offset is from zero) and
 * magnitude parameters. An example is y=m*x+c can be scaled to a*y=a*m*x+a*c. Interface methods that compute the
 * function value are rescaled after computation.
 */
public class MLEScaledFunctionSolver extends ScaledFunctionSolver implements MLEFunctionSolver
{
	protected final MLEFunctionSolver mleSolver;

	/**
	 * Instantiates a new MLE scaled function solver.
	 * <p>
	 * Indexed parameters are up-scaled prior to calling the inner function solver. Output parameters, deviations and
	 * the function value are are down-scaled upon completion.
	 *
	 * @param solver
	 *            the solver
	 * @param scale
	 *            the scale
	 * @param indices
	 *            the indices of the parameters to scale
	 */
	public MLEScaledFunctionSolver(MLEFunctionSolver solver, double scale, int[] indices)
	{
		super(solver, scale, indices);
		mleSolver = solver;
	}

	@Override
	public double getLogLikelihood()
	{
		return mleSolver.getLogLikelihood();
	}

	@Override
	public double getLogLikelihoodRatio()
	{
		return mleSolver.getLogLikelihoodRatio();
	}

	@Override
	public double getQ()
	{
		return mleSolver.getQ();
	}
}
