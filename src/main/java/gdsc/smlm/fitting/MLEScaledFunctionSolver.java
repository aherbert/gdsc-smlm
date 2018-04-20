package gdsc.smlm.fitting;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

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

	public double getLogLikelihood()
	{
		return mleSolver.getLogLikelihood();
	}

	public double getLogLikelihoodRatio()
	{
		return mleSolver.getLogLikelihoodRatio();
	}

	public double getQ()
	{
		return mleSolver.getQ();
	}
}
