package gdsc.smlm.fitting;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Defines methods to fit a function with coefficients (a) for a
 * set of data points (x, y).
 */
public interface FunctionSolver
{
	/**
	 * Fit a function with coefficients (a) for a set of data points (x, y).
	 * <p>
	 * It is assumed that the data points x[i] corresponding to y[i] are consecutive integers from zero.
	 * 
	 * @param n
	 *            The number of points to fit, n <= y.length (allows input data y to be used as a buffer)
	 * @param y
	 *            Set of n data points to fit (input)
	 * @param y_fit
	 *            Fitted data points (output)
	 * @param a
	 *            Set of m coefficients (input/output)
	 * @param a_dev
	 *            Standard deviation of the set of m coefficients (output)
	 * @param error
	 *            Output parameter. The Mean-Squared Error (MSE) for the fit if noise is 0. If noise is provided then
	 *            this will be applied to create a reduced chi-square measure.
	 * @param noise
	 *            Estimate of the noise in the individual measurements
	 * @return The fit status
	 */
	public FitStatus fit(final int n, final float[] y, final float[] y_fit, final float[] a, final float[] a_dev,
			final double[] error, final double noise);

	/**
	 * @return the total Sum Of Squares of the input data points
	 */
	public double getTotalSumOfSquares();

	/**
	 * @return the residual Sum Of Squares
	 */
	public double getResidualSumOfSquares();

	/**
	 * @return the number Of Fitted Parameters
	 */
	public int getNumberOfFittedParameters();

	/**
	 * @return the number Of Fitted Points
	 */
	public int getNumberOfFittedPoints();

	/**
	 * @return The number of iterations used to solve the function
	 */
	public int getIterations();

	/**
	 * @return The number of function evaluations used to solve the function
	 */
	public int getEvaluations();

	/**
	 * Specifies if the function solver supports a bounded search (i.e. a search of parameter space within the total
	 * allowed space of valid parameters, or the parameter constraints). If true then the bounds can be set before a
	 * call to the {@link #fit(int, float[], float[], float[], float[], double[], double)} method.
	 * 
	 * @return True if the function solver supports a bounded search
	 */
	public boolean isBounded();

	/**
	 * Specifies if the function solver supports constraints on the parameters. If true then the constraints can be set
	 * before a call to the {@link #fit(int, float[], float[], float[], float[], double[], double)} method.
	 * <p>
	 * Note that constraints are to be used to specify the values that are absolutely not allowed. They are not meant to
	 * be as restrictive as the bounds for a solver that supports a bounded search. For example the constraints on a
	 * parameter may be 0 - Infinity but the bounds may be 5 - 15. A bounded solver can be used to search within the
	 * expected range for a parameter.
	 * 
	 * @return True if the function solver supports a constrained search
	 */
	public boolean isConstrained();

	/**
	 * Set the bounds for each of the parameters. If a subset of the parameters are fitted then the bounds can be
	 * ignored for the fixed parameters.
	 * <p>
	 * The bounds can be used to set the expected range for a parameter.
	 * 
	 * @param lower
	 * @param upper
	 */
	public void setBounds(float[] lower, float[] upper);

	/**
	 * Set the constraints for each of the parameters. If a subset of the parameters are fitted then the bounds can be
	 * ignored for the fixed parameters.
	 * 
	 * @param lower
	 * @param upper
	 */
	public void setConstraints(float[] lower, float[] upper);
}
