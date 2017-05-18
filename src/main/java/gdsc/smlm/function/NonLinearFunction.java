package gdsc.smlm.function;

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
 * Defines the non-linear fitting function
 */
public interface NonLinearFunction extends GradientFunction
{
	/**
	 * The non-linear fitting function. Produce an output predicted value for a given input
	 * predictor (x) and partial gradient for each of the coefficients (a).
	 * 
	 * @param x
	 *            Predictor
	 * @param dyda
	 *            Partial gradient of function with respect to each coefficient identified by {@link #gradientIndices()}
	 *            . Note: dyda.length must be >= to gradientIndices().length
	 * @return The predicted value y
	 */
	double eval(final int x, final double[] dyda);

	/**
	 * The non-linear fitting function. Produce an output predicted value for a given input
	 * predictor (x).
	 * 
	 * @param x
	 *            Predictor
	 * @return The predicted value y
	 */
	double eval(final int x);

	/**
	 * The non-linear fitting function. Produce an output predicted value for a given input
	 * predictor (x) and partial gradient for each of the coefficients (a).
	 * 
	 * @param x
	 *            Predictor
	 * @param dyda
	 *            Partial gradient of function with respect to each coefficient identified by {@link #gradientIndices()}
	 *            . Note: dyda.length must be >= to gradientIndices().length
	 * @param w
	 *            The output weight. Equivalent to the expected variance of the predicted value. This should not be zero
	 *            to avoid divide by zero error.
	 * @return The predicted value y
	 * @throws NullPointerException
	 *             If the output weight argument is null
	 * @throws ArrayIndexOutOfBoundsException
	 *             If the output weight argument is length 0
	 */
	double eval(final int x, final double[] dyda, final double[] w);

	/**
	 * The non-linear fitting function. Produce an output predicted value for a given input
	 * predictor (x).
	 * 
	 * @param x
	 *            Predictor
	 * @param w
	 *            The output weight. Equivalent to the expected variance of the predicted value. This should not be zero
	 *            to avoid divide by zero error.
	 * @return The predicted value y
	 */
	double evalw(final int x, final double[] w);

	/**
	 * @return True if the {@link #eval(int, double[], double[])} can compute weights other than 1
	 */
	boolean canComputeWeights();
}
