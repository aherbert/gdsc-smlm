package gdsc.smlm.fitting.function;

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
public interface NonLinearFunction
{
	/**
	 * Set the predictor coefficients (a) that will be used to predict each value. Allows the function to perform
	 * initialisation.
	 * 
	 * @param a
	 *            An array of coefficients
	 */
	void initialise(final float[] a);

	/**
	 * The function will evaluate the gradient for up to n parameters where n <= a.length. This method
	 * returns the indices that are evaluated.
	 * 
	 * @return The gradient indices
	 */
	int[] gradientIndices();
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
	float eval(final int x, final float[] dyda);

	/**
	 * The non-linear fitting function. Produce an output predicted value for a given input
	 * predictor (x).
	 * 
	 * @param x
	 *            Predictor
	 * @return The predicted value y
	 */
	float eval(final int x);

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
	 *            The output weight. Equivalent to the expected variance of the predicted value
	 * @return The predicted value y
	 * @throws NullPointerException If the output weight argument is null
	 * @throws ArrayIndexOutOfBoundsException If the output weight argument is length 0
	 */
	float eval(final int x, final float[] dyda, final float[] w);
	
	/**
	 * @return True if the {@link #eval(int, float[], float[])} can compute weights other than 1
	 */
	boolean canComputeWeights();
}
