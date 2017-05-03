package gdsc.smlm.function;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Defines function that can produce gradients
 */
public interface GradientFunction
{
	/**
	 * Returns the size of the valid range of the function. Procedures passed to the forEach methods will be expected to
	 * be called this number of times.
	 *
	 * @return the size
	 */
	int size();

	/**
	 * Set the predictor coefficients (a) that will be used to predict each value. Allows the function to perform
	 * initialisation.
	 * 
	 * @param a
	 *            An array of coefficients
	 */
	void initialise(final double[] a);

	/**
	 * The function will evaluate the gradient for up to n parameters where n <= a.length. This method
	 * returns the indices that are evaluated.
	 * 
	 * @return The gradient indices
	 */
	int[] gradientIndices();
	
	/**
	 * Gets the number of gradients. The function will evaluate this many partial derivatives.
	 *
	 * @return the number of gradients
	 */
	int getNumberOfGradients();

	/**
	 * Applies the procedure for the valid range of the function.
	 *
	 * @param procedure
	 *            the procedure
	 */
	public void forEach(ValueProcedure procedure);

	/**
	 * Applies the procedure for the valid range of the function.
	 *
	 * @param procedure
	 *            the procedure
	 */
	public void forEach(Gradient1Procedure procedure);
}
