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
 * Defines function that can produce second order gradients
 */
public interface Gradient2Function extends Gradient1Function
{
	/**
	 * Set the predictor coefficients (a) that will be used to predict each value. Allows the function to perform
	 * initialisation.
	 * 
	 * @param a
	 *            An array of coefficients
	 */
	void initialise2(final double[] a);

	/**
	 * Applies the procedure for the valid range of the function.
	 *
	 * @param procedure
	 *            the procedure
	 */
	public void forEach(Gradient2Procedure procedure);
}
