package gdsc.smlm.function;

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
 * Wraps a set of function values to implement the forEach procedure
 */
public class PrecomputedGradient2Function extends PrecomputedGradient1Function implements Gradient2Function
{
	protected final double[][] g2;

	/**
	 * Instantiates a new pre-computed value function.
	 *
	 * @param values
	 *            the pre-computed values
	 * @param g1
	 *            the first order gradient
	 * @param g2
	 *            the second order gradient
	 * @throws IllegalArgumentException
	 *             if the values length does not match the function size
	 */
	public PrecomputedGradient2Function(double[] values, double[][] g1, double[][] g2)
	{
		super(values, g1);
		checkGradient(g2);
		this.g2 = g2;
	}

	public double[][] getGradient2Ref()
	{
		return g2;
	}

	public void initialise2(double[] a)
	{
		// Ignore
	}

	public void forEach(Gradient2Procedure procedure)
	{
		for (int i = 0; i < values.length; i++)
			procedure.execute(values[i], g1[i], g2[i]);
	}
}
