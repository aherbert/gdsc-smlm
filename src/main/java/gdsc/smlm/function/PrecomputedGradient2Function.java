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
 * Wraps a value function to add pre-computed values to the forEach procedure
 */
public class PrecomputedGradient2Function extends PrecomputedGradient1Function
		implements Gradient2Function, Gradient2Procedure
{
	protected final Gradient2Function f2;
	protected Gradient2Procedure procedure;

	/**
	 * Instantiates a new precomputed gradient2 function.
	 *
	 * @param f
	 *            the function
	 * @param values
	 *            the precomputed values
	 * @throws IllegalArgumentException
	 *             if the values length does not match the function size
	 */
	public PrecomputedGradient2Function(Gradient2Function f, double[] values)
	{
		super(f, values);
		f2 = f;
	}

	public int[] gradientIndices()
	{
		return f2.gradientIndices();
	}

	public int getNumberOfGradients()
	{
		return f2.getNumberOfGradients();
	}

	public void initialise2(double[] a)
	{
		f2.initialise2(a);
		i = 0;
	}

	public void forEach(Gradient2Procedure procedure)
	{
		this.procedure = procedure;
		f2.forEach((Gradient2Procedure)this);
	}

	public void execute(double value, double[] dy_da, double[] d2y_da2)
	{
		procedure.execute(value + values[i++], dy_da, d2y_da2);
	}
}
