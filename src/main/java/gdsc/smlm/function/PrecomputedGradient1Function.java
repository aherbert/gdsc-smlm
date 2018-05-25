package gdsc.smlm.function;

import gdsc.core.utils.SimpleArrayUtils;

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
public class PrecomputedGradient1Function extends PrecomputedValueFunction implements Gradient1Function
{
	protected final int[] gradientIndices;
	protected final double[][] g1;

	/**
	 * Instantiates a new pre-computed value function.
	 *
	 * @param values
	 *            the pre-computed values
	 * @param g1
	 *            the first order gradient
	 * @throws IllegalArgumentException
	 *             if the values length does not match the function size
	 */
	public PrecomputedGradient1Function(double[] values, double[][] g1)
	{
		super(values);
		int numberOfGradients = checkGradient(g1);
		gradientIndices = SimpleArrayUtils.newArray(numberOfGradients, 0, 1);
		this.g1 = g1;
	}

	protected int checkGradient(double[][] g)
	{
		if (g == null)
			throw new IllegalArgumentException("Gradient is null");
		if (g.length != values.length)
			throw new IllegalArgumentException("Gradient is not same size as values");
		if (g.length == 0)
			return 0;
		if (g[0] == null)
			throw new IllegalArgumentException("Gradient[0][] is null");
		int n = g[0].length;
		for (int i = 1; i < g.length; i++)
			if (g[i] == null || g[i].length != n)
				throw new IllegalArgumentException("Gradient[" + i + "][] is incorrect size");
		return n;
	}

	public double[][] getGradient1Ref()
	{
		return g1;
	}

	public void initialise(double[] a)
	{
		// Ignore
	}

	public int[] gradientIndices()
	{
		return gradientIndices;
	}

	public int getNumberOfGradients()
	{
		return gradientIndices.length;
	}

	public void initialise1(double[] a)
	{
		// Ignore
	}

	public void forEach(Gradient1Procedure procedure)
	{
		for (int i = 0; i < values.length; i++)
			procedure.execute(values[i], g1[i]);
	}
}
