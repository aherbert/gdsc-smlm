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
 * Implement a fixed value non-linear fitting function
 */
public class FixedNonLinearFunction implements NonLinearFunction
{
	final double[] values;

	public FixedNonLinearFunction(double[] values)
	{
		this.values = values;
	}

	public void initialise(double[] a)
	{

	}

	public int[] gradientIndices()
	{
		return new int[0];
	}

	public int getNumberOfGradients()
	{
		return 0;
	}

	public double eval(int x, double[] dyda)
	{
		return values[x];
	}

	public double eval(int x)
	{
		return values[x];
	}

	public double eval(int x, double[] dyda, double[] w)
	{
		return values[x];
	}

	public double evalw(int x, double[] w)
	{
		return values[x];
	}

	public boolean canComputeWeights()
	{
		return false;
	}
}
