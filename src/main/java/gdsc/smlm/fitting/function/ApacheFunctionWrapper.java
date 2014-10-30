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
 * Wrap the NonLinearFunction with classes that implement the required Apache Commons Math interfaces
 */
public abstract class ApacheFunctionWrapper
{
	NonLinearFunction fun;
	double[] a;
	int n;

	public ApacheFunctionWrapper(NonLinearFunction fun, double[] a, int n)
	{
		this.fun = fun;
		this.a = a.clone();
		this.n = n;
	}

	void initialiseFunction(double[] variables)
	{
		int[] gradientIndices = fun.gradientIndices();
		for (int i = 0; i < gradientIndices.length; i++)
			a[gradientIndices[i]] = variables[i];
		fun.initialise(a);
	}
}