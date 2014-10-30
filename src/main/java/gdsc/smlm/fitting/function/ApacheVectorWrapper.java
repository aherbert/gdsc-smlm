package gdsc.smlm.fitting.function;

import org.apache.commons.math3.analysis.MultivariateVectorFunction;

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

public class ApacheVectorWrapper extends ApacheFunctionWrapper implements MultivariateVectorFunction
{
	public ApacheVectorWrapper(NonLinearFunction fun, double[] a, int n)
	{
		super(fun, a, n);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.apache.commons.math3.analysis.MultivariateFunction#value(double[])
	 */
	@Override
	public double[] value(double[] point)
	{
		initialiseFunction(point);
		double[] values = new double[n];
		for (int i = 0; i < values.length; i++)
		{
			// Assume linear X from 0..N
			values[i] = fun.eval(i);
		}
		return values;
	}
}