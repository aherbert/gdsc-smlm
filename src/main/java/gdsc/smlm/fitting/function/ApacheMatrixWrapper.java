package gdsc.smlm.fitting.function;

import org.apache.commons.math3.analysis.MultivariateMatrixFunction;

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

public class ApacheMatrixWrapper extends ApacheFunctionWrapper implements MultivariateMatrixFunction
{
	public ApacheMatrixWrapper(NonLinearFunction fun, double[] a, int n)
	{
		super(fun, a, n);
	}

	@Override
	public double[][] value(double[] point) throws IllegalArgumentException
	{
		initialiseFunction(point);
		
		double[][] jacobian = new double[n][point.length];
		double[] dyda = new double[point.length];

		for (int i = 0; i < jacobian.length; ++i)
		{
			//float y = gf.eval(x.get(i).intValue());
			// Assume linear X from 0..N
			fun.eval(i, dyda);

			// Differentiate with respect to each parameter:
			for (int j = 0; j < dyda.length; j++)
				jacobian[i][j] = dyda[j];
		}

		return jacobian;
	}
}