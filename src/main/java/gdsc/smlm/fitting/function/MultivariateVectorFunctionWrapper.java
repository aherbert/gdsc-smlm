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

/**
 * Wrap the NonLinearFunction to allow use with the Apache Commons Math library
 */
public class MultivariateVectorFunctionWrapper extends NonLinearFunctionWrapper implements MultivariateVectorFunction
{
	public MultivariateVectorFunctionWrapper(NonLinearFunction fun, double[] a, int n)
	{
		super(fun, a, n);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.apache.commons.math3.analysis.MultivariateFunction#value(double[])
	 */
	public double[] value(double[] point)
	{
		return computeValue(point);
	}
}