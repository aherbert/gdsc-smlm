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
 * Class for evaluating a function
 */
public class StandardGradient1Procedure implements Gradient1Procedure
{
	private int i;

	/** The values from the last call to {@link #getValues(Gradient1Function, double[])}. */
	public double[] values;
	/** The gradients from the last call to {@link #getValues(Gradient1Function, double[])}. */
	public double[][] dyda;

	/**
	 * Gets the values.
	 *
	 * @param f
	 *            the function
	 * @param a
	 *            the function coefficients
	 * @return the values
	 */
	public double[] getValues(Gradient1Function f, double[] a)
	{
		values = new double[f.size()];
		dyda = new double[values.length][];
		i = 0;
		f.initialise1(a);
		f.forEach(this);
		return values;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Procedure#execute(double, double[])
	 */
	public void execute(double value, double[] dy_da)
	{
		values[i] = value;
		dyda[i] = dy_da.clone();
		i++;
	}
}