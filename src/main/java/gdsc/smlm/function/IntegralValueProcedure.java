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
 * Class for evaluating the integral (sum) of a function
 */
public class IntegralValueProcedure implements ValueProcedure
{
	/** The integral (sum) or the values from the last call to {@link #getIntegral(ValueFunction, double[])} */
	public double integral;

	/**
	 * Gets the integral.
	 *
	 * @param f
	 *            the function
	 * @param a
	 *            the function coefficients
	 * @return the integral
	 */
	public double getIntegral(ValueFunction f, double[] a)
	{
		integral = 0;
		f.initialise0(a);
		f.forEach(this);
		return integral;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.ValueProcedure#execute(double)
	 */
	public void execute(double value)
	{
		integral += value;
	}
}