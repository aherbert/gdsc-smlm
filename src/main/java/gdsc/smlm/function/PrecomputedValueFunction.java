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
public class PrecomputedValueFunction implements ValueFunction
{
	protected final double[] values;

	/**
	 * Instantiates a new pre-computed value function.
	 *
	 * @param values
	 *            the pre-computed values
	 * @throws IllegalArgumentException
	 *             if the values length does not match the function size
	 */
	public PrecomputedValueFunction(double[] values)
	{
		if (values == null)
			throw new IllegalArgumentException("Value is null");
		this.values = values;
	}

	public double[] getValuesRef()
	{
		return values;
	}

	public int size()
	{
		return values.length;
	}

	public void initialise0(double[] a)
	{
		// Ignore
	}

	public void forEach(ValueProcedure procedure)
	{
		for (int i = 0; i < values.length; i++)
			procedure.execute(values[i]);
	}
}
