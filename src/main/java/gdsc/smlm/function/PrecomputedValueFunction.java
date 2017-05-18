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
public class PrecomputedValueFunction implements ValueFunction, ValueProcedure
{
	protected final ValueFunction f;
	protected final double[] values;
	protected int i;
	protected ValueProcedure procedure;

	/**
	 * Instantiates a new precomputed value function.
	 *
	 * @param f
	 *            the function
	 * @param values
	 *            the precomputed values
	 * @throws IllegalArgumentException
	 *             if the values length does not match the function size
	 */
	public PrecomputedValueFunction(ValueFunction f, double[] values)
	{
		if (f.size() != values.length)
			throw new IllegalArgumentException("Length of precomputed values must match function size");
		this.f = f;
		this.values = values;
	}

	public int size()
	{
		return f.size();
	}

	public void initialise0(double[] a)
	{
		f.initialise0(a);
		i = 0;
	}

	public void forEach(ValueProcedure procedure)
	{
		this.procedure = procedure;
		f.forEach(this);
	}

	public void execute(double value)
	{
		procedure.execute(value + values[i++]);
	}
}
