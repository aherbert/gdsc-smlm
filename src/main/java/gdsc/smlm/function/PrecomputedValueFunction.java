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
	protected PrecomputedValueFunction(ValueFunction f, double[] values)
	{
		if (f.size() != values.length)
			throw new IllegalArgumentException("Length of precomputed values must match function size");
		this.f = f;
		this.values = values;
	}

	/**
	 * Instantiates a new precomputed value function by combining the current precomputed values with more precomputed
	 * values. This is used internally and so no checks are made on the size of values arrays (which must match).
	 *
	 * @param pre
	 *            the pre-computed function
	 * @param values2
	 *            the second set of precomputed values
	 */
	protected PrecomputedValueFunction(PrecomputedValueFunction pre, double[] values2)
	{
		this.f = pre.f;
		final int n = f.size();
		final double[] values1 = pre.values;
		values = new double[n];
		for (int i = 0; i < n; i++)
			values[i] = values1[i] + values2[i];
	}

	public ValueFunction getValueFunction()
	{
		return f;
	}

	public int size()
	{
		return f.size();
	}

	public void initialise0(double[] a)
	{
		f.initialise0(a);
	}

	public void forEach(ValueProcedure procedure)
	{
		this.procedure = procedure;
		i = 0;
		f.forEach(this);
	}

	public void execute(double value)
	{
		procedure.execute(value + values[i++]);
	}

	/**
	 * Wrap a function with pre-computed values.
	 *
	 * @param func
	 *            the function
	 * @param b
	 *            Baseline pre-computed y-values
	 * @return the wrapped function (or the original if pre-computed values are null or wrong length)
	 */
	public static ValueFunction wrapValueFunction(final ValueFunction func, final double[] b)
	{
		if (b != null && b.length == func.size())
		{
			// Avoid multiple wrapping
			if (func instanceof PrecomputedValueFunction)
			{
				return new PrecomputedValueFunction((PrecomputedValueFunction) func, b);
			}
			return new PrecomputedValueFunction(func, b);
		}
		return func;
	}
}
