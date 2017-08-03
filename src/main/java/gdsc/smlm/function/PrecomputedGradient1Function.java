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
public class PrecomputedGradient1Function extends PrecomputedValueFunction
		implements Gradient1Function, Gradient1Procedure
{
	protected final Gradient1Function f1;
	protected Gradient1Procedure procedure;

	/**
	 * Instantiates a new precomputed gradient1 function.
	 *
	 * @param f
	 *            the function
	 * @param values
	 *            the precomputed values
	 * @throws IllegalArgumentException
	 *             if the values length does not match the function size
	 */
	protected PrecomputedGradient1Function(Gradient1Function f, double[] values)
	{
		super(f, values);
		f1 = f;
	}

	protected PrecomputedGradient1Function(PrecomputedGradient1Function pre, double[] values2)
	{
		super(pre, values2);
		f1 = (Gradient1Function) f;
	}
	
	public Gradient1Function getGradient1Function()
	{
		return f1;
	}
	
	public void initialise(double[] a)
	{
		f1.initialise(a);
		i = 0;
	}

	public int[] gradientIndices()
	{
		return f1.gradientIndices();
	}

	public int getNumberOfGradients()
	{
		return f1.getNumberOfGradients();
	}

	public void initialise1(double[] a)
	{
		f1.initialise1(a);
	}

	public void forEach(Gradient1Procedure procedure)
	{
		this.procedure = procedure;
		i = 0;
		f1.forEach((Gradient1Procedure) this);
	}

	public void execute(double value, double[] dy_da)
	{
		procedure.execute(value + values[i++], dy_da);
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
	public static Gradient1Function wrapGradient1Function(final Gradient1Function func, final double[] b)
	{
		if (b != null && b.length == func.size())
		{
			// Avoid multiple wrapping
			if (func instanceof PrecomputedGradient1Function)
			{
				return new PrecomputedGradient1Function((PrecomputedGradient1Function)func, b);
			}
			return new PrecomputedGradient1Function(func, b);
		}
		return func;
	}
}
