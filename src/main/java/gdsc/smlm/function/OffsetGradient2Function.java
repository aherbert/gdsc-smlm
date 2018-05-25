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
 * Wraps a value function to add a pre-computed offset to the value during the forEach procedure
 */
public class OffsetGradient2Function extends OffsetGradient1Function
		implements Gradient2Function, Gradient2Procedure
{
	protected final Gradient2Function f2;
	protected Gradient2Procedure procedure;

	/**
	 * Instantiates a new precomputed gradient2 function.
	 *
	 * @param f
	 *            the function
	 * @param values
	 *            the precomputed values
	 * @throws IllegalArgumentException
	 *             if the values length does not match the function size
	 */
	protected OffsetGradient2Function(Gradient2Function f, double[] values)
	{
		super(f, values);
		f2 = f;
	}
	
	protected OffsetGradient2Function(OffsetGradient2Function pre, double[] values2)
	{
		super(pre, values2);
		f2 = (Gradient2Function) f;
	}
	
	public Gradient2Function getGradient2Function()
	{
		return f2;
	}

	public void initialise2(double[] a)
	{
		f2.initialise2(a);
	}

	public void forEach(Gradient2Procedure procedure)
	{
		this.procedure = procedure;
		i = 0;
		f2.forEach((Gradient2Procedure) this);
	}

	public void execute(double value, double[] dy_da, double[] d2y_da2)
	{
		procedure.execute(value + offset[i++], dy_da, d2y_da2);
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
	public static Gradient2Function wrapGradient2Function(final Gradient2Function func, final double[] b)
	{
		if (b != null && b.length == func.size())
		{
			// Avoid multiple wrapping
			if (func instanceof OffsetGradient2Function)
			{
				return new OffsetGradient2Function((OffsetGradient2Function)func, b);
			}
			return new OffsetGradient2Function(func, b);
		}
		return func;
	}
}
