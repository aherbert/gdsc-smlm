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
public class PrecomputedExtendedGradient2Function extends PrecomputedGradient2Function
		implements ExtendedGradient2Function, ExtendedGradient2Procedure
{
	protected final ExtendedGradient2Function ef2;
	protected ExtendedGradient2Procedure procedure;

	/**
	 * Instantiates a new precomputed extended gradient2 function.
	 *
	 * @param f
	 *            the function
	 * @param values
	 *            the precomputed values
	 * @throws IllegalArgumentException
	 *             if the values length does not match the function size
	 */
	protected PrecomputedExtendedGradient2Function(ExtendedGradient2Function f, double[] values)
	{
		super(f, values);
		ef2 = f;
	}

	protected PrecomputedExtendedGradient2Function(PrecomputedExtendedGradient2Function pre, double[] values2)
	{
		super(pre, values2);
		ef2 = (ExtendedGradient2Function) f;
	}

	public ExtendedGradient2Function getExtendedGradient2Function()
	{
		return ef2;
	}

	public void initialiseExtended2(double[] a)
	{
		ef2.initialiseExtended2(a);
	}

	public void forEach(ExtendedGradient2Procedure procedure)
	{
		this.procedure = procedure;
		i = 0;
		ef2.forEach((ExtendedGradient2Procedure) this);
	}

	public void executeExtended(double value, double[] dy_da, double[] d2y_dadb)
	{
		procedure.executeExtended(value + values[i++], dy_da, d2y_dadb);
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
	public static ExtendedGradient2Function wrapExtendedGradient2Function(final ExtendedGradient2Function func,
			final double[] b)
	{
		if (b != null && b.length == func.size())
		{
			// Avoid multiple wrapping
			if (func instanceof PrecomputedExtendedGradient2Function)
			{
				return new PrecomputedExtendedGradient2Function((PrecomputedExtendedGradient2Function) func, b);
			}
			return new PrecomputedExtendedGradient2Function(func, b);
		}
		return func;
	}
}
