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
public class OffsetFunctionFactory
{
	/**
	 * Wrap a function with a pre-computed offset.
	 *
	 * @param func
	 *            the function
	 * @param b
	 *            Baseline pre-computed offset for the values
	 * @return the wrapped function (or the original if pre-computed values are null or wrong length)
	 */
	public static ValueFunction wrapFunction(final ValueFunction func, final double[] b)
	{
		if (b != null && b.length == func.size())
		{
			// Wrap appropriately
			if (func instanceof ExtendedGradient2Function)
			{
				return OffsetExtendedGradient2Function
						.wrapExtendedGradient2Function((ExtendedGradient2Function) func, b);
			}
			if (func instanceof Gradient2Function)
			{
				return OffsetGradient2Function.wrapGradient2Function((Gradient2Function) func, b);
			}
			if (func instanceof Gradient1Function)
			{
				return OffsetGradient1Function.wrapGradient1Function((Gradient1Function) func, b);
			}
			return OffsetValueFunction.wrapValueFunction(func, b);
		}
		return func;
	}
}
