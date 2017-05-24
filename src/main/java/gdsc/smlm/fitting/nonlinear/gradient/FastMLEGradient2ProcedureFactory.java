package gdsc.smlm.fitting.nonlinear.gradient;

import gdsc.smlm.function.Gradient2Function;
import gdsc.smlm.function.PrecomputedGradient2Function;

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
 * Create a Fast MLE gradient procedure
 */
public class FastMLEGradient2ProcedureFactory
{
	/**
	 * Create a new gradient procedure.
	 *
	 * @param x
	 *            Data to fit (must be positive, i.e. the value of a Poisson process)
	 * @param func
	 *            Gradient function
	 * @return the gradient procedure
	 */
	public static FastMLEGradient2Procedure create(final double[] x, final Gradient2Function func)
	{
		return new FastMLEGradient2Procedure(x, func);
		// Note:
		// JUnit speed tests show the unrolled version are slower, i.e. the JVM is able to 
		// efficiently optimise the single for loops in the procedure. So just return the 
		// default implementation.
		//return createUnrolled(x, func);
	}

	/**
	 * Create a new gradient procedure.
	 *
	 * @param x
	 *            Data to fit (must be positive, i.e. the value of a Poisson process)
	 * @param b
	 *            Baseline pre-computed y-values
	 * @param func
	 *            Gradient function
	 * @return the gradient procedure
	 */
	public static FastMLEGradient2Procedure create(final double[] x, final double[] b,
			final Gradient2Function func)
	{
		return create(x, PrecomputedGradient2Function.wrapGradient2Function(func, b));
	}

	/**
	 * Create a new gradient procedure that has the loops unrolled.
	 *
	 * @param func
	 *            Gradient function
	 * @return the gradient procedure
	 */
	static FastMLEGradient2Procedure createUnrolled(final double[] x, final Gradient2Function func)
	{
		switch (func.getNumberOfGradients())
		{
			case 5:
				return new FastMLEGradient2Procedure5(x, func);
			case 4:
				return new FastMLEGradient2Procedure4(x, func);
			case 6:
				return new FastMLEGradient2Procedure6(x, func);
			default:
				return new FastMLEGradient2Procedure(x, func);
		}
	}
}
