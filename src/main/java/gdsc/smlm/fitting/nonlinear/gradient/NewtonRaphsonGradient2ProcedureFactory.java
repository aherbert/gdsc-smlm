package gdsc.smlm.fitting.nonlinear.gradient;

import gdsc.smlm.function.Gradient2Function;

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
 * Create a Newton-Raphson gradient procedure
 */
public class NewtonRaphsonGradient2ProcedureFactory
{
	/**
	 * Create a new gradient procedure.
	 *
	 * @param func
	 *            Gradient function
	 * @return the gradient procedure
	 */
	public static NewtonRaphsonGradient2Procedure create(final double[] x, final Gradient2Function func)
	{
		switch (func.getNumberOfGradients())
		{
			case 5:
				return new NewtonRaphsonGradient2Procedure5(x, func);
			case 4:
				return new NewtonRaphsonGradient2Procedure4(x, func);
			case 6:
				return new NewtonRaphsonGradient2Procedure6(x, func);
			default:
				return new NewtonRaphsonGradient2Procedure(x, func);
		}
	}
}
