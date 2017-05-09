package gdsc.smlm.fitting.nonlinear.gradient;

import gdsc.smlm.function.Gradient1Function;

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
 * Create a Poisson gradient procedure
 */
public class PoissonGradientProcedureFactory
{
	/**
	 * Create a new gradient procedure.
	 *
	 * @param func
	 *            Gradient function
	 * @return the gradient procedure
	 */
	public static PoissonGradientProcedure create(final Gradient1Function func)
	{
		switch (func.getNumberOfGradients())
		{
			case 5:
				return new PoissonGradientProcedure5(func);
			case 4:
				return new PoissonGradientProcedure4(func);
			case 6:
				return new PoissonGradientProcedure6(func);

			default:
				return new PoissonGradientProcedure(func);
		}
	}
}
