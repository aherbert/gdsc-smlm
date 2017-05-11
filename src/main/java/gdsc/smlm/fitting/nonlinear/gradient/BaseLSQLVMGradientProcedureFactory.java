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
 * Create a gradient procedure.
 */
abstract class BaseLSQLVMGradientProcedureFactory
{
	// Instance methods for testing
	BaseLSQLVMGradientProcedure createProcedure(final double[] y, final Gradient1Function func)
	{
		return createProcedure(y, null, func);
	}

	abstract BaseLSQLVMGradientProcedure createProcedure(final double[] y, final double[] b,
			final Gradient1Function func);
}
