package gdsc.smlm.fitting.nonlinear.gradient;

import gdsc.smlm.function.DFastLog;
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
 * Create a gradient procedure for use in the Levenberg–Marquardt (LVM) algorithm
 */
public class LVMGradientProcedureFactory
{
	public enum Type
	{
		LSQ, MLE, WLSQ, FastLogMLE
	}

	/**
	 * Create a new gradient calculator.
	 *
	 * @param y
	 *            Data to fit
	 * @param func
	 *            Gradient function
	 * @param type
	 *            the type
	 * @param fastLog
	 *            the fast log
	 * @return the gradient procedure
	 */
	public static LVMGradientProcedure create(final double[] y, final Gradient1Function func, Type type,
			DFastLog fastLog)
	{
		switch (type)
		{
			case WLSQ:
				// Do not support per observation weights
				return WLSQLVMGradientProcedureFactory.create(y, null, func);
			case MLE:
				return MLELVMGradientProcedureFactory.create(y, func);
			case LSQ:
				return LSQLVMGradientProcedureFactory.create(y, func);
			case FastLogMLE:
				return MLELVMGradientProcedureFactory.create(y, func, fastLog);
			default:
				throw new IllegalStateException("Unknown type: " + type);
		}
	}
}
