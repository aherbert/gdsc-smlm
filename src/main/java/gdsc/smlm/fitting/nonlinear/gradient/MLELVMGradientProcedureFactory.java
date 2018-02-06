package gdsc.smlm.fitting.nonlinear.gradient;

import gdsc.smlm.function.FastLog;
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
public class MLELVMGradientProcedureFactory
{
	/**
	 * Create a new gradient procedure
	 * 
	 * @param y
	 *            Data to fit
	 * @param func
	 *            Gradient function
	 * @return the gradient procedure
	 */
	public static MLELVMGradientProcedure create(final double[] y, final Gradient1Function func)
	{
		if (isStrictlyPositive(y))
		{
			switch (func.getNumberOfGradients())
			{
				case 5:
					return new MLELVMGradientProcedureX5(y, func);
				case 4:
					return new MLELVMGradientProcedureX4(y, func);
				case 6:
					return new MLELVMGradientProcedureX6(y, func);

				default:
					return new MLELVMGradientProcedureX(y, func);
			}
		}

		switch (func.getNumberOfGradients())
		{
			case 5:
				return new MLELVMGradientProcedure5(y, func);
			case 4:
				return new MLELVMGradientProcedure4(y, func);
			case 6:
				return new MLELVMGradientProcedure6(y, func);

			default:
				return new MLELVMGradientProcedure(y, func);
		}
	}

	/**
	 * Create a new gradient procedure using a fast log function that discards insignificant bits from floating-point
	 * values.
	 *
	 * @param y
	 *            Data to fit
	 * @param func
	 *            Gradient function
	 * @param fastLog
	 *            the fast log instance
	 * @return the gradient procedure
	 */
	public static MLELVMGradientProcedure create(final double[] y, final Gradient1Function func, FastLog fastLog)
	{
		if (isStrictlyPositive(y))
		{
			switch (func.getNumberOfGradients())
			{
				case 5:
					return new FastLogMLELVMGradientProcedureX5(y, func, fastLog);
				case 4:
					return new FastLogMLELVMGradientProcedureX4(y, func, fastLog);
				case 6:
					return new FastLogMLELVMGradientProcedureX6(y, func, fastLog);

				default:
					return new FastLogMLELVMGradientProcedureX(y, func, fastLog);
			}
		}

		switch (func.getNumberOfGradients())
		{
			case 5:
				return new FastLogMLELVMGradientProcedure5(y, func, fastLog);
			case 4:
				return new FastLogMLELVMGradientProcedure4(y, func, fastLog);
			case 6:
				return new FastLogMLELVMGradientProcedure6(y, func, fastLog);

			default:
				return new FastLogMLELVMGradientProcedure(y, func, fastLog);
		}
	}

	/**
	 * Checks if is strictly positive (above zero). Ignores NaN or infinity checks.
	 *
	 * @param y
	 *            the y
	 * @return true, if is strictly positive
	 */
	private static boolean isStrictlyPositive(double[] y)
	{
		for (int i = 0; i < y.length; i++)
			if (y[i] <= 0)
				return false;
		return true;
	}
}
