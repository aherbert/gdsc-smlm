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
 * Create a gradient procedure for use in the Levenberg–Marquardt (LVM) algorithm
 */
public class LVMGradientProcedureFactory
{
	/**
	 * Create a new gradient calculator.
	 *
	 * @param y
	 *            Data to fit
	 * @param b
	 *            Baseline pre-computed y-values
	 * @param func
	 *            Gradient function
	 * @param mle
	 *            Set to true to create a Maximum Likelihood Estimator for Poisson data (default is Least Squares)
	 * @return the gradient procedure
	 */
	public static LVMGradientProcedure create(final double[] y, final double[] b, final Gradient1Function func,
			boolean mle)
	{
		return create(y, GradientProcedureHelper.wrapGradient1Function(func, b), mle);

		//if (mle)
		//{
		//	// Use baseline version if appropriate
		//	if (b != null && b.length == y.length)
		//	{
		//		switch (func.getNumberOfGradients())
		//		{
		//			case 5:
		//				return new MLELVMGradientProcedureB5(y, b, func);
		//			case 4:
		//				return new MLELVMGradientProcedureB4(y, b, func);
		//			case 6:
		//				return new MLELVMGradientProcedureB6(y, b, func);
		//			default:
		//				return new MLELVMGradientProcedureB(y, b, func);
		//		}
		//	}
		//	else
		//	{
		//		switch (func.getNumberOfGradients())
		//		{
		//			case 5:
		//				return new MLELVMGradientProcedure5(y, func);
		//			case 4:
		//				return new MLELVMGradientProcedure4(y, func);
		//			case 6:
		//				return new MLELVMGradientProcedure6(y, func);
		//			default:
		//				return new MLELVMGradientProcedure(y, func);
		//		}
		//	}
		//}
		//
		//switch (func.getNumberOfGradients())
		//{
		//	case 5:
		//		return new LSQLVMGradientProcedure5(y, b, func);
		//	case 4:
		//		return new LSQLVMGradientProcedure4(y, b, func);
		//	case 6:
		//		return new LSQLVMGradientProcedure6(y, b, func);
		//
		//	default:
		//		return new LSQLVMGradientProcedure(y, b, func);
		//}
	}

	/**
	 * Create a new gradient calculator.
	 *
	 * @param y
	 *            Data to fit
	 * @param func
	 *            Gradient function
	 * @param mle
	 *            Set to true to create a Maximum Likelihood Estimator for Poisson data (default is Least Squares)
	 * @return the gradient procedure
	 */
	public static LVMGradientProcedure create(final double[] y, final Gradient1Function func, boolean mle)
	{
		if (mle)
		{
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

		switch (func.getNumberOfGradients())
		{
			case 5:
				return new LSQLVMGradientProcedure5(y, func);
			case 4:
				return new LSQLVMGradientProcedure4(y, func);
			case 6:
				return new LSQLVMGradientProcedure6(y, func);

			default:
				return new LSQLVMGradientProcedure(y, func);
		}
	}
}
