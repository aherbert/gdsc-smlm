package gdsc.smlm.fitting.nonlinear.gradient;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Create a gradient calculator.
 */
public class GradientCalculatorFactory
{
	/**
	 * Create a new gradient calculator
	 * 
	 * @param nparams
	 *            the number of gradient parameters
	 * @return the calculator
	 */
	public static GradientCalculator newCalculator(int nparams)
	{
		return newCalculator(nparams, false);
	}

	/**
	 * Create a new gradient calculator.
	 *
	 * @param nparams
	 *            the number of gradient parameters
	 * @param mle
	 *            true to compute for Maximum Likelihood Estimation
	 * @return the calculator
	 */
	public static GradientCalculator newCalculator(int nparams, boolean mle)
	{
		if (mle)
		{
			switch (nparams)
			{
				case 4:
					// fixed width single Gaussian
					// circular single Gaussian, no background
					return new MLEGradientCalculator4();

				case 5:
					// circular single Gaussian
					// free circular single Gaussian, no background
					return new MLEGradientCalculator5();

				case 6:
					// free circular single Gaussian
					// elliptical single Gaussian, no background
					return new MLEGradientCalculator6();

				case 7:
					// elliptical single Gaussian
					return new MLEGradientCalculator7();

				case 3:
					// fixed width single Gaussian, no background
					return new MLEGradientCalculator3();

				default:
					return new MLEGradientCalculator(nparams);
			}
		}
		else
		{
			switch (nparams)
			{
				case 4:
					// fixed width single Gaussian
					// circular single Gaussian, no background
					return new GradientCalculator4();

				case 5:
					// circular single Gaussian
					// free circular single Gaussian, no background
					return new GradientCalculator5();

				case 6:
					// free circular single Gaussian
					// elliptical single Gaussian, no background
					return new GradientCalculator6();

				case 7:
					// elliptical single Gaussian
					return new GradientCalculator7();

				case 3:
					// fixed width single Gaussian, no background
					return new GradientCalculator3();

				default:
					return new GradientCalculator(nparams);
			}
		}
	}
}
