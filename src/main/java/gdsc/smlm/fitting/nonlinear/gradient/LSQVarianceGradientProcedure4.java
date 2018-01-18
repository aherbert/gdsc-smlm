package gdsc.smlm.fitting.nonlinear.gradient;

import java.util.Arrays;

import gdsc.smlm.fitting.linear.EJMLLinearSolver;
import gdsc.smlm.function.Gradient1Function;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Compute the variance of the parameters of the function assuming a least squares fit of a Poisson process.
 * <p>
 * Uses the Mortensen formula (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 25.
 */
public class LSQVarianceGradientProcedure4 extends LSQVarianceGradientProcedure
{
	/**
	 * Instantiates a new LSQ variance gradient procedure.
	 *
	 * @param func
	 *            Gradient function
	 */
	public LSQVarianceGradientProcedure4(final Gradient1Function func)
	{
		super(func);
		if (n != 4)
			throw new IllegalArgumentException("Function must compute 4 gradients");
	}

	/**
	 * Instantiates a new LSQ variance gradient procedure.
	 *
	 * @param func
	 *            Gradient function
	 * @param solver
	 *            the solver
	 * @throws IllegalArgumentException
	 *             if the solver is null
	 */
	public LSQVarianceGradientProcedure4(final Gradient1Function func, EJMLLinearSolver solver)
			throws IllegalArgumentException
	{
		super(func, solver);
		if (n != 4)
			throw new IllegalArgumentException("Function must compute 4 gradients");
	}

	@Override
	protected void initialise()
	{
		I[0] = 0;
		E[0] = 0;
		I[4] = 0;
		E[4] = 0;
		I[5] = 0;
		E[5] = 0;
		I[8] = 0;
		E[8] = 0;
		I[9] = 0;
		E[9] = 0;
		I[10] = 0;
		E[10] = 0;
		I[12] = 0;
		E[12] = 0;
		I[13] = 0;
		E[13] = 0;
		I[14] = 0;
		E[14] = 0;
		I[15] = 0;
		E[15] = 0;
		Arrays.fill(variance, 0);
	}

	@Override
	protected boolean finish()
	{
		if (I[0] != I[0] || I[4] != I[4] || I[5] != I[5] || I[8] != I[8] || I[9] != I[9] || I[10] != I[10] ||
				I[12] != I[12] || I[13] != I[13] || I[14] != I[14] || I[15] != I[15])
			return true;

		I[1] = I[4];
		E[1] = E[4];
		I[2] = I[8];
		E[2] = E[8];
		I[6] = I[9];
		E[6] = E[9];
		I[3] = I[12];
		E[3] = E[12];
		I[7] = I[13];
		E[7] = E[13];
		I[11] = I[14];
		E[11] = E[14];
		return false;
	}

	@Override
	protected void computeVariance()
	{
		variance[0] = I[0] * E[0] * I[0] + I[0] * E[1] * I[4] + I[0] * E[2] * I[8] + I[0] * E[3] * I[12] +
				I[1] * E[4] * I[0] + I[1] * E[5] * I[4] + I[1] * E[6] * I[8] + I[1] * E[7] * I[12] +
				I[2] * E[8] * I[0] + I[2] * E[9] * I[4] + I[2] * E[10] * I[8] + I[2] * E[11] * I[12] +
				I[3] * E[12] * I[0] + I[3] * E[13] * I[4] + I[3] * E[14] * I[8] + I[3] * E[15] * I[12];
		variance[1] = I[4] * E[0] * I[1] + I[4] * E[1] * I[5] + I[4] * E[2] * I[9] + I[4] * E[3] * I[13] +
				I[5] * E[4] * I[1] + I[5] * E[5] * I[5] + I[5] * E[6] * I[9] + I[5] * E[7] * I[13] +
				I[6] * E[8] * I[1] + I[6] * E[9] * I[5] + I[6] * E[10] * I[9] + I[6] * E[11] * I[13] +
				I[7] * E[12] * I[1] + I[7] * E[13] * I[5] + I[7] * E[14] * I[9] + I[7] * E[15] * I[13];
		variance[2] = I[8] * E[0] * I[2] + I[8] * E[1] * I[6] + I[8] * E[2] * I[10] + I[8] * E[3] * I[14] +
				I[9] * E[4] * I[2] + I[9] * E[5] * I[6] + I[9] * E[6] * I[10] + I[9] * E[7] * I[14] +
				I[10] * E[8] * I[2] + I[10] * E[9] * I[6] + I[10] * E[10] * I[10] + I[10] * E[11] * I[14] +
				I[11] * E[12] * I[2] + I[11] * E[13] * I[6] + I[11] * E[14] * I[10] + I[11] * E[15] * I[14];
		variance[3] = I[12] * E[0] * I[3] + I[12] * E[1] * I[7] + I[12] * E[2] * I[11] + I[12] * E[3] * I[15] +
				I[13] * E[4] * I[3] + I[13] * E[5] * I[7] + I[13] * E[6] * I[11] + I[13] * E[7] * I[15] +
				I[14] * E[8] * I[3] + I[14] * E[9] * I[7] + I[14] * E[10] * I[11] + I[14] * E[11] * I[15] +
				I[15] * E[12] * I[3] + I[15] * E[13] * I[7] + I[15] * E[14] * I[11] + I[15] * E[15] * I[15];
	}
}
