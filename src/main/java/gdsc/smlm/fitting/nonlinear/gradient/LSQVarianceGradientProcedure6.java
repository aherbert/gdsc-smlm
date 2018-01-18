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
public class LSQVarianceGradientProcedure6 extends LSQVarianceGradientProcedure
{
	/**
	 * Instantiates a new LSQ variance gradient procedure.
	 *
	 * @param func
	 *            Gradient function
	 */
	public LSQVarianceGradientProcedure6(final Gradient1Function func)
	{
		super(func);
		if (n != 6)
			throw new IllegalArgumentException("Function must compute 6 gradients");
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
	public LSQVarianceGradientProcedure6(final Gradient1Function func, EJMLLinearSolver solver)
			throws IllegalArgumentException
	{
		super(func, solver);
		if (n != 6)
			throw new IllegalArgumentException("Function must compute 6 gradients");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Procedure#execute(double, double[])
	 */
	public void execute(final double Ei, double[] Eix)
	{
		for (int a = 0; a < n; a++)
		{
			for (int b = 0, j = a * n; b <= a; b++, j++)
			{
				double v = Eix[a] * Eix[b];
				I[j] += v;
				E[j] += Ei * v;
			}
		}
	}

	@Override
	protected void initialise()
	{
		I[0] = 0;
		E[0] = 0;
		I[6] = 0;
		E[6] = 0;
		I[7] = 0;
		E[7] = 0;
		I[12] = 0;
		E[12] = 0;
		I[13] = 0;
		E[13] = 0;
		I[14] = 0;
		E[14] = 0;
		I[18] = 0;
		E[18] = 0;
		I[19] = 0;
		E[19] = 0;
		I[20] = 0;
		E[20] = 0;
		I[21] = 0;
		E[21] = 0;
		I[24] = 0;
		E[24] = 0;
		I[25] = 0;
		E[25] = 0;
		I[26] = 0;
		E[26] = 0;
		I[27] = 0;
		E[27] = 0;
		I[28] = 0;
		E[28] = 0;
		I[30] = 0;
		E[30] = 0;
		I[31] = 0;
		E[31] = 0;
		I[32] = 0;
		E[32] = 0;
		I[33] = 0;
		E[33] = 0;
		I[34] = 0;
		E[34] = 0;
		I[35] = 0;
		E[35] = 0;
		Arrays.fill(variance, 0);
	}

	@Override
	protected boolean finish()
	{
		if (I[0] != I[0] || I[6] != I[6] || I[7] != I[7] || I[12] != I[12] || I[13] != I[13] || I[14] != I[14] ||
				I[18] != I[18] || I[19] != I[19] || I[20] != I[20] || I[21] != I[21] || I[24] != I[24] ||
				I[25] != I[25] || I[26] != I[26] || I[27] != I[27] || I[28] != I[28] || I[30] != I[30] ||
				I[31] != I[31] || I[32] != I[32] || I[33] != I[33] || I[34] != I[34] || I[35] != I[35])
			return true;

		I[1] = I[6];
		E[1] = E[6];
		I[2] = I[12];
		E[2] = E[12];
		I[8] = I[13];
		E[8] = E[13];
		I[3] = I[18];
		E[3] = E[18];
		I[9] = I[19];
		E[9] = E[19];
		I[15] = I[20];
		E[15] = E[20];
		I[4] = I[24];
		E[4] = E[24];
		I[10] = I[25];
		E[10] = E[25];
		I[16] = I[26];
		E[16] = E[26];
		I[22] = I[27];
		E[22] = E[27];
		I[5] = I[30];
		E[5] = E[30];
		I[11] = I[31];
		E[11] = E[31];
		I[17] = I[32];
		E[17] = E[32];
		I[23] = I[33];
		E[23] = E[33];
		I[29] = I[34];
		E[29] = E[34];
		return false;
	}

	@Override
	protected void computeVariance()
	{
		variance[0] = I[0] * E[0] * I[0] + I[0] * E[1] * I[6] + I[0] * E[2] * I[12] + I[0] * E[3] * I[18] +
				I[0] * E[4] * I[24] + I[0] * E[5] * I[30] + I[1] * E[6] * I[0] + I[1] * E[7] * I[6] +
				I[1] * E[8] * I[12] + I[1] * E[9] * I[18] + I[1] * E[10] * I[24] + I[1] * E[11] * I[30] +
				I[2] * E[12] * I[0] + I[2] * E[13] * I[6] + I[2] * E[14] * I[12] + I[2] * E[15] * I[18] +
				I[2] * E[16] * I[24] + I[2] * E[17] * I[30] + I[3] * E[18] * I[0] + I[3] * E[19] * I[6] +
				I[3] * E[20] * I[12] + I[3] * E[21] * I[18] + I[3] * E[22] * I[24] + I[3] * E[23] * I[30] +
				I[4] * E[24] * I[0] + I[4] * E[25] * I[6] + I[4] * E[26] * I[12] + I[4] * E[27] * I[18] +
				I[4] * E[28] * I[24] + I[4] * E[29] * I[30] + I[5] * E[30] * I[0] + I[5] * E[31] * I[6] +
				I[5] * E[32] * I[12] + I[5] * E[33] * I[18] + I[5] * E[34] * I[24] + I[5] * E[35] * I[30];
		variance[1] = I[6] * E[0] * I[1] + I[6] * E[1] * I[7] + I[6] * E[2] * I[13] + I[6] * E[3] * I[19] +
				I[6] * E[4] * I[25] + I[6] * E[5] * I[31] + I[7] * E[6] * I[1] + I[7] * E[7] * I[7] +
				I[7] * E[8] * I[13] + I[7] * E[9] * I[19] + I[7] * E[10] * I[25] + I[7] * E[11] * I[31] +
				I[8] * E[12] * I[1] + I[8] * E[13] * I[7] + I[8] * E[14] * I[13] + I[8] * E[15] * I[19] +
				I[8] * E[16] * I[25] + I[8] * E[17] * I[31] + I[9] * E[18] * I[1] + I[9] * E[19] * I[7] +
				I[9] * E[20] * I[13] + I[9] * E[21] * I[19] + I[9] * E[22] * I[25] + I[9] * E[23] * I[31] +
				I[10] * E[24] * I[1] + I[10] * E[25] * I[7] + I[10] * E[26] * I[13] + I[10] * E[27] * I[19] +
				I[10] * E[28] * I[25] + I[10] * E[29] * I[31] + I[11] * E[30] * I[1] + I[11] * E[31] * I[7] +
				I[11] * E[32] * I[13] + I[11] * E[33] * I[19] + I[11] * E[34] * I[25] + I[11] * E[35] * I[31];
		variance[2] = I[12] * E[0] * I[2] + I[12] * E[1] * I[8] + I[12] * E[2] * I[14] + I[12] * E[3] * I[20] +
				I[12] * E[4] * I[26] + I[12] * E[5] * I[32] + I[13] * E[6] * I[2] + I[13] * E[7] * I[8] +
				I[13] * E[8] * I[14] + I[13] * E[9] * I[20] + I[13] * E[10] * I[26] + I[13] * E[11] * I[32] +
				I[14] * E[12] * I[2] + I[14] * E[13] * I[8] + I[14] * E[14] * I[14] + I[14] * E[15] * I[20] +
				I[14] * E[16] * I[26] + I[14] * E[17] * I[32] + I[15] * E[18] * I[2] + I[15] * E[19] * I[8] +
				I[15] * E[20] * I[14] + I[15] * E[21] * I[20] + I[15] * E[22] * I[26] + I[15] * E[23] * I[32] +
				I[16] * E[24] * I[2] + I[16] * E[25] * I[8] + I[16] * E[26] * I[14] + I[16] * E[27] * I[20] +
				I[16] * E[28] * I[26] + I[16] * E[29] * I[32] + I[17] * E[30] * I[2] + I[17] * E[31] * I[8] +
				I[17] * E[32] * I[14] + I[17] * E[33] * I[20] + I[17] * E[34] * I[26] + I[17] * E[35] * I[32];
		variance[3] = I[18] * E[0] * I[3] + I[18] * E[1] * I[9] + I[18] * E[2] * I[15] + I[18] * E[3] * I[21] +
				I[18] * E[4] * I[27] + I[18] * E[5] * I[33] + I[19] * E[6] * I[3] + I[19] * E[7] * I[9] +
				I[19] * E[8] * I[15] + I[19] * E[9] * I[21] + I[19] * E[10] * I[27] + I[19] * E[11] * I[33] +
				I[20] * E[12] * I[3] + I[20] * E[13] * I[9] + I[20] * E[14] * I[15] + I[20] * E[15] * I[21] +
				I[20] * E[16] * I[27] + I[20] * E[17] * I[33] + I[21] * E[18] * I[3] + I[21] * E[19] * I[9] +
				I[21] * E[20] * I[15] + I[21] * E[21] * I[21] + I[21] * E[22] * I[27] + I[21] * E[23] * I[33] +
				I[22] * E[24] * I[3] + I[22] * E[25] * I[9] + I[22] * E[26] * I[15] + I[22] * E[27] * I[21] +
				I[22] * E[28] * I[27] + I[22] * E[29] * I[33] + I[23] * E[30] * I[3] + I[23] * E[31] * I[9] +
				I[23] * E[32] * I[15] + I[23] * E[33] * I[21] + I[23] * E[34] * I[27] + I[23] * E[35] * I[33];
		variance[4] = I[24] * E[0] * I[4] + I[24] * E[1] * I[10] + I[24] * E[2] * I[16] + I[24] * E[3] * I[22] +
				I[24] * E[4] * I[28] + I[24] * E[5] * I[34] + I[25] * E[6] * I[4] + I[25] * E[7] * I[10] +
				I[25] * E[8] * I[16] + I[25] * E[9] * I[22] + I[25] * E[10] * I[28] + I[25] * E[11] * I[34] +
				I[26] * E[12] * I[4] + I[26] * E[13] * I[10] + I[26] * E[14] * I[16] + I[26] * E[15] * I[22] +
				I[26] * E[16] * I[28] + I[26] * E[17] * I[34] + I[27] * E[18] * I[4] + I[27] * E[19] * I[10] +
				I[27] * E[20] * I[16] + I[27] * E[21] * I[22] + I[27] * E[22] * I[28] + I[27] * E[23] * I[34] +
				I[28] * E[24] * I[4] + I[28] * E[25] * I[10] + I[28] * E[26] * I[16] + I[28] * E[27] * I[22] +
				I[28] * E[28] * I[28] + I[28] * E[29] * I[34] + I[29] * E[30] * I[4] + I[29] * E[31] * I[10] +
				I[29] * E[32] * I[16] + I[29] * E[33] * I[22] + I[29] * E[34] * I[28] + I[29] * E[35] * I[34];
		variance[5] = I[30] * E[0] * I[5] + I[30] * E[1] * I[11] + I[30] * E[2] * I[17] + I[30] * E[3] * I[23] +
				I[30] * E[4] * I[29] + I[30] * E[5] * I[35] + I[31] * E[6] * I[5] + I[31] * E[7] * I[11] +
				I[31] * E[8] * I[17] + I[31] * E[9] * I[23] + I[31] * E[10] * I[29] + I[31] * E[11] * I[35] +
				I[32] * E[12] * I[5] + I[32] * E[13] * I[11] + I[32] * E[14] * I[17] + I[32] * E[15] * I[23] +
				I[32] * E[16] * I[29] + I[32] * E[17] * I[35] + I[33] * E[18] * I[5] + I[33] * E[19] * I[11] +
				I[33] * E[20] * I[17] + I[33] * E[21] * I[23] + I[33] * E[22] * I[29] + I[33] * E[23] * I[35] +
				I[34] * E[24] * I[5] + I[34] * E[25] * I[11] + I[34] * E[26] * I[17] + I[34] * E[27] * I[23] +
				I[34] * E[28] * I[29] + I[34] * E[29] * I[35] + I[35] * E[30] * I[5] + I[35] * E[31] * I[11] +
				I[35] * E[32] * I[17] + I[35] * E[33] * I[23] + I[35] * E[34] * I[29] + I[35] * E[35] * I[35];
	}
}
