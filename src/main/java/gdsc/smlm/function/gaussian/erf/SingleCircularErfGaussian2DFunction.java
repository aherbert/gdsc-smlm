package gdsc.smlm.function.gaussian.erf;

import org.apache.commons.math3.util.Pair;

import gdsc.smlm.function.gaussian.Gaussian2DFunction;

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
 * Evaluates a 2-dimensional Gaussian function for a single peak.
 */
public class SingleCircularErfGaussian2DFunction extends SingleFreeCircularErfGaussian2DFunction
{
	private static final int[] gradientIndices;
	static
	{
		gradientIndices = createGradientIndices(1, new SingleCircularErfGaussian2DFunction(1, 1, 0));
	}

	/**
	 * Constructor.
	 *
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param derivativeOrder
	 *            Set to the order of partial derivatives required
	 */
	public SingleCircularErfGaussian2DFunction(int maxx, int maxy, int derivativeOrder)
	{
		super(maxx, maxy, derivativeOrder);
	}

	@Override
	public Gaussian2DFunction create(int derivativeOrder)
	{
		if (derivativeOrder == this.derivativeOrder)
			return this;
		return new SingleCircularErfGaussian2DFunction(maxx, maxy, derivativeOrder);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.fitting.function.NonLinearFunction#initialise(double[])
	 */
	public void initialise(double[] a)
	{
		tI = a[Gaussian2DFunction.SIGNAL];
		tB = a[Gaussian2DFunction.BACKGROUND];
		// Pre-compute the offset by 0.5
		final double tx = a[Gaussian2DFunction.X_POSITION] + 0.5;
		final double ty = a[Gaussian2DFunction.Y_POSITION] + 0.5;
		final double s = a[Gaussian2DFunction.X_SD];

		// We can pre-compute part of the derivatives for position and sd in arrays 
		// since the Gaussian is XY separable

		if (derivativeOrder == (byte) 2)
		{
			final double one_sSqrt2pi = ONE_OVER_ROOT2PI / s;
			final double ss = s * s;
			final double one_sSqrt2 = ONE_OVER_ROOT2 / s;
			final double one_2ss = 0.5 / ss;
			final double I_sSqrt2pi = tI * ONE_OVER_ROOT2PI / s;
			final double I_ssSqrt2pi = tI * ONE_OVER_ROOT2PI / ss;
			final double I_sssSqrt2pi = I_sSqrt2pi / ss;
			final double one_sssSqrt2pi = one_sSqrt2pi / ss;
			final double one_sssssSqrt2pi = one_sssSqrt2pi / ss;
			createSecondOrderTables(tI, one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi, I_sssSqrt2pi, ss, one_sssSqrt2pi,
					one_sssssSqrt2pi, deltaEx, du_dtx, du_dtsx, d2u_dtx2, d2u_dtsx2, tx);
			createSecondOrderTables(tI, one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi, I_sssSqrt2pi, ss, one_sssSqrt2pi,
					one_sssssSqrt2pi, deltaEy, du_dty, du_dtsy, d2u_dty2, d2u_dtsy2, ty);
		}
		else if (derivativeOrder == (byte) 1)
		{
			final double one_sSqrt2 = ONE_OVER_ROOT2 / s;
			final double one_2ss = 0.5 / (s * s);
			final double I_sSqrt2pi = tI * ONE_OVER_ROOT2PI / s;
			final double I_ssSqrt2pi = tI * ONE_OVER_ROOT2PI / (s * s);
			createFirstOrderTables(one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi, deltaEx, du_dtx, du_dtsx, tx);
			createFirstOrderTables(one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi, deltaEy, du_dty, du_dtsy, ty);
		}
		else
		{
			final double one_sSqrt2 = ONE_OVER_ROOT2 / s;
			createDeltaETable(one_sSqrt2, deltaEx, tx);
			createDeltaETable(one_sSqrt2, deltaEy, ty);
		}
	}

	/**
	 * Evaluates an 2-dimensional Gaussian function for a single peak.
	 * 
	 * @param i
	 *            Input predictor
	 * @param duda
	 *            Partial gradient of function with respect to each coefficient
	 * @return The predicted value
	 * 
	 * @see gdsc.smlm.function.NonLinearFunction#eval(int, double[])
	 */
	public double eval(final int i, final double[] duda)
	{
		// Unpack the predictor into the dimensions
		final int y = i / maxx;
		final int x = i % maxx;

		// Return in order of Gaussian2DFunction.createGradientIndices().
		// Use pre-computed gradients
		duda[0] = 1.0;
		duda[1] = deltaEx[x] * deltaEy[y];
		duda[2] = du_dtx[x] * deltaEy[y];
		duda[3] = du_dty[y] * deltaEx[x];
		duda[4] = du_dtsx[x] * deltaEy[y] + du_dtsy[y] * deltaEx[x];

		return tB + tI * duda[1];
	}

	/**
	 * Evaluates an 2-dimensional Gaussian function for a single peak.
	 * 
	 * @param i
	 *            Input predictor
	 * @param duda
	 *            Partial first gradient of function with respect to each coefficient
	 * @param d2uda2
	 *            Partial second gradient of function with respect to each coefficient
	 * @return The predicted value
	 */
	public double eval(final int i, final double[] duda, final double[] d2uda2)
	{
		// Unpack the predictor into the dimensions
		final int y = i / maxx;
		final int x = i % maxx;

		// Return in order of Gaussian2DFunction.createGradientIndices().
		// Use pre-computed gradients
		duda[0] = 1.0;
		duda[1] = deltaEx[x] * deltaEy[y];
		duda[2] = du_dtx[x] * deltaEy[y];
		duda[3] = du_dty[y] * deltaEx[x];
		duda[4] = du_dtsx[x] * deltaEy[y] + du_dtsy[y] * deltaEx[x];
		d2uda2[0] = 0;
		d2uda2[1] = 0;
		d2uda2[2] = d2u_dtx2[x] * deltaEy[y];
		d2uda2[3] = d2u_dty2[y] * deltaEx[x];
		// Working  example of this in GraspJ source code:
		// https://github.com/isman7/graspj/blob/master/graspj/src/main/java/eu/brede/graspj/opencl/src/functions/psfmodel_derivatives_sigma.cl
		//@formatter:off
		d2uda2[4] = d2u_dtsx2[x] * deltaEy[y] + 
				    d2u_dtsy2[y] * deltaEx[x] + 
				    2 * du_dtsx[x] * du_dtsy[y] / tI;
		//@formatter:on

		return tB + tI * duda[1];
	}

	// Support for ExtendedNonLinear Function. This can take advantage of x/y iteration.
	@Override
	public double[][] computeJacobian(double[] variables)
	{
		initialise(variables);

		final int n = maxx * maxy;
		final int nParams = 5;
		final double[][] jacobian = new double[n][];

		for (int y = 0, i = 0; y < maxy; y++)
		{
			final double dudy = this.du_dty[y];
			final double deltaEy = this.deltaEy[y];
			final double dudsy = this.du_dtsy[y];
			for (int x = 0; x < maxx; x++, i++)
			{
				final double[] duda = new double[nParams];
				duda[0] = 1.0;
				duda[1] = deltaEx[x] * deltaEy;
				duda[2] = du_dtx[x] * deltaEy;
				duda[3] = dudy * deltaEx[x];
				duda[4] = du_dtsx[x] * deltaEy + dudsy * deltaEx[x];
				jacobian[i] = duda;
			}
		}

		return jacobian;
	}

	// Support for ExtendedNonLinear Function. This can take advantage of x/y iteration.
	@Override
	public Pair<double[], double[][]> computeValuesAndJacobian(double[] variables)
	{
		initialise(variables);

		final int n = maxx * maxy;
		final int nParams = 5;
		final double[][] jacobian = new double[n][];
		final double[] values = new double[n];

		for (int y = 0, i = 0; y < maxy; y++)
		{
			final double dudy = this.du_dty[y];
			final double deltaEy = this.deltaEy[y];
			final double dudsy = this.du_dtsy[y];
			for (int x = 0; x < maxx; x++, i++)
			{
				final double[] duda = new double[nParams];
				duda[0] = 1.0;
				duda[1] = deltaEx[x] * deltaEy;
				duda[2] = du_dtx[x] * deltaEy;
				duda[3] = dudy * deltaEx[x];
				duda[4] = du_dtsx[x] * deltaEy + dudsy * deltaEx[x];
				values[i] = tB + tI * duda[1];
				jacobian[i] = duda;
			}
		}

		return new Pair<double[], double[][]>(values, jacobian);
	}

	@Override
	public int getNPeaks()
	{
		return 1;
	}

	@Override
	public boolean evaluatesBackground()
	{
		return true;
	}

	@Override
	public boolean evaluatesSignal()
	{
		return true;
	}

	@Override
	public boolean evaluatesShape()
	{
		return false;
	}

	@Override
	public boolean evaluatesPosition()
	{
		return true;
	}

	@Override
	public boolean evaluatesSD0()
	{
		return true;
	}

	@Override
	public boolean evaluatesSD1()
	{
		return false;
	}

	@Override
	public int getParametersPerPeak()
	{
		return 4;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.fitting.function.NonLinearFunction#gradientIndices()
	 */
	public int[] gradientIndices()
	{
		return gradientIndices;
	}
}
