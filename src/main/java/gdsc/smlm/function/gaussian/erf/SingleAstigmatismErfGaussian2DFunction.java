package gdsc.smlm.function.gaussian.erf;

import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.Gradient2Procedure;
import gdsc.smlm.function.gaussian.AstimatismZModel;
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
public class SingleAstigmatismErfGaussian2DFunction extends SingleFreeCircularErfGaussian2DFunction
{
	private static final int[] gradientIndices;
	static
	{
		gradientIndices = createGradientIndices(1, new SingleAstigmatismErfGaussian2DFunction(1, 1, 0, null));
	}

	protected final AstimatismZModel zModel;

	// Required for the z-depth gradients
	protected double dtsx_dtz, d2tsx_dtz2, dtsy_dtz, d2tsy_dtz2;

	/**
	 * Constructor.
	 *
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param derivativeOrder
	 *            Set to the order of partial derivatives required
	 * @param zModel
	 *            the z model
	 */
	public SingleAstigmatismErfGaussian2DFunction(int maxx, int maxy, int derivativeOrder, AstimatismZModel zModel)
	{
		super(maxx, maxy, derivativeOrder);
		this.zModel = zModel;
	}

	@Override
	public ErfGaussian2DFunction create(int derivativeOrder)
	{
		if (derivativeOrder == this.derivativeOrder)
			return this;
		return new SingleAstigmatismErfGaussian2DFunction(maxx, maxy, derivativeOrder, zModel);
	}

	@Override
	public ErfGaussian2DFunction copy()
	{
		return new SingleAstigmatismErfGaussian2DFunction(maxx, maxy, derivativeOrder, zModel);
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
		final double tsx = a[Gaussian2DFunction.X_SD];
		final double tsy = a[Gaussian2DFunction.Y_SD];
		final double tz = a[ErfGaussian2DFunction.Z_POSITION];

		// We can pre-compute part of the derivatives for position and sd in arrays 
		// since the Gaussian is XY separable

		if (derivativeOrder == (byte) 2)
		{
			final double[] ds_dz = new double[2];
			final double sx = tsx * zModel.getSx2(tz, ds_dz);
			dtsx_dtz = tsx * ds_dz[0];
			d2tsx_dtz2 = tsx * ds_dz[1];
			final double sy = tsy * zModel.getSy2(tz, ds_dz);
			dtsy_dtz = tsy * ds_dz[0];
			d2tsy_dtz2 = tsy * ds_dz[1];
			createSecondOrderTables(tI, deltaEx, du_dtx, du_dtsx, d2u_dtx2, d2u_dtsx2, tx, sx);
			createSecondOrderTables(tI, deltaEy, du_dty, du_dtsy, d2u_dty2, d2u_dtsy2, ty, sy);
		}
		else if (derivativeOrder == (byte) 1)
		{
			final double[] ds_dz = new double[2];
			final double sx = tsx * zModel.getSx(tz, ds_dz);
			dtsx_dtz = tsx * ds_dz[0];
			final double sy = tsy * zModel.getSy(tz, ds_dz);
			dtsy_dtz = tsy * ds_dz[0];
			createFirstOrderTables(tI, deltaEx, du_dtx, du_dtsx, tx, sx);
			createFirstOrderTables(tI, deltaEy, du_dty, du_dtsy, ty, sy);
		}
		else
		{
			final double sx = tsx * zModel.getSx(tz);
			final double sy = tsy * zModel.getSy(tz);
			createDeltaETable(ONE_OVER_ROOT2 / sx, deltaEx, tx);
			createDeltaETable(ONE_OVER_ROOT2 / sy, deltaEy, ty);
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
		duda[2] = du_dtsx[x] * deltaEy[y] * dtsx_dtz + du_dtsy[y] * deltaEx[x] * dtsy_dtz;
		duda[3] = du_dtx[x] * deltaEy[y];
		duda[4] = du_dty[y] * deltaEx[x];

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
		final double du_dsx = du_dtsx[x] * deltaEy[y];
		final double du_dsy = du_dtsy[y] * deltaEx[x];

		duda[0] = 1.0;
		duda[1] = deltaEx[x] * deltaEy[y];
		duda[2] = du_dsx * dtsx_dtz + du_dsy * dtsy_dtz;
		duda[3] = du_dtx[x] * deltaEy[y];
		duda[4] = du_dty[y] * deltaEx[x];
		d2uda2[0] = 0;
		d2uda2[1] = 0;
		//@formatter:off
		d2uda2[2] =
				d2u_dtsx2[x] * deltaEy[y] * dtsx_dtz * dtsx_dtz +
				du_dsx * d2tsx_dtz2 +
				d2u_dtsy2[y] * deltaEx[x] * dtsy_dtz * dtsy_dtz + 
				du_dsy * d2tsy_dtz2 +
				// Add the equivalent term we add in the circular version.
				// Note: this is not in the Smith, et al (2010) paper but is 
				// in the GraspJ source code and it works in JUnit tests.
				+ 2 * du_dtsx[x] * dtsx_dtz * du_dtsy[y] * dtsy_dtz / tI;		
		//@formatter:on
		d2uda2[3] = d2u_dtx2[x] * deltaEy[y];
		d2uda2[4] = d2u_dty2[y] * deltaEx[x];

		return tB + tI * duda[1];
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
		return true;
	}

	@Override
	public boolean evaluatesPosition()
	{
		return true;
	}

	@Override
	public boolean evaluatesSD0()
	{
		return false;
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.GradientFunction#getNumberOfGradients()
	 */
	public int getNumberOfGradients()
	{
		return 5;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.GradientFunction#forEach(gdsc.smlm.function.GradientFunction.Gradient1Procedure)
	 */
	public void forEach(Gradient1Procedure procedure)
	{
		final double[] duda = new double[getNumberOfGradients()];
		duda[0] = 1.0;
		for (int y = 0; y < maxy; y++)
		{
			final double du_dty = this.du_dty[y];
			final double deltaEy = this.deltaEy[y];
			final double deltaEy_by_dtsx_dtz = deltaEy * dtsx_dtz;
			final double du_dtsy_by_dtsy_dtz = this.du_dtsy[y] * dtsy_dtz;
			for (int x = 0; x < maxx; x++)
			{
				duda[1] = deltaEx[x] * deltaEy;
				duda[2] = du_dtsx[x] * deltaEy_by_dtsx_dtz + du_dtsy_by_dtsy_dtz * deltaEx[x];
				duda[3] = du_dtx[x] * deltaEy;
				duda[4] = du_dty * deltaEx[x];
				procedure.execute(tB + tI * duda[1], duda);
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient2Function#forEach(gdsc.smlm.function.Gradient2Procedure)
	 */
	public void forEach(Gradient2Procedure procedure)
	{
		final double[] duda = new double[getNumberOfGradients()];
		final double[] d2uda2 = new double[getNumberOfGradients()];
		duda[0] = 1.0;
		final double dtsx_dtz_2 = dtsx_dtz * dtsx_dtz;
		final double dtsy_dtz_2 = dtsy_dtz * dtsy_dtz;
		final double two_dtsx_dtz_by_dtsy_dtz_tI = 2 * dtsx_dtz * dtsy_dtz / tI;
		for (int y = 0; y < maxy; y++)
		{
			final double du_dty = this.du_dty[y];
			final double deltaEy = this.deltaEy[y];
			final double du_dtsy = this.du_dtsy[y];
			final double d2u_dty2 = this.d2u_dty2[y];
			final double deltaEy_by_dtsx_dtz_2 = deltaEy * dtsx_dtz_2;
			final double d2u_dtsy2_by_dtsy_dtz_2 = this.d2u_dtsy2[y] * dtsy_dtz_2;
			final double two_dtsx_dtz_by_du_dtsy_by_dtsy_dtz_tI = two_dtsx_dtz_by_dtsy_dtz_tI * du_dtsy;
			for (int x = 0; x < maxx; x++)
			{
				final double du_dsx = du_dtsx[x] * deltaEy;
				final double du_dsy = du_dtsy * deltaEx[x];

				duda[1] = deltaEx[x] * deltaEy;
				duda[2] = du_dsx * dtsx_dtz + du_dsy * dtsy_dtz;
				duda[3] = du_dtx[x] * deltaEy;
				duda[4] = du_dty * deltaEx[x];
				//@formatter:off
				d2uda2[2] =
						d2u_dtsx2[x] * deltaEy_by_dtsx_dtz_2 +
						du_dsx * d2tsx_dtz2 +
						//d2u_dtsy2[y] * deltaEx[x] * dtsy_dtz_2 +
						d2u_dtsy2_by_dtsy_dtz_2 * deltaEx[x] +
						du_dsy * d2tsy_dtz2 +
						// Add the equivalent term we add in the circular version.
						// Note: this is not in the Smith, et al (2010) paper but is 
						// in the GraspJ source code and it works in JUnit tests.
						//2 * du_dtsx[x] * dtsx_dtz * du_dtsy * dtsy_dtz / tI;
						two_dtsx_dtz_by_du_dtsy_by_dtsy_dtz_tI * du_dtsx[x]; 
				//@formatter:on
				d2uda2[3] = d2u_dtx2[x] * deltaEy;
				d2uda2[4] = d2u_dty2 * deltaEx[x];

				procedure.execute(tB + tI * duda[1], duda, d2uda2);
			}
		}
	}
}
