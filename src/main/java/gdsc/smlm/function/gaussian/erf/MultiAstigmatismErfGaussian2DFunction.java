package gdsc.smlm.function.gaussian.erf;

import gdsc.smlm.function.ExtendedGradient2Procedure;
import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.Gradient2Procedure;
import gdsc.smlm.function.gaussian.AstigmatismZModel;
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
public class MultiAstigmatismErfGaussian2DFunction extends MultiFreeCircularErfGaussian2DFunction
{
	protected final AstigmatismZModel zModel;

	// Required for the z-depth gradients
	protected double[] dtsx_dtz, d2tsx_dtz2, dtsy_dtz, d2tsy_dtz2;

	/**
	 * Constructor.
	 *
	 * @param nPeaks
	 *            The number of peaks
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param zModel
	 *            the z model
	 */
	public MultiAstigmatismErfGaussian2DFunction(int nPeaks, int maxx, int maxy, AstigmatismZModel zModel)
	{
		super(nPeaks, maxx, maxy);
		this.zModel = zModel;
	}

	@Override
	protected void create1Arrays()
	{
		if (du_dtx != null)
			return;
		du_dtx = new double[deltaEx.length];
		du_dty = new double[deltaEy.length];
		du_dtsx = new double[deltaEx.length];
		du_dtsy = new double[deltaEy.length];
		dtsx_dtz = new double[deltaEx.length];
		dtsy_dtz = new double[deltaEy.length];
	}

	@Override
	protected void create2Arrays()
	{
		if (d2u_dtx2 != null)
			return;
		d2u_dtx2 = new double[deltaEx.length];
		d2u_dty2 = new double[deltaEy.length];
		d2u_dtsx2 = new double[deltaEx.length];
		d2u_dtsy2 = new double[deltaEy.length];
		d2tsx_dtz2 = new double[deltaEx.length];
		d2tsy_dtz2 = new double[deltaEy.length];
		create1Arrays();
	}

	@Override
	protected int[] createGradientIndices()
	{
		return replicateGradientIndices(SingleAstigmatismErfGaussian2DFunction.gradientIndices);
	}

	@Override
	public ErfGaussian2DFunction copy()
	{
		return new MultiAstigmatismErfGaussian2DFunction(nPeaks, maxx, maxy, zModel);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction#initialise0(double[])
	 */
	public void initialise0(double[] a)
	{
		tB = a[Gaussian2DFunction.BACKGROUND];
		for (int n = 0, i = 0; n < nPeaks; n++, i += PARAMETERS_PER_PEAK)
		{
			tI[n] = a[i + Gaussian2DFunction.SIGNAL];
			// Pre-compute the offset by 0.5
			final double tx = a[i + Gaussian2DFunction.X_POSITION] + 0.5;
			final double ty = a[i + Gaussian2DFunction.Y_POSITION] + 0.5;
			final double tsx = abs(a[i + Gaussian2DFunction.X_SD]);
			final double tsy = abs(a[i + Gaussian2DFunction.Y_SD]);
			final double tz = a[i + Gaussian2DFunction.Z_POSITION];

			final double sx = tsx * zModel.getSx(tz);
			final double sy = tsy * zModel.getSy(tz);
			createDeltaETable(n, maxx, ONE_OVER_ROOT2 / sx, deltaEx, tx);
			createDeltaETable(n, maxy, ONE_OVER_ROOT2 / sy, deltaEy, ty);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction#initialise1(double[])
	 */
	public void initialise1(double[] a)
	{
		create1Arrays();
		final double[] ds_dz = new double[1];
		tB = a[Gaussian2DFunction.BACKGROUND];
		for (int n = 0, i = 0; n < nPeaks; n++, i += PARAMETERS_PER_PEAK)
		{
			tI[n] = a[i + Gaussian2DFunction.SIGNAL];
			// Pre-compute the offset by 0.5
			final double tx = a[i + Gaussian2DFunction.X_POSITION] + 0.5;
			final double ty = a[i + Gaussian2DFunction.Y_POSITION] + 0.5;
			final double tsx = abs(a[i + Gaussian2DFunction.X_SD]);
			final double tsy = abs(a[i + Gaussian2DFunction.Y_SD]);
			final double tz = a[i + Gaussian2DFunction.Z_POSITION];

			// We can pre-compute part of the derivatives for position and sd in arrays 
			// since the Gaussian is XY separable
			final double sx = tsx * zModel.getSx(tz, ds_dz);
			dtsx_dtz[n] = tsx * ds_dz[0];
			final double sy = tsy * zModel.getSy(tz, ds_dz);
			dtsy_dtz[n] = tsy * ds_dz[0];
			createFirstOrderTables(n, maxx, tI[n], deltaEx, du_dtx, du_dtsx, tx, sx);
			createFirstOrderTables(n, maxy, tI[n], deltaEy, du_dty, du_dtsy, ty, sy);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction#initialise2(double[])
	 */
	public void initialise2(double[] a)
	{
		create2Arrays();
		final double[] ds_dz = new double[2];
		tB = a[Gaussian2DFunction.BACKGROUND];
		for (int n = 0, i = 0; n < nPeaks; n++, i += PARAMETERS_PER_PEAK)
		{
			tI[n] = a[i + Gaussian2DFunction.SIGNAL];
			// Pre-compute the offset by 0.5
			final double tx = a[i + Gaussian2DFunction.X_POSITION] + 0.5;
			final double ty = a[i + Gaussian2DFunction.Y_POSITION] + 0.5;
			final double tsx = abs(a[i + Gaussian2DFunction.X_SD]);
			final double tsy = abs(a[i + Gaussian2DFunction.Y_SD]);
			final double tz = a[i + Gaussian2DFunction.Z_POSITION];

			// We can pre-compute part of the derivatives for position and sd in arrays 
			// since the Gaussian is XY separable
			final double sx = tsx * zModel.getSx2(tz, ds_dz);
			dtsx_dtz[n] = tsx * ds_dz[0];
			d2tsx_dtz2[n] = tsx * ds_dz[1];
			final double sy = tsy * zModel.getSy2(tz, ds_dz);
			dtsy_dtz[n] = tsy * ds_dz[0];
			d2tsy_dtz2[n] = tsy * ds_dz[1];
			createSecondOrderTables(n, maxx, tI[n], deltaEx, du_dtx, du_dtsx, d2u_dtx2, d2u_dtsx2, tx, sx);
			createSecondOrderTables(n, maxy, tI[n], deltaEy, du_dty, du_dtsy, d2u_dty2, d2u_dtsy2, ty, sy);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.gaussian.erf.MultiFreeCircularErfGaussian2DFunction#initialiseExtended2(double[])
	 */
	public void initialiseExtended2(double[] a)
	{
		createEx2Arrays();
		final double[] ds_dz = new double[2];
		tB = a[Gaussian2DFunction.BACKGROUND];
		for (int n = 0, i = 0; n < nPeaks; n++, i += PARAMETERS_PER_PEAK)
		{
			tI[n] = a[i + Gaussian2DFunction.SIGNAL];
			// Pre-compute the offset by 0.5
			final double tx = a[i + Gaussian2DFunction.X_POSITION] + 0.5;
			final double ty = a[i + Gaussian2DFunction.Y_POSITION] + 0.5;
			final double tsx = abs(a[i + Gaussian2DFunction.X_SD]);
			final double tsy = abs(a[i + Gaussian2DFunction.Y_SD]);
			final double tz = a[i + ErfGaussian2DFunction.Z_POSITION];

			// We can pre-compute part of the derivatives for position and sd in arrays 
			// since the Gaussian is XY separable
			final double sx = tsx * zModel.getSx2(tz, ds_dz);
			dtsx_dtz[n] = tsx * ds_dz[0];
			d2tsx_dtz2[n] = tsx * ds_dz[1];
			final double sy = tsy * zModel.getSy2(tz, ds_dz);
			dtsy_dtz[n] = tsy * ds_dz[0];
			d2tsy_dtz2[n] = tsy * ds_dz[1];
			createExSecondOrderTables(n, maxx, tI[n], deltaEx, du_dtx, du_dtsx, d2u_dtx2, d2u_dtsx2, d2deltaEx_dtsxdx,
					tx, sx);
			createExSecondOrderTables(n, maxy, tI[n], deltaEy, du_dty, du_dtsy, d2u_dty2, d2u_dtsy2, d2deltaEy_dtsydy,
					ty, sy);
		}
		// Pre-apply the gradient mapping from width to z
		for (int x = 0; x < maxx; x++)
			for (int n = 0, xx = x; n < nPeaks; n++, xx += maxx)
				d2deltaEx_dtsxdx[xx] *= dtsx_dtz[n];
		for (int y = 0; y < maxy; y++)
			for (int n = 0, yy = y; n < nPeaks; n++, yy += maxy)
				d2deltaEy_dtsydy[yy] *= dtsy_dtz[n];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.gaussian.erf.MultiErfGaussian2DFunction#eval(int, double[])
	 */
	public double eval(final int i, final double[] duda)
	{
		// Unpack the predictor into the dimensions
		int yy = i / maxx;
		int xx = i % maxx;

		// Return in order of Gaussian2DFunction.createGradientIndices().
		// Use pre-computed gradients
		duda[0] = 1.0;
		double I = tB;
		for (int n = 0, a = 1; n < nPeaks; n++, xx += maxx, yy += maxy)
		{
			duda[a] = deltaEx[xx] * deltaEy[yy];
			I += tI[n] * duda[a++];
			duda[a++] = du_dtx[xx] * deltaEy[yy];
			duda[a++] = du_dty[yy] * deltaEx[xx];
			duda[a++] = du_dtsx[xx] * deltaEy[yy] * dtsx_dtz[n] + du_dtsy[yy] * deltaEx[xx] * dtsy_dtz[n];
		}
		return I;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.gaussian.erf.MultiErfGaussian2DFunction#eval(int, double[], double[])
	 */
	public double eval(final int i, final double[] duda, final double[] d2uda2)
	{
		// Unpack the predictor into the dimensions
		int yy = i / maxx;
		int xx = i % maxx;

		// Return in order of Gaussian2DFunction.createGradientIndices().
		// Use pre-computed gradients
		duda[0] = 1.0;
		d2uda2[0] = 0;
		double I = tB;
		for (int n = 0, a = 1; n < nPeaks; n++, xx += maxx, yy += maxy)
		{
			final double du_dsx = du_dtsx[xx] * deltaEy[yy];
			final double du_dsy = du_dtsy[yy] * deltaEx[xx];

			duda[a] = deltaEx[xx] * deltaEy[yy];
			I += tI[n] * duda[a];
			d2uda2[a++] = 0;
			duda[a] = du_dtx[xx] * deltaEy[yy];
			d2uda2[a++] = d2u_dtx2[xx] * deltaEy[yy];
			duda[a] = du_dty[yy] * deltaEx[xx];
			d2uda2[a++] = d2u_dty2[yy] * deltaEx[xx];
			duda[a] = du_dsx * dtsx_dtz[n] + du_dsy * dtsy_dtz[n];
			//@formatter:off
			d2uda2[a++] =
					d2u_dtsx2[xx] * deltaEy[yy] * dtsx_dtz[n] * dtsx_dtz[n] +
					du_dsx * d2tsx_dtz2[n] +
					d2u_dtsy2[yy] * deltaEx[xx] * dtsy_dtz[n] * dtsy_dtz[n] + 
					du_dsy * d2tsy_dtz2[n] +
					// Add the equivalent term we add in the circular version.
					// Note: this is not in the Smith, et al (2010) paper but is 
					// in the GraspJ source code and it works in JUnit tests.
					2 * du_dtsx[xx] * dtsx_dtz[n] * du_dtsy[yy] * dtsy_dtz[n] / tI[n];		
			//@formatter:on
		}
		return I;
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
	public boolean evaluatesPosition()
	{
		return true;
	}

	@Override
	public boolean evaluatesZ()
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
	public int getGradientParametersPerPeak()
	{
		return 4;
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
		final double[] deltaEy_by_dtsx_dtz = new double[nPeaks];
		final double[] du_dtsy_by_dtsy_dtz = new double[nPeaks];
		for (int y = 0; y < maxy; y++)
		{
			for (int n = 0, yy = y; n < nPeaks; n++, yy += maxy)
			{
				deltaEy_by_dtsx_dtz[n] = deltaEy[yy] * dtsx_dtz[n];
				du_dtsy_by_dtsy_dtz[n] = du_dtsy[yy] * dtsy_dtz[n];
			}

			for (int x = 0; x < maxx; x++)
			{
				double I = tB;
				for (int n = 0, xx = x, yy = y, a = 1; n < nPeaks; n++, xx += maxx, yy += maxy)
				{
					duda[a] = deltaEx[xx] * deltaEy[yy];
					I += tI[n] * duda[a++];
					duda[a++] = du_dtx[xx] * deltaEy[yy];
					duda[a++] = du_dty[yy] * deltaEx[xx];
					duda[a++] = du_dtsx[xx] * deltaEy_by_dtsx_dtz[n] + du_dtsy_by_dtsy_dtz[n] * deltaEx[xx];
				}
				procedure.execute(I, duda);
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
		final double[] dtsx_dtz_2 = new double[nPeaks];
		final double[] dtsy_dtz_2 = new double[nPeaks];
		final double[] two_dtsx_dtz_by_dtsy_dtz_tI = new double[nPeaks];
		for (int n = 0; n < nPeaks; n++)
		{
			dtsx_dtz_2[n] = dtsx_dtz[n] * dtsx_dtz[n];
			dtsy_dtz_2[n] = dtsy_dtz[n] * dtsy_dtz[n];
			two_dtsx_dtz_by_dtsy_dtz_tI[n] = 2 * dtsx_dtz[n] * dtsy_dtz[n] / tI[n];
		}
		final double[] deltaEy_by_dtsx_dtz_2 = new double[nPeaks];
		final double[] d2u_dtsy2_by_dtsy_dtz_2 = new double[nPeaks];
		final double[] two_dtsx_dtz_by_du_dtsy_by_dtsy_dtz_tI = new double[nPeaks];
		for (int y = 0; y < maxy; y++)
		{
			for (int n = 0, yy = y; n < nPeaks; n++, yy += maxy)
			{
				deltaEy_by_dtsx_dtz_2[n] = deltaEy[yy] * dtsx_dtz_2[n];
				d2u_dtsy2_by_dtsy_dtz_2[n] = d2u_dtsy2[yy] * dtsy_dtz_2[n];
				two_dtsx_dtz_by_du_dtsy_by_dtsy_dtz_tI[n] = two_dtsx_dtz_by_dtsy_dtz_tI[n] * du_dtsy[yy];
			}

			for (int x = 0; x < maxx; x++)
			{
				double I = tB;
				for (int n = 0, xx = x, yy = y, a = 1; n < nPeaks; n++, xx += maxx, yy += maxy)
				{
					final double du_dsx = du_dtsx[xx] * deltaEy[yy];
					final double du_dsy = du_dtsy[yy] * deltaEx[xx];

					duda[a] = deltaEx[xx] * deltaEy[yy];
					I += tI[n] * duda[a++];
					duda[a] = du_dtx[xx] * deltaEy[yy];
					d2uda2[a++] = d2u_dtx2[xx] * deltaEy[yy];
					duda[a] = du_dty[yy] * deltaEx[xx];
					d2uda2[a++] = d2u_dty2[yy] * deltaEx[xx];
					duda[a] = du_dsx * dtsx_dtz[n] + du_dsy * dtsy_dtz[n];
					//@formatter:off
					d2uda2[a++] = 
							d2u_dtsx2[xx] * deltaEy_by_dtsx_dtz_2[n] + 
							du_dsx * d2tsx_dtz2[n] +
							d2u_dtsy2_by_dtsy_dtz_2[n] * deltaEx[xx] + 
							du_dsy * d2tsy_dtz2[n] +
							// Add the equivalent term we add in the circular version.
							// Note: this is not in the Smith, et al (2010) paper but is 
							// in the GraspJ source code and it works in JUnit tests.
							//2 * du_dtsx[x] * dtsx_dtz * du_dtsy * dtsy_dtz / tI;
							two_dtsx_dtz_by_du_dtsy_by_dtsy_dtz_tI[n] * du_dtsx[xx];
					//@formatter:on
				}
				procedure.execute(I, duda, d2uda2);
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.ExtendedGradient2Function#forEach(gdsc.smlm.function.ExtendedGradient2Procedure)
	 */
	public void forEach(ExtendedGradient2Procedure procedure)
	{
		final int ng = getNumberOfGradients();
		final double[] duda = new double[ng];
		final double[] d2udadb = new double[ng * ng];
		duda[0] = 1.0;
		final double[] du_dtsx_tI = new double[du_dtsx.length];
		for (int x = 0; x < maxx; x++)
			for (int n = 0, xx = x; n < nPeaks; n++, xx += maxx)
				du_dtsx_tI[xx] = du_dtsx[xx] / tI[n];
		final double[] du_dty_tI = new double[nPeaks];
		final double[] du_dtsy_by_dtsy_dtz_tI = new double[nPeaks];
		final double[] du_dty_by_dtsx_dtz_tI = new double[nPeaks];
		final double[] deltaEy_by_dtsx_dtz_2 = new double[nPeaks];
		final double[] d2u_dtsy2_by_dtsy_dtz_2 = new double[nPeaks];
		final double[] two_dtsx_dtz_by_du_dtsy_by_dtsy_dtz_tI = new double[nPeaks];

		final double[] dtsx_dtz_2 = new double[nPeaks];
		final double[] dtsy_dtz_2 = new double[nPeaks];
		final double[] two_dtsx_dtz_by_dtsy_dtz_tI = new double[nPeaks];
		final double[] dtsx_dtz_tI = new double[nPeaks];
		final double[] dtsy_dtz_tI = new double[nPeaks];
		for (int n = 0; n < nPeaks; n++)
		{
			dtsx_dtz_2[n] = dtsx_dtz[n] * dtsx_dtz[n];
			dtsy_dtz_2[n] = dtsy_dtz[n] * dtsy_dtz[n];
			two_dtsx_dtz_by_dtsy_dtz_tI[n] = 2 * dtsx_dtz[n] * dtsy_dtz[n] / tI[n];
			dtsx_dtz_tI[n] = dtsx_dtz[n] / tI[n];
			dtsy_dtz_tI[n] = dtsy_dtz[n] / tI[n];
		}

		for (int y = 0; y < maxy; y++)
		{
			for (int n = 0, yy = y; n < nPeaks; n++, yy += maxy)
			{
				du_dty_tI[n] = du_dty[yy] / tI[n];
				du_dtsy_by_dtsy_dtz_tI[n] = du_dtsy[yy] * dtsy_dtz_tI[n];
				du_dty_by_dtsx_dtz_tI[n] = du_dty[yy] * dtsx_dtz_tI[n];
				deltaEy_by_dtsx_dtz_2[n] = deltaEy[yy] * dtsx_dtz_2[n];
				d2u_dtsy2_by_dtsy_dtz_2[n] = d2u_dtsy2[yy] * dtsy_dtz_2[n];
				two_dtsx_dtz_by_du_dtsy_by_dtsy_dtz_tI[n] = two_dtsx_dtz_by_dtsy_dtz_tI[n] * du_dtsy[yy];
			}
			for (int x = 0; x < maxx; x++)
			{
				double I = tB;
				for (int n = 0, xx = x, yy = y, a = 1; n < nPeaks; n++, xx += maxx, yy += maxy)
				{
					final double du_dsx = du_dtsx[xx] * deltaEy[yy];
					final double du_dsy = du_dtsy[yy] * deltaEx[xx];

					duda[a] = deltaEx[xx] * deltaEy[yy];
					I += tI[n] * duda[a];
					duda[a + 1] = du_dtx[xx] * deltaEy[yy];
					duda[a + 2] = du_dty[yy] * deltaEx[xx];
					duda[a + 3] = du_dsx * dtsx_dtz[n] + du_dsy * dtsy_dtz[n];

					// Compute all the partial second order derivatives
					final double tI = this.tI[n];

					// Background are all 0

					int k = a * ng + a;
					// Signal,X
					d2udadb[k + 1] = duda[a + 1] / tI;
					// Signal,Y
					d2udadb[k + 2] = duda[a + 2] / tI;
					// Signal,Z
					d2udadb[k + 3] = duda[a + 3] / tI;

					a += 4;

					int kk = k + ng;
					// X,Signal
					d2udadb[kk] = d2udadb[k + 1];
					// X,X
					d2udadb[kk + 1] = d2u_dtx2[xx] * deltaEy[yy];
					// X,Y
					d2udadb[kk + 2] = du_dtx[xx] * du_dty_tI[n];
					// X,Z
					d2udadb[kk + 3] = deltaEy[yy] * d2deltaEx_dtsxdx[xx] + du_dtx[xx] * du_dtsy_by_dtsy_dtz_tI[n];

					int kkk = kk + ng;
					// Y,Signal
					d2udadb[kkk] = d2udadb[k + 2];
					// Y,X
					d2udadb[kkk + 1] = d2udadb[kk + 2];
					// Y,Y
					d2udadb[kkk + 2] = d2u_dty2[yy] * deltaEx[xx];
					// X,Z
					d2udadb[kkk + 3] = du_dtsx[xx] * du_dty_by_dtsx_dtz_tI[n] + deltaEx[xx] * d2deltaEy_dtsydy[yy];

					int kkkk = kkk + ng;
					// Z,Signal
					d2udadb[kkkk] = d2udadb[k + 3];
					// Z,X
					d2udadb[kkkk + 1] = d2udadb[kk + 3];
					// Z,Y
					d2udadb[kkkk + 2] = d2udadb[kkk + 3];
					// Z,Z
					//@formatter:off
					d2udadb[kkkk + 3] =
							d2u_dtsx2[xx] * deltaEy_by_dtsx_dtz_2[n] + 
							du_dsx * d2tsx_dtz2[n] +
							d2u_dtsy2_by_dtsy_dtz_2[n] * deltaEx[xx] + 
							du_dsy * d2tsy_dtz2[n] +
							// Add the equivalent term we add in the circular version.
							// Note: this is not in the Smith, et al (2010) paper but is 
							// in the GraspJ source code and it works in JUnit tests.
							//2 * du_dtsx[x] * dtsx_dtz * du_dtsy * dtsy_dtz / tI;
							two_dtsx_dtz_by_du_dtsy_by_dtsy_dtz_tI[n] * du_dtsx[xx];
					//@formatter:on
				}
				procedure.executeExtended(I, duda, d2udadb);
			}
		}
	}
}
