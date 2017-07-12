package gdsc.smlm.function.gaussian.erf;

import gdsc.smlm.function.ExtendedGradient2Procedure;
import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.Gradient2Procedure;
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
public class MultiCircularErfGaussian2DFunction extends MultiFreeCircularErfGaussian2DFunction
{
	/**
	 * Constructor.
	 *
	 * @param nPeaks
	 *            The number of peaks
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 */
	public MultiCircularErfGaussian2DFunction(int nPeaks, int maxx, int maxy)
	{
		super(nPeaks, maxx, maxy);
	}

	@Override
	protected int[] createGradientIndices()
	{
		return replicateGradientIndices(SingleCircularErfGaussian2DFunction.gradientIndices);
	}

	@Override
	public ErfGaussian2DFunction copy()
	{
		return new MultiCircularErfGaussian2DFunction(nPeaks, maxx, maxy);
	}

	public void initialise0(double[] a)
	{
		tB = a[Gaussian2DFunction.BACKGROUND];
		for (int n = 0, i = 0; n < nPeaks; n++, i += PARAMETERS_PER_PEAK)
		{
			tI[n] = a[i + Gaussian2DFunction.SIGNAL];
			// Pre-compute the offset by 0.5
			final double tx = a[i + Gaussian2DFunction.X_POSITION] + 0.5;
			final double ty = a[i + Gaussian2DFunction.Y_POSITION] + 0.5;
			final double s = abs(a[i + Gaussian2DFunction.X_SD]);
			final double one_sSqrt2 = ONE_OVER_ROOT2 / s;

			createDeltaETable(n, maxx, one_sSqrt2, deltaEx, tx);
			createDeltaETable(n, maxy, one_sSqrt2, deltaEy, ty);
		}
	}

	public void initialise1(double[] a)
	{
		create1Arrays();
		tB = a[Gaussian2DFunction.BACKGROUND];
		for (int n = 0, i = 0; n < nPeaks; n++, i += PARAMETERS_PER_PEAK)
		{
			tI[n] = a[i + Gaussian2DFunction.SIGNAL];
			// Pre-compute the offset by 0.5
			final double tx = a[i + Gaussian2DFunction.X_POSITION] + 0.5;
			final double ty = a[i + Gaussian2DFunction.Y_POSITION] + 0.5;
			final double s = abs(a[i + Gaussian2DFunction.X_SD]);

			// We can pre-compute part of the derivatives for position and sd in arrays 
			// since the Gaussian is XY separable
			final double one_sSqrt2 = ONE_OVER_ROOT2 / s;
			final double one_2ss = 0.5 / (s * s);
			final double I_sSqrt2pi = tI[n] * ONE_OVER_ROOT2PI / s;
			final double I_ssSqrt2pi = tI[n] * ONE_OVER_ROOT2PI / (s * s);

			// We can pre-compute part of the derivatives for position and sd in arrays 
			// since the Gaussian is XY separable
			createFirstOrderTables(n, maxx, one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi, deltaEx, du_dtx, du_dtsx, tx);
			createFirstOrderTables(n, maxy, one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi, deltaEy, du_dty, du_dtsy, ty);
		}
	}

	public void initialise2(double[] a)
	{
		create2Arrays();
		tB = a[Gaussian2DFunction.BACKGROUND];
		for (int n = 0, i = 0; n < nPeaks; n++, i += PARAMETERS_PER_PEAK)
		{
			tI[n] = a[i + Gaussian2DFunction.SIGNAL];
			// Pre-compute the offset by 0.5
			final double tx = a[i + Gaussian2DFunction.X_POSITION] + 0.5;
			final double ty = a[i + Gaussian2DFunction.Y_POSITION] + 0.5;
			final double s = abs(a[i + Gaussian2DFunction.X_SD]);

			// We can pre-compute part of the derivatives for position and sd in arrays 
			// since the Gaussian is XY separable
			final double one_sSqrt2pi = ONE_OVER_ROOT2PI / s;
			final double ss = s * s;
			final double one_sSqrt2 = ONE_OVER_ROOT2 / s;
			final double one_2ss = 0.5 / ss;
			final double I_sSqrt2pi = tI[n] * ONE_OVER_ROOT2PI / s;
			final double I_ssSqrt2pi = tI[n] * ONE_OVER_ROOT2PI / ss;
			final double I_sssSqrt2pi = I_sSqrt2pi / ss;
			final double one_sssSqrt2pi = one_sSqrt2pi / ss;
			final double one_sssssSqrt2pi = one_sssSqrt2pi / ss;

			// We can pre-compute part of the derivatives for position and sd in arrays 
			// since the Gaussian is XY separable
			createSecondOrderTables(n, maxx, tI[n], one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi, I_sssSqrt2pi, ss,
					one_sssSqrt2pi, one_sssssSqrt2pi, deltaEx, du_dtx, du_dtsx, d2u_dtx2, d2u_dtsx2, tx);
			createSecondOrderTables(n, maxy, tI[n], one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi, I_sssSqrt2pi, ss,
					one_sssSqrt2pi, one_sssssSqrt2pi, deltaEy, du_dty, du_dtsy, d2u_dty2, d2u_dtsy2, ty);
		}
	}

	public void initialiseExtended2(double[] a)
	{
		createEx2Arrays();
		tB = a[Gaussian2DFunction.BACKGROUND];
		for (int n = 0, i = 0; n < nPeaks; n++, i += PARAMETERS_PER_PEAK)
		{
			tI[n] = a[i + Gaussian2DFunction.SIGNAL];
			// Pre-compute the offset by 0.5
			final double tx = a[i + Gaussian2DFunction.X_POSITION] + 0.5;
			final double ty = a[i + Gaussian2DFunction.Y_POSITION] + 0.5;
			final double s = abs(a[i + Gaussian2DFunction.X_SD]);

			// We can pre-compute part of the derivatives for position and sd in arrays 
			// since the Gaussian is XY separable
			final double one_sSqrt2pi = ONE_OVER_ROOT2PI / s;
			final double ss = s * s;
			final double one_sSqrt2 = ONE_OVER_ROOT2 / s;
			final double one_2ss = 0.5 / ss;
			final double I_sSqrt2pi = tI[n] * ONE_OVER_ROOT2PI / s;
			final double I_ssSqrt2pi = tI[n] * ONE_OVER_ROOT2PI / ss;
			final double I_sssSqrt2pi = I_sSqrt2pi / ss;
			final double one_sssSqrt2pi = one_sSqrt2pi / ss;
			final double one_sssssSqrt2pi = one_sssSqrt2pi / ss;

			// We can pre-compute part of the derivatives for position and sd in arrays 
			// since the Gaussian is XY separable
			createExSecondOrderTables(n, maxx, tI[n], one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi, I_sssSqrt2pi, ss,
					one_sssSqrt2pi, one_sssssSqrt2pi, deltaEx, du_dtx, du_dtsx, d2u_dtx2, d2u_dtsx2, d2deltaEx_dtsxdx,
					tx);
			createExSecondOrderTables(n, maxy, tI[n], one_sSqrt2, one_2ss, I_sSqrt2pi, I_ssSqrt2pi, I_sssSqrt2pi, ss,
					one_sssSqrt2pi, one_sssssSqrt2pi, deltaEy, du_dty, du_dtsy, d2u_dty2, d2u_dtsy2, d2deltaEy_dtsydy,
					ty);
		}
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
			duda[a++] = du_dtsx[xx] * deltaEy[yy] + du_dtsy[yy] * deltaEx[xx];
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
			duda[a] = deltaEx[xx] * deltaEy[yy];
			I += tI[n] * duda[a];
			d2uda2[a++] = 0;
			duda[a] = du_dtx[xx] * deltaEy[yy];
			d2uda2[a++] = d2u_dtx2[xx] * deltaEy[yy];
			duda[a] = du_dty[yy] * deltaEx[xx];
			d2uda2[a++] = d2u_dty2[yy] * deltaEx[xx];
			duda[a] = du_dtsx[xx] * deltaEy[yy] + du_dtsy[yy] * deltaEx[xx];
			//@formatter:off
			d2uda2[a++] = d2u_dtsx2[xx] * deltaEy[yy] + 
					      d2u_dtsy2[yy] * deltaEx[xx] +
					      2 * du_dtsx[xx] * du_dtsy[yy] / tI[n];
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
	public boolean evaluatesAngle()
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
		for (int y = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++)
			{
				double I = tB;
				for (int n = 0, xx = x, yy = y, a = 1; n < nPeaks; n++, xx += maxx, yy += maxy)
				{
					duda[a] = deltaEx[xx] * deltaEy[yy];
					I += tI[n] * duda[a++];
					duda[a++] = du_dtx[xx] * deltaEy[yy];
					duda[a++] = du_dty[yy] * deltaEx[xx];
					duda[a++] = du_dtsx[xx] * deltaEy[yy] + du_dtsy[yy] * deltaEx[xx];
				}
				//invalidGradients(duda);
				procedure.execute(I, duda);
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.GradientFunction#forEach(gdsc.smlm.function.GradientFunction.Gradient2Procedure)
	 */
	public void forEach(Gradient2Procedure procedure)
	{
		final double[] duda = new double[getNumberOfGradients()];
		final double[] d2uda2 = new double[getNumberOfGradients()];
		final double[] two_du_dtsy_tI = new double[nPeaks];
		duda[0] = 1.0;
		for (int y = 0; y < maxy; y++)
		{
			for (int n = 0, yy = y; n < nPeaks; n++, yy += maxy)
				two_du_dtsy_tI[n] = 2 * this.du_dtsy[yy] / tI[n];
			for (int x = 0; x < maxx; x++)
			{
				double I = tB;
				for (int n = 0, xx = x, yy = y, a = 1; n < nPeaks; n++, xx += maxx, yy += maxy)
				{
					duda[a] = deltaEx[xx] * deltaEy[yy];
					I += tI[n] * duda[a++];
					duda[a] = du_dtx[xx] * deltaEy[yy];
					d2uda2[a++] = d2u_dtx2[xx] * deltaEy[yy];
					duda[a] = du_dty[yy] * deltaEx[xx];
					d2uda2[a++] = d2u_dty2[yy] * deltaEx[xx];
					duda[a] = du_dtsx[xx] * deltaEy[yy] + du_dtsy[yy] * deltaEx[xx];
					//@formatter:off
    				d2uda2[a++] = d2u_dtsx2[xx] * deltaEy[yy] + 
    							  d2u_dtsy2[yy] * deltaEx[xx] + 
    						      du_dtsx[xx] * two_du_dtsy_tI[n];
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
		final double[] du_dtsy_tI = new double[nPeaks];
		final double[] two_du_dtsy_tI = new double[nPeaks];
		for (int y = 0; y < maxy; y++)
		{
			for (int n = 0, yy = y; n < nPeaks; n++, yy += maxy)
			{
				du_dty_tI[n] = du_dty[yy] / tI[n];
				du_dtsy_tI[n] = du_dtsy[yy] / tI[n];
				two_du_dtsy_tI[n] = 2 * du_dtsy[yy] / tI[n];
			}
			for (int x = 0; x < maxx; x++)
			{
				double I = tB;
				for (int n = 0, xx = x, yy = y, a = 1; n < nPeaks; n++, xx += maxx, yy += maxy)
				{
					duda[a] = deltaEx[xx] * deltaEy[yy];
					I += tI[n] * duda[a];
					duda[a + 1] = du_dtx[xx] * deltaEy[yy];
					duda[a + 2] = du_dty[yy] * deltaEx[xx];
					duda[a + 3] = du_dtsx[xx] * deltaEy[yy] + du_dtsy[yy] * deltaEx[xx];

					// Compute all the partial second order derivatives
					final double tI = this.tI[n];

					// Background are all 0

					int k = a * ng + a;
					// Signal,X
					d2udadb[k + 1] = duda[a + 1] / tI;
					// Signal,Y
					d2udadb[k + 2] = duda[a + 2] / tI;
					// Signal,X SD
					d2udadb[k + 3] = duda[a + 3] / tI;

					a += 4;

					int kk = k + ng;
					// X,Signal
					d2udadb[kk] = d2udadb[k + 1];
					// X,X
					d2udadb[kk + 1] = d2u_dtx2[xx] * deltaEy[yy];
					// X,Y
					d2udadb[kk + 2] = du_dtx[xx] * du_dty_tI[n];
					// X,X SD
					d2udadb[kk + 3] = deltaEy[yy] * d2deltaEx_dtsxdx[xx] + du_dtx[xx] * du_dtsy_tI[n];

					int kkk = kk + ng;
					// Y,Signal
					d2udadb[kkk] = d2udadb[k + 2];
					// Y,X
					d2udadb[kkk + 1] = d2udadb[kk + 2];
					// Y,Y
					d2udadb[kkk + 2] = d2u_dty2[yy] * deltaEx[xx];
					// Y,X SD
					d2udadb[kkk + 3] = du_dty[yy] * du_dtsx_tI[xx] + deltaEx[xx] * d2deltaEy_dtsydy[yy];

					int kkkk = kkk + ng;
					// X SD,Signal
					d2udadb[kkkk] = d2udadb[k + 3];
					// X SD,X
					d2udadb[kkkk + 1] = d2udadb[kk + 3];
					// X SD,Y
					d2udadb[kkkk + 2] = d2udadb[kkk + 3];
					// X SD,X SD
					//@formatter:off
					d2udadb[kkkk + 3] = d2u_dtsx2[xx] * deltaEy[yy] + 
    						            d2u_dtsy2[yy] * deltaEx[xx] + 
    						            du_dtsx[xx] * two_du_dtsy_tI[n];
    				//@formatter:on
				}
				procedure.executeExtended(I, duda, d2udadb);
			}
		}
	}
}
