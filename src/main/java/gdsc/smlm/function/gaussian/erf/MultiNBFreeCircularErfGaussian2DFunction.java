package gdsc.smlm.function.gaussian.erf;

import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.Gradient2Procedure;
import gdsc.smlm.function.ValueProcedure;

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
public class MultiNBFreeCircularErfGaussian2DFunction extends MultiFreeCircularErfGaussian2DFunction
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
	public MultiNBFreeCircularErfGaussian2DFunction(int nPeaks, int maxx, int maxy)
	{
		super(nPeaks, maxx, maxy);
	}

	@Override
	protected int[] createGradientIndices()
	{
		return replicateGradientIndices(SingleNBFreeCircularErfGaussian2DFunction.gradientIndices);
	}

	@Override
	public ErfGaussian2DFunction copy()
	{
		return new MultiNBFreeCircularErfGaussian2DFunction(nPeaks, maxx, maxy);
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
		int yy = i / maxx;
		int xx = i % maxx;

		// Return in order of Gaussian2DFunction.createGradientIndices().
		// Use pre-computed gradients
		double I = tB;
		for (int n = 0, a = 0; n < nPeaks; n++, xx += maxx, yy += maxy)
		{
			duda[a] = deltaEx[xx] * deltaEy[yy];
			I += tI[n] * duda[a++];
			duda[a++] = du_dtx[xx] * deltaEy[yy];
			duda[a++] = du_dty[yy] * deltaEx[xx];
			duda[a++] = du_dtsx[xx] * deltaEy[yy];
			duda[a++] = du_dtsy[yy] * deltaEx[xx];
		}
		return I;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction#eval(int, double[], double[])
	 */
	public double eval(final int i, final double[] duda, final double[] d2uda2)
	{
		// Unpack the predictor into the dimensions
		int yy = i / maxx;
		int xx = i % maxx;

		// Return in order of Gaussian2DFunction.createGradientIndices().
		// Use pre-computed gradients
		double I = tB;
		for (int n = 0, a = 0; n < nPeaks; n++, xx += maxx, yy += maxy)
		{
			duda[a] = deltaEx[xx] * deltaEy[yy];
			I += tI[n] * duda[a];
			d2uda2[a++] = 0;
			duda[a] = du_dtx[xx] * deltaEy[yy];
			d2uda2[a++] = d2u_dtx2[xx] * deltaEy[yy];
			duda[a] = du_dty[yy] * deltaEx[xx];
			d2uda2[a++] = d2u_dty2[yy] * deltaEx[xx];
			duda[a] = du_dtsx[xx] * deltaEy[yy];
			d2uda2[a++] = d2u_dtsx2[xx] * deltaEy[yy];
			duda[a] = du_dtsy[yy] * deltaEx[xx];
			d2uda2[a++] = d2u_dtsy2[yy] * deltaEx[xx];
		}
		return I;
	}

	@Override
	public boolean evaluatesBackground()
	{
		return false;
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
		return true;
	}

	@Override
	public int getParametersPerPeak()
	{
		return 5;
	}

	@Override
	public void forEach(ValueProcedure procedure)
	{
		if (tB == 0 && nPeaks == 2)
		{
			// Specialised implementation without a background.
			// (This function is likely to be used to compute the Gaussian integral
			// without a background.)
			for (int y = 0; y < maxy; y++)
			{
				// Pre-compute
				final double tI_deltaEy0 = tI[0] * deltaEy[y];
				final double tI_deltaEy1 = tI[1] * deltaEy[y + maxy];

				for (int x = 0; x < maxx; x++)
				{
					procedure.execute(tI_deltaEy0 * deltaEx[x] + tI_deltaEy1 * deltaEx[x + maxx]);
				}
			}
		}
		else
		{
			super.forEach(procedure);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.GradientFunction#forEach(gdsc.smlm.function.GradientFunction.Gradient1Procedure)
	 */
	public void forEach(Gradient1Procedure procedure)
	{
		final double[] duda = new double[getNumberOfGradients()];
		for (int y = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++)
			{
				double I = tB;
				for (int n = 0, xx = x, yy = y, a = 0; n < nPeaks; n++, xx += maxx, yy += maxy)
				{
					duda[a] = deltaEx[xx] * deltaEy[yy];
					I += tI[n] * duda[a++];
					duda[a++] = du_dtx[xx] * deltaEy[yy];
					duda[a++] = du_dty[yy] * deltaEx[xx];
					duda[a++] = du_dtsx[xx] * deltaEy[yy];
					duda[a++] = du_dtsy[yy] * deltaEx[xx];
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
		for (int y = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++)
			{
				double I = tB;
				for (int n = 0, xx = x, yy = y, a = 0; n < nPeaks; n++, xx += maxx, yy += maxy)
				{
					duda[a] = deltaEx[xx] * deltaEy[yy];
					I += tI[n] * duda[a++];
					duda[a] = du_dtx[xx] * deltaEy[yy];
					d2uda2[a++] = d2u_dtx2[xx] * deltaEy[yy];
					duda[a] = du_dty[yy] * deltaEx[xx];
					d2uda2[a++] = d2u_dty2[yy] * deltaEx[xx];
					duda[a] = du_dtsx[xx] * deltaEy[yy];
					d2uda2[a++] = d2u_dtsx2[xx] * deltaEy[yy];
					duda[a] = du_dtsy[yy] * deltaEx[xx];
					d2uda2[a++] = d2u_dtsy2[yy] * deltaEx[xx];
				}
				procedure.execute(I, duda, d2uda2);
			}
		}
	}
}
