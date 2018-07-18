/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.function.gaussian.erf;

import gdsc.smlm.function.ExtendedGradient2Procedure;
import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.Gradient2Procedure;
import gdsc.smlm.function.ValueProcedure;

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

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.gaussian.erf.MultiErfGaussian2DFunction#eval(int, double[])
	 */
	@Override
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
	 * @see gdsc.smlm.function.gaussian.erf.MultiErfGaussian2DFunction#eval(int, double[], double[])
	 */
	@Override
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
		return true;
	}

	@Override
	public int getGradientParametersPerPeak()
	{
		return 5;
	}

	@Override
	public void forEach(ValueProcedure procedure)
	{
		if (tB == 0 && nPeaks == 2)
			// Specialised implementation without a background.
			// (This function is likely to be used to compute the Gaussian integral
			// without a background.)
			for (int y = 0; y < maxy; y++)
			{
				// Pre-compute
				final double tI_deltaEy0 = tI[0] * deltaEy[y];
				final double tI_deltaEy1 = tI[1] * deltaEy[y + maxy];

				for (int x = 0; x < maxx; x++)
					procedure.execute(tI_deltaEy0 * deltaEx[x] + tI_deltaEy1 * deltaEx[x + maxx]);
			}
		else
			super.forEach(procedure);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.GradientFunction#forEach(gdsc.smlm.function.GradientFunction.Gradient1Procedure)
	 */
	@Override
	public void forEach(Gradient1Procedure procedure)
	{
		final double[] duda = new double[getNumberOfGradients()];
		for (int y = 0; y < maxy; y++)
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

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.Gradient2Function#forEach(gdsc.smlm.function.Gradient2Procedure)
	 */
	@Override
	public void forEach(Gradient2Procedure procedure)
	{
		final double[] duda = new double[getNumberOfGradients()];
		final double[] d2uda2 = new double[getNumberOfGradients()];
		for (int y = 0; y < maxy; y++)
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

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.ExtendedGradient2Function#forEach(gdsc.smlm.function.ExtendedGradient2Procedure)
	 */
	@Override
	public void forEach(ExtendedGradient2Procedure procedure)
	{
		final int ng = getNumberOfGradients();
		final double[] duda = new double[ng];
		final double[] d2udadb = new double[ng * ng];
		final double[] du_dtsx_tI = new double[du_dtsx.length];
		for (int x = 0; x < maxx; x++)
			for (int n = 0, xx = x; n < nPeaks; n++, xx += maxx)
				du_dtsx_tI[xx] = du_dtsx[xx] / tI[n];
		final double[] du_dty_tI = new double[nPeaks];
		final double[] du_dtsy_tI = new double[nPeaks];
		for (int y = 0; y < maxy; y++)
		{
			for (int n = 0, yy = y; n < nPeaks; n++, yy += maxy)
			{
				du_dty_tI[n] = du_dty[yy] / tI[n];
				du_dtsy_tI[n] = du_dtsy[yy] / tI[n];
			}
			for (int x = 0; x < maxx; x++)
			{
				double I = tB;
				for (int n = 0, xx = x, yy = y, a = 0; n < nPeaks; n++, xx += maxx, yy += maxy)
				{
					duda[a] = deltaEx[xx] * deltaEy[yy];
					I += tI[n] * duda[a];
					duda[a + 1] = du_dtx[xx] * deltaEy[yy];
					duda[a + 2] = du_dty[yy] * deltaEx[xx];
					duda[a + 3] = du_dtsx[xx] * deltaEy[yy];
					duda[a + 4] = du_dtsy[yy] * deltaEx[xx];

					// Compute all the partial second order derivatives
					final double tI = this.tI[n];

					final int k = a * ng + a;
					// Signal,X
					d2udadb[k + 1] = duda[a + 1] / tI;
					// Signal,Y
					d2udadb[k + 2] = duda[a + 2] / tI;
					// Signal,X SD
					d2udadb[k + 3] = duda[a + 3] / tI;
					// Signal,Y SD
					d2udadb[k + 4] = duda[a + 4] / tI;

					a += 5;

					final int kk = k + ng;
					// X,Signal
					d2udadb[kk] = d2udadb[k + 1];
					// X,X
					d2udadb[kk + 1] = d2u_dtx2[xx] * deltaEy[yy];
					// X,Y
					d2udadb[kk + 2] = du_dtx[xx] * du_dty_tI[n];
					// X,X SD
					d2udadb[kk + 3] = deltaEy[yy] * d2deltaEx_dtsxdx[xx];
					// X,Y SD
					d2udadb[kk + 4] = du_dtx[xx] * du_dtsy_tI[n];

					final int kkk = kk + ng;
					// Y,Signal
					d2udadb[kkk] = d2udadb[k + 2];
					// Y,X
					d2udadb[kkk + 1] = d2udadb[kk + 2];
					// Y,Y
					d2udadb[kkk + 2] = d2u_dty2[yy] * deltaEx[xx];
					// Y,X SD
					d2udadb[kkk + 3] = du_dty[yy] * du_dtsx_tI[xx];
					// Y,Y SD
					d2udadb[kkk + 4] = deltaEx[xx] * d2deltaEy_dtsydy[yy];

					final int kkkk = kkk + ng;
					// X SD,Signal
					d2udadb[kkkk] = d2udadb[k + 3];
					// X SD,X
					d2udadb[kkkk + 1] = d2udadb[kk + 3];
					// X SD,Y
					d2udadb[kkkk + 2] = d2udadb[kkk + 3];
					// X SD,X SD
					d2udadb[kkkk + 3] = d2u_dtsx2[xx] * deltaEy[yy];
					// X SD,Y SD
					d2udadb[kkkk + 4] = du_dtsy[yy] * du_dtsx_tI[xx];

					final int kkkkk = kkkk + ng;
					// Y SD,Signal
					d2udadb[kkkkk] = d2udadb[k + 4];
					// Y SD,X
					d2udadb[kkkkk + 1] = d2udadb[kk + 4];
					// Y SD,Y
					d2udadb[kkkkk + 2] = d2udadb[kkk + 4];
					// Y SD,X SD
					d2udadb[kkkkk + 3] = d2udadb[kkkk + 4];
					// Y SD,Y SD
					d2udadb[kkkkk + 4] = d2u_dtsy2[yy] * deltaEx[xx];

				}
				procedure.executeExtended(I, duda, d2udadb);
			}
		}
	}
}
