package gdsc.smlm.function.gaussian.erf;

import gdsc.smlm.function.ExtendedGradient2Procedure;
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
public class SingleNBFreeCircularErfGaussian2DFunction extends SingleFreeCircularErfGaussian2DFunction
{
	static final int[] gradientIndices;
	static
	{
		gradientIndices = createGradientIndices(1, new SingleNBFreeCircularErfGaussian2DFunction(1, 1));
	}

	/**
	 * Constructor.
	 *
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 */
	public SingleNBFreeCircularErfGaussian2DFunction(int maxx, int maxy)
	{
		super(maxx, maxy);
	}

	@Override
	public ErfGaussian2DFunction copy()
	{
		return new SingleNBFreeCircularErfGaussian2DFunction(maxx, maxy);
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
		duda[0] = deltaEx[x] * deltaEy[y];
		duda[1] = du_dtx[x] * deltaEy[y];
		duda[2] = du_dty[y] * deltaEx[x];
		duda[3] = du_dtsx[x] * deltaEy[y];
		duda[4] = du_dtsy[y] * deltaEx[x];

		return tB + tI * duda[0];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction#eval(int, double[], double[])
	 */
	public double eval(final int i, final double[] duda, final double[] d2uda2)
	{
		// Unpack the predictor into the dimensions
		final int y = i / maxx;
		final int x = i % maxx;

		// Return in order of Gaussian2DFunction.createGradientIndices().
		// Use pre-computed gradients
		duda[0] = deltaEx[x] * deltaEy[y];
		duda[1] = du_dtx[x] * deltaEy[y];
		duda[2] = du_dty[y] * deltaEx[x];
		duda[3] = du_dtsx[x] * deltaEy[y];
		duda[4] = du_dtsy[y] * deltaEx[x];
		d2uda2[0] = 0;
		d2uda2[1] = d2u_dtx2[x] * deltaEy[y];
		d2uda2[2] = d2u_dty2[y] * deltaEx[x];
		d2uda2[3] = d2u_dtsx2[x] * deltaEy[y];
		d2uda2[4] = d2u_dtsy2[y] * deltaEx[x];

		return tB + tI * duda[0];
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
	 * @see gdsc.smlm.function.GradientFunction#forEach(gdsc.smlm.function.GradientFunction.ValueProcedure)
	 */
	public void forEach(ValueProcedure procedure)
	{
		if (tB == 0)
		{
			// Specialised implementation without a background.
			// (This function is likely to be used to compute the Gaussian integral
			// without a background.)
			for (int y = 0; y < maxy; y++)
			{
				final double tI_deltaEy = tI * deltaEy[y];
				for (int x = 0; x < maxx; x++)
				{
					procedure.execute(tI_deltaEy * deltaEx[x]);
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
			final double du_dty = this.du_dty[y];
			final double deltaEy = this.deltaEy[y];
			final double du_dtsy = this.du_dtsy[y];
			for (int x = 0; x < maxx; x++)
			{
				duda[0] = deltaEx[x] * deltaEy;
				duda[1] = du_dtx[x] * deltaEy;
				duda[2] = du_dty * deltaEx[x];
				duda[3] = du_dtsx[x] * deltaEy;
				duda[4] = du_dtsy * deltaEx[x];
				procedure.execute(tB + tI * duda[0], duda);
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
			final double du_dty = this.du_dty[y];
			final double deltaEy = this.deltaEy[y];
			final double du_dtsy = this.du_dtsy[y];
			final double d2u_dty2 = this.d2u_dty2[y];
			final double d2u_dtsy2 = this.d2u_dtsy2[y];
			for (int x = 0; x < maxx; x++)
			{
				duda[0] = deltaEx[x] * deltaEy;
				duda[1] = du_dtx[x] * deltaEy;
				duda[2] = du_dty * deltaEx[x];
				duda[3] = du_dtsx[x] * deltaEy;
				duda[4] = du_dtsy * deltaEx[x];
				d2uda2[1] = d2u_dtx2[x] * deltaEy;
				d2uda2[2] = d2u_dty2 * deltaEx[x];
				d2uda2[3] = d2u_dtsx2[x] * deltaEy;
				d2uda2[4] = d2u_dtsy2 * deltaEx[x];
				procedure.execute(tB + tI * duda[0], duda, d2uda2);
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
		final int n = getNumberOfGradients();
		final double[] duda = new double[n];
		final double[] d2udadb = new double[n * n];
		final double[] du_dtsx_tI = new double[maxx];
		for (int x = 0; x < maxx; x++)
			du_dtsx_tI[x] = du_dtsx[x] / tI;
		for (int y = 0; y < maxy; y++)
		{
			final double du_dty = this.du_dty[y];
			final double du_dty_tI = du_dty / tI;
			final double deltaEy = this.deltaEy[y];
			final double du_dtsy = this.du_dtsy[y];
			final double du_dtsy_tI = du_dtsy / tI;
			final double d2u_dty2 = this.d2u_dty2[y];
			final double d2u_dtsy2 = this.d2u_dtsy2[y];
			final double d2deltaEy_dtsydy = this.d2deltaEy_dtsydy[y];
			for (int x = 0; x < maxx; x++)
			{
				duda[0] = deltaEx[x] * deltaEy;
				duda[1] = du_dtx[x] * deltaEy;
				duda[2] = du_dty * deltaEx[x];
				duda[3] = du_dtsx[x] * deltaEy;
				duda[4] = du_dtsy * deltaEx[x];

				// Compute all the partial second order derivatives

				// Signal,X
				d2udadb[1] = duda[1] / tI;
				// Signal,Y
				d2udadb[2] = duda[2] / tI;
				// Signal,X SD
				d2udadb[3] = duda[3] / tI;
				// Signal,Y SD
				d2udadb[4] = duda[4] / tI;

				// X,Signal
				d2udadb[5] = d2udadb[1];
				// X,X
				d2udadb[6] = d2u_dtx2[x] * deltaEy;
				// X,Y
				d2udadb[7] = du_dtx[x] * du_dty_tI;
				// X,X SD
				d2udadb[8] = deltaEy * d2deltaEx_dtsxdx[x];
				// X,Y SD
				d2udadb[9] = du_dtx[x] * du_dtsy_tI;

				// Y,Signal
				d2udadb[10] = d2udadb[2];
				// Y,X
				d2udadb[11] = d2udadb[7];
				// Y,Y
				d2udadb[12] = d2u_dty2 * deltaEx[x];
				// Y,X SD
				d2udadb[13] = du_dty * du_dtsx_tI[x];
				// Y,Y SD
				d2udadb[14] = deltaEx[x] * d2deltaEy_dtsydy;

				// X SD,Signal
				d2udadb[15] = d2udadb[3];
				// X SD,X
				d2udadb[16] = d2udadb[8];
				// X SD,Y
				d2udadb[17] = d2udadb[13];
				// X SD,X SD
				d2udadb[18] = d2u_dtsx2[x] * deltaEy;
				// X SD,Y SD
				d2udadb[19] = du_dtsy * du_dtsx_tI[x];

				// Y SD,Signal
				d2udadb[20] = d2udadb[4];
				// Y SD,X
				d2udadb[21] = d2udadb[9];
				// Y SD,Y
				d2udadb[22] = d2udadb[14];
				// Y SD,X SD
				d2udadb[23] = d2udadb[19];
				// Y SD,Y SD
				d2udadb[24] = d2u_dtsy2 * deltaEx[x];

				procedure.executeExtended(tB + tI * duda[0], duda, d2udadb);
			}
		}
	}	
}
