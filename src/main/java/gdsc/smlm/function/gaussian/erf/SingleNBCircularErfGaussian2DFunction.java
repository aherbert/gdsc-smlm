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
public class SingleNBCircularErfGaussian2DFunction extends SingleCircularErfGaussian2DFunction
{
	static final int[] gradientIndices;
	static
	{
		gradientIndices = createGradientIndices(1, new SingleNBCircularErfGaussian2DFunction(1, 1));
	}

	/**
	 * Constructor.
	 *
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 */
	public SingleNBCircularErfGaussian2DFunction(int maxx, int maxy)
	{
		super(maxx, maxy);
	}

	@Override
	public ErfGaussian2DFunction copy()
	{
		return new SingleNBCircularErfGaussian2DFunction(maxx, maxy);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.gaussian.erf.SingleErfGaussian2DFunction#eval(int, double[])
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
		duda[3] = du_dtsx[x] * deltaEy[y] + du_dtsy[y] * deltaEx[x];

		return tB + tI * duda[0];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.gaussian.erf.SingleErfGaussian2DFunction#eval(int, double[], double[])
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
		duda[3] = du_dtsx[x] * deltaEy[y] + du_dtsy[y] * deltaEx[x];
		d2uda2[0] = 0;
		d2uda2[1] = d2u_dtx2[x] * deltaEy[y];
		d2uda2[2] = d2u_dty2[y] * deltaEx[x];
		// Working  example of this in GraspJ source code:
		// https://github.com/isman7/graspj/blob/master/graspj/src/main/java/eu/brede/graspj/opencl/src/functions/psfmodel_derivatives_sigma.cl
		//@formatter:off
		d2uda2[3] = d2u_dtsx2[x] * deltaEy[y] + 
				    d2u_dtsy2[y] * deltaEx[x] + 
				    2 * du_dtsx[x] * du_dtsy[y] / tI;
		//@formatter:on

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
		return 4;
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
				duda[3] = du_dtsx[x] * deltaEy + du_dtsy * deltaEx[x];
				//invalidGradients(duda);
				procedure.execute(tB + tI * duda[0], duda);
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
		for (int y = 0; y < maxy; y++)
		{
			final double du_dty = this.du_dty[y];
			final double deltaEy = this.deltaEy[y];
			final double du_dtsy = this.du_dtsy[y];
			final double two_du_dtsy_tI = 2 * this.du_dtsy[y] / tI;
			final double d2u_dty2 = this.d2u_dty2[y];
			final double d2u_dtsy2 = this.d2u_dtsy2[y];
			for (int x = 0; x < maxx; x++)
			{
				duda[0] = deltaEx[x] * deltaEy;
				duda[1] = du_dtx[x] * deltaEy;
				duda[2] = du_dty * deltaEx[x];
				duda[3] = du_dtsx[x] * deltaEy + du_dtsy * deltaEx[x];
				d2uda2[1] = d2u_dtx2[x] * deltaEy;
				d2uda2[2] = d2u_dty2 * deltaEx[x];
				//@formatter:off
				d2uda2[3] = d2u_dtsx2[x] * deltaEy + 
					        d2u_dtsy2 * deltaEx[x] + 
					        du_dtsx[x] * two_du_dtsy_tI;
				//@formatter:on
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
			final double two_du_dtsy_tI = 2 * this.du_dtsy[y] / tI;
			final double d2u_dty2 = this.d2u_dty2[y];
			final double d2u_dtsy2 = this.d2u_dtsy2[y];
			final double d2deltaEy_dtsydy = this.d2deltaEy_dtsydy[y];
			for (int x = 0; x < maxx; x++)
			{
				duda[0] = deltaEx[x] * deltaEy;
				duda[1] = du_dtx[x] * deltaEy;
				duda[2] = du_dty * deltaEx[x];
				duda[3] = du_dtsx[x] * deltaEy + du_dtsy * deltaEx[x];

				// Compute all the partial second order derivatives

				// Signal,X
				d2udadb[1] = duda[1] / tI;
				// Signal,Y
				d2udadb[2] = duda[2] / tI;
				// Signal,X SD
				d2udadb[3] = duda[3] / tI;

				// X,Signal
				d2udadb[4] = d2udadb[1];
				// X,X
				d2udadb[5] = d2u_dtx2[x] * deltaEy;
				// X,Y
				d2udadb[6] = du_dtx[x] * du_dty_tI;
				// X,X SD
				d2udadb[7] = deltaEy * d2deltaEx_dtsxdx[x] + du_dtx[x] * du_dtsy_tI;

				// Y,Signal
				d2udadb[8] = d2udadb[2];
				// Y,X
				d2udadb[9] = d2udadb[6];
				// Y,Y
				d2udadb[10] = d2u_dty2 * deltaEx[x];
				// Y,X SD
				d2udadb[11] = du_dty * du_dtsx_tI[x] + deltaEx[x] * d2deltaEy_dtsydy;

				// X SD,Signal
				d2udadb[12] = d2udadb[3];
				// X SD,X
				d2udadb[13] = d2udadb[7];
				// X SD,Y
				d2udadb[14] = d2udadb[11];
				// X SD,X SD
				//@formatter:off
				d2udadb[15] = d2u_dtsx2[x] * deltaEy + 
         				      d2u_dtsy2 * deltaEx[x] + 
         				      du_dtsx[x] * two_du_dtsy_tI;
				//@formatter:on

				procedure.executeExtended(tB + tI * duda[0], duda, d2udadb);
			}
		}
	}
}
