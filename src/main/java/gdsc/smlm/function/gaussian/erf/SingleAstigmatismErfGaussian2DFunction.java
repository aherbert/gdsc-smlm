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
import gdsc.smlm.function.gaussian.AstigmatismZModel;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;

/**
 * Evaluates a 2-dimensional Gaussian function for a single peak.
 */
public class SingleAstigmatismErfGaussian2DFunction extends SingleFreeCircularErfGaussian2DFunction
{
	/** The gradient indices. */
	static final int[] gradientIndices;
	static
	{
		gradientIndices = createGradientIndices(1, new SingleAstigmatismErfGaussian2DFunction(1, 1, null));
	}

	/** The z model. */
	protected final AstigmatismZModel zModel;

	// Required for the z-depth gradients

	/** The x|z pre-factors for first-order partial derivatives. */
	protected double dtsx_dtz;
	/** The x|z pre-factors for first-order partial derivatives. */
	protected double d2tsx_dtz2;
	/** The y|z pre-factors for second-order partial derivatives. */
	protected double dtsy_dtz;
	/** The y|z pre-factors for second-order partial derivatives. */
	protected double d2tsy_dtz2;

	/**
	 * Constructor.
	 *
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param zModel
	 *            the z model
	 */
	public SingleAstigmatismErfGaussian2DFunction(int maxx, int maxy, AstigmatismZModel zModel)
	{
		super(maxx, maxy);
		this.zModel = zModel;
	}

	@Override
	public ErfGaussian2DFunction copy()
	{
		return new SingleAstigmatismErfGaussian2DFunction(maxx, maxy, zModel);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction#initialise0(double[])
	 */
	@Override
	public void initialise0(double[] a)
	{
		tB = a[Gaussian2DFunction.BACKGROUND];
		tI = a[Gaussian2DFunction.SIGNAL];
		// Pre-compute the offset by 0.5
		final double tx = a[Gaussian2DFunction.X_POSITION] + 0.5;
		final double ty = a[Gaussian2DFunction.Y_POSITION] + 0.5;
		final double tz = a[Gaussian2DFunction.Z_POSITION];

		final double sx = zModel.getSx(tz);
		final double sy = zModel.getSy(tz);
		createDeltaETable(ONE_OVER_ROOT2 / sx, deltaEx, tx);
		createDeltaETable(ONE_OVER_ROOT2 / sy, deltaEy, ty);
	}

	@Override
	public double integral(double[] a)
	{
		final double tB = a[Gaussian2DFunction.BACKGROUND];
		final double tI = a[Gaussian2DFunction.SIGNAL];
		// Pre-compute the offset by 0.5
		final double tx = a[Gaussian2DFunction.X_POSITION] + 0.5;
		final double ty = a[Gaussian2DFunction.Y_POSITION] + 0.5;
		final double tz = a[Gaussian2DFunction.Z_POSITION];

		final double sx = zModel.getSx(tz);
		final double sy = zModel.getSy(tz);

		return tB * size() + tI * compute1DIntegral(ONE_OVER_ROOT2 / sx, maxx, tx) *
				compute1DIntegral(ONE_OVER_ROOT2 / sy, maxy, ty);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction#initialise1(double[])
	 */
	@Override
	public void initialise1(double[] a)
	{
		create1Arrays();
		tB = a[Gaussian2DFunction.BACKGROUND];
		tI = a[Gaussian2DFunction.SIGNAL];
		// Pre-compute the offset by 0.5
		final double tx = a[Gaussian2DFunction.X_POSITION] + 0.5;
		final double ty = a[Gaussian2DFunction.Y_POSITION] + 0.5;
		final double tz = a[Gaussian2DFunction.Z_POSITION];

		// We can pre-compute part of the derivatives for position and sd in arrays
		// since the Gaussian is XY separable
		final double[] ds_dz = new double[1];
		final double sx = zModel.getSx(tz, ds_dz);
		dtsx_dtz = ds_dz[0];
		final double sy = zModel.getSy(tz, ds_dz);
		dtsy_dtz = ds_dz[0];
		createFirstOrderTables(tI, deltaEx, du_dtx, du_dtsx, tx, sx);
		createFirstOrderTables(tI, deltaEy, du_dty, du_dtsy, ty, sy);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction#initialise2(double[])
	 */
	@Override
	public void initialise2(double[] a)
	{
		create2Arrays();
		tB = a[Gaussian2DFunction.BACKGROUND];
		tI = a[Gaussian2DFunction.SIGNAL];
		// Pre-compute the offset by 0.5
		final double tx = a[Gaussian2DFunction.X_POSITION] + 0.5;
		final double ty = a[Gaussian2DFunction.Y_POSITION] + 0.5;
		final double tz = a[Gaussian2DFunction.Z_POSITION];

		// We can pre-compute part of the derivatives for position and sd in arrays
		// since the Gaussian is XY separable
		final double[] ds_dz = new double[2];
		final double sx = zModel.getSx2(tz, ds_dz);
		dtsx_dtz = ds_dz[0];
		d2tsx_dtz2 = ds_dz[1];
		final double sy = zModel.getSy2(tz, ds_dz);
		dtsy_dtz = ds_dz[0];
		d2tsy_dtz2 = ds_dz[1];
		createSecondOrderTables(tI, deltaEx, du_dtx, du_dtsx, d2u_dtx2, d2u_dtsx2, tx, sx);
		createSecondOrderTables(tI, deltaEy, du_dty, du_dtsy, d2u_dty2, d2u_dtsy2, ty, sy);
	}

	@Override
	public void initialiseExtended2(double[] a)
	{
		createEx2Arrays();
		tB = a[Gaussian2DFunction.BACKGROUND];
		tI = a[Gaussian2DFunction.SIGNAL];
		// Pre-compute the offset by 0.5
		final double tx = a[Gaussian2DFunction.X_POSITION] + 0.5;
		final double ty = a[Gaussian2DFunction.Y_POSITION] + 0.5;
		final double tz = a[Gaussian2DFunction.Z_POSITION];

		// We can pre-compute part of the derivatives for position and sd in arrays
		// since the Gaussian is XY separable
		final double[] ds_dz = new double[2];
		final double sx = zModel.getSx2(tz, ds_dz);
		dtsx_dtz = ds_dz[0];
		d2tsx_dtz2 = ds_dz[1];
		final double sy = zModel.getSy2(tz, ds_dz);
		dtsy_dtz = ds_dz[0];
		d2tsy_dtz2 = ds_dz[1];
		createExSecondOrderTables(tI, deltaEx, du_dtx, du_dtsx, d2u_dtx2, d2u_dtsx2, d2deltaEx_dtsxdx, tx, sx);
		createExSecondOrderTables(tI, deltaEy, du_dty, du_dtsy, d2u_dty2, d2u_dtsy2, d2deltaEy_dtsydy, ty, sy);
		// Pre-apply the gradient mapping from width to z
		for (int i = 0; i < d2deltaEx_dtsxdx.length; i++)
			d2deltaEx_dtsxdx[i] *= dtsx_dtz;
		for (int i = 0; i < d2deltaEy_dtsydy.length; i++)
			d2deltaEy_dtsydy[i] *= dtsy_dtz;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.gaussian.erf.SingleErfGaussian2DFunction#eval(int, double[])
	 */
	@Override
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
		duda[4] = du_dtsx[x] * deltaEy[y] * dtsx_dtz + du_dtsy[y] * deltaEx[x] * dtsy_dtz;

		return tB + tI * duda[1];
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.gaussian.erf.SingleErfGaussian2DFunction#eval(int, double[], double[])
	 */
	@Override
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
		duda[2] = du_dtx[x] * deltaEy[y];
		duda[3] = du_dty[y] * deltaEx[x];
		duda[4] = du_dtsx[x] * deltaEy[y] * dtsx_dtz + du_dtsy[y] * deltaEx[x] * dtsy_dtz;
		d2uda2[0] = 0;
		d2uda2[1] = 0;
		d2uda2[2] = d2u_dtx2[x] * deltaEy[y];
		d2uda2[3] = d2u_dty2[y] * deltaEx[x];
		//@formatter:off
		d2uda2[4] =
				d2u_dtsx2[x] * deltaEy[y] * dtsx_dtz * dtsx_dtz +
				du_dsx * d2tsx_dtz2 +
				d2u_dtsy2[y] * deltaEx[x] * dtsy_dtz * dtsy_dtz +
				du_dsy * d2tsy_dtz2 +
				// Add the equivalent term we add in the circular version.
				// Note: this is not in the Smith, et al (2010) paper but is
				// in the GraspJ source code and it works in JUnit tests.
				2 * du_dtsx[x] * dtsx_dtz * du_dtsy[y] * dtsy_dtz / tI;
		//@formatter:on

		return tB + tI * duda[1];
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
	 * @see gdsc.fitting.function.NonLinearFunction#gradientIndices()
	 */
	@Override
	public int[] gradientIndices()
	{
		return gradientIndices;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.GradientFunction#getNumberOfGradients()
	 */
	@Override
	public int getNumberOfGradients()
	{
		return 5;
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
				duda[2] = du_dtx[x] * deltaEy;
				duda[3] = du_dty * deltaEx[x];
				duda[4] = du_dtsx[x] * deltaEy_by_dtsx_dtz + du_dtsy_by_dtsy_dtz * deltaEx[x];
				procedure.execute(tB + tI * duda[1], duda);
			}
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
				duda[2] = du_dtx[x] * deltaEy;
				duda[3] = du_dty * deltaEx[x];
				duda[4] = du_dsx * dtsx_dtz + du_dsy * dtsy_dtz;
				d2uda2[2] = d2u_dtx2[x] * deltaEy;
				d2uda2[3] = d2u_dty2 * deltaEx[x];
				//@formatter:off
				d2uda2[4] =
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

				procedure.execute(tB + tI * duda[1], duda, d2uda2);
			}
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
		final int n = getNumberOfGradients();
		final double[] duda = new double[n];
		final double[] d2udadb = new double[n * n];
		duda[0] = 1.0;
		final double dtsx_dtz_2 = dtsx_dtz * dtsx_dtz;
		final double dtsy_dtz_2 = dtsy_dtz * dtsy_dtz;
		final double two_dtsx_dtz_by_dtsy_dtz_tI = 2 * dtsx_dtz * dtsy_dtz / tI;
		final double dtsx_dtz_tI = dtsx_dtz / tI;
		final double dtsy_dtz_tI = dtsy_dtz / tI;
		for (int y = 0; y < maxy; y++)
		{
			final double du_dty = this.du_dty[y];
			final double du_dty_tI = du_dty / tI;
			final double deltaEy = this.deltaEy[y];
			final double du_dtsy = this.du_dtsy[y];
			final double d2u_dty2 = this.d2u_dty2[y];
			final double du_dtsy_by_dtsy_dtz_tI = du_dtsy * dtsy_dtz_tI;
			final double du_dty_by_dtsx_dtz_tI = du_dty * dtsx_dtz_tI;
			final double d2deltaEy_dtsydy = this.d2deltaEy_dtsydy[y];
			final double deltaEy_by_dtsx_dtz_2 = deltaEy * dtsx_dtz_2;
			final double d2u_dtsy2_by_dtsy_dtz_2 = this.d2u_dtsy2[y] * dtsy_dtz_2;
			final double two_dtsx_dtz_by_du_dtsy_by_dtsy_dtz_tI = two_dtsx_dtz_by_dtsy_dtz_tI * du_dtsy;
			for (int x = 0; x < maxx; x++)
			{
				final double du_dsx = du_dtsx[x] * deltaEy;
				final double du_dsy = du_dtsy * deltaEx[x];

				duda[1] = deltaEx[x] * deltaEy;
				duda[2] = du_dtx[x] * deltaEy;
				duda[3] = du_dty * deltaEx[x];
				duda[4] = du_dsx * dtsx_dtz + du_dsy * dtsy_dtz;

				// Compute all the partial second order derivatives

				// Background are all 0

				// Signal,X
				d2udadb[7] = duda[2] / tI;
				// Signal,Y
				d2udadb[8] = duda[3] / tI;
				// Signal,Z
				d2udadb[9] = duda[4] / tI;

				// X,Signal
				d2udadb[11] = d2udadb[7];
				// X,X
				d2udadb[12] = d2u_dtx2[x] * deltaEy;
				// X,Y
				d2udadb[13] = du_dtx[x] * du_dty_tI;
				// X,Z
				d2udadb[14] = deltaEy * d2deltaEx_dtsxdx[x] + du_dtx[x] * du_dtsy_by_dtsy_dtz_tI;

				// Y,Signal
				d2udadb[16] = d2udadb[8];
				// Y,X
				d2udadb[17] = d2udadb[13];
				// Y,Y
				d2udadb[18] = d2u_dty2 * deltaEx[x];
				// Y,Z
				d2udadb[19] = du_dtsx[x] * du_dty_by_dtsx_dtz_tI + deltaEx[x] * d2deltaEy_dtsydy;

				// Z,Signal
				d2udadb[21] = d2udadb[9];
				// Z,X
				d2udadb[22] = d2udadb[14];
				// Z,Y
				d2udadb[23] = d2udadb[19];
				// Z,Z
				//@formatter:off
				d2udadb[24] =
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

				procedure.executeExtended(tB + tI * duda[1], duda, d2udadb);
			}
		}
	}
}
