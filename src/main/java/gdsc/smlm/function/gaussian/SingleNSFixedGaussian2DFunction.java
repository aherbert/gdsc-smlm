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
package gdsc.smlm.function.gaussian;

import org.apache.commons.math3.util.FastMath;

/**
 * Evaluates an 2-dimensional Gaussian function for a single peak.
 * <p>
 * The single parameter x in the {@link #eval(int, double[])} function is assumed to be a linear index into
 * 2-dimensional
 * data. The dimensions of the data must be specified to allow unpacking to coordinates.
 * <p>
 * Data should be packed in descending dimension order, e.g. Y,X : Index for [x,y] = MaxX*y + x.
 */
public class SingleNSFixedGaussian2DFunction extends Gaussian2DFunction
{
	private static final int[] gradientIndices;
	static
	{
		gradientIndices = createGradientIndices(1, new SingleNSFixedGaussian2DFunction(1, 1));
	}

	protected double width;

	protected double background;
	protected double x0pos;
	protected double x1pos;

	protected double n;
	protected double height;
	protected double aa;
	protected double aa2;

	/**
	 * Constructor
	 *
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 */
	public SingleNSFixedGaussian2DFunction(int maxx, int maxy)
	{
		super(maxx, maxy);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.gaussian.Gaussian2DFunction#copy()
	 */
	@Override
	public Gaussian2DFunction copy()
	{
		return new SingleNSFixedGaussian2DFunction(maxx, maxy);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.fitting.function.NonLinearFunction#initialise(double[])
	 */
	@Override
	public void initialise(double[] a)
	{
		background = a[BACKGROUND];
		x0pos = a[X_POSITION];
		x1pos = a[Y_POSITION];
		width = a[X_SD];

		final double sx = a[X_SD];
		final double sx2 = sx * sx;

		n = ONE_OVER_TWO_PI / sx2;
		height = a[SIGNAL] * n;

		// All prefactors are negated since the Gaussian uses the exponential to the negative:
		// A * exp( -( a(x-x0)^2 + 2b(x-x0)(y-y0) + c(y-y0)^2 ) )

		aa = -0.5 / sx2;
		aa2 = -2.0 * aa;
	}

	/**
	 * Evaluates an 2-dimensional fixed circular Gaussian function for a single peak.
	 * <p>
	 * {@inheritDoc}
	 *
	 * @see gdsc.smlm.function.gaussian.Gaussian2DFunction#eval(int, double[])
	 */
	@Override
	public double eval(final int x, final double[] dyda)
	{
		// First parameter is the background level
		dyda[0] = 1.0; // Gradient for a constant background is 1

		// Unpack the predictor into the dimensions
		final int x1 = x / maxx;
		final int x0 = x % maxx;

		return background + gaussian(x0, x1, dyda);
	}

	private double gaussian(final int x0, final int x1, final double[] dy_da)
	{
		final double dx = x0 - x0pos;
		final double dy = x1 - x1pos;

		// Calculate gradients
		final double y = height * FastMath.exp(aa * (dx * dx + dy * dy));
		final double yaa2 = y * aa2;
		dy_da[1] = yaa2 * dx;
		dy_da[2] = yaa2 * dy;

		return y;
	}

	/**
	 * Evaluates an 2-dimensional fixed circular Gaussian function for a single peak.
	 * <p>
	 * {@inheritDoc}
	 *
	 * @see gdsc.smlm.function.gaussian.Gaussian2DFunction#eval(int)
	 */
	@Override
	public double eval(final int x)
	{
		// Unpack the predictor into the dimensions
		final int x1 = x / maxx;
		final int x0 = x % maxx;

		final double dx = x0 - x0pos;
		final double dy = x1 - x1pos;

		return background + height * FastMath.exp(aa * (dx * dx + dy * dy));
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
		return false;
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
		return 2;
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
}
