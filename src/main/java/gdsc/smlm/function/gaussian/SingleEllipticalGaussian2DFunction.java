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
 * Evaluates an 2-dimensional elliptical Gaussian function for a single peak.
 * <p>
 * The single parameter x in the {@link #eval(int, double[])} function is assumed to be a linear index into
 * 2-dimensional data. The dimensions of the data must be specified to allow unpacking to coordinates.
 * <p>
 * Data should be packed in descending dimension order, e.g. Y,X : Index for [x,y] = MaxX*y + x.
 */
public class SingleEllipticalGaussian2DFunction extends Gaussian2DFunction
{
	private static final int[] gradientIndices;
	static
	{
		gradientIndices = createGradientIndices(1, new SingleEllipticalGaussian2DFunction(1, 1));
	}

	protected double background;
	protected double x0pos;
	protected double x1pos;

	protected boolean zeroAngle;
	protected double n;
	protected double height;
	protected double aa;
	protected double bb;
	protected double cc;
	protected double aa2;
	protected double bb2;
	protected double cc2;
	protected double nx;
	protected double ax;
	protected double bx;
	protected double cx;
	protected double ny;
	protected double ay;
	protected double by;
	protected double cy;

	/**
	 * Constructor
	 * 
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 */
	public SingleEllipticalGaussian2DFunction(int maxx, int maxy)
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
		return new SingleEllipticalGaussian2DFunction(maxx, maxy);
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

		// Precalculate multiplication factors
		final double theta = a[ANGLE];
		final double sx = a[X_SD];
		final double sy = a[Y_SD];
		final double sx2 = sx * sx;
		final double sy2 = sy * sy;
		final double sx3 = sx2 * sx;
		final double sy3 = sy2 * sy;

		n = ONE_OVER_TWO_PI / (sx * sy);
		height = a[SIGNAL] * n;

		// All prefactors are negated since the Gaussian uses the exponential to the negative:
		// (A/2*pi*sx*sy) * exp( -( a(x-x0)^2 + 2b(x-x0)(y-y0) + c(y-y0)^2 ) )

		if (theta == 0)
		{
			zeroAngle = true;

			// cosSqt = 1
			// sinSqt = 0
			// sincost = 0
			// sin2t = 0
			// cos2t = 1

			aa = -0.5 / sx2;
			bb = 0;
			cc = -0.5 / sy2;

			// For the angle gradient
			aa2 = 0;
			bb2 = -0.5 * (-1.0 / sx2 + 1.0 / sy2);
			cc2 = 0;

			// For the x-width gradient
			nx = -1 / sx;
			ax = 1.0 / sx3;
			bx = 0;
			cx = 0;

			// For the y-width gradient
			ny = -1 / sy;
			ay = 0;
			by = 0;
			cy = 1.0 / sy3;
		}
		else
		{
			zeroAngle = false;

			final double cosSqt = Math.cos(theta) * Math.cos(theta);
			final double sinSqt = Math.sin(theta) * Math.sin(theta);
			final double sincost = Math.sin(theta) * Math.cos(theta);
			final double sin2t = Math.sin(2 * theta);
			final double cos2t = Math.cos(2 * theta);

			aa = -0.5 * (cosSqt / sx2 + sinSqt / sy2);
			bb = -0.25 * (-sin2t / sx2 + sin2t / sy2);
			cc = -0.5 * (sinSqt / sx2 + cosSqt / sy2);

			// For the angle gradient
			aa2 = -(-sincost / sx2 + sincost / sy2);
			bb2 = -0.5 * (-cos2t / sx2 + cos2t / sy2);
			cc2 = -(sincost / sx2 - sincost / sy2);

			// For the x-width gradient
			nx = -1.0 / sx;
			ax = cosSqt / sx3;
			bx = -0.5 * sin2t / sx3;
			cx = sinSqt / sx3;

			// For the y-width gradient
			ny = -1.0 / sy;
			ay = sinSqt / sy3;
			by = 0.5 * sin2t / sy3;
			cy = cosSqt / sy3;
		}
	}

	/**
	 * Evaluates an 2-dimensional elliptical Gaussian function for a single peak.
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
		final double dx2 = dx * dx;
		final double dxy = dx * dy;
		final double dy2 = dy * dy;

		// Calculate gradients
		if (zeroAngle)
		{
			final double exp = FastMath.exp(aa * dx2 + cc * dy2);
			dy_da[1] = n * exp;
			final double y = height * exp;

			dy_da[2] = y * (-2.0 * aa * dx);
			dy_da[3] = y * (-2.0 * cc * dy);

			dy_da[4] = y * (nx + ax * dx2);
			dy_da[5] = y * (ny + cy * dy2);

			dy_da[6] = y * (bb2 * dxy);

			return y;
		}
		else
		{
			final double exp = FastMath.exp(aa * dx2 + bb * dxy + cc * dy2);
			dy_da[1] = n * exp;
			final double y = height * exp;

			dy_da[2] = y * (-2.0 * aa * dx - bb * dy);
			dy_da[3] = y * (-2.0 * cc * dy - bb * dx);

			dy_da[4] = y * (nx + ax * dx2 + bx * dxy + cx * dy2);
			dy_da[5] = y * (ny + ay * dx2 + by * dxy + cy * dy2);

			dy_da[6] = y * (aa2 * dx2 + bb2 * dxy + cc2 * dy2);

			return y;
		}
	}

	/**
	 * Evaluates an 2-dimensional elliptical Gaussian function for a single peak.
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

		if (zeroAngle)
			return background + height * FastMath.exp(aa * dx * dx + cc * dy * dy);
		else
			return background + height * FastMath.exp(aa * dx * dx + bb * dx * dy + cc * dy * dy);
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
	public boolean evaluatesAngle()
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
		return 6;
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
