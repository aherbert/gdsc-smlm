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
package gdsc.smlm.function.cspline;

import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.Gradient2Procedure;
import gdsc.smlm.function.ValueProcedure;
import gdsc.smlm.results.PeakResult;

/**
 * Represent a cubic spline function for a single points.
 */
public class SingleCubicSplineFunction extends CubicSplineFunction
{
	/** The gradient indices for a single point [Background, Intensity, X, Y, Z] */
	private static final int[] gradientIndices = new int[] { 0, 1, 2, 3, 4 };

	/** The single target spline. */
	private final TargetSpline t;

	/**
	 * The working spline for the current evaluation.
	 * This is null if the point is outside the target range
	 */
	private TargetSpline working;

	/**
	 * Instantiates a new cubic spline function.
	 *
	 * @param splineData
	 *            the spline data
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data
	 * @throws IllegalArgumentException
	 *             If the function does not have an integer grid spacing from the origin
	 */
	public SingleCubicSplineFunction(CubicSplineData splineData, int maxx, int maxy) throws IllegalArgumentException
	{
		super(splineData, maxx, maxy);
		t = (splines[0][0].isSinglePrecision()) ? new FloatTargetSpline() : new DoubleTargetSpline();
	}

	/**
	 * Instantiates a new cubic spline function.
	 *
	 * @param splineData
	 *            the spline data
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data
	 * @param cx
	 *            the x centre of the spline data
	 * @param cy
	 *            the y centre of the spline data
	 * @param cz
	 *            the z centre of the spline data
	 * @param scale
	 *            the scale of the spline data
	 * @throws IllegalArgumentException
	 *             If the function does not have an integer grid spacing from the origin
	 */
	public SingleCubicSplineFunction(CubicSplineData splineData, int maxx, int maxy, double cx, double cy, double cz,
			int scale) throws IllegalArgumentException
	{
		super(splineData, maxx, maxy, cx, cy, cz, scale);
		t = (splines[0][0].isSinglePrecision()) ? new FloatTargetSpline() : new DoubleTargetSpline();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.cspline.CubicSplineFunction#getN()
	 */
	@Override
	public int getN()
	{
		return 1;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.GradientFunction#gradientIndices()
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
	 * @see gdsc.smlm.function.cspline.CubicSplineFunction#initialise(double[], int)
	 */
	@Override
	protected void initialise(double[] a, int order)
	{
		tB = a[PeakResult.BACKGROUND];
		working = (t.initialise(0, a[PeakResult.INTENSITY], a[PeakResult.X], a[PeakResult.Y], a[PeakResult.Z], order))
				? t : null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.ValueFunction#forEach(gdsc.smlm.function.ValueProcedure)
	 */
	@Override
	public void forEach(ValueProcedure procedure)
	{
		if (working == null)
		{
			// Special case as the spline does not overlap the target region
			for (int i = maxx * maxy; i-- > 0;)
				procedure.execute(tB);
		}
		else
		{
			working.reset();
			for (int y = 0; y < maxy; y++)
			{
				if (working.isNextYActive())
				{
					for (int x = 0; x < maxx; x++)
					{
						procedure.execute(tB + working.value(x));
					}
				}
				else
				{
					for (int x = 0; x < maxx; x++)
						procedure.execute(tB);
				}
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Function#forEach(gdsc.smlm.function.Gradient1Procedure)
	 */
	@Override
	public void forEach(Gradient1Procedure procedure)
	{
		final double[] duda = new double[getNumberOfGradients()];
		duda[0] = 1.0;

		if (working == null)
		{
			// Special case as the spline does not overlap the target region
			for (int i = maxx * maxy; i-- > 0;)
				procedure.execute(tB, duda);
		}
		else
		{
			working.reset();
			for (int y = 0; y < maxy; y++)
			{
				if (working.isNextYActive(duda))
				{
					for (int x = 0; x < maxx; x++)
					{
						// Because the call to value(...) occurs before passing the arguments
						procedure.execute(tB + working.value(x, duda), duda);
					}
				}
				else
				{
					for (int x = 0; x < maxx; x++)
						procedure.execute(tB, duda);
				}
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

		if (working == null)
		{
			// Special case as the spline does not overlap the target region
			for (int i = maxx * maxy; i-- > 0;)
				procedure.execute(tB, duda, d2uda2);
		}
		else
		{
			working.reset();
			for (int y = 0; y < maxy; y++)
			{
				if (working.isNextYActive(duda, d2uda2))
				{
					for (int x = 0; x < maxx; x++)
					{
						// Because the call to value(...) occurs before passing the arguments
						procedure.execute(tB + working.value(x, duda, d2uda2), duda, d2uda2);
					}
				}
				else
				{
					for (int x = 0; x < maxx; x++)
						procedure.execute(tB, duda, d2uda2);
				}
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.cspline.CubicSplineFunction#isNodeBoundary(int)
	 */
	@Override
	public boolean isNodeBoundary(int gradientIndex)
	{
		int parameterIndex = gradientIndices[gradientIndex];
		if (parameterIndex == BACKGROUND)
			return false;

		int dimension = (parameterIndex - 1) % PARAMETERS_PER_PEAK;
		if (dimension == 0)
			return false; // Signal

		int peak = getPeak(parameterIndex);
		if (peak == 0 && working != null)
			return working.isNodeBoundary(dimension - 1);
		return false;
	}
}
