package gdsc.smlm.function.cspline;

import java.util.Arrays;

import gdsc.core.math.interpolation.CubicSplinePosition;
import gdsc.core.math.interpolation.CustomTricubicFunction;
import gdsc.core.math.interpolation.CustomTricubicInterpolatingFunction;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.Gradient2Function;
import gdsc.smlm.function.Gradient2Procedure;
import gdsc.smlm.function.ValueProcedure;
import gdsc.smlm.results.PeakResult;

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
 * Represent a cubic spline function.
 */
public class CubicSplineFunction implements Gradient2Function
{
	private class TargetSpline
	{
		double tI, tX, tY, tZ;
		int ix1, ix2, iy1, iy2;
		int ix0, iy0;
		int xindex, yindex, zindex;
		double[] table1 = new double[64];
		double[] table2;
		double[] table3;
		double[] table6;
		boolean[] activeX = new boolean[maxx];

		public boolean initialise(double tI, double tX, double tY, double tZ, int order)
		{
			this.tI = tI;

			// Map z to a position in the spline
			// We want 0 to be in the centre. 
			// Note: Scale up the input parameter to the spline scale.
			tZ = cz + scale * tZ;

			if (tZ < 0 || tZ > uz)
			{
				return false;
			}

			// Map the scaled XY spline bounds to the target region (maxx x maxy)
			double x1 = lx + tX;
			double x2 = ux + tX;
			double y1 = ly + tY;
			double y2 = uy + tY;

			// Check if it is within the region 
			if (!(x2 > 0 && y2 > 0 && x1 < maxx && y1 < maxy))
			{
				return false;
			}

			// Convert to integer grid in the target region
			// We want the first integer that the function overlaps, i.e. can interpolate a value for
			ix1 = (int) Math.ceil(x1);
			iy1 = (int) Math.ceil(y1);
			// We want the last integer that the function overlaps
			ix2 = (int) Math.floor(x2);
			iy2 = (int) Math.floor(y2);

			// How far into the unscaled function is the first point
			double x = scale * (ix1 - x1);
			double y = scale * (iy1 - y1);

			// Get the spline index position for 0,0
			ix0 = (int) x - scale * ix1;
			iy0 = (int) y - scale * iy1;
			// Store the z spline node
			zindex = (int) tZ;

			// Create the interpolation table
			CubicSplinePosition sz = new CubicSplinePosition(tZ - zindex);
			CubicSplinePosition sx = new CubicSplinePosition(x - (int) x);
			CubicSplinePosition sy = new CubicSplinePosition(y - (int) y);

			table1 = CustomTricubicFunction.computePowerTable(sx, sy, sz);

			// Set the working flag for all x
			xindex = ix0;
			for (int i = 0; i < maxx; i++)
			{
				activeX[i] = xindex >= 0 && xindex <= maxSx;
				xindex += scale;
			}

			// Create tables for derivatives
			if (order == 1)
			{
				table2 = multiply(2, table2);
				table3 = multiply(3, table3);
			}
			else if (order == 2)
			{
				table2 = multiply(2, table2);
				table3 = multiply(3, table3);
				table6 = multiply(6, table6);
			}

			return true;
		}

		private double[] multiply(int n, double[] table)
		{
			if (table == null)
				table = new double[64];
			for (int i = 0; i < 64; i++)
				table[i] = table1[i] * n;
			return table;
		}

		public boolean isNextYActive()
		{
			boolean active = yindex >= 0 && yindex <= maxSy;
			// Increment yindex
			yindex += scale;
			// Reset xindex
			xindex = ix0;
			return active;
		}

		public double value(int x)
		{
			if (activeX[x])
			{
				double v = function.value(xindex, yindex, zindex, table1);
				xindex += scale;
				return v;
			}
			xindex += scale;
			return 0;
		}

		public void reset()
		{
			xindex = ix0;
			yindex = iy0;
		}
	}

	/** The cubic spline function. */
	private final CustomTricubicInterpolatingFunction function;

	private final int maxSx, maxSy;

	/** The scale. */
	private int scale = 1;

	private double cx, cy, cz;
	// The scaled bounds of the function with the centre at 0,0,0
	private double lx, ly, ux, uy;
	private final double uz;

	private final int maxx, maxy;

	private int n = 1;

	private int[] gradientIndices;

	private double tB;

	private TargetSpline[] t = new TargetSpline[0];

	private TargetSpline[] working;
	private int w;

	private TargetSpline[] workingY;
	private int wy;

	/**
	 * Instantiates a new cubic spline function.
	 *
	 * @param function
	 *            the function
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data
	 * @throws IllegalArgumentException
	 *             If the function does not have an integer grid spacing from the origin
	 */
	public CubicSplineFunction(CustomTricubicInterpolatingFunction function, int maxx, int maxy)
			throws IllegalArgumentException
	{
		checkFunction(function);
		this.function = function;
		maxSx = function.getMaxXSplinePosition();
		maxSy = function.getMaxYSplinePosition();
		this.maxx = (maxx < 1) ? 1 : maxx;
		this.maxy = (maxy < 1) ? 1 : maxy;
		// Centre in the middle, assuming the min is zero
		uz = function.getMaxZ();
		cx = (function.getMaxX() / 2.0);
		cy = (function.getMaxY() / 2.0);
		cz = (uz / 2.0);
		setScale(1);
	}

	private void checkFunction(CustomTricubicInterpolatingFunction function)
	{
		if (function == null || !function.isInteger())
			throw new IllegalArgumentException("Require integer grid spacing for the cubic spline");
		if (function.getMinX() != 0 || function.getMinY() != 0 || function.getMinZ() != 0)
			throw new IllegalArgumentException("Require integer grid spacing at the origin (0,0,0)");
	}

	/**
	 * Instantiates a new cubic spline function.
	 *
	 * @param function
	 *            the function
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data
	 * @param cx
	 *            the centre of the cubic spline x dimension
	 * @param cy
	 *            the centre of the cubic spline y dimension
	 * @param cz
	 *            the centre of the cubic spline z dimension
	 * @param scale
	 *            the scale
	 * @throws IllegalArgumentException
	 *             If the function does not have an integer grid spacing from the origin
	 */
	public CubicSplineFunction(CustomTricubicInterpolatingFunction function, int maxx, int maxy, double cx, double cy,
			double cz, int scale) throws IllegalArgumentException
	{
		checkFunction(function);
		this.function = function;
		maxSx = function.getMaxXSplinePosition();
		maxSy = function.getMaxYSplinePosition();
		this.maxx = (maxx < 1) ? 1 : maxx;
		this.maxy = (maxy < 1) ? 1 : maxy;
		this.cx = cx;
		this.cy = cy;
		this.cz = cz;
		uz = function.getMaxZ();
		setScale(scale);
	}

	private void updateFunctionBounds()
	{
		// Store the bounds of the cubic spline if it were positioned at 0,0
		lx = -cx / scale;
		ux = (function.getMaxX() - cx) / scale;
		ly = -cy / scale;
		uy = (function.getMaxY() - cy) / scale;
		// Check if the centre was within the function
		if (lx > 0 || ly > 0 || ux < 0 || uy < 0 || cz < 0 || cz > uz)
			throw new IllegalArgumentException("Require the centre within the cubic spline");
	}

	/**
	 * Gets the scale to map the cubic spline function to the integer grid. E.g. set a scale of 2 to render the spline
	 * at half its size.
	 *
	 * @return the scale
	 */
	public int getScale()
	{
		return scale;
	}

	/**
	 * Sets the scale to map the cubic spline function to the integer grid. E.g. set a scale of 2 to render the spline
	 * at half its size.
	 *
	 * @param scale
	 *            the new scale
	 * @throws IllegalArgumentException
	 *             If the scale is not strictly positive
	 */
	public void setScale(int scale) throws IllegalArgumentException
	{
		if (scale < 1)
			throw new IllegalArgumentException();
		this.scale = scale;
		updateFunctionBounds();
	}

	/**
	 * Gets the number of splines to draw.
	 *
	 * @return the number of splines to draw
	 */
	public int getN()
	{
		return n;
	}

	/**
	 * Sets the number of splines to draw.
	 *
	 * @param n
	 *            the new number of splines to draw
	 * @throws IllegalArgumentException
	 *             If the number is not strictly positive
	 */
	public void setN(int n) throws IllegalArgumentException
	{
		if (n < 1)
			throw new IllegalArgumentException();
		if (n != this.n)
			gradientIndices = null;
		this.n = n;
	}

	/**
	 * Gets the centre X.
	 *
	 * @return the centre X
	 */
	public double getCentreX()
	{
		return cx;
	}

	/**
	 * Sets the centre X.
	 *
	 * @param cx
	 *            the new centre X
	 */
	public void setCentreX(double cx)
	{
		this.cx = cx;
		updateFunctionBounds();
	}

	/**
	 * Gets the centre Y.
	 *
	 * @return the centre Y
	 */
	public double getCentreY()
	{
		return cy;
	}

	/**
	 * Sets the centre Y.
	 *
	 * @param cy
	 *            the new centre Y
	 */
	public void setCentreY(double cy)
	{
		this.cy = cy;
		updateFunctionBounds();
	}

	/**
	 * Gets the centre Z.
	 *
	 * @return the centre Z
	 */
	public double getCentreZ()
	{
		return cz;
	}

	/**
	 * Sets the centre Z.
	 *
	 * @param cz
	 *            the new centre Z
	 */
	public void setCentreZ(double cz)
	{
		this.cz = cz;
		updateFunctionBounds();
	}

	/**
	 * Gets the function.
	 *
	 * @return the function
	 */
	public CustomTricubicInterpolatingFunction getFunction()
	{
		return function;
	}

	public int size()
	{
		return maxx * maxy;
	}

	public int[] gradientIndices()
	{
		if (gradientIndices == null)
			gradientIndices = SimpleArrayUtils.newIntArray(getNumberOfGradients(), 0);
		return gradientIndices;
	}

	public int getNumberOfGradients()
	{
		return 1 + 4 * n;
	}

	public void initialise(double[] a)
	{
		initialise(a, 0);
	}

	public void initialise0(double[] a)
	{
		initialise(a, 0);
	}

	public void initialise1(double[] a)
	{
		initialise(a, 1);
	}

	public void initialise2(double[] a)
	{
		initialise(a, 2);
	}

	public void initialise(double[] a, int order)
	{
		tB = a[PeakResult.BACKGROUND];
		// Ensure we have enough room
		if (t.length < n)
		{
			int m = t.length;
			t = Arrays.copyOf(t, n); // Preserve memory space
			while (m < n)
				t[m++] = new TargetSpline();
			working = new TargetSpline[n];
			workingY = new TargetSpline[n];
		}
		// Convert the target parameters to spline offset tables
		w = 0;
		for (int i = 0, j = PeakResult.INTENSITY; i < n; i++)
		{
			double tI = a[j++];
			double tX = a[j++];
			double tY = a[j++];
			double tZ = a[j++];
			if (t[i].initialise(tI, tX, tY, tZ, order))
				working[w++] = t[i];
		}
	}

	public void forEach(ValueProcedure procedure)
	{
		for (int n = 0; n < w; n++)
			working[w].reset();

		for (int y = 0; y < maxy; y++)
		{
			// Get the working targets for this Y 
			wy = 0;
			for (int n = 0; n < w; n++)
			{
				if (working[w].isNextYActive())
					workingY[wy++] = working[w];
			}

			if (wy == 0)
			{
				for (int x = 0; x < maxx; x++)
					procedure.execute(tB);
			}
			else
			{
				for (int x = 0; x < maxx; x++)
				{
					double I = tB;
					for (int n = 0; n < wy; n++)
					{
						I += workingY[wy].value(x);
					}
					procedure.execute(I);
				}
			}
		}
	}

	public void forEach(Gradient1Procedure procedure)
	{

	}

	public void forEach(Gradient2Procedure procedure)
	{

	}
}