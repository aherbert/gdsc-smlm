package gdsc.smlm.function.cspline;

import gdsc.core.math.interpolation.CustomTricubicFunction;
import gdsc.smlm.function.Gradient2Function;

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
public abstract class CubicSplineFunction implements Gradient2Function
{
	/**
	 * Internal class to control visiting the correct cubic spline node for each [x][y] index in the
	 * target range [0 <= x < maxx], [0 <= y < maxy]
	 */
	protected abstract class TargetSpline
	{
		int ix1, iy1;
		int ix0, iy0;
		int index;
		int yindex;
		CustomTricubicFunction[] xySplines;
		boolean[] activeX = new boolean[maxx];

		public boolean initialise(double tI, double tX, double tY, double tZ, int order)
		{
			// Map z to a position in the spline
			// We want 0 to be in the centre. 
			// Note: Scale up the input parameter to the spline scale.
			double z = cz + scale * tZ;

			if (z < 0 || z > maxSz)
			{
				return false;
			}

			// Shift the scaled XY spline bounds by the target centre 
			double x1 = lx + tX;
			double x2 = ux + tX;
			double y1 = ly + tY;
			double y2 = uy + tY;

			// Check if it is within the region 
			if (!(x2 > 0 && y2 > 0 && x1 < maxx && y1 < maxy))
			{
				return false;
			}

			// Convert the lower bounds to integer grid in the target region, 
			// i.e. we sample the region at x=0,1,2,...
			// We want the first integer that the function overlaps, 
			// i.e. can interpolate a value for so it must be above the lower bounds
			ix1 = (int) Math.ceil(x1);
			iy1 = (int) Math.ceil(y1);

			// How far into the unscaled function is the first point.
			// i.e. x=1 may be 0.6 above the scaled lower bound (0.4) but that would require
			// the first sample to be taken at spline[1] @ 0.2 if the scale is 2.
			double x = scale * (ix1 - x1);
			double y = scale * (iy1 - y1);

			// This is the first index for the spline sample
			int ix = (int) x;
			int iy = (int) y;
			int iz = (int) z;

			// Get the spline index position for 0,0 offset by the scale (for pre-increment loops)
			ix0 = ix - scale * ix1 - scale;
			iy0 = iy - scale * iy1 - scale;
			// Store the xy splines for the z position
			xySplines = splines[iz];

			// Set the working flag for all x
			for (int i = 0, xindex = ix0; i < maxx; i++)
			{
				xindex += scale;
				activeX[i] = xindex >= 0 && xindex < maxSx;
			}

			computePowerTable(x - ix, y - iy, z - iz, order);

			return true;
		}

		public void reset()
		{
			yindex = iy0;
		}

		public boolean isNextYActive()
		{
			// pre-increment yindex
			yindex += scale;
			if (yindex >= 0 && yindex < maxSy)
			{
				// The y-index is inside the XY spline data
				// Reset the index
				index = yindex * maxSx + ix0;
			}
			return false;
		}

		abstract public void computePowerTable(double x, double y, double z, int order);

		public double value(int x)
		{
			index += scale; // pre-increment
			return (activeX[x]) ? computeValue(xySplines[index]) : 0;
		}

		abstract public double computeValue(CustomTricubicFunction customTricubicFunction);
	}

	/**
	 * Double precision computation of the target spline
	 */
	protected class DoubleTargetSpline extends TargetSpline
	{
		double[] table1 = new double[64];
		double[] table2;
		double[] table3;
		double[] table6;

		@Override
		public void computePowerTable(double x, double y, double z, int order)
		{
			table1 = CustomTricubicFunction.computePowerTable(x, y, z);

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
		}

		private double[] multiply(int n, double[] table)
		{
			if (table == null)
				table = new double[64];
			for (int i = 0; i < 64; i++)
				table[i] = table1[i] * n;
			return table;
		}

		public double computeValue(CustomTricubicFunction customTricubicFunction)
		{
			return customTricubicFunction.value(table1);
		}
	}

	/**
	 * Double precision computation of the target spline
	 */
	protected class FloatTargetSpline extends TargetSpline
	{
		float[] table1 = new float[64];
		float[] table2;
		float[] table3;
		float[] table6;

		@Override
		public void computePowerTable(double x, double y, double z, int order)
		{
			table1 = CustomTricubicFunction.computeFloatPowerTable(x, y, z);

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
		}

		private float[] multiply(int n, float[] table)
		{
			if (table == null)
				table = new float[64];
			for (int i = 0; i < 64; i++)
				table[i] = table1[i] * n;
			return table;
		}

		public double computeValue(CustomTricubicFunction customTricubicFunction)
		{
			return customTricubicFunction.value(table1);
		}
	}

	/** The scale to map the target range (maxx * maxy) to the spline. */
	protected int scale = 1;

	// The centre of the spline (unscaled)
	protected double cx, cy, cz;

	// The scaled bounds of the function with the centre at 0,0,0
	protected double lx, ly, ux, uy;

	// Max size of spline data
	protected final int maxSx, maxSy, maxSz;

	// The tricubic spline packed as Z * YX arrays 
	protected final CustomTricubicFunction[][] splines;

	// The target range
	protected final int maxx, maxy;

	protected double tB;

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
	public CubicSplineFunction(CubicSplineData splineData, int maxx, int maxy) throws IllegalArgumentException
	{
		this.splines = splineData.splines;
		this.maxx = (maxx < 1) ? 1 : maxx;
		this.maxy = (maxy < 1) ? 1 : maxy;
		maxSx = splineData.maxx;
		maxSy = splineData.maxy;
		maxSz = splines.length;
		// Centre in the middle, assuming the min is zero
		cx = (maxSx / 2.0);
		cy = (maxSy / 2.0);
		cz = (maxSz / 2.0);
		setScale(1);
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
	public CubicSplineFunction(CubicSplineData splineData, int maxx, int maxy, double cx, double cy, double cz,
			int scale) throws IllegalArgumentException
	{
		this.splines = splineData.splines;
		this.maxx = (maxx < 1) ? 1 : maxx;
		this.maxy = (maxy < 1) ? 1 : maxy;
		maxSx = splineData.maxx;
		maxSy = splineData.maxy;
		maxSz = splines.length;
		this.cx = cx;
		this.cy = cy;
		this.cz = cz;
		setScale(scale);
	}

	private void updateFunctionBounds()
	{
		// Store the bounds of the cubic spline if it were positioned at 0,0
		lx = -cx / scale;
		ux = (maxSx - cx) / scale;
		ly = -cy / scale;
		uy = (maxSy - cy) / scale;
		// Check if the centre was within the function
		if (lx > 0 || ly > 0 || ux < 0 || uy < 0 || cz < 0 || cz > maxSz)
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
	abstract public int getN();

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

	public int size()
	{
		return maxx * maxy;
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

	abstract protected void initialise(double[] a, int order);
}