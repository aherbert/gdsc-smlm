package gdsc.smlm.utils;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2015 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Provides sampling from a 2D histogram
 */
public class PDF2D
{
	private final PDF[] rows;
	private final PDF sum;
	public final int nx, ny;
	/**
	 * The cumulative sum of the original input data
	 */
	public final double cumulative;

	/**
	 * Default constructor. Assumes the x-range and y-range increment from zero in integers.
	 * 
	 * @param data
	 *            The data (packed in XY order, i = nx*y + x)
	 * @param nx
	 *            The X-dimension size
	 * @param ny
	 *            The y-dimension size
	 * @throws InvalidArgumentException
	 *             if the dimensions are not above zero
	 * @throws InvalidArgumentException
	 *             if the input data length is not at least nx * ny
	 * @throws InvalidArgumentException
	 *             if the input data contains negatives
	 */
	public PDF2D(double[] data, int nx, int ny)
	{
		if (nx < 1 || ny < 1)
			throw new IllegalArgumentException("Dimensions must be above zero");
		this.nx = nx;
		this.ny = ny;

		if (data == null || data.length < nx * ny)
			throw new IllegalArgumentException("Input data must be at least equal to nx * ny");
		this.rows = new PDF[ny];
		double[] sum = new double[ny];

		// Build a PDF for each row of data
		double[] tmp = new double[nx];
		for (int y = 0, i = 0; y < ny; y++, i += nx)
		{
			System.arraycopy(data, i, tmp, 0, nx);
			rows[y] = new PDF(tmp);
			sum[y] = rows[y].cumulative;
		}

		// Build a PDF for the sum of the rows
		this.sum = new PDF(sum);
		cumulative = this.sum.cumulative;
	}

	/**
	 * Sample from the histogram using two uniform random numbers
	 * 
	 * @param r1
	 * @param r2
	 * @param point
	 *            The output coordinates buffer
	 * @return true if a sample was produced
	 */
	public boolean sample(double r1, double r2, double[] point)
	{
		if (point == null || point.length < 2)
			return false;

		// Sample within the sum of the rows to find the y-coordinate
		point[1] = sum.sample(r1);
		if (point[1] < 0)
			return false;
		
		// Sample within the specific row to find the x-coordinate
		point[0] = rows[(int) point[1]].sample(r2);
		if (point[0] < 0)
			return false;

		return true;
	}
}
