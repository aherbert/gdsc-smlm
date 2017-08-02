package gdsc.smlm.results.filter;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Stores a set of results within a grid arrangement at a given resolution. Allows checking for duplicates.
 * <p>
 * Uses a block resolution of 1.
 */
public class GridCoordinateStore1 extends GridCoordinateStore
{
	// Note: We have package level constructors so that the factory must be used to create an instance.

	/**
	 * Create a grid for coordinates.
	 *
	 * @param minx
	 *            the min x coordinate value
	 * @param miny
	 *            the min y coordinate value
	 * @param width
	 *            the width
	 * @param height
	 *            the height
	 * @param xyResolution
	 *            the xy resolution
	 * @param zResolution
	 *            the z resolution
	 */
	GridCoordinateStore1(int minx, int miny, int width, int height, double xyResolution, double zResolution)
	{
		super(minx, miny, width, height, xyResolution, zResolution);
		checkResolution(xyResolution);
	}

	private void checkResolution(double xyResolution)
	{
		if (xyResolution > 1)
			throw new IllegalArgumentException("XY Resolution must be <= 1");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.GridCoordinateStore#newInstance(int, int, int, int)
	 */
	@Override
	public GridCoordinateStore newInstance(int minx, int miny, int width, int height)
	{
		return new GridCoordinateStore1(minx, miny, width, height, getXYResolution(), getZResolution());
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.GridCoordinateStore#getBlock(double)
	 */
	@Override
	protected int getBlock(final double x)
	{
		// blockResolution is always 1
		return (int) x;
	}

	@Override
	public void changeXYResolution(double xyResolution)
	{
		checkResolution(xyResolution);
		super.changeXYResolution(xyResolution);
	}
}