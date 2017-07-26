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
	 * @param maxx
	 *            the max x coordinate value
	 * @param maxy
	 *            the max y coordinate value
	 * @param xyResolution
	 *            the xy resolution
	 * @param zResolution
	 *            the z resolution
	 */
	GridCoordinateStore1(int maxx, int maxy, double xyResolution, double zResolution)
	{
		super(maxx, maxy, xyResolution, zResolution);
		checkResolution(xyResolution);
	}

	private void checkResolution(double xyResolution)
	{
		if (xyResolution > 1)
			throw new IllegalArgumentException("XY Resolution must be <= 1");
	}

	/**
	 * Create a grid for coordinates.
	 *
	 * @param xyResolution
	 *            the xy resolution
	 * @param zResolution
	 *            the z resolution
	 * @param xBlocks
	 *            the x blocks
	 * @param yBlocks
	 *            the y blocks
	 */
	private GridCoordinateStore1(double xyResolution, double zResolution, int xBlocks, int yBlocks)
	{
		super(xyResolution, zResolution, xBlocks, yBlocks);
		checkResolution(xyResolution);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.GridCoordinateStore#newInstance(int, int)
	 */
	@Override
	protected GridCoordinateStore newInstance(int xBlocks, int yBlocks)
	{
		return new GridCoordinateStore1(getXYResolution(), getZResolution(), xBlocks, yBlocks);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.GridCoordinateStore#getBlock(double)
	 */
	@Override
	protected int getBlock(final double x)
	{
		return (int) x;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.GridCoordinateStore#getCoordinate(int)
	 */
	@Override
	protected int getCoordinate(int blocks)
	{
		return blocks;
	}

	@Override
	public void changeXYResolution(double xyResolution)
	{
		checkResolution(xyResolution);
		super.changeXYResolution(xyResolution);
	}
}