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
	 * Create an empty grid for coordinates. The grid should be resized to the max dimensions of the data using
	 * {@link #resize(int, int)}
	 *
	 * @param resolution
	 *            the resolution
	 */
	GridCoordinateStore1(double resolution)
	{
		this(Math.max(0, resolution) * Math.max(0, resolution), 1, 1);
	}

	/**
	 * Create a grid for coordinates.
	 *
	 * @param maxx
	 *            the max x coordinate value
	 * @param maxy
	 *            the max y coordinate value
	 * @param resolution
	 *            the resolution
	 */
	GridCoordinateStore1(int maxx, int maxy, double resolution)
	{
		// Number of blocks is just the max dimensions since we have a block resolution of 1
		this(Math.max(0, resolution) * Math.max(0, resolution), Math.max(0, maxx) + 1, Math.max(0, maxy) + 1);
	}

	/**
	 * Create a grid for coordinates.
	 *
	 * @param d2
	 *            the d 2
	 * @param xBlocks
	 *            the x blocks
	 * @param yBlocks
	 *            the y blocks
	 */
	private GridCoordinateStore1(double d2, int xBlocks, int yBlocks)
	{
		// Use a block resolution of 1
		super(1, d2, xBlocks, yBlocks);
		if (d2 > 1)
			throw new IllegalArgumentException("Resolution must be less than 1");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.GridCoordinateStore#newInstance(int, int)
	 */
	@Override
	protected GridCoordinateStore newInstance(int xBlocks, int yBlocks)
	{
		return new GridCoordinateStore1(getSquaredDistance(), xBlocks, yBlocks);
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
		return (int) Math.round(blocks);
	}

	@Override
	public void changeResolution(double resolution)
	{
		if (resolution > 1)
			throw new IllegalArgumentException("Resolution must be less than 1");
		super.changeResolution(resolution);
	}
}