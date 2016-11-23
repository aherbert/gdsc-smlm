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

import java.util.Arrays;

/**
 * Stores a set of results within a grid arrangement at a given resolution. Allows checking for duplicates.
 */
public class GridCoordinateStore implements CoordinateStore
{
	private class CoordinateList
	{
		int size = 0;
		double[] list = new double[4];

		void add(double x, double y)
		{
			if (list.length == size)
				list = Arrays.copyOf(list, size * 2);
			list[size++] = x;
			list[size++] = y;
		}

		void clear()
		{
			size = 0;
		}
	}

	private final CoordinateList[][] grid;
	private final CoordinateList queue = new CoordinateList();
	private final double blockResolution;
	private final double d2;
	private final int xBlocks, yBlocks;

	/**
	 * Create an empty grid for coordinates. The grid should be resized to the max dimensions of the data using
	 * {@link #resize(int, int)}
	 *
	 * @param resolution
	 *            the resolution
	 */
	public GridCoordinateStore(double resolution)
	{
		this(0, 0, resolution);
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
	public GridCoordinateStore(int maxx, int maxy, double resolution)
	{
		if (maxx < 0)
			maxx = 0;
		if (maxy < 0)
			maxy = 0;
		if (resolution < 0)
			resolution = 0;
		this.blockResolution = Math.max(0.5, resolution);
		this.d2 = resolution * resolution;
		this.xBlocks = getBlock(maxx) + 1;
		this.yBlocks = getBlock(maxy) + 1;

		grid = new CoordinateList[xBlocks][yBlocks];
		for (int x = 0; x < xBlocks; x++)
		{
			final CoordinateList[] list = grid[x];
			for (int y = 0; y < yBlocks; y++)
			{
				list[y] = new CoordinateList();
			}
		}
	}

	/**
	 * Create a grid for coordinates.
	 *
	 * @param blockResolution
	 *            the block resolution
	 * @param d2
	 *            the d 2
	 * @param xBlocks
	 *            the x blocks
	 * @param yBlocks
	 *            the y blocks
	 */
	private GridCoordinateStore(double blockResolution, double d2, int xBlocks, int yBlocks)
	{
		this.blockResolution = blockResolution;
		this.d2 = d2;
		this.xBlocks = xBlocks;
		this.yBlocks = yBlocks;

		grid = new CoordinateList[xBlocks][yBlocks];
		for (int x = 0; x < xBlocks; x++)
		{
			final CoordinateList[] list = grid[x];
			for (int y = 0; y < yBlocks; y++)
			{
				list[y] = new CoordinateList();
			}
		}
	}

	/**
	 * Create a new instance with the same settings
	 *
	 * @return the grid coordinate store
	 */
	public GridCoordinateStore newInstance()
	{
		return new GridCoordinateStore(blockResolution, d2, xBlocks, yBlocks);
	}

	/**
	 * Resize to the given dimensions. If these match the existing dimensions the current store is returned. Otherwise a
	 * new store is returned.
	 *
	 * @param maxx
	 *            the max x coordinate value
	 * @param maxy
	 *            the max y coordinate value
	 * @return the grid coordinate store
	 */
	public GridCoordinateStore resize(int maxx, int maxy)
	{
		if (maxx < 0)
			maxx = 0;
		if (maxy < 0)
			maxy = 0;
		int xBlocks = getBlock(maxx) + 1;
		int yBlocks = getBlock(maxy) + 1;
		if (this.xBlocks == xBlocks && this.yBlocks == yBlocks)
			return this;
		return new GridCoordinateStore(blockResolution, d2, xBlocks, yBlocks);
	}

	/**
	 * Gets the block for the coordinate.
	 *
	 * @param x
	 *            the coordinate
	 * @return the block
	 */
	private int getBlock(final double x)
	{
		return (int) (x / blockResolution);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.CoordinateStore#getResolution()
	 */
	public double getResolution()
	{
		return Math.sqrt(d2);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.CoordinateStore#queue(double, double)
	 */
	public boolean queue(double x, double y)
	{
		if (contains(x, y))
			return false;
		queue.add(x, y);
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.CoordinateStore#flush()
	 */
	public void flush()
	{
		for (int i = 0; i < queue.size; i += 2)
			add(queue.list[i], queue.list[i + 1]);
		queue.clear();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.CoordinateStore#add(double, double)
	 */
	public void add(double x, double y)
	{
		grid[getBlock(x)][getBlock(y)].add(x, y);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.CoordinateStore#clear()
	 */
	public void clear()
	{
		for (int x = 0; x < xBlocks; x++)
		{
			final CoordinateList[] list = grid[x];
			for (int y = 0; y < yBlocks; y++)
			{
				list[y].clear();
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.CoordinateStore#contains(double, double)
	 */
	public boolean contains(double x, double y)
	{
		final int xBlock = getBlock(x);
		final int yBlock = getBlock(y);

		for (int xx = xBlock - 1; xx <= xBlock + 1; xx++)
		{
			if (xx < 0 || xx >= xBlocks)
				continue;
			for (int yy = yBlock - 1; yy <= yBlock + 1; yy++)
			{
				if (yy < 0 || yy >= yBlocks)
					continue;
				final int size = grid[xx][yy].size;
				if (size == 0)
					continue;
				final double[] list = grid[xx][yy].list;
				for (int i = 0; i < size; i += 2)
				{
					if (distance2(x, y, list[i], list[i + 1]) < d2)
						return true;
				}
			}
		}

		return false;
	}

	/**
	 * Get the squared distance
	 *
	 * @param x
	 *            the x coordinate
	 * @param y
	 *            the y coordinate
	 * @param x2
	 *            the x2 coordinate
	 * @param y2
	 *            the y2 coordinate
	 * @return the squared distance
	 */
	private double distance2(double x, double y, double x2, double y2)
	{
		final double dx = x - x2;
		final double dy = y - y2;
		return dx * dx + dy * dy;
	}

	public double[] find(double x, double y)
	{
		final int xBlock = getBlock(x);
		final int yBlock = getBlock(y);

		double[] match = new double[2];
		double min = d2;
		for (int xx = xBlock - 1; xx <= xBlock + 1; xx++)
		{
			if (xx < 0 || xx >= xBlocks)
				continue;
			for (int yy = yBlock - 1; yy <= yBlock + 1; yy++)
			{
				if (yy < 0 || yy >= yBlocks)
					continue;
				final int size = grid[xx][yy].size;
				if (size == 0)
					continue;
				final double[] list = grid[xx][yy].list;
				for (int i = 0; i < size; i += 2)
				{
					final double d = distance2(x, y, list[i], list[i + 1]);
					if (d < min)
					{
						min = d;
						match[0] = list[i];
						match[1] = list[i + 1];
					}
				}
			}
		}

		return (min < d2) ? match : null;
	}
}