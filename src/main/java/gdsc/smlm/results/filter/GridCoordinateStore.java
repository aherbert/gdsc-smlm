package gdsc.smlm.results.filter;

import gdsc.core.utils.Maths;

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
 */
public class GridCoordinateStore implements CoordinateStore
{
	/**
	 * The minimum size for the block. This prevents excess use of memory when the resolution is smaller than this
	 * value.
	 */
	public final static double MINIMUM_BLOCK_SIZE = 0.5;

	private class CoordinateList
	{
		int timestamp = 0;
		int size = 0;
		double[] list = new double[3];

		// Not needed as we only create lists in the constructor when the timestamp is zero
		//CoordinateList()
		//{
		//	refresh();
		//}

		void add(double x, double y, double z)
		{
			if (list.length == size)
			{
				final double[] list2 = new double[size * 2];
				System.arraycopy(list, 0, list2, 0, size);
				list = list2;
			}
			list[size] = x;
			list[size + 1] = y;
			list[size + 2] = z;
			size += 3;
		}

		void refresh()
		{
			// If this is out-of-date then clear the contents 
			if (this.timestamp != GridCoordinateStore.this.timestamp)
				clear();
		}

		void clear()
		{
			size = 0;
			this.timestamp = GridCoordinateStore.this.timestamp;
		}
	}

	private int timestamp = 0;

	private CoordinateList[][] grid;
	private final CoordinateList queue = new CoordinateList();
	private double blockResolution;
	private double xyResolution;
	/**
	 * The squared XY resolution. If the XY resolution is negative then this should never be used as the store is not
	 * active.
	 */
	private double xy2;
	/**
	 * The z resolution. If this is negative then this is ignored and the store behaves as if processing 2D coordinates.
	 */
	private double zResolution;
	private int xBlocks, yBlocks;
	/**
	 * Flag to indicate that the store is active (i.e. storing coordinates). It is not active if the XY resolution is
	 * negative.
	 */
	private boolean isActive;
	/** Flag to indicate that the store is ignoring the z coordinate. This is true if the z resolution is negative. */
	private boolean is2D;

	// Note: This is a public constructor so this can be used without the factory

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
	public GridCoordinateStore(int maxx, int maxy, double xyResolution, double zResolution)
	{
		if (maxx < 0)
			maxx = 0;
		if (maxy < 0)
			maxy = 0;

		setXYResolution(xyResolution);
		setZResolution(zResolution);
		this.xBlocks = getBlock(maxx) + 1;
		this.yBlocks = getBlock(maxy) + 1;

		createGrid();
	}

	private void setXYResolution(double xyResolution)
	{
		this.xyResolution = xyResolution;
		this.blockResolution = Math.max(MINIMUM_BLOCK_SIZE, xyResolution);
		this.xy2 = xyResolution * xyResolution;
		isActive = xyResolution >= 0;
	}

	private void setZResolution(double zResolution)
	{
		this.zResolution = zResolution;
		is2D = zResolution < 0;
	}

	/**
	 * Gets the block for the coordinate.
	 *
	 * @param x
	 *            the coordinate
	 * @return the block
	 */
	protected int getBlock(final double x)
	{
		return (int) (x / blockResolution);
	}

	private void createGrid()
	{
		if (isActive)
		{
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
	}

	/**
	 * Create a grid for coordinates.
	 *
	 * @param blockResolution
	 *            the block resolution
	 * @param xy2
	 *            the xy resolution squared
	 * @param zResolution
	 *            the z resolution
	 * @param xBlocks
	 *            the x blocks
	 * @param yBlocks
	 *            the y blocks
	 */
	protected GridCoordinateStore(double xyResolution, double zResolution, int xBlocks, int yBlocks)
	{
		setXYResolution(xyResolution);
		setZResolution(zResolution);
		this.xBlocks = xBlocks;
		this.yBlocks = yBlocks;

		createGrid();
	}

	/**
	 * Create a new instance with the same settings
	 *
	 * @return the grid coordinate store
	 * @see gdsc.smlm.results.filter.CoordinateStore#newInstance()
	 */
	public GridCoordinateStore newInstance()
	{
		return newInstance(xBlocks, yBlocks);
	}

	/**
	 * Create a new instance with the same settings but different number of blocks
	 *
	 * @return the grid coordinate store
	 */
	protected GridCoordinateStore newInstance(int xBlocks, int yBlocks)
	{
		return new GridCoordinateStore(xyResolution, zResolution, xBlocks, yBlocks);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.CoordinateStore#resize(int, int)
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
		{
			clear(); // For consistency with a new instance also being empty.
			return this;
		}
		return newInstance(xBlocks, yBlocks);
	}

	/**
	 * Change the XY resolution of the store. The max dimensions are unchanged. Changing the resolution clears the
	 * store.
	 * 
	 * @param xyResolution
	 *            The new XY resolution
	 */
	public void changeXYResolution(double xyResolution)
	{
		clear();

		if (xyResolution == this.xyResolution || xyResolution < 0 && this.xyResolution < 0)
			// No size change
			return;

		int maxx = getCoordinate(xBlocks - 1);
		int maxy = getCoordinate(yBlocks - 1);

		setXYResolution(xyResolution);

		this.xBlocks = getBlock(maxx) + 1;
		this.yBlocks = getBlock(maxy) + 1;

		if (grid.length < xBlocks || grid[0].length < yBlocks)
			// The grid is larger so recreate
			createGrid();
	}

	/**
	 * Change Z resolution. Changing the resolution clears the store.
	 *
	 * @param zResolution
	 *            the z resolution
	 */
	public void changeZResolution(double zResolution)
	{
		// This is not strictly necessary since we use a 2D grid but is 
		// done for consistency with changeXYResolution(...).
		clear();
		setZResolution(zResolution);
	}

	/**
	 * Gets the coordinate for the block.
	 *
	 * @param blocks
	 *            the blocks
	 * @return the coordinate
	 */
	protected int getCoordinate(final int blocks)
	{
		return (int) Math.round(blocks * this.blockResolution);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.CoordinateStore#getXYResolution()
	 */
	public double getXYResolution()
	{
		return xyResolution;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.CoordinateStore#getZResolution()
	 */
	public double getZResolution()
	{
		return zResolution;
	}

	/**
	 * Note: This does not check that the x,y coordinates are within the correct bounds.
	 * 
	 * @see gdsc.smlm.results.filter.CoordinateStore#addToQueue(double, double, double)
	 */
	public void addToQueue(double x, double y, double z)
	{
		if (isActive)
		{
			queue.add(x, y, z);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.CoordinateStore#flush()
	 */
	public void flush()
	{
		for (int i = 0; i < queue.size; i += 3)
			addf(queue.list[i], queue.list[i + 1], queue.list[i + 2]);
		//queue.clear(); // Avoid the timestamp refresh
		queue.size = 0;
	}

	/**
	 * Note: This does not check that the x,y coordinates are within the correct bounds. Use
	 * {@link #safeAdd(double, double)} to do a bounds check.
	 * 
	 * @see gdsc.smlm.results.filter.CoordinateStore#add(double, double, double)
	 */
	public void add(final double x, final double y, final double z)
	{
		if (isActive)
		{
			addf(x, y, z);
		}
	}

	private void addf(final double x, final double y, final double z)
	{
		// Note: A bounds check is not currently necessary since the this code is only used
		// by the MultiPathFilter on results that are inside the bounds
		getList(getBlock(x), getBlock(y)).add(x, y, z);
	}

	/**
	 * Add a coordinate to the store. Assumes that the coordinates are within the size of the grid otherwise they will
	 * be ignored.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param z
	 *            the z
	 */
	public void safeAdd(final double x, final double y, final double z)
	{
		if (isActive)
		{
			// Check bounds 
			final int xBlock = getBlock(x);
			if (xBlock < 0 || xBlock >= xBlocks)
				return;
			final int yBlock = getBlock(y);
			if (yBlock < 0 || yBlock >= yBlocks)
				return;
			getList(xBlock, yBlock).add(x, y, z);
		}
	}

	/**
	 * Gets the list for the grid point. Refresh the list using the current timestamp (i.e. clear any old data).
	 *
	 * @param xBlock
	 *            the x grid point
	 * @param yBlock
	 *            the y grid point
	 * @return the list
	 */
	private CoordinateList getList(final int xBlock, final int yBlock)
	{
		final CoordinateList l = grid[xBlock][yBlock];
		l.refresh();
		return l;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.CoordinateStore#clear()
	 */
	public void clear()
	{
		// Clearing each item in the grid is a big overhead when the grid is large and the number of additions to the grid is small.
		// So store a timestamp for the clear and we refresh each list when we next use it.
		timestamp++;
		queue.size = 0;

		// Reset after an entire cycle of timestamps
		if (timestamp == 0)
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
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.CoordinateStore#contains(double, double, double)
	 */
	public boolean contains(final double x, final double y, final double z)
	{
		// If not active then nothing could have been added
		if (!isActive)
			return false;

		final int xBlock = getBlock(x);
		final int yBlock = getBlock(y);

		final int xmin = Math.max(0, xBlock - 1);
		final int ymin = Math.max(0, yBlock - 1);
		final int xmax = Math.min(xBlocks, xBlock + 2);
		final int ymax = Math.min(yBlocks, yBlock + 2);

		for (int xx = xmin; xx < xmax; xx++)
		{
			for (int yy = ymin; yy < ymax; yy++)
			{
				final CoordinateList l = getList(xx, yy);
				final int size = l.size;
				if (size == 0)
					continue;
				final double[] list = l.list;
				for (int i = 0; i < size; i += 3)
				{
					if (distance2(x, y, list[i], list[i + 1]) <= xy2)
					{
						if (is2D)
							return true;
						// Otherwise z resolution is not negative and we check that
						if (distance(z, list[i + 2]) <= zResolution)
							return true;
					}
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
	private static double distance2(final double x, final double y, final double x2, final double y2)
	{
		final double dx = x - x2;
		final double dy = y - y2;
		return dx * dx + dy * dy;
	}

	/**
	 * Get the absolute distance.
	 *
	 * @param x
	 *            the x coordinate
	 * @param x2
	 *            the x2 coordinate
	 * @return the absolute distance
	 */
	private static double distance(final double x, final double x2)
	{
		return (x > x2) ? x - x2 : x2 - x;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.CoordinateStore#find(double, double, double)
	 */
	public double[] find(final double x, final double y, final double z)
	{
		if (!isActive)
			return null;

		final int xBlock = getBlock(x);
		final int yBlock = getBlock(y);

		double[] match = new double[3];
		double min = Double.POSITIVE_INFINITY;

		final int xmin = Math.max(0, xBlock - 1);
		final int ymin = Math.max(0, yBlock - 1);
		final int xmax = Math.min(xBlocks, xBlock + 2);
		final int ymax = Math.min(yBlocks, yBlock + 2);

		for (int xx = xmin; xx < xmax; xx++)
		{
			for (int yy = ymin; yy < ymax; yy++)
			{
				final CoordinateList l = getList(xx, yy);
				final int size = l.size;
				if (size == 0)
					continue;
				final double[] list = l.list;
				for (int i = 0; i < size; i += 3)
				{
					double d = distance2(x, y, list[i], list[i + 1]);
					if (d <= xy2)
					{
						if (!is2D)
						{
							// z resolution is not negative and we check that
							final double dd = distance(z, list[i + 2]);
							if (dd > zResolution)
								continue;
							// Get a combined Euclidean squared distance
							d += Maths.pow2(dd);
						}

						if (d < min)
						{
							min = d;
							match[0] = list[i];
							match[1] = list[i + 1];
							match[2] = list[i + 2];
						}
					}
				}
			}
		}

		return (min < Double.POSITIVE_INFINITY) ? match : null;
	}
}