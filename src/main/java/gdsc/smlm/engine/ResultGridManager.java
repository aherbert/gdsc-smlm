package gdsc.smlm.engine;

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

import gdsc.core.utils.Maths;
import gdsc.smlm.filters.Spot;
import gdsc.smlm.results.PeakResult;

/**
 * Stores a set of results within a grid arrangement at a given resolution. Allows listing neighbours of a given
 * position.
 */
public class ResultGridManager
{
	private class PeakList
	{
		int c = 0;
		PeakResult[] list = null;

		void add(PeakResult peak)
		{
			if (list == null)
				list = new PeakResult[5];
			else if (list.length == c)
				list = Arrays.copyOf(list, (int) (c * 1.5));
			list[c++] = peak;
		}
	}

	private class SpotList
	{
		int c = 0;
		Spot[] list = null;

		void add(Spot spot)
		{
			if (list == null)
				list = new Spot[5];
			else if (list.length == c)
				list = Arrays.copyOf(list, (int) (c * 1.5));
			list[c++] = spot;
		}
	}

	private SpotList[][] spotGrid;
	private PeakList[][] peakGrid;
	private final int resolution, xBlocks, yBlocks;

	private PeakResult[] peakCache = null;
	private int peakCacheX = -1, peakCacheY = -1;
	private Spot[] spotCache = null;
	private int spotCacheX = -1, spotCacheY = -1;

	/**
	 * Clear the cache. This should be called when more data has been added to the grid.
	 */
	public void clearCache()
	{
		peakCache = null;
		peakCacheX = -1;
		peakCacheY = -1;
		spotCache = null;
		spotCacheX = -1;
		spotCacheY = -1;
	}

	/**
	 * Create a grid for spots or peak results
	 * 
	 * @param maxx
	 * @param maxy
	 * @param resolution
	 */
	public ResultGridManager(int maxx, int maxy, int resolution)
	{
		this.resolution = resolution;
		xBlocks = getBlock(maxx) + 1;
		yBlocks = getBlock(maxy) + 1;

		spotGrid = new SpotList[xBlocks][yBlocks];
		peakGrid = new PeakList[xBlocks][yBlocks];
	}

	/**
	 * Create a grid only of peak results. No spots can be added to the grid.
	 * 
	 * @param results
	 * @param resolution
	 */
	public ResultGridManager(PeakResult[] results, double resolution)
	{
		this.resolution = Maths.max(1, (int) Math.ceil(resolution));
		double maxx = 0, maxy = 0;
		for (PeakResult p : results)
		{
			if (maxx < p.getXPosition())
				maxx = p.getXPosition();
			if (maxy < p.getYPosition())
				maxy = p.getYPosition();
		}
		xBlocks = getBlock((int) maxx) + 1;
		yBlocks = getBlock((int) maxy) + 1;

		peakGrid = new PeakList[xBlocks][yBlocks];
		for (PeakResult p : results)
			putOnGrid(p);
	}

	private int getBlock(final int x)
	{
		return x / resolution;
	}

	/**
	 * Add a peak to the grid. Assumes that the coordinates are within the size of the grid.
	 * 
	 * @param peak
	 */
	public void addToGrid(PeakResult peak)
	{
		final int xBlock = getBlock((int) peak.getXPosition());
		final int yBlock = getBlock((int) peak.getYPosition());
		peakGrid[xBlock][yBlock].add(peak);
		clearCache();
	}
	
	/**
	 * Add a peak to the grid. Assumes that the coordinates are within the size of the grid.
	 * <p>
	 * This method does not clear the cache and should be called only when initialising the grid.
	 *
	 * @param peak the peak
	 */
	public void putOnGrid(PeakResult peak)
	{
		final int xBlock = getBlock((int) peak.getXPosition());
		final int yBlock = getBlock((int) peak.getYPosition());
		peakGrid[xBlock][yBlock].add(peak);
	}

	/**
	 * Add a spot to the grid. Assumes that the coordinates are within the size of the grid.
	 * 
	 * @param spot
	 */
	public void addToGrid(Spot spot)
	{
		final int xBlock = getBlock(spot.x);
		final int yBlock = getBlock(spot.y);
		spotGrid[xBlock][yBlock].add(spot);
		clearCache();
	}

	/**
	 * Add a spot to the grid. Assumes that the coordinates are within the size of the grid.
	 * <p>
	 * This method does not clear the cache and should be called only when initialising the grid.
	 *
	 * @param spot the spot
	 */
	public void putOnGrid(Spot spot)
	{
		final int xBlock = getBlock(spot.x);
		final int yBlock = getBlock(spot.y);
		spotGrid[xBlock][yBlock].add(spot);
	}
	
	/**
	 * Get the neighbours in the local region (defined by the input resolution). All neighbours within the resolution
	 * distance will be returned, plus there may be others and so distances should be checked.
	 * 
	 * @param x
	 * @param y
	 * @return the neighbours
	 */
	public PeakResult[] getPeakResultNeighbours(final int x, final int y)
	{
		final int xBlock = getBlock(x);
		final int yBlock = getBlock(y);

		if (peakCache != null && peakCacheX == xBlock && peakCacheY == yBlock)
			return peakCache;

		int size = 0;
		for (int xx = xBlock - 1; xx <= xBlock + 1; xx++)
		{
			if (xx < 0 || xx >= xBlocks)
				continue;
			for (int yy = yBlock - 1; yy <= yBlock + 1; yy++)
			{
				if (yy < 0 || yy >= yBlocks)
					continue;
				size += peakGrid[xx][yy].c;
			}
		}
		final PeakResult[] list = new PeakResult[size];
		size = 0;
		for (int xx = xBlock - 1; xx <= xBlock + 1; xx++)
		{
			if (xx < 0 || xx >= xBlocks)
				continue;
			for (int yy = yBlock - 1; yy <= yBlock + 1; yy++)
			{
				if (yy < 0 || yy >= yBlocks)
					continue;
				System.arraycopy(peakGrid[xx][yy].list, 0, list, size, peakGrid[xx][yy].c);
				size += peakGrid[xx][yy].c;
			}
		}

		peakCache = list;
		peakCacheX = xBlock;
		peakCacheY = yBlock;

		return list;
	}

	/**
	 * Get the neighbours in the local region (defined by the input resolution). All neighbours within the resolution
	 * distance will be returned, plus there may be others and so distances should be checked.
	 * 
	 * @param x
	 * @param y
	 * @param minScore
	 *            The minimum score that the the spot must have
	 * @return the neighbours
	 */
	public Spot[] getSpotNeighbours(final int x, final int y, float minScore)
	{
		final int xBlock = getBlock(x);
		final int yBlock = getBlock(y);

		int size = 0;
		for (int xx = xBlock - 1; xx <= xBlock + 1; xx++)
		{
			if (xx < 0 || xx >= xBlocks)
				continue;
			for (int yy = yBlock - 1; yy <= yBlock + 1; yy++)
			{
				if (yy < 0 || yy >= yBlocks)
					continue;
				size += spotGrid[xx][yy].c;
			}
		}
		final Spot[] list = new Spot[size];
		size = 0;
		for (int xx = xBlock - 1; xx <= xBlock + 1; xx++)
		{
			if (xx < 0 || xx >= xBlocks)
				continue;
			for (int yy = yBlock - 1; yy <= yBlock + 1; yy++)
			{
				if (yy < 0 || yy >= yBlocks)
					continue;
				final SpotList spotList = spotGrid[xx][yy];
				for (int i = 0; i < spotList.c; i++)
					if (spotList.list[i].getScore() >= minScore)
						list[size++] = spotList.list[i];
			}
		}

		return (size < list.length) ? Arrays.copyOf(list, size) : list;
	}

	/**
	 * Get the neighbours in the local region (defined by the input resolution). All neighbours within the resolution
	 * distance will be returned, plus there may be others and so distances should be checked.
	 * 
	 * @param x
	 * @param y
	 * @return the neighbours
	 */
	public Spot[] getSpotNeighbours(final int x, final int y)
	{
		final int xBlock = getBlock(x);
		final int yBlock = getBlock(y);

		if (spotCache != null && spotCacheX == xBlock && spotCacheY == yBlock)
			return spotCache;

		int size = 0;
		for (int xx = xBlock - 1; xx <= xBlock + 1; xx++)
		{
			if (xx < 0 || xx >= xBlocks)
				continue;
			for (int yy = yBlock - 1; yy <= yBlock + 1; yy++)
			{
				if (yy < 0 || yy >= yBlocks)
					continue;
				size += spotGrid[xx][yy].c;
			}
		}
		final Spot[] list = new Spot[size];
		size = 0;
		for (int xx = xBlock - 1; xx <= xBlock + 1; xx++)
		{
			if (xx < 0 || xx >= xBlocks)
				continue;
			for (int yy = yBlock - 1; yy <= yBlock + 1; yy++)
			{
				if (yy < 0 || yy >= yBlocks)
					continue;
				System.arraycopy(spotGrid[xx][yy].list, 0, list, size, spotGrid[xx][yy].c);
				size += spotGrid[xx][yy].c;
			}
		}

		spotCache = list;
		spotCacheX = xBlock;
		spotCacheY = yBlock;

		return list;
	}
}