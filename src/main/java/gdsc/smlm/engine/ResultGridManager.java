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
import gdsc.smlm.results.PeakResult;

/**
 * Stores a set of results within a grid arrangement at a given resolution. Allows listing neighbours of a given
 * position.
 */
public class ResultGridManager
{
	private class PeakList
	{
		int size = 0;
		PeakResult[] list = null;

		void add(PeakResult peak)
		{
			if (list == null)
				list = new PeakResult[5];
			else if (list.length == size)
				list = Arrays.copyOf(list, (int) (size * 1.5));
			list[size++] = peak;
		}
	}

	public class CandidateList
	{
		int size = 0;
		Candidate[] list = null;

		private CandidateList(int size, Candidate[] list)
		{
			this.size = size;
			this.list = list;
		}

		public void add(Candidate spot)
		{
			if (list == null)
				list = new Candidate[5];
			else if (list.length == size)
				list = Arrays.copyOf(list, (int) (size * 1.5));
			list[size++] = spot;
		}
	}

	private CandidateList[][] spotGrid;
	private PeakList[][] peakGrid;
	private final int resolution, xBlocks, yBlocks;

	private PeakResult[] peakCache = null;
	private int peakCacheX = -1, peakCacheY = -1;
	private CandidateList neighbourCache = null;
	private Candidate neighbourCacheCandidate = null;

	/**
	 * Clear the cache. This should be called when more data has been added to the grid.
	 */
	public void clearCache()
	{
		peakCache = null;
		peakCacheX = -1;
		peakCacheY = -1;
		neighbourCache = null;
		neighbourCacheCandidate = null;
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

		spotGrid = new CandidateList[xBlocks][yBlocks];
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
	 * @param peak
	 *            the peak
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
	public void addToGrid(Candidate spot)
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
	 * @param spot
	 *            the spot
	 */
	public void putOnGrid(Candidate spot)
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
				size += peakGrid[xx][yy].size;
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
				System.arraycopy(peakGrid[xx][yy].list, 0, list, size, peakGrid[xx][yy].size);
				size += peakGrid[xx][yy].size;
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
	 * @param candidate
	 *            the candidate
	 * @return the neighbours
	 */
	public CandidateList getNeighbours(final Candidate candidate)
	{
		if (neighbourCache != null && neighbourCacheCandidate != null &&
				neighbourCacheCandidate.index == candidate.index)
			return neighbourCache;

		// Get all
		final Candidate[] list = getCandidates(candidate.x, candidate.y);

		// Remove the candidate. 

		//		// Assumes non-unique candidates
		//		int size = 0;
		//		for (int i = 0; i < list.length; i++)
		//		{
		//			if (list[i].index == c.index)
		//				continue;
		//			list[size++] = list[i];
		//		}
		//		list = (size < list.length) ? Arrays.copyOf(list, size) : list;

		int size = list.length;
		for (int i = 0; i < size; i++)
		{
			if (list[i].index == candidate.index)
			{
				int remaining = list.length - i - 1;
				if (remaining != 0)
					System.arraycopy(list, i + 1, list, i, remaining);
				size--;
				// Assume a unique candidate index
				break;
			}
		}

		neighbourCache = new CandidateList(size, list);
		neighbourCacheCandidate = candidate;

		return neighbourCache;
	}

	/**
	 * Get the neighbours in the local region (defined by the input resolution). All neighbours within the resolution
	 * distance will be returned, plus there may be others and so distances should be checked.
	 * 
	 * @param x
	 * @param y
	 * @return the neighbours
	 */
	private Candidate[] getCandidates(final int x, final int y)
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
				size += spotGrid[xx][yy].size;
			}
		}
		final Candidate[] list = new Candidate[size];
		size = 0;
		for (int xx = xBlock - 1; xx <= xBlock + 1; xx++)
		{
			if (xx < 0 || xx >= xBlocks)
				continue;
			for (int yy = yBlock - 1; yy <= yBlock + 1; yy++)
			{
				if (yy < 0 || yy >= yBlocks)
					continue;
				System.arraycopy(spotGrid[xx][yy].list, 0, list, size, spotGrid[xx][yy].size);
				size += spotGrid[xx][yy].size;
			}
		}

		return list;
	}
}