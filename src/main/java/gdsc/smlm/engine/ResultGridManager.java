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
import java.util.Comparator;

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
				list = new PeakResult[4];
			else if (list.length == size)
				list = Arrays.copyOf(list, size * 2);
			list[size++] = peak;
		}
	}

	private class CandidateComparator implements Comparator<Candidate>
	{
		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
		 */
		public int compare(Candidate o1, Candidate o2)
		{
			return o1.index - o2.index;
		}
	}

	private final static CandidateComparator comp;
	static
	{
		comp = new ResultGridManager().new CandidateComparator();
	}

	public class CandidateList
	{
		int size = 0;
		Candidate[] list = null;

		/**
		 * Instantiates a new candidate list.
		 */
		private CandidateList()
		{
		}

		/**
		 * Instantiates a new candidate list.
		 *
		 * @param size
		 *            the size
		 * @param list
		 *            the list
		 */
		private CandidateList(int size, Candidate[] list)
		{
			this.size = size;
			this.list = list;
		}

		/**
		 * Add a candidate
		 * 
		 * @param candidate
		 */
		public void add(Candidate candidate)
		{
			if (list == null)
				list = new Candidate[4];
			else if (list.length == size)
				list = Arrays.copyOf(list, size * 2);
			list[size++] = candidate;
		}

		/**
		 * Sort in ascending order of Id
		 */
		public void sort()
		{
			if (size != 0)
				Arrays.sort(list, 0, size, comp);
		}

		/**
		 * Gets the size.
		 *
		 * @return the size
		 */
		public int getSize()
		{
			return size;
		}

		/**
		 * Gets the candidate
		 *
		 * @param index
		 *            the index
		 * @return the candidate
		 */
		public Candidate get(int index)
		{
			return list[index];
		}

		/**
		 * Copy this list.
		 *
		 * @return the new candidate list
		 */
		public CandidateList copy()
		{
			return new CandidateList(size, Arrays.copyOf(list, size));
		}
	}

	private CandidateList[][] candidateGrid;
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

	private ResultGridManager()
	{
		resolution = xBlocks = yBlocks = 0;
	}

	/**
	 * Create a grid for candidates or peak results
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

		createCandidateGrid();
		createPeakGrid();
	}

	private void createPeakGrid()
	{
		peakGrid = new PeakList[xBlocks][yBlocks];
		for (int x = 0; x < xBlocks; x++)
		{
			final PeakList[] list = peakGrid[x];
			for (int y = 0; y < yBlocks; y++)
			{
				list[y] = new PeakList();
			}
		}
	}

	private void createCandidateGrid()
	{
		candidateGrid = new CandidateList[xBlocks][yBlocks];
		for (int x = 0; x < xBlocks; x++)
		{
			final CandidateList[] list = candidateGrid[x];
			for (int y = 0; y < yBlocks; y++)
			{
				list[y] = new CandidateList();
			}
		}
	}

	/**
	 * Create a grid only of peak results. No candidates can be added to the grid.
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

		createPeakGrid();
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
	 * Add a candidate to the grid. Assumes that the coordinates are within the size of the grid.
	 * 
	 * @param candidate
	 */
	public void addToGrid(Candidate candidate)
	{
		final int xBlock = getBlock(candidate.x);
		final int yBlock = getBlock(candidate.y);
		candidateGrid[xBlock][yBlock].add(candidate);
		clearCache();
	}

	/**
	 * Add a candidate to the grid. Assumes that the coordinates are within the size of the grid.
	 * <p>
	 * This method does not clear the cache and should be called only when initialising the grid.
	 *
	 * @param candidate
	 *            the candidate
	 */
	public void putOnGrid(Candidate candidate)
	{
		final int xBlock = getBlock(candidate.x);
		final int yBlock = getBlock(candidate.y);
		candidateGrid[xBlock][yBlock].add(candidate);
	}

	/**
	 * Get the neighbours in the local region (defined by the input resolution). All neighbours within the
	 * resolution
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
		if (size != 0)
		{
			size = 0;
			for (int xx = xBlock - 1; xx <= xBlock + 1; xx++)
			{
				if (xx < 0 || xx >= xBlocks)
					continue;
				for (int yy = yBlock - 1; yy <= yBlock + 1; yy++)
				{
					if (yy < 0 || yy >= yBlocks || peakGrid[xx][yy].size == 0)
						continue;
					System.arraycopy(peakGrid[xx][yy].list, 0, list, size, peakGrid[xx][yy].size);
					size += peakGrid[xx][yy].size;
				}
			}
		}

		peakCache = list;
		peakCacheX = xBlock;
		peakCacheY = yBlock;

		return list;
	}

	/**
	 * Get the neighbours in the local region (defined by the input resolution). All neighbours within the
	 * resolution distance will be returned, plus there may be others and so distances should be checked.
	 * <p>
	 * Note: Assumes candidate indices are unique.
	 *
	 * @param candidate
	 *            the candidate
	 * @return the neighbours
	 */
	public CandidateList getNeighbours(final Candidate candidate)
	{
		return getNeighbours(candidate, false);
	}

	/**
	 * Get the neighbours in the local region (defined by the input resolution). All neighbours within the
	 * resolution distance will be returned, plus there may be others and so distances should be checked.
	 *
	 * @param candidate
	 *            the candidate
	 * @param nonUnique
	 *            True if candidate indices are not unique unique
	 * @return the neighbours
	 */
	public CandidateList getNeighbours(final Candidate candidate, boolean nonUnique)
	{
		if (neighbourCache != null && neighbourCacheCandidate != null &&
				neighbourCacheCandidate.index == candidate.index)
			return neighbourCache;

		// Get all
		final Candidate[] list = getCandidates(candidate.x, candidate.y);

		// Remove the candidate.
		int size;
		if (nonUnique)
		{
			// Assumes non-unique candidates
			size = 0;
			for (int i = 0; i < list.length; i++)
			{
				if (list[i].index == candidate.index)
					continue;
				list[size++] = list[i];
			}
		}
		else
		{
			size = list.length;
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
		}

		neighbourCache = new CandidateList(size, list);
		neighbourCacheCandidate = candidate;

		return neighbourCache;
	}

	/**
	 * Get the neighbours in the local region (defined by the input resolution). All neighbours within the
	 * resolution
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
				size += candidateGrid[xx][yy].size;
			}
		}
		final Candidate[] list = new Candidate[size];
		if (size != 0)
		{
			size = 0;
			for (int xx = xBlock - 1; xx <= xBlock + 1; xx++)
			{
				if (xx < 0 || xx >= xBlocks)
					continue;
				for (int yy = yBlock - 1; yy <= yBlock + 1; yy++)
				{
					if (yy < 0 || yy >= yBlocks || candidateGrid[xx][yy].size == 0)
						continue;
					System.arraycopy(candidateGrid[xx][yy].list, 0, list, size, candidateGrid[xx][yy].size);
					size += candidateGrid[xx][yy].size;
				}
			}
		}

		return list;
	}
}