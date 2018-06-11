/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 * 
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.engine;


/**
 * Stores a set of results within a grid arrangement at a given resolution. Allows listing neighbours of a given
 * position.
 */
class CandidateGridManager
{
	private CandidateList[][] candidateGrid;
	private CandidateList[][] fittedCandidateGrid;
	private final int resolution, xBlocks, yBlocks;

	private CandidateList fittedCandidateCache = null;
	private int fittedCandidateCacheX = -1, fittedCandidateCacheY = -1;
	private CandidateList neighbourCache = null;
	private Candidate neighbourCacheCandidate = null;
	
	private CandidateList fitted = new CandidateList();

	/**
	 * Clear the cache. This should be called when more data has been added to the grid.
	 */
	public void clearCache()
	{
		fittedCandidateCache = null;
		fittedCandidateCacheX = -1;
		fittedCandidateCacheY = -1;
		neighbourCache = null;
		neighbourCacheCandidate = null;
	}

	/**
	 * Create a grid for candidates or fittedCandidate results
	 * 
	 * @param maxx
	 * @param maxy
	 * @param resolution
	 */
	public CandidateGridManager(int maxx, int maxy, int resolution)
	{
		this.resolution = resolution;
		xBlocks = getBlock(maxx) + 1;
		yBlocks = getBlock(maxy) + 1;

		createGrids();
	}

	private void createGrids()
	{
		candidateGrid = new CandidateList[xBlocks][yBlocks];
		fittedCandidateGrid = new CandidateList[xBlocks][yBlocks];
		for (int x = 0; x < xBlocks; x++)
		{
			final CandidateList[] list = candidateGrid[x];
			final CandidateList[] list2 = fittedCandidateGrid[x];
			for (int y = 0; y < yBlocks; y++)
			{
				list[y] = new CandidateList();
				list2[y] = new CandidateList();
			}
		}
	}

	private int getBlock(final int x)
	{
		return x / resolution;
	}

	/**
	 * Add a fittedCandidate to the grid. Assumes that the coordinates are within the size of the grid.
	 *
	 * @param candidate
	 *            the candidate
	 */
	public void addFittedToGrid(Candidate candidate)
	{
		putFittedOnGrid(candidate);
		clearCache();
	}

	/**
	 * Add a fittedCandidate to the grid. Does not assume that the coordinates are within the size of the grid.
	 * <p>
	 * This method does not clear the cache and should be called only when initialising the grid.
	 *
	 * @param candidate
	 *            the candidate
	 */
	public void putFittedOnGrid(Candidate candidate)
	{
		int xBlock = getBlock(candidate.x);
		int yBlock = getBlock(candidate.y);
		// The fitted position may be off the grid
		if (xBlock < 0)
			xBlock = 0;
		else if (xBlock >= xBlocks)
			xBlock = xBlocks - 1;
		if (yBlock < 0)
			yBlock = 0;
		else if (yBlock >= yBlocks)
			yBlock = yBlocks - 1;
		fittedCandidateGrid[xBlock][yBlock].add(candidate);
		fitted.add(candidate);
	}
	
	/**
	 * Gets the fitted candidates.
	 *
	 * @return the fitted candidates
	 */
	public CandidateList getFittedCandidates()
	{
		return fitted;
	}
	
	/**
	 * Gets the fitted candidates size.
	 *
	 * @return the fitted candidates size
	 */
	public int getFittedCandidatesSize()
	{
		return fitted.getSize();
	}

	/**
	 * Add a candidate to the grid. Assumes that the coordinates are within the size of the grid.
	 * 
	 * @param candidate
	 */
	public void addCandidateToGrid(Candidate candidate)
	{
		putCandidateOnGrid(candidate);
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
	public void putCandidateOnGrid(Candidate candidate)
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
	public CandidateList getFittedNeighbours(final int x, final int y)
	{
		final int xBlock = getBlock(x);
		final int yBlock = getBlock(y);

		if (fittedCandidateCache != null && fittedCandidateCacheX == xBlock && fittedCandidateCacheY == yBlock)
			return fittedCandidateCache;

		final Candidate[] list = getCandidates(fittedCandidateGrid, x, y);

		fittedCandidateCache = new CandidateList(list);
		fittedCandidateCacheX = xBlock;
		fittedCandidateCacheY = yBlock;

		return fittedCandidateCache;
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
	public CandidateList getCandidateNeighbours(final Candidate candidate)
	{
		return getCandidateNeighbours(candidate, false);
	}

	/**
	 * Get the neighbours in the local region (defined by the input resolution). All neighbours within the
	 * resolution distance will be returned, plus there may be others and so distances should be checked.
	 *
	 * @param candidate
	 *            the candidate
	 * @param nonUnique
	 *            True if candidate indices are not unique
	 * @return the neighbours
	 */
	public CandidateList getCandidateNeighbours(final Candidate candidate, boolean nonUnique)
	{
		if (neighbourCache != null && neighbourCacheCandidate != null &&
				neighbourCacheCandidate.index == candidate.index)
			return neighbourCache;

		// Get all
		final Candidate[] list = getCandidates(candidateGrid, candidate.x, candidate.y);

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
	 * @param grid
	 *            the candidate grid
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @return the neighbours
	 */
	private Candidate[] getCandidates(CandidateList[][] grid, final int x, final int y)
	{
		final int xBlock = getBlock(x);
		final int yBlock = getBlock(y);

		int size = 0;

		final int xmin = Math.max(0, xBlock - 1);
		final int ymin = Math.max(0, yBlock - 1);
		final int xmax = Math.min(xBlocks, xBlock + 2);
		final int ymax = Math.min(yBlocks, yBlock + 2);

		for (int xx = xmin; xx < xmax; xx++)
		{
			for (int yy = ymin; yy < ymax; yy++)
			{
				size += grid[xx][yy].getSize();
			}
		}
		final Candidate[] list = new Candidate[size];
		if (size != 0)
		{
			size = 0;
			for (int xx = xmin; xx < xmax; xx++)
			{
				for (int yy = ymin; yy < ymax; yy++)
				{
					if (grid[xx][yy].getSize() == 0)
						continue;
					grid[xx][yy].copyTo(list, size);
					size += grid[xx][yy].getSize();
				}
			}
		}

		return list;
	}
}
