package gdsc.smlm.results;

import gdsc.core.utils.Maths;
import gdsc.smlm.results.PeakResult;

/**
 * Stores a set of results within a grid arrangement at a given resolution. Allows listing neighbours of a given
 * position.
 */
public class PeakResultGridManager
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
			{
				final PeakResult[] list2 = new PeakResult[size * 2];
				System.arraycopy(list, 0, list2, 0, size);
				list = list2;
			}
			list[size++] = peak;
		}
	}

	private PeakList[][] peakGrid;
	private final int resolution, xBlocks, yBlocks;

	private PeakResult[] peakCache = null;
	private int peakCacheX = -1, peakCacheY = -1;

	/**
	 * Clear the cache. This should be called when more data has been added to the grid.
	 */
	public void clearCache()
	{
		peakCache = null;
		peakCacheX = -1;
		peakCacheY = -1;
	}

	/**
	 * Create a grid for candidates or peak results
	 * 
	 * @param maxx
	 * @param maxy
	 * @param resolution
	 */
	public PeakResultGridManager(int maxx, int maxy, int resolution)
	{
		this.resolution = resolution;
		xBlocks = getBlock(maxx) + 1;
		yBlocks = getBlock(maxy) + 1;

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

	/**
	 * Create a grid only of peak results. No candidates can be added to the grid.
	 * 
	 * @param results
	 * @param resolution
	 */
	public PeakResultGridManager(PeakResult[] results, double resolution)
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
		putOnGrid(peak);
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
	private void putOnGrid(PeakResult peak)
	{
		final int xBlock = getBlock((int) peak.getXPosition());
		final int yBlock = getBlock((int) peak.getYPosition());
		peakGrid[xBlock][yBlock].add(peak);
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
		
		final int xmin = Math.max(0, xBlock - 1);
		final int ymin = Math.max(0, yBlock - 1);
		final int xmax = Math.min(xBlocks, xBlock + 2);
		final int ymax = Math.min(yBlocks, yBlock + 2);
		
		for (int xx = xmin; xx < xmax; xx++)
		{
			for (int yy = ymin; yy < ymax; yy++)
			{
				size += peakGrid[xx][yy].size;
			}
		}
		final PeakResult[] list = new PeakResult[size];
		if (size != 0)
		{
			size = 0;
			for (int xx = xmin; xx < xmax; xx++)
			{
				for (int yy = ymin; yy < ymax; yy++)
				{
					if (peakGrid[xx][yy].size == 0)
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
}