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
package gdsc.smlm.filters;

import java.util.Arrays;

/**
 * Helper for the spot filter
 */
public class SpotFilterHelper
{
	private IntBlockSumFilter sumFilter = null;
	private int[] data = null;

	/**
	 * Count neighbours within a 2n+1 region around each spot.
	 * <p>
	 * This is performed using a block sum filter which may sub-optimal for small lists of spots.
	 * <p>
	 * The dimensions of the data will be extracted from the spot x/y coordinates.
	 *
	 * @param spots
	 *            the spots
	 * @param n
	 *            The block size
	 * @return the neighbour count for each spot
	 */
	public int[] countNeighbours(Spot[] spots, int n)
	{
		if (spots.length <= 1 || n <= 0)
			// No neighbours are possible
			return new int[spots.length];

		// Get the range for the sum filter using the limits.
		// This prevents filtering too large an image.
		int minx = spots[0].x;
		int maxx = minx;
		int miny = spots[0].y;
		int maxy = miny;
		for (int i = 1; i < spots.length; i++)
		{
			if (maxx < spots[i].x)
				maxx = spots[i].x;
			else if (minx > spots[i].x)
				minx = spots[i].x;
			if (maxy < spots[i].y)
				maxy = spots[i].y;
			else if (miny > spots[i].y)
				miny = spots[i].y;
		}

		return countNeighbours(spots, minx, miny, maxx, maxy, n);
	}

	/**
	 * Count neighbours within a 2n+1 region around each spot.
	 * <p>
	 * This is performed using a block sum filter which may sub-optimal for small lists of spots.
	 *
	 * @param spots
	 *            the spots
	 * @param width
	 *            The width of the data
	 * @param height
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @return the neighbour count for each spot
	 */
	public int[] countNeighbours(Spot[] spots, int width, int height, int n)
	{
		if (spots.length <= 1 || n <= 0)
			// No neighbours are possible
			return new int[spots.length];

		return countNeighbours(spots, 0, 0, width - 1, height - 1, n);
	}

	/**
	 * Count neighbours within a 2n+1 region around each spot.
	 * <p>
	 * This is performed using a block sum filter which may sub-optimal for small lists of spots.
	 *
	 * @param spots
	 *            the spots
	 * @param minx
	 *            the minx
	 * @param miny
	 *            the miny
	 * @param maxx
	 *            the maxx
	 * @param maxy
	 *            the maxy
	 * @param n
	 *            The block size
	 * @return the neighbour count for each spot
	 */
	private int[] countNeighbours(Spot[] spots, int minx, int miny, int maxx, int maxy, int n)
	{
		if (spots.length <= 1 || n <= 0)
			// No neighbours are possible
			return new int[spots.length];

		// Initialise
		if (sumFilter == null)
			sumFilter = new IntBlockSumFilter();
		int width = maxx - minx + 1;
		int height = maxy - miny + 1;
		int size = width * height;
		if (data == null || data.length < size)
			data = new int[size];
		else
			Arrays.fill(data, 0, size, 0);

		// Add the spots
		for (int i = 0; i < spots.length; i++)
			data[(spots[i].x - minx) + (spots[i].y - miny) * width] = 1;

		sumFilter.rollingBlockFilter(data, width, height, n);

		int[] count = new int[spots.length];
		for (int i = 0; i < spots.length; i++)
			// Subtract the actual spot from the count
			count[i] = data[(spots[i].x - minx) + (spots[i].y - miny) * width] - 1;

		return count;
	}

	/**
	 * Count neighbours within a 2n+1 region around each spot.
	 * <p>
	 * This is performed using a block sum filter which may sub-optimal for small lists of spots.
	 * <p>
	 * The dimensions of the data will be extracted from the spot x/y coordinates.
	 *
	 * @param spots
	 *            the spots
	 * @param n
	 *            The block size
	 * @return the neighbour count for each spot
	 */
	public static int[] runCountNeighbours(Spot[] spots, int n)
	{
		return new SpotFilterHelper().countNeighbours(spots, n);
	}

	/**
	 * Count neighbours within a 2n+1 region around each spot.
	 * <p>
	 * This is performed using a block sum filter which may sub-optimal for small lists of spots.
	 *
	 * @param spots
	 *            the spots
	 * @param width
	 *            The width of the data
	 * @param height
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @return the neighbour count for each spot
	 */
	public static int[] runCountNeighbours(Spot[] spots, int width, int height, int n)
	{
		return new SpotFilterHelper().countNeighbours(spots, width, height, n);
	}
}
