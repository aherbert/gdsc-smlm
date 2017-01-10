package gdsc.smlm.filters;

import java.util.Arrays;

import org.apache.commons.math3.util.FastMath;

import gdsc.core.utils.FixedIntList;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Computes the local maxima.
 */
public class NonMaximumSuppression implements Cloneable
{
	private float background = 0;
	private float fractionAboveBackground = 0;
	private float minimumHeight = 0;
	private float minimumWidth = 0;
	private boolean neighbourCheck = false;
	private boolean dataBuffer = true;

	private float[] newDataFloat = null;
	private int[] newDataInt = null;
	private FixedIntList resultsBuffer = null;
	private boolean[] maximaFlagBuffer = null;

	/**
	 * Compute the local-maxima within a 2n+1 block
	 * <p>
	 * Any maxima below the configured fraction above background are ignored. Fraction = (Max-background)/Max within the
	 * 2n+1 neighbourhood. Maxima below the minimum height (above background) or the minimum peak-width at half maximum
	 * in any dimension are ignored.
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @return Indices to the local maxima (index = maxx * y + x)
	 */
	public int[] maxFind(float[] data, int maxx, int maxy, int n)
	{
		FixedIntList results = getResultsBuffer(data.length / 4);
		boolean[] maximaFlag = getFlagBuffer(data.length);

		// Boundary control
		int xwidth = FastMath.min(n, maxx - 1);
		int ywidth = FastMath.min(n, maxy - 1);
		int xlimit = maxx - xwidth;
		int ylimit = maxy - ywidth;

		int[] offset = new int[(2 * xwidth + 1) * (2 * ywidth + 1) - 1];
		int[] xoffset = new int[offset.length];
		int[] yoffset = new int[offset.length];
		int d = 0;
		for (int y = -ywidth; y <= ywidth; y++)
			for (int x = -xwidth; x <= xwidth; x++)
				if (x != 0 || y != 0)
				{
					offset[d] = maxx * y + x;
					xoffset[d] = x;
					yoffset[d] = y;
					d++;
				}

		// Compare all points
		float heightThreshold = (float) getHeightThreshold();
		float floatBackground = (float) background;

		int index = 0;
		for (int y = 0; y < maxy; y++)
		{
			FIND_MAXIMUM: for (int x = 0; x < maxx; x++, index++)
			{
				float v = data[index];
				if (v < heightThreshold)
					continue;

				// Flag to indicate this pixels has a complete (2n+1) neighbourhood 
				boolean isInnerXY = (y >= ywidth && y < ylimit) && (x >= xwidth && x < xlimit);

				// Sweep neighbourhood
				if (isInnerXY)
				{
					for (int i = 0; i < offset.length; i++)
					{
						if (maximaFlag[index + offset[i]])
							continue FIND_MAXIMUM;
						if (data[index + offset[i]] > v)
							continue FIND_MAXIMUM;
					}
				}
				else
				{
					for (d = offset.length; d-- > 0;)
					{
						// Get the coords and check if it is within the data
						int yy = y + yoffset[d];
						int xx = x + xoffset[d];
						boolean isWithin = (yy >= 0 && yy < maxy) && (xx >= 0 && xx < maxx);

						if (isWithin)
						{
							if (maximaFlag[index + offset[d]])
								continue FIND_MAXIMUM;
							if (data[index + offset[d]] > v)
								continue FIND_MAXIMUM;
						}
					}
				}

				// Check the maximum width
				if (minimumWidth > 0)
				{
					// Get the width at half maximum.
					float v_half = floatHalfMaximum(floatBackground, v);
					int index2;
					// Scan right
					int x1 = x + 1;
					index2 = index + 1;
					while (x1 < maxx && data[index2] > v_half)
					{
						x1++;
						index2++;
					}
					// Scan left
					int x2 = x - 1;
					index2 = index - 1;
					while (x2 >= 0 && data[index2] > v_half)
					{
						x2--;
						index2--;
					}
					if (x1 - x2 < minimumWidth)
						continue;
					// Scan up
					int y1 = y + 1;
					index2 = index + maxx;
					while (y1 < maxy && data[index2] > v_half)
					{
						y1++;
						index2 += maxx;
					}
					// Scan down
					int y2 = y - 1;
					index2 = index - 1;
					while (y2 >= 0 && data[index2] > v_half)
					{
						y2--;
						index2 -= maxx;
					}
					if (y1 - y2 < minimumWidth)
						continue;
				}

				results.add(index);
				maximaFlag[index] = true;
			} // end FIND_MAXIMA
		}

		return results.toArray();
	}

	/**
	 * Compute the local-maxima within a 2n+1 block
	 * An inner boundary of N is ignored as potential maxima.
	 * <p>
	 * Any maxima below the configured fraction above background are ignored. Fraction = (Max-background)/Max within the
	 * 2n+1 neighbourhood. Maxima below the minimum height (above background) or the minimum peak-width at half maximum
	 * in any dimension are ignored.
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @return Indices to the local maxima (index = maxx * y + x)
	 */
	public int[] maxFindInternal(float[] data, int maxx, int maxy, int n)
	{
		FixedIntList results = getResultsBuffer(data.length / 4);
		boolean[] maximaFlag = getFlagBuffer(data.length);

		// Boundary control
		int xwidth = FastMath.min(n, maxx - 1);
		int ywidth = FastMath.min(n, maxy - 1);

		int[] offset = new int[(2 * xwidth + 1) * (2 * ywidth + 1) - 1];
		int d = 0;
		for (int y = -ywidth; y <= ywidth; y++)
			for (int x = -xwidth; x <= xwidth; x++)
				if (x != 0 || y != 0)
				{
					offset[d] = maxx * y + x;
					d++;
				}

		// Compare all points
		float heightThreshold = (float) getHeightThreshold();
		float floatBackground = (float) background;

		for (int y = n; y < maxy - n; y++)
		{
			int index = y * maxx + n;
			FIND_MAXIMUM: for (int x = n; x < maxx - n; x++, index++)
			{
				float v = data[index];
				if (v < heightThreshold)
					continue;

				// Sweep neighbourhood - 
				// No check for boundaries as this should be an internal sweep.
				for (int i = 0; i < offset.length; i++)
				{
					if (maximaFlag[index + offset[i]])
						continue FIND_MAXIMUM;
					if (data[index + offset[i]] > v)
						continue FIND_MAXIMUM;
				}

				// Check the maximum width
				if (minimumWidth > 0)
				{
					float v_half = floatHalfMaximum(floatBackground, v);
					int index2;
					// Scan right
					int x1 = x + 1;
					index2 = index + 1;
					while (x1 < maxx && data[index2] > v_half)
					{
						x1++;
						index2++;
					}
					// Scan left
					int x2 = x - 1;
					index2 = index - 1;
					while (x2 >= 0 && data[index2] > v_half)
					{
						x2--;
						index2--;
					}
					if (x1 - x2 < minimumWidth)
						continue;
					// Scan up
					int y1 = y + 1;
					index2 = index + maxx;
					while (y1 < maxy && data[index2] > v_half)
					{
						y1++;
						index2 += maxx;
					}
					// Scan down
					int y2 = y - 1;
					index2 = index - 1;
					while (y2 >= 0 && data[index2] > v_half)
					{
						y2--;
						index2 -= maxx;
					}
					if (y1 - y2 < minimumWidth)
						continue;
				}

				results.add(index);
				maximaFlag[index] = true;
			} // end FIND_MAXIMA
		}

		return results.toArray();
	}

	/**
	 * Compute the local-maxima within a 2n+1 block
	 * An inner boundary is ignored as potential maxima.
	 * <p>
	 * Any maxima below the configured fraction above background are ignored. Fraction = (Max-background)/Max within the
	 * 2n+1 neighbourhood. Maxima below the minimum height (above background) or the minimum peak-width at half maximum
	 * in any dimension are ignored.
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @param border
	 *            The internal border
	 * @return Indices to the local maxima (index = maxx * y + x)
	 */
	public int[] maxFindInternal(float[] data, int maxx, int maxy, int n, int border)
	{
		if (n == border)
			// Faster algorithm as there is no requirement for bounds checking.
			return maxFindInternal(data, maxx, maxy, n);

		FixedIntList results = getResultsBuffer(data.length / 4);
		boolean[] maximaFlag = getFlagBuffer(data.length);

		// Boundary control
		int xwidth = FastMath.min(n, maxx - 1);
		int ywidth = FastMath.min(n, maxy - 1);
		int xlimit = maxx - xwidth;
		int ylimit = maxy - ywidth;

		int[] offset = new int[(2 * xwidth + 1) * (2 * ywidth + 1) - 1];
		int[] xoffset = new int[offset.length];
		int[] yoffset = new int[offset.length];
		int d = 0;
		for (int y = -ywidth; y <= ywidth; y++)
			for (int x = -xwidth; x <= xwidth; x++)
				if (x != 0 || y != 0)
				{
					offset[d] = maxx * y + x;
					xoffset[d] = x;
					yoffset[d] = y;
					d++;
				}

		// All blocks fit within the border
		boolean inner = (n < border);

		// Compare all points
		float heightThreshold = (float) getHeightThreshold();
		float floatBackground = (float) background;

		for (int y = border; y < maxy - border; y++)
		{
			int index = y * maxx + border;
			FIND_MAXIMUM: for (int x = border; x < maxx - border; x++, index++)
			{

				float v = data[index];
				if (v < heightThreshold)
					continue;

				if (inner)
				{
					for (int i = 0; i < offset.length; i++)
					{
						if (maximaFlag[index + offset[i]])
							continue FIND_MAXIMUM;
						if (data[index + offset[i]] > v)
							continue FIND_MAXIMUM;
					}
				}
				else
				{
					// Flag to indicate this pixels has a complete (2n+1) neighbourhood 
					boolean isInnerXY = (y >= ywidth && y < ylimit) && (x >= xwidth && x < xlimit);

					// Sweep neighbourhood
					for (d = offset.length; d-- > 0;)
					{
						boolean isWithin = isInnerXY;
						if (!isWithin)
						{
							// Get the coords and check if it is within the data
							int yy = y + yoffset[d];
							int xx = x + xoffset[d];
							isWithin = (yy >= 0 && yy < maxy) && (xx >= 0 && xx < maxx);
						}

						if (isWithin)
						{
							if (maximaFlag[index + offset[d]])
								continue FIND_MAXIMUM;
							if (data[index + offset[d]] > v)
								continue FIND_MAXIMUM;
						}
					}
				}

				// Check the maximum width
				if (minimumWidth > 0)
				{
					float v_half = floatHalfMaximum(floatBackground, v);
					int index2;
					// Scan right
					int x1 = x + 1;
					index2 = index + 1;
					while (x1 < maxx && data[index2] > v_half)
					{
						x1++;
						index2++;
					}
					// Scan left
					int x2 = x - 1;
					index2 = index - 1;
					while (x2 >= 0 && data[index2] > v_half)
					{
						x2--;
						index2--;
					}
					if (x1 - x2 < minimumWidth)
						continue;
					// Scan up
					int y1 = y + 1;
					index2 = index + maxx;
					while (y1 < maxy && data[index2] > v_half)
					{
						y1++;
						index2 += maxx;
					}
					// Scan down
					int y2 = y - 1;
					index2 = index - 1;
					while (y2 >= 0 && data[index2] > v_half)
					{
						y2--;
						index2 -= maxx;
					}
					if (y1 - y2 < minimumWidth)
						continue;
				}

				results.add(index);
				maximaFlag[index] = true;
			} // end FIND_MAXIMA
		}

		return results.toArray();
	}

	/**
	 * Compute the local-maxima within a 2n+1 block.
	 * <p>
	 * Any maxima below the configured fraction above background are ignored. Fraction = (Max-background)/Max within the
	 * 2n+1 neighbourhood. Maxima below the minimum height (above background) or the minimum peak-width at half maximum
	 * in any dimension are ignored.
	 * <p>
	 * Uses the 2D block algorithm of Neubeck and Van Gool (2006).
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @return Indices to the local maxima (index = maxx * y + x)
	 */
	public int[] blockFind(float[] data, int maxx, int maxy, int n)
	{
		if (n == 1)
			// optimised version for the special case 
			return blockFind3x3(data, maxx, maxy);

		return blockFindNxN(data, maxx, maxy, n);
	}

	/**
	 * Compute the local-maxima within a 2n+1 block.
	 * An inner boundary of N is ignored as potential maxima.
	 * <p>
	 * Any maxima below the configured fraction above background are ignored. Fraction = (Max-background)/Max within the
	 * 2n+1 neighbourhood. Maxima below the minimum height (above background) or the minimum peak-width at half maximum
	 * in any dimension are ignored.
	 * <p>
	 * Uses the 2D block algorithm of Neubeck and Van Gool (2006).
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @param border
	 *            The internal border
	 * @return Indices to the local maxima (index = maxx * y + x)
	 */
	public int[] blockFindInternal(float[] data, int maxx, int maxy, int n, int border)
	{
		if (n == 1)
			// optimised version for the special case 
			return blockFind3x3Internal(data, maxx, maxy, border);

		return blockFindNxNInternal(data, maxx, maxy, n, border);
	}

	/**
	 * Compute the local-maxima within a 2n+1 block.
	 * <p>
	 * Any maxima below the configured fraction above background are ignored. Fraction = (Max-background)/Max within the
	 * 2n+1 neighbourhood. Maxima below the minimum height (above background) or the minimum peak-width at half maximum
	 * in any dimension are ignored.
	 * <p>
	 * Uses the 2D block algorithm of Neubeck and Van Gool (2006).
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @return Indices to the local maxima (index = maxx * y + x)
	 */
	public int[] blockFindNxN(float[] data, int maxx, int maxy, int n)
	{
		int[] blockMaxima = findBlockMaximaNxN(data, maxx, maxy, n);
		int nMaxima = 0;
		float heightThreshold = (float) getHeightThreshold();
		float floatBackground = (float) background;

		boolean[] maximaFlag = null;
		if (isNeighbourCheck())
			maximaFlag = getFlagBuffer(data.length);

		FIND_MAXIMUM: for (int index : blockMaxima)
		{

			float v = data[index];

			if (v < heightThreshold)
				continue;

			int x = index % maxx;
			int y = index / maxx;

			// Compare the maxima to the surroundings. Ignore the block region already processed.
			//
			//.......(mi-n,mj-n)----------------------------------(mi+n,mj-n)
			//............|aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa|
			//............|aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa|
			//............|aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa|
			//............|aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa|
			//............|aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa|
			//............|bbbbbbbb.(i,j)-------------(i+n,j).ccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbbbbb|.....(mi,mj)......|ccccccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbb.(i,j+n)-----------(i+n,j+n).ccccccc|
			//............|dddddddddddddddddddddddddddddddddddddddddddd|
			//............|dddddddddddddddddddddddddddddddddddddddddddd|
			//............|dddddddddddddddddddddddddddddddddddddddddddd|
			//......(mi-n,mj+n)----------------------------------(mi+n,mj+n)
			//
			// This must be done without over-running boundaries
			//int steps = 0;
			int mi = x;
			int mj = y;
			int i = (n + 1) * (mi / (n + 1));
			int j = (n + 1) * (mj / (n + 1));
			int i_plus_n = FastMath.min(i + n + 1, maxx - 1);
			int j_plus_n = FastMath.min(j + n + 1, maxy - 1);
			int mi_minus_n = FastMath.max(mi - n, 0);
			int mi_plus_n = FastMath.min(mi + n, maxx - 1);
			int mj_minus_n = FastMath.max(mj - n, 0);
			int mj_plus_n = FastMath.min(mj + n, maxy - 1);

			// A
			for (int jj = mj_minus_n; jj < j; jj++)
			{
				int indexStart = jj * maxx + mi_minus_n;
				int indexEnd = jj * maxx + mi_plus_n;
				for (; indexStart <= indexEnd; indexStart++)
				{
					//steps++;
					if (v < data[indexStart])
						continue FIND_MAXIMUM;
				}
			}
			for (int jj = j; jj < j_plus_n; jj++)
			{
				// B
				{
					int indexStart = jj * maxx + mi_minus_n;
					int indexEnd = jj * maxx + i;
					for (; indexStart < indexEnd; indexStart++)
					{
						//steps++;
						if (v < data[indexStart])
							continue FIND_MAXIMUM;
					}
				}

				// C
				{
					int indexStart = jj * maxx + i_plus_n;
					int indexEnd = jj * maxx + mi_plus_n;
					for (; indexStart <= indexEnd; indexStart++)
					{
						//steps++;
						if (v < data[indexStart])
							continue FIND_MAXIMUM;
					}
				}
			}
			// D
			for (int jj = j_plus_n; jj <= mj_plus_n; jj++)
			{
				int indexStart = jj * maxx + mi_minus_n;
				int indexEnd = jj * maxx + mi_plus_n;
				for (; indexStart <= indexEnd; indexStart++)
				{
					//steps++;
					if (v < data[indexStart])
						continue FIND_MAXIMUM;
				}
			}
			//System.out.printf("%.2f @ %d,%d. Steps = %d, max = %b\n", v, x, y, steps, isMax);

			if (isNeighbourCheck())
			{
				// A
				for (int jj = mj_minus_n; jj < j; jj++)
				{
					int indexStart = jj * maxx + mi_minus_n;
					int indexEnd = jj * maxx + mi_plus_n;
					for (; indexStart <= indexEnd; indexStart++)
					{
						//steps++;
						if (maximaFlag[indexStart])
							continue FIND_MAXIMUM;
					}
				}
				for (int jj = j; jj < j_plus_n; jj++)
				{
					// B
					{
						int indexStart = jj * maxx + mi_minus_n;
						int indexEnd = jj * maxx + i;
						for (; indexStart < indexEnd; indexStart++)
						{
							//steps++;
							if (maximaFlag[indexStart])
								continue FIND_MAXIMUM;
						}
					}

					// C
					{
						int indexStart = jj * maxx + i_plus_n;
						int indexEnd = jj * maxx + mi_plus_n;
						for (; indexStart <= indexEnd; indexStart++)
						{
							//steps++;
							if (maximaFlag[indexStart])
								continue FIND_MAXIMUM;
						}
					}
				}
				// D
				for (int jj = j_plus_n; jj <= mj_plus_n; jj++)
				{
					int indexStart = jj * maxx + mi_minus_n;
					int indexEnd = jj * maxx + mi_plus_n;
					for (; indexStart <= indexEnd; indexStart++)
					{
						//steps++;
						if (maximaFlag[indexStart])
							continue FIND_MAXIMUM;
					}
				}
			}

			// Check the maximum width
			if (minimumWidth > 0)
			{
				float v_half = floatHalfMaximum(floatBackground, v);
				int index2;
				// Scan right
				int x1 = x + 1;
				index2 = index + 1;
				while (x1 < maxx && data[index2] > v_half)
				{
					x1++;
					index2++;
				}
				// Scan left
				int x2 = x - 1;
				index2 = index - 1;
				while (x2 >= 0 && data[index2] > v_half)
				{
					x2--;
					index2--;
				}
				if (x1 - x2 < minimumWidth)
					continue;
				// Scan up
				int y1 = y + 1;
				index2 = index + maxx;
				while (y1 < maxy && data[index2] > v_half)
				{
					y1++;
					index2 += maxx;
				}
				// Scan down
				int y2 = y - 1;
				index2 = index - 1;
				while (y2 >= 0 && data[index2] > v_half)
				{
					y2--;
					index2 -= maxx;
				}
				if (y1 - y2 < minimumWidth)
					continue;
			}

			//System.out.printf("blockFind [%d,%d]\n", mi, mj);

			// Re-use storage space
			blockMaxima[nMaxima++] = index;
			if (isNeighbourCheck())
				maximaFlag[index] = true;
		} // end FIND_MAXIMA
		  //System.out.printf("---\n");

		return truncate(blockMaxima, nMaxima);
	}

	/**
	 * Compute the local-maxima within a 2n+1 block.
	 * An inner boundary of N is ignored as potential maxima.
	 * <p>
	 * Any maxima below the configured fraction above background are ignored. Fraction = (Max-background)/Max within the
	 * 2n+1 neighbourhood. Maxima below the minimum height (above background) or the minimum peak-width at half maximum
	 * in any dimension are ignored.
	 * <p>
	 * Uses the 2D block algorithm of Neubeck and Van Gool (2006).
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @param border
	 *            The internal border
	 * @return Indices to the local maxima (index = maxx * y + x)
	 */
	public int[] blockFindNxNInternal(float[] data, int maxx, int maxy, int n, int border)
	{
		int[] blockMaxima = findBlockMaximaNxNInternal(data, maxx, maxy, n, border);
		int nMaxima = 0;
		float heightThreshold = (float) getHeightThreshold();
		float floatBackground = (float) background;

		boolean[] maximaFlag = null;
		if (isNeighbourCheck())
			maximaFlag = getFlagBuffer(data.length);

		FIND_MAXIMUM: for (int index : blockMaxima)
		{
			float v = data[index];

			if (v < heightThreshold)
				continue;

			int x = index % maxx;
			int y = index / maxx;

			// Compare the maxima to the surroundings. Ignore the block region already processed.
			//
			//.......(mi-n,mj-n)----------------------------------(mi+n,mj-n)
			//............|aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa|
			//............|aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa|
			//............|aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa|
			//............|aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa|
			//............|aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa|
			//............|bbbbbbbb.(i,j)-------------(i+n,j).ccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbbbbb|.....(mi,mj)......|ccccccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbb.(i,j+n)-----------(i+n,j+n).ccccccc|
			//............|dddddddddddddddddddddddddddddddddddddddddddd|
			//............|dddddddddddddddddddddddddddddddddddddddddddd|
			//............|dddddddddddddddddddddddddddddddddddddddddddd|
			//......(mi-n,mj+n)----------------------------------(mi+n,mj+n)
			//
			// No check for over-running boundaries since this is the internal version
			//int steps = 0;
			int mi = x;
			int mj = y;
			int i = (n + 1) * ((mi - border) / (n + 1)) + border; // Blocks n+1 wide
			int j = (n + 1) * ((mj - border) / (n + 1)) + border; // Blocks n+1 wide
			// The block boundaries will have been truncated on the final block. Ensure this is swept
			int i_plus_n = FastMath.min(i + n + 1, maxx - border - 1);
			int j_plus_n = FastMath.min(j + n + 1, maxy - border - 1);
			int mi_minus_n = FastMath.max(mi - n, 0);
			int mi_plus_n = FastMath.min(mi + n, maxx - 1);
			int mj_minus_n = FastMath.max(mj - n, 0);
			int mj_plus_n = FastMath.min(mj + n, maxy - 1);

			//System.out.printf("Block [%d,%d] => [%d,%d]\n", x, y, i, j);

			// A
			for (int jj = mj_minus_n; jj < j; jj++)
			{
				int indexStart = jj * maxx + mi_minus_n;
				int indexEnd = jj * maxx + mi_plus_n;
				for (; indexStart <= indexEnd; indexStart++)
				{
					//steps++;
					if (v < data[indexStart])
						continue FIND_MAXIMUM;
				}
			}
			for (int jj = j; jj < j_plus_n; jj++)
			{
				// B
				{
					int indexStart = jj * maxx + mi_minus_n;
					int indexEnd = jj * maxx + i;
					for (; indexStart < indexEnd; indexStart++)
					{
						//steps++;
						if (v < data[indexStart])
							continue FIND_MAXIMUM;
					}
				}

				// C
				{
					int indexStart = jj * maxx + i_plus_n;
					int indexEnd = jj * maxx + mi_plus_n;
					for (; indexStart <= indexEnd; indexStart++)
					{
						//steps++;
						if (v < data[indexStart])
							continue FIND_MAXIMUM;
					}
				}
			}
			// D
			for (int jj = j_plus_n; jj <= mj_plus_n; jj++)
			{
				int indexStart = jj * maxx + mi_minus_n;
				int indexEnd = jj * maxx + mi_plus_n;
				for (; indexStart <= indexEnd; indexStart++)
				{
					//steps++;
					if (v < data[indexStart])
						continue FIND_MAXIMUM;
				}
			}
			//System.out.printf("%.2f @ %d,%d. Steps = %d, max = %b\n", v, x, y, steps, isMax);

			if (isNeighbourCheck())
			{
				// A
				for (int jj = mj_minus_n; jj < j; jj++)
				{
					int indexStart = jj * maxx + mi_minus_n;
					int indexEnd = jj * maxx + mi_plus_n;
					for (; indexStart <= indexEnd; indexStart++)
					{
						//steps++;
						if (maximaFlag[indexStart])
							continue FIND_MAXIMUM;
					}
				}
				for (int jj = j; jj < j_plus_n; jj++)
				{
					// B
					{
						int indexStart = jj * maxx + mi_minus_n;
						int indexEnd = jj * maxx + i;
						for (; indexStart < indexEnd; indexStart++)
						{
							//steps++;
							if (maximaFlag[indexStart])
								continue FIND_MAXIMUM;
						}
					}

					// C
					{
						int indexStart = jj * maxx + i_plus_n;
						int indexEnd = jj * maxx + mi_plus_n;
						for (; indexStart <= indexEnd; indexStart++)
						{
							//steps++;
							if (maximaFlag[indexStart])
								continue FIND_MAXIMUM;
						}
					}
				}
				// D
				for (int jj = j_plus_n; jj <= mj_plus_n; jj++)
				{
					int indexStart = jj * maxx + mi_minus_n;
					int indexEnd = jj * maxx + mi_plus_n;
					for (; indexStart <= indexEnd; indexStart++)
					{
						//steps++;
						if (maximaFlag[indexStart])
							continue FIND_MAXIMUM;
					}
				}
			}

			// Check the maximum width
			if (minimumWidth > 0)
			{
				float v_half = floatHalfMaximum(floatBackground, v);
				int index2;
				// Scan right
				int x1 = x + 1;
				index2 = index + 1;
				while (x1 < maxx && data[index2] > v_half)
				{
					x1++;
					index2++;
				}
				// Scan left
				int x2 = x - 1;
				index2 = index - 1;
				while (x2 >= 0 && data[index2] > v_half)
				{
					x2--;
					index2--;
				}
				if (x1 - x2 < minimumWidth)
					continue;
				// Scan up
				int y1 = y + 1;
				index2 = index + maxx;
				while (y1 < maxy && data[index2] > v_half)
				{
					y1++;
					index2 += maxx;
				}
				// Scan down
				int y2 = y - 1;
				index2 = index - 1;
				while (y2 >= 0 && data[index2] > v_half)
				{
					y2--;
					index2 -= maxx;
				}
				if (y1 - y2 < minimumWidth)
					continue;
			}

			//System.out.printf("blockFind [%d,%d]\n", mi, mj);

			// Re-use storage space
			blockMaxima[nMaxima++] = index;
			if (isNeighbourCheck())
				maximaFlag[index] = true;
		} // end FIND_MAXIMA
		  //System.out.printf("---\n");

		return truncate(blockMaxima, nMaxima);
	}

	/**
	 * Search the data for the index of the maximum in each block of size n*n.
	 * <p>
	 * E.g. Max [ (i,i+n) x (i,j+n) ] for (i=0; i&lt;maxx; i+=n) x (j=0, j&lt;maxy; j+=n)
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @return The maxima indices
	 */
	public int[] findBlockMaximaNxN(float[] data, int maxx, int maxy, int n)
	{
		// Include the actual pixel in the block
		// This makes the block search inclusive: i=0; i<=n; i++
		n++;

		// The number of blocks in x and y
		int xblocks = getBlocks(maxx, n);
		int yblocks = getBlocks(maxy, n);

		// The final index in each dimension (where an incomplete block is found).
		// This equals maxx/maxy if the number of blocks fits exactly.
		int xfinal = n * (maxx / n);
		int yfinal = n * (maxy / n);

		int[] maxima = new int[xblocks * yblocks];

		int block = 0;
		for (int y = 0; y < maxy; y += n)
		{
			for (int x = 0; x < maxx; x += n)
			{
				// Find the sweep size in each direction
				int xsize = (x != xfinal) ? n : maxx - xfinal;
				int ysize = (y != yfinal) ? n : maxy - yfinal;

				int index = y * maxx + x;
				int maxIndex = index;
				float max = data[maxIndex];

				while (ysize-- > 0)
				{
					for (int x2 = xsize; x2-- > 0;)
					{
						if (max < data[index])
						{
							max = data[index];
							maxIndex = index;
						}
						index++;
					}
					index += maxx - xsize;
				}

				maxima[block++] = maxIndex;
			}
		}

		return maxima;
	}

	/**
	 * Search the data for the index of the maximum in each block of size n*n.
	 * An inner boundary of N is ignored as potential maxima.
	 * <p>
	 * E.g. Max [ (i,i+n) x (i,j+n) ] for (i=n; i&lt;maxx-n; i+=n) x (j=n, j&lt;maxy-n; j+=n)
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @param border
	 *            The internal border
	 * @return The maxima indices
	 */
	public int[] findBlockMaximaNxNInternal(float[] data, int maxx, int maxy, int n, int border)
	{
		// Include the actual pixel in the block
		// This makes the block search inclusive: i=0; i<=n; i++
		n++;

		// The number of blocks in x and y. Subtract the boundary blocks.
		int xblocks = getBlocks(maxx - border, n);
		int yblocks = getBlocks(maxy - border, n);

		if (xblocks < 1 || yblocks < 1)
			return new int[0];

		// The final index in each dimension (where an incomplete block is found).
		// This equals maxx/maxy if the number of blocks fits exactly.
		int xfinal = n * ((maxx - 2 * border) / n) + border;
		int yfinal = n * ((maxy - 2 * border) / n) + border;

		int[] maxima = new int[xblocks * yblocks];

		int block = 0;
		for (int y = border; y < maxy - border; y += n)
		{
			for (int x = border; x < maxx - border; x += n)
			{
				//System.out.printf("Block [%d,%d]\n", x, y);

				// Find the sweep size in each direction
				int xsize = (x != xfinal) ? n : maxx - xfinal - border;
				int ysize = (y != yfinal) ? n : maxy - yfinal - border;

				int index = y * maxx + x;
				int maxIndex = index;
				float max = data[maxIndex];

				while (ysize-- > 0)
				{
					for (int x2 = xsize; x2-- > 0;)
					{
						if (max < data[index])
						{
							max = data[index];
							maxIndex = index;
						}
						index++;
					}
					index += maxx - xsize;
				}

				maxima[block++] = maxIndex;
			}
		}

		return truncate(maxima, block);
	}

	/**
	 * Search the data for the index of the maximum in each block of size n*n.
	 * <p>
	 * E.g. Max [ (i,i+n) x (i,j+n) ] for (i=0; i&lt;maxx; i+=n) x (j=0, j&lt;maxy; j+=n)
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @return The maxima indices
	 */
	public int[] findBlockMaximaOptimised(float[] data, int maxx, int maxy, int n)
	{
		// =-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=
		// The following code DOES NOT increase speed. 
		// I believe the Java optimiser cannot unroll all the nested loops.
		// =-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=

		// Include the actual pixel in the block
		// This makes the block search inclusive: i=0; i<=n; i++
		n++;

		// The number of blocks in x and y
		int xblocks = getBlocks(maxx, n);
		int yblocks = getBlocks(maxy, n);

		// The final index in each dimension (where an incomplete block is found).
		// This equals maxx/maxy if the number of blocks fits exactly.
		int xfinal = n * (maxx / n);
		int yfinal = n * (maxy / n);

		int[] maxima = new int[xblocks * yblocks];

		// Sweep 4 regions:
		// x.............xfinal  
		// |.............|..maxx
		// |.............|..|  
		// aaaaaaaaaaaaaaabbb-- y
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb--yfinal
		// cccccccccccccccddd
		// cccccccccccccccddd--maxy

		int block = 0;
		for (int y = 0; y < yfinal; y += n)
		{
			// A
			for (int x = 0; x < xfinal; x += n)
			{
				// Find the sweep size in each direction
				int xsize = n;
				int ysize = n;

				int index = y * maxx + x;
				int maxIndex = index;
				float max = data[maxIndex];

				while (ysize-- > 0)
				{
					for (int x2 = xsize; x2-- > 0;)
					{
						if (max < data[index])
						{
							max = data[index];
							maxIndex = index;
						}
						index++;
					}
					index += maxx - xsize;
				}

				maxima[block++] = maxIndex;
			}

			// B
			if (xfinal != maxx)
			{
				// Find the sweep size in each direction
				int xsize = maxx - xfinal;
				int ysize = n;

				int index = y * maxx + xfinal;
				int maxIndex = index;
				float max = data[maxIndex];

				while (ysize-- > 0)
				{
					for (int x2 = xsize; x2-- > 0;)
					{
						if (max < data[index])
						{
							max = data[index];
							maxIndex = index;
						}
						index++;
					}
					index += maxx - xsize;
				}

				maxima[block++] = maxIndex;
			}
		}
		if (yfinal != maxy)
		{
			// C
			for (int x = 0; x < xfinal; x += n)
			{
				// Find the sweep size in each direction
				int xsize = n;
				int ysize = maxy - yfinal;

				int index = yfinal * maxx + x;
				int maxIndex = index;
				float max = data[maxIndex];

				while (ysize-- > 0)
				{
					for (int x2 = xsize; x2-- > 0;)
					{
						if (max < data[index])
						{
							max = data[index];
							maxIndex = index;
						}
						index++;
					}
					index += maxx - xsize;
				}

				maxima[block++] = maxIndex;
			}

			// D
			if (xfinal != maxx)
			{
				// Find the sweep size in each direction
				int xsize = maxx - xfinal;
				int ysize = maxy - yfinal;

				int index = yfinal * maxx + xfinal;
				int maxIndex = index;
				float max = data[maxIndex];

				while (ysize-- > 0)
				{
					for (int x2 = xsize; x2-- > 0;)
					{
						if (max < data[index])
						{
							max = data[index];
							maxIndex = index;
						}
						index++;
					}
					index += maxx - xsize;
				}

				maxima[block++] = maxIndex;
			}
		}

		return maxima;
	}

	/**
	 * Search the data for the index of the maximum in each block of size 2*2.
	 * <p>
	 * E.g. Max [ (i,i+n) x (i,j+n) ] for (i=0; i&lt;maxx; i+=n) x (j=0, j&lt;maxy; j+=n)
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @return The maxima indices
	 */
	public int[] findBlockMaxima2x2(float[] data, int maxx, int maxy)
	{
		// Optimised for 2x2 block
		//final int n = 2;

		// The number of blocks in x and y
		int xblocks = getBlocks(maxx, 2);
		int yblocks = getBlocks(maxy, 2);

		// The final index in each dimension (where an incomplete block is found).
		// This equals maxx/maxy if the number of blocks fits exactly.
		int xfinal = 2 * (maxx / 2);
		int yfinal = 2 * (maxy / 2);

		// TODO - Try this by expanding the data if xfinal != maxx || yfinal != maxy
		// This will allow less management of boundaries.

		int[] maxima = new int[xblocks * yblocks];

		// Sweep 4 regions:
		// x.............xfinal  
		// |.............|..maxx
		// |.............|..|  
		// aaaaaaaaaaaaaaabbb-- y
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb--yfinal
		// cccccccccccccccddd
		// cccccccccccccccddd--maxy

		int block = 0;
		for (int y = 0; y < yfinal; y += 2)
		{
			// A
			int xindex = y * maxx;
			for (int x = 0; x < xfinal; x += 2)
			{
				// Compare 2x2 block
				int max1 = (data[xindex] > data[xindex + 1]) ? xindex : xindex + 1;
				int max2 = (data[xindex + maxx] > data[xindex + maxx + 1]) ? xindex + maxx : xindex + maxx + 1;

				maxima[block++] = (data[max1] > data[max2]) ? max1 : max2;

				xindex += 2;
			}
			// B
			if (xfinal != maxx)
			{
				// Compare 1x2 block
				int index = y * maxx + xfinal;
				maxima[block++] = (data[index] > data[index + maxx]) ? index : index + maxx;
			}
		}
		if (yfinal != maxy)
		{
			// C
			int xindex = yfinal * maxx;
			for (int x = 0; x < xfinal; x += 2)
			{
				// Compare 2x1 block
				maxima[block++] = (data[xindex] > data[xindex + 1]) ? xindex : xindex + 1;
				xindex += 2;
			}
			// D
			if (xfinal != maxx)
			{
				// Compare 1x1 block
				maxima[block++] = yfinal * maxx + xfinal;
			}
		}

		return maxima;
	}

	/**
	 * Compute the local-maxima within a 3x3 block
	 * <p>
	 * Any maxima below the configured fraction above background are ignored. Fraction = (Max-background)/Max within the
	 * 3x3 neighbourhood. Maxima below the minimum height (above background) or the minimum peak-width at half maximum
	 * in any dimension are ignored.
	 * <p>
	 * Note: The height thresholding step is ignored if the height threshold is zero.
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @return Indices to the local maxima (index = maxx * y + x)
	 */
	public int[] blockFind3x3(float[] data, int maxx, int maxy)
	{
		// The number of blocks in x and y
		int xblocks = getBlocks(maxx, 2);
		int yblocks = getBlocks(maxy, 2);

		// The final index in each dimension (where an incomplete block is found).
		// This equals maxx/maxy if the number of blocks fits exactly.
		int xfinal = 2 * (maxx / 2);
		int yfinal = 2 * (maxy / 2);

		// Expand the canvas
		int newx = maxx + (maxx - xfinal) + 2;
		int newy = maxy + (maxy - yfinal) + 2;
		data = expand(data, maxx, maxy, newx, newy);

		int[] maxima = new int[xblocks * yblocks];

		float heightThreshold = (float) getHeightThreshold();
		float floatBackground = (float) background;
		boolean validations = (heightThreshold > 0 || minimumWidth > 0 || isNeighbourCheck());

		boolean[] maximaFlag = null;
		if (isNeighbourCheck())
			maximaFlag = getFlagBuffer(data.length);

		// Compare 2x2 block
		// ....
		// .AB.
		// .CD.
		// ....

		// Create the scan arrays to search the remaining 5 locations around the 2x2 block maxima
		// ....  aaa.  .bbb  ....  .... 
		// .AB.  aA..  ..Bb  c...  ...d 
		// .CD.  a...  ...b  cC..  ..Dd 
		// ....  ....  ....  ccc.  .ddd 
		int[] a = new int[] { -newx - 1, -newx, -newx + 1, -1, +newx - 1 };
		int[] b = new int[] { -newx - 1, -newx, -newx + 1, +1, +newx + 1 };
		int[] c = new int[] { -newx - 1, -1, +newx - 1, +newx, +newx + 1 };
		int[] d = new int[] { -newx + 1, +1, +newx - 1, +newx, +newx + 1 };

		int nMaxima = 0;
		for (int y = 0; y < maxy; y += 2)
		{
			int xindex = (y + 1) * newx + 1;
			FIND_MAXIMUM: for (int x = 0; x < maxx; x += 2, xindex += 2)
			{
				int[] scan = a;
				int maxIndex = xindex;
				if (data[maxIndex] < data[xindex + 1])
				{
					scan = b;
					maxIndex = xindex + 1;
				}
				if (data[maxIndex] < data[xindex + newx])
				{
					scan = c;
					maxIndex = xindex + newx;
				}
				if (data[maxIndex] < data[xindex + newx + 1])
				{
					scan = d;
					maxIndex = xindex + newx + 1;
				}

				// Check the remaining region
				for (int offset : scan)
				{
					if (data[maxIndex] < data[maxIndex + offset])
						continue FIND_MAXIMUM;
				}

				if (validations)
				{
					if (data[maxIndex] < heightThreshold)
						continue;

					if (isNeighbourCheck())
					{
						for (int offset : scan)
						{
							if (maximaFlag[maxIndex + offset])
								continue FIND_MAXIMUM;
						}
					}

					// Check the maximum width
					if (minimumWidth > 0)
					{
						int x0 = maxIndex % maxx;
						int y0 = maxIndex / maxx;

						// Get the width at half maximum.
						float v_half = floatHalfMaximum(floatBackground, data[maxIndex]);
						int index2;
						// Scan right
						int x1 = x0 + 1;
						index2 = maxIndex + 1;
						while (x1 < maxx && data[index2] > v_half)
						{
							x1++;
							index2++;
						}
						// Scan left
						int x2 = x0 - 1;
						index2 = maxIndex - 1;
						while (x2 >= 0 && data[index2] > v_half)
						{
							x2--;
							index2--;
						}
						if (x1 - x2 < minimumWidth)
							continue;
						// Scan up
						int y1 = y0 + 1;
						index2 = maxIndex + newx;
						while (y1 < maxy && data[index2] > v_half)
						{
							y1++;
							index2 += newx;
						}
						// Scan down
						int y2 = y0 - 1;
						index2 = maxIndex - newx;
						while (y2 >= 0 && data[index2] > v_half)
						{
							y2--;
							index2 -= newx;
						}
						if (y1 - y2 < minimumWidth)
							continue;
					}

					if (isNeighbourCheck())
						maximaFlag[maxIndex] = true;
				}

				// Remap the maxima
				int xx = maxIndex % newx;
				int yy = maxIndex / newx;

				//System.out.printf("blockFind3x3 [%d,%d]\n", xx-1, yy-1);
				maxima[nMaxima++] = (yy - 1) * maxx + xx - 1;
			} // end FIND_MAXIMA
		}

		return truncate(maxima, nMaxima);
	}

	/**
	 * Compute the local-maxima within a 3x3 block.
	 * An inner boundary of 1 is ignored as potential maxima on the top and left, and a boundary of 1 or 2 on the right
	 * or bottom (depending if the image is even/odd dimensions).
	 * <p>
	 * Any maxima below the configured fraction above background are ignored. Fraction = (Max-background)/Max within the
	 * 3x3 neighbourhood. Maxima below the minimum height (above background) or the minimum peak-width at half maximum
	 * in any dimension are ignored.
	 * <p>
	 * Note: The height thresholding step is ignored if the height threshold is zero.
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param border
	 *            The internal border
	 * @return Indices to the local maxima (index = maxx * y + x)
	 */
	public int[] blockFind3x3Internal(float[] data, int maxx, int maxy, int border)
	{
		if (border < 1)
			return blockFind3x3(data, maxx, maxy);

		// The number of blocks in x and y
		int xblocks = getBlocks(maxx - border, 2);
		int yblocks = getBlocks(maxy - border, 2);

		int[] maxima = new int[xblocks * yblocks];

		float heightThreshold = (float) getHeightThreshold();
		float floatBackground = (float) background;
		boolean validations = (heightThreshold > 0 || minimumWidth > 0 || isNeighbourCheck());

		boolean[] maximaFlag = null;
		if (isNeighbourCheck())
			maximaFlag = getFlagBuffer(data.length);

		// Compare 2x2 block
		// ....
		// .AB.
		// .CD.
		// ....

		// Create the scan arrays to search the remaining 5 locations around the 2x2 block maxima
		// ....  aaa.  .bbb  ....  .... 
		// .AB.  aA..  ..Bb  c...  ...d 
		// .CD.  a...  ...b  cC..  ..Dd 
		// ....  ....  ....  ccc.  .ddd 
		int[] a = new int[] { -maxx - 1, -maxx, -maxx + 1, -1, +maxx - 1 };
		int[] b = new int[] { -maxx - 1, -maxx, -maxx + 1, +1, +maxx + 1 };
		int[] c = new int[] { -maxx - 1, -1, +maxx - 1, +maxx, +maxx + 1 };
		int[] d = new int[] { -maxx + 1, +1, +maxx - 1, +maxx, +maxx + 1 };

		int nMaxima = 0;
		for (int y = border; y < maxy - border - 1; y += 2)
		{
			int xindex = y * maxx + border;
			FIND_MAXIMUM: for (int x = border; x < maxx - border - 1; x += 2, xindex += 2)
			{
				int[] scan = a;
				int maxIndex = xindex;
				if (data[maxIndex] < data[xindex + 1])
				{
					scan = b;
					maxIndex = xindex + 1;
				}
				if (data[maxIndex] < data[xindex + maxx])
				{
					scan = c;
					maxIndex = xindex + maxx;
				}
				if (data[maxIndex] < data[xindex + maxx + 1])
				{
					scan = d;
					maxIndex = xindex + maxx + 1;
				}

				// Check the remaining region
				for (int offset : scan)
				{
					if (data[maxIndex] < data[maxIndex + offset])
						continue FIND_MAXIMUM;
				}

				if (validations)
				{
					if (data[maxIndex] < heightThreshold)
						continue;

					if (isNeighbourCheck())
					{
						for (int offset : scan)
						{
							if (maximaFlag[maxIndex + offset])
								continue FIND_MAXIMUM;
						}
					}

					// Check the maximum width
					if (minimumWidth > 0)
					{
						int x0 = maxIndex % maxx;
						int y0 = maxIndex / maxx;

						// Get the width at half maximum.
						float v_half = floatHalfMaximum(floatBackground, data[maxIndex]);
						int index2;
						// Scan right
						int x1 = x0 + 1;
						index2 = maxIndex + 1;
						while (x1 < maxx && data[index2] > v_half)
						{
							x1++;
							index2++;
						}
						// Scan left
						int x2 = x0 - 1;
						index2 = maxIndex - 1;
						while (x2 >= 0 && data[index2] > v_half)
						{
							x2--;
							index2--;
						}
						if (x1 - x2 < minimumWidth)
							continue;
						// Scan up
						int y1 = y0 + 1;
						index2 = maxIndex + maxx;
						while (y1 < maxy && data[index2] > v_half)
						{
							y1++;
							index2 += maxx;
						}
						// Scan down
						int y2 = y0 - 1;
						index2 = maxIndex - maxx;
						while (y2 >= 0 && data[index2] > v_half)
						{
							y2--;
							index2 -= maxx;
						}
						if (y1 - y2 < minimumWidth)
							continue;
					}

					if (isNeighbourCheck())
						maximaFlag[maxIndex] = true;
				}

				maxima[nMaxima++] = maxIndex;
			} // end FIND_MAXIMA
		}

		return truncate(maxima, nMaxima);
	}

	/**
	 * Search the data for the index of the maximum in each block of size n*n.
	 * <p>
	 * E.g. Max [ (i,i+n) x (i,j+n) ] for (i=0; i&lt;maxx; i+=n) x (j=0, j&lt;maxy; j+=n)
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @return The maxima indices
	 */
	public int[] findBlockMaxima(float[] data, int maxx, int maxy, int n)
	{
		if (n == 1)
			return findBlockMaxima2x2(data, maxx, maxy);

		return findBlockMaximaNxN(data, maxx, maxy, n);
	}

	/**
	 * Expand the image to the new dimensions with a 1-pixel border
	 */
	private float[] expand(float[] data, int maxx, int maxy, int newx, int newy)
	{
		int size = newx * newy;

		if (!dataBuffer || newDataFloat == null || newDataFloat.length < size)
			newDataFloat = new float[size];

		// Zero first row
		for (int x = 0; x < newx; x++)
		{
			newDataFloat[x] = Float.NEGATIVE_INFINITY;
		}
		// Zero last rows
		for (int y = maxy + 1; y < newy; y++)
		{
			int newIndex = y * newx;
			for (int x = 0; x < newx; x++)
			{
				newDataFloat[newIndex++] = Float.NEGATIVE_INFINITY;
			}
		}

		int index = 0;
		for (int y = 0; y < maxy; y++)
		{
			int newIndex = (y + 1) * newx;

			// Zero first column
			newDataFloat[newIndex++] = Float.NEGATIVE_INFINITY;

			// Copy data
			for (int x = 0; x < maxx; x++)
			{
				newDataFloat[newIndex++] = data[index++];
			}

			// Zero remaining columns
			for (int x = maxx + 1; x < newx; x++)
				newDataFloat[newIndex++] = Float.NEGATIVE_INFINITY;
		}

		return newDataFloat;
	}

	/**
	 * Expand the image to the new dimensions with a 1-pixel border
	 */
	private int[] expand(int[] data, int maxx, int maxy, int newx, int newy)
	{
		int size = newx * newy;

		if (!dataBuffer || newDataInt == null || newDataInt.length < size)
			newDataInt = new int[size];

		// Zero first row
		for (int x = 0; x < newx; x++)
		{
			newDataInt[x] = Integer.MIN_VALUE;
		}
		// Zero last rows
		for (int y = maxy + 1; y < newy; y++)
		{
			int newIndex = y * newx;
			for (int x = 0; x < newx; x++)
			{
				newDataInt[newIndex++] = Integer.MIN_VALUE;
			}
		}

		int index = 0;
		for (int y = 0; y < maxy; y++)
		{
			int newIndex = (y + 1) * newx;

			// Zero first column
			newDataInt[newIndex++] = Integer.MIN_VALUE;

			// Copy data
			for (int x = 0; x < maxx; x++)
			{
				newDataInt[newIndex++] = data[index++];
			}

			// Zero remaining columns
			for (int x = maxx + 1; x < newx; x++)
				newDataInt[newIndex++] = Integer.MIN_VALUE;
		}

		return newDataInt;
	}

	/**
	 * Get a buffer for storing result indices.
	 *
	 * @param size
	 *            the size
	 * @return the results buffer
	 */
	private FixedIntList getResultsBuffer(int size)
	{
		if (!dataBuffer)
			return new FixedIntList(size);

		if (resultsBuffer == null || resultsBuffer.capacity() < size)
		{
			resultsBuffer = new FixedIntList(size);
		}
		else
		{
			resultsBuffer.clear();
		}

		return resultsBuffer;
	}

	/**
	 * Get a buffer for flagging maxima.
	 *
	 * @param size
	 *            the size
	 * @return the flag buffer
	 */
	private boolean[] getFlagBuffer(int size)
	{
		if (!dataBuffer)
			return new boolean[size];

		if (maximaFlagBuffer == null || maximaFlagBuffer.length < size)
		{
			maximaFlagBuffer = new boolean[size];
		}
		else
		{
			// Reset flags
			for (int x = size; x-- > 0;)
				maximaFlagBuffer[x] = false;
		}

		return maximaFlagBuffer;
	}

	private int getBlocks(int max, int n)
	{
		int blocks = (int) Math.ceil((1.0 * max) / n);
		return blocks;
	}

	/**
	 * Get the height threshold for peaks using the current minimum height and fraction above background.
	 * 
	 * @return the height threshold
	 */
	public float getHeightThreshold()
	{
		float heightThreshold = minimumHeight + background;
		if (fractionAboveBackground != 1)
		{
			float heightThreshold2 = background / (1 - fractionAboveBackground);
			heightThreshold = FastMath.max(heightThreshold, heightThreshold2);
		}
		return heightThreshold;
	}

	/**
	 * Truncate the array to the specified size
	 * 
	 * @param array
	 * @param size
	 * @return The truncated array
	 */
	private int[] truncate(int[] array, int size)
	{
		if (array.length == size)
			return array;
		if (size == 0)
			return new int[0];
		return Arrays.copyOf(array, size);
	}

	/**
	 * @param background
	 * @param maximum
	 * @return The half-maximum value (accounting for background)
	 */
	private float floatHalfMaximum(float background, float maximum)
	{
		return (maximum - background) * 0.5f;
	}

	/**
	 * @param background
	 * @param maximum
	 * @return The half-maximum value (accounting for background)
	 */
	private int intHalfMaximum(int background, int maximum)
	{
		return (int) (0.5 + (maximum - background) / 2);
	}

	/**
	 * @param fractionAboveBackground
	 *            the fraction above background
	 */
	public void setFractionAboveBackground(float fractionAboveBackground)
	{
		this.fractionAboveBackground = fractionAboveBackground;
	}

	/**
	 * @return the fraction above background
	 */
	public float getFractionAboveBackground()
	{
		return fractionAboveBackground;
	}

	/**
	 * @param minimumHeight
	 *            the minimumHeight to set
	 */
	public void setMinimumHeight(float minimumHeight)
	{
		this.minimumHeight = minimumHeight;
	}

	/**
	 * @return the minimumHeight
	 */
	public float getMinimumHeight()
	{
		return minimumHeight;
	}

	/**
	 * @param minimumWidth
	 *            the minimumWidth to set
	 */
	public void setMinimumWidth(float minimumWidth)
	{
		this.minimumWidth = minimumWidth;
	}

	/**
	 * @return the minimumWidth
	 */
	public float getMinimumWidth()
	{
		return minimumWidth;
	}

	/**
	 * @param background
	 *            the background to set
	 */
	public void setBackground(float background)
	{
		this.background = background;
	}

	/**
	 * @return the background
	 */
	public float getBackground()
	{
		return background;
	}

	/**
	 * Neighbour checking performs an additional comparison between the local maxima within the block (size n)
	 * and the neighbours within the 2n+1 perimeter. If any neighbour is already a maxima then the local maxima within
	 * the block is elminiated. This step is only relevant when neighbour data points have equal values since the
	 * search for maxima uses the &lt; operator.
	 * <p>
	 * Applies to the blockFind algorithms.
	 * 
	 * @param neighbourCheck
	 *            Enable neighbour checking
	 */
	public void setNeighbourCheck(boolean neighbourCheck)
	{
		this.neighbourCheck = neighbourCheck;
	}

	/**
	 * @return True if neighbour checking is enabled
	 */
	public boolean isNeighbourCheck()
	{
		return neighbourCheck;
	}

	/**
	 * Allow the class to keep a data buffer for processing images with the blockFind3x3 algorithm
	 * 
	 * @param dataBuffer
	 *            Enable the data buffer
	 */
	public void setDataBuffer(boolean dataBuffer)
	{
		this.dataBuffer = dataBuffer;
		if (!dataBuffer)
			newDataFloat = null;
	}

	/**
	 * @return True if the data buffer is enabled
	 */
	public boolean isBufferData()
	{
		return dataBuffer;
	}

	// ----------------------------------------------------
	// NOTE:
	// The following code is copied directly from above. 
	// All 'float' have been replaced with 'int'
	// ----------------------------------------------------

	/**
	 * Compute the local-maxima within a 2n+1 block
	 * <p>
	 * Any maxima below the configured fraction above background are ignored. Fraction = (Max-background)/Max within the
	 * 2n+1 neighbourhood. Maxima below the minimum height (above background) or the minimum peak-width at half maximum
	 * in any dimension are ignored.
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @return Indices to the local maxima (index = maxx * y + x)
	 */
	public int[] maxFind(int[] data, int maxx, int maxy, int n)
	{
		FixedIntList results = getResultsBuffer(data.length / 4);
		boolean[] maximaFlag = getFlagBuffer(data.length);

		// Boundary control
		int xwidth = FastMath.min(n, maxx - 1);
		int ywidth = FastMath.min(n, maxy - 1);
		int xlimit = maxx - xwidth;
		int ylimit = maxy - ywidth;

		int[] offset = new int[(2 * xwidth + 1) * (2 * ywidth + 1) - 1];
		int[] xoffset = new int[offset.length];
		int[] yoffset = new int[offset.length];
		int d = 0;
		for (int y = -ywidth; y <= ywidth; y++)
			for (int x = -xwidth; x <= xwidth; x++)
				if (x != 0 || y != 0)
				{
					offset[d] = maxx * y + x;
					xoffset[d] = x;
					yoffset[d] = y;
					d++;
				}

		// Compare all points
		int heightThreshold = (int) getHeightThreshold();
		int floatBackground = (int) background;

		int index = 0;
		for (int y = 0; y < maxy; y++)
		{
			FIND_MAXIMUM: for (int x = 0; x < maxx; x++, index++)
			{
				int v = data[index];
				if (v < heightThreshold)
					continue;

				// Flag to indicate this pixels has a complete (2n+1) neighbourhood 
				boolean isInnerXY = (y >= ywidth && y < ylimit) && (x >= xwidth && x < xlimit);

				// Sweep neighbourhood
				if (isInnerXY)
				{
					for (int i = 0; i < offset.length; i++)
					{
						if (maximaFlag[index + offset[i]])
							continue FIND_MAXIMUM;
						if (data[index + offset[i]] > v)
							continue FIND_MAXIMUM;
					}
				}
				else
				{
					for (d = offset.length; d-- > 0;)
					{
						// Get the coords and check if it is within the data
						int yy = y + yoffset[d];
						int xx = x + xoffset[d];
						boolean isWithin = (yy >= 0 && yy < maxy) && (xx >= 0 && xx < maxx);

						if (isWithin)
						{
							if (maximaFlag[index + offset[d]])
								continue FIND_MAXIMUM;
							if (data[index + offset[d]] > v)
								continue FIND_MAXIMUM;
						}
					}
				}

				// Check the maximum width
				if (minimumWidth > 0)
				{
					// Get the width at half maximum.
					int v_half = intHalfMaximum(floatBackground, v);
					int index2;
					// Scan right
					int x1 = x + 1;
					index2 = index + 1;
					while (x1 < maxx && data[index2] > v_half)
					{
						x1++;
						index2++;
					}
					// Scan left
					int x2 = x - 1;
					index2 = index - 1;
					while (x2 >= 0 && data[index2] > v_half)
					{
						x2--;
						index2--;
					}
					if (x1 - x2 < minimumWidth)
						continue;
					// Scan up
					int y1 = y + 1;
					index2 = index + maxx;
					while (y1 < maxy && data[index2] > v_half)
					{
						y1++;
						index2 += maxx;
					}
					// Scan down
					int y2 = y - 1;
					index2 = index - 1;
					while (y2 >= 0 && data[index2] > v_half)
					{
						y2--;
						index2 -= maxx;
					}
					if (y1 - y2 < minimumWidth)
						continue;
				}

				results.add(index);
				maximaFlag[index] = true;
			} // end FIND_MAXIMA
		}

		return results.toArray();
	}

	/**
	 * Compute the local-maxima within a 2n+1 block
	 * An inner boundary of N is ignored as potential maxima.
	 * <p>
	 * Any maxima below the configured fraction above background are ignored. Fraction = (Max-background)/Max within the
	 * 2n+1 neighbourhood. Maxima below the minimum height (above background) or the minimum peak-width at half maximum
	 * in any dimension are ignored.
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @return Indices to the local maxima (index = maxx * y + x)
	 */
	public int[] maxFindInternal(int[] data, int maxx, int maxy, int n)
	{
		FixedIntList results = getResultsBuffer(data.length / 4);
		boolean[] maximaFlag = getFlagBuffer(data.length);

		// Boundary control
		int xwidth = FastMath.min(n, maxx - 1);
		int ywidth = FastMath.min(n, maxy - 1);

		int[] offset = new int[(2 * xwidth + 1) * (2 * ywidth + 1) - 1];
		int d = 0;
		for (int y = -ywidth; y <= ywidth; y++)
			for (int x = -xwidth; x <= xwidth; x++)
				if (x != 0 || y != 0)
				{
					offset[d] = maxx * y + x;
					d++;
				}

		// Compare all points
		int heightThreshold = (int) getHeightThreshold();
		int floatBackground = (int) background;

		for (int y = n; y < maxy - n; y++)
		{
			int index = y * maxx + n;
			FIND_MAXIMUM: for (int x = n; x < maxx - n; x++, index++)
			{
				int v = data[index];
				if (v < heightThreshold)
					continue;

				// Sweep neighbourhood - 
				// No check for boundaries as this should be an internal sweep.
				for (int i = 0; i < offset.length; i++)
				{
					if (maximaFlag[index + offset[i]])
						continue FIND_MAXIMUM;
					if (data[index + offset[i]] > v)
						continue FIND_MAXIMUM;
				}

				// Check the maximum width
				if (minimumWidth > 0)
				{
					int v_half = intHalfMaximum(floatBackground, v);
					int index2;
					// Scan right
					int x1 = x + 1;
					index2 = index + 1;
					while (x1 < maxx && data[index2] > v_half)
					{
						x1++;
						index2++;
					}
					// Scan left
					int x2 = x - 1;
					index2 = index - 1;
					while (x2 >= 0 && data[index2] > v_half)
					{
						x2--;
						index2--;
					}
					if (x1 - x2 < minimumWidth)
						continue;
					// Scan up
					int y1 = y + 1;
					index2 = index + maxx;
					while (y1 < maxy && data[index2] > v_half)
					{
						y1++;
						index2 += maxx;
					}
					// Scan down
					int y2 = y - 1;
					index2 = index - 1;
					while (y2 >= 0 && data[index2] > v_half)
					{
						y2--;
						index2 -= maxx;
					}
					if (y1 - y2 < minimumWidth)
						continue;
				}

				results.add(index);
				maximaFlag[index] = true;
			} // end FIND_MAXIMA
		}

		return results.toArray();
	}

	/**
	 * Compute the local-maxima within a 2n+1 block
	 * An inner boundary is ignored as potential maxima.
	 * <p>
	 * Any maxima below the configured fraction above background are ignored. Fraction = (Max-background)/Max within the
	 * 2n+1 neighbourhood. Maxima below the minimum height (above background) or the minimum peak-width at half maximum
	 * in any dimension are ignored.
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @param border
	 *            The internal border
	 * @return Indices to the local maxima (index = maxx * y + x)
	 */
	public int[] maxFindInternal(int[] data, int maxx, int maxy, int n, int border)
	{
		if (n == border)
			// Faster algorithm as there is no requirement for bounds checking.
			return maxFindInternal(data, maxx, maxy, n);

		FixedIntList results = getResultsBuffer(data.length / 4);
		boolean[] maximaFlag = getFlagBuffer(data.length);

		// Boundary control
		int xwidth = FastMath.min(n, maxx - 1);
		int ywidth = FastMath.min(n, maxy - 1);
		int xlimit = maxx - xwidth;
		int ylimit = maxy - ywidth;

		int[] offset = new int[(2 * xwidth + 1) * (2 * ywidth + 1) - 1];
		int[] xoffset = new int[offset.length];
		int[] yoffset = new int[offset.length];
		int d = 0;
		for (int y = -ywidth; y <= ywidth; y++)
			for (int x = -xwidth; x <= xwidth; x++)
				if (x != 0 || y != 0)
				{
					offset[d] = maxx * y + x;
					xoffset[d] = x;
					yoffset[d] = y;
					d++;
				}

		// All blocks fit within the border
		boolean inner = (n < border);

		// Compare all points
		int heightThreshold = (int) getHeightThreshold();
		int floatBackground = (int) background;

		for (int y = border; y < maxy - border; y++)
		{
			int index = y * maxx + border;
			FIND_MAXIMUM: for (int x = border; x < maxx - border; x++, index++)
			{
				int v = data[index];
				if (v < heightThreshold)
					continue;

				if (inner)
				{
					for (int i = 0; i < offset.length; i++)
					{
						if (maximaFlag[index + offset[i]])
							continue FIND_MAXIMUM;
						if (data[index + offset[i]] > v)
							continue FIND_MAXIMUM;
					}
				}
				else
				{
					// Flag to indicate this pixels has a complete (2n+1) neighbourhood 
					boolean isInnerXY = (y >= ywidth && y < ylimit) && (x >= xwidth && x < xlimit);

					// Sweep neighbourhood
					for (d = offset.length; d-- > 0;)
					{
						boolean isWithin = isInnerXY;
						if (!isWithin)
						{
							// Get the coords and check if it is within the data
							int yy = y + yoffset[d];
							int xx = x + xoffset[d];
							isWithin = (yy >= 0 && yy < maxy) && (xx >= 0 && xx < maxx);
						}

						if (isWithin)
						{
							if (maximaFlag[index + offset[d]])
								continue FIND_MAXIMUM;
							if (data[index + offset[d]] > v)
								continue FIND_MAXIMUM;
						}
					}
				}

				// Check the maximum width
				if (minimumWidth > 0)
				{
					int v_half = intHalfMaximum(floatBackground, v);
					int index2;
					// Scan right
					int x1 = x + 1;
					index2 = index + 1;
					while (x1 < maxx && data[index2] > v_half)
					{
						x1++;
						index2++;
					}
					// Scan left
					int x2 = x - 1;
					index2 = index - 1;
					while (x2 >= 0 && data[index2] > v_half)
					{
						x2--;
						index2--;
					}
					if (x1 - x2 < minimumWidth)
						continue;
					// Scan up
					int y1 = y + 1;
					index2 = index + maxx;
					while (y1 < maxy && data[index2] > v_half)
					{
						y1++;
						index2 += maxx;
					}
					// Scan down
					int y2 = y - 1;
					index2 = index - 1;
					while (y2 >= 0 && data[index2] > v_half)
					{
						y2--;
						index2 -= maxx;
					}
					if (y1 - y2 < minimumWidth)
						continue;
				}

				results.add(index);
				maximaFlag[index] = true;
			} // end FIND_MAXIMA
		}

		return results.toArray();
	}

	/**
	 * Compute the local-maxima within a 2n+1 block.
	 * <p>
	 * Any maxima below the configured fraction above background are ignored. Fraction = (Max-background)/Max within the
	 * 2n+1 neighbourhood. Maxima below the minimum height (above background) or the minimum peak-width at half maximum
	 * in any dimension are ignored.
	 * <p>
	 * Uses the 2D block algorithm of Neubeck and Van Gool (2006).
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @return Indices to the local maxima (index = maxx * y + x)
	 */
	public int[] blockFind(int[] data, int maxx, int maxy, int n)
	{
		if (n == 1)
			// optimised version for the special case 
			return blockFind3x3(data, maxx, maxy);

		return blockFindNxN(data, maxx, maxy, n);
	}

	/**
	 * Compute the local-maxima within a 2n+1 block.
	 * An inner boundary of N is ignored as potential maxima.
	 * <p>
	 * Any maxima below the configured fraction above background are ignored. Fraction = (Max-background)/Max within the
	 * 2n+1 neighbourhood. Maxima below the minimum height (above background) or the minimum peak-width at half maximum
	 * in any dimension are ignored.
	 * <p>
	 * Uses the 2D block algorithm of Neubeck and Van Gool (2006).
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @param border
	 *            The internal border
	 * @return Indices to the local maxima (index = maxx * y + x)
	 */
	public int[] blockFindInternal(int[] data, int maxx, int maxy, int n, int border)
	{
		if (n == 1)
			// optimised version for the special case 
			return blockFind3x3Internal(data, maxx, maxy, border);

		return blockFindNxNInternal(data, maxx, maxy, n, border);
	}

	/**
	 * Compute the local-maxima within a 2n+1 block.
	 * <p>
	 * Any maxima below the configured fraction above background are ignored. Fraction = (Max-background)/Max within the
	 * 2n+1 neighbourhood. Maxima below the minimum height (above background) or the minimum peak-width at half maximum
	 * in any dimension are ignored.
	 * <p>
	 * Uses the 2D block algorithm of Neubeck and Van Gool (2006).
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @return Indices to the local maxima (index = maxx * y + x)
	 */
	public int[] blockFindNxN(int[] data, int maxx, int maxy, int n)
	{
		int[] blockMaxima = findBlockMaximaNxN(data, maxx, maxy, n);
		int nMaxima = 0;
		int heightThreshold = (int) getHeightThreshold();
		int floatBackground = (int) background;

		boolean[] maximaFlag = null;
		if (isNeighbourCheck())
			maximaFlag = getFlagBuffer(data.length);

		FIND_MAXIMUM: for (int index : blockMaxima)
		{
			int v = data[index];

			if (v < heightThreshold)
				continue;

			int x = index % maxx;
			int y = index / maxx;

			// Compare the maxima to the surroundings. Ignore the block region already processed.
			//
			//.......(mi-n,mj-n)----------------------------------(mi+n,mj-n)
			//............|aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa|
			//............|aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa|
			//............|aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa|
			//............|aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa|
			//............|aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa|
			//............|bbbbbbbb.(i,j)-------------(i+n,j).ccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbbbbb|.....(mi,mj)......|ccccccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbb.(i,j+n)-----------(i+n,j+n).ccccccc|
			//............|dddddddddddddddddddddddddddddddddddddddddddd|
			//............|dddddddddddddddddddddddddddddddddddddddddddd|
			//............|dddddddddddddddddddddddddddddddddddddddddddd|
			//......(mi-n,mj+n)----------------------------------(mi+n,mj+n)
			//
			// This must be done without over-running boundaries
			//int steps = 0;
			int mi = x;
			int mj = y;
			int i = (n + 1) * (mi / (n + 1));
			int j = (n + 1) * (mj / (n + 1));
			int i_plus_n = FastMath.min(i + n + 1, maxx - 1);
			int j_plus_n = FastMath.min(j + n + 1, maxy - 1);
			int mi_minus_n = FastMath.max(mi - n, 0);
			int mi_plus_n = FastMath.min(mi + n, maxx - 1);
			int mj_minus_n = FastMath.max(mj - n, 0);
			int mj_plus_n = FastMath.min(mj + n, maxy - 1);

			// A
			for (int jj = mj_minus_n; jj < j; jj++)
			{
				int indexStart = jj * maxx + mi_minus_n;
				int indexEnd = jj * maxx + mi_plus_n;
				for (; indexStart <= indexEnd; indexStart++)
				{
					//steps++;
					if (v < data[indexStart])
						continue FIND_MAXIMUM;
				}
			}
			for (int jj = j; jj < j_plus_n; jj++)
			{
				// B
				{
					int indexStart = jj * maxx + mi_minus_n;
					int indexEnd = jj * maxx + i;
					for (; indexStart < indexEnd; indexStart++)
					{
						//steps++;
						if (v < data[indexStart])
							continue FIND_MAXIMUM;
					}
				}

				// C
				{
					int indexStart = jj * maxx + i_plus_n;
					int indexEnd = jj * maxx + mi_plus_n;
					for (; indexStart <= indexEnd; indexStart++)
					{
						//steps++;
						if (v < data[indexStart])
							continue FIND_MAXIMUM;
					}
				}
			}
			// D
			for (int jj = j_plus_n; jj <= mj_plus_n; jj++)
			{
				int indexStart = jj * maxx + mi_minus_n;
				int indexEnd = jj * maxx + mi_plus_n;
				for (; indexStart <= indexEnd; indexStart++)
				{
					//steps++;
					if (v < data[indexStart])
						continue FIND_MAXIMUM;
				}
			}
			//System.out.printf("%.2f @ %d,%d. Steps = %d, max = %b\n", v, x, y, steps, isMax);

			if (isNeighbourCheck())
			{
				// A
				for (int jj = mj_minus_n; jj < j; jj++)
				{
					int indexStart = jj * maxx + mi_minus_n;
					int indexEnd = jj * maxx + mi_plus_n;
					for (; indexStart <= indexEnd; indexStart++)
					{
						//steps++;
						if (maximaFlag[indexStart])
							continue FIND_MAXIMUM;
					}
				}
				for (int jj = j; jj < j_plus_n; jj++)
				{
					// B
					{
						int indexStart = jj * maxx + mi_minus_n;
						int indexEnd = jj * maxx + i;
						for (; indexStart < indexEnd; indexStart++)
						{
							//steps++;
							if (maximaFlag[indexStart])
								continue FIND_MAXIMUM;
						}
					}

					// C
					{
						int indexStart = jj * maxx + i_plus_n;
						int indexEnd = jj * maxx + mi_plus_n;
						for (; indexStart <= indexEnd; indexStart++)
						{
							//steps++;
							if (maximaFlag[indexStart])
								continue FIND_MAXIMUM;
						}
					}
				}
				// D
				for (int jj = j_plus_n; jj <= mj_plus_n; jj++)
				{
					int indexStart = jj * maxx + mi_minus_n;
					int indexEnd = jj * maxx + mi_plus_n;
					for (; indexStart <= indexEnd; indexStart++)
					{
						//steps++;
						if (maximaFlag[indexStart])
							continue FIND_MAXIMUM;
					}
				}
			}

			// Check the maximum width
			if (minimumWidth > 0)
			{
				int v_half = intHalfMaximum(floatBackground, v);
				int index2;
				// Scan right
				int x1 = x + 1;
				index2 = index + 1;
				while (x1 < maxx && data[index2] > v_half)
				{
					x1++;
					index2++;
				}
				// Scan left
				int x2 = x - 1;
				index2 = index - 1;
				while (x2 >= 0 && data[index2] > v_half)
				{
					x2--;
					index2--;
				}
				if (x1 - x2 < minimumWidth)
					continue;
				// Scan up
				int y1 = y + 1;
				index2 = index + maxx;
				while (y1 < maxy && data[index2] > v_half)
				{
					y1++;
					index2 += maxx;
				}
				// Scan down
				int y2 = y - 1;
				index2 = index - 1;
				while (y2 >= 0 && data[index2] > v_half)
				{
					y2--;
					index2 -= maxx;
				}
				if (y1 - y2 < minimumWidth)
					continue;
			}

			//System.out.printf("blockFind [%d,%d]\n", mi, mj);

			// Re-use storage space
			blockMaxima[nMaxima++] = index;
			if (isNeighbourCheck())
				maximaFlag[index] = true;
		} // end FIND_MAXIMA
		  //System.out.printf("---\n");

		return truncate(blockMaxima, nMaxima);
	}

	/**
	 * Compute the local-maxima within a 2n+1 block.
	 * An inner boundary of N is ignored as potential maxima.
	 * <p>
	 * Any maxima below the configured fraction above background are ignored. Fraction = (Max-background)/Max within the
	 * 2n+1 neighbourhood. Maxima below the minimum height (above background) or the minimum peak-width at half maximum
	 * in any dimension are ignored.
	 * <p>
	 * Uses the 2D block algorithm of Neubeck and Van Gool (2006).
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @param border
	 *            The internal border
	 * @return Indices to the local maxima (index = maxx * y + x)
	 */
	public int[] blockFindNxNInternal(int[] data, int maxx, int maxy, int n, int border)
	{
		int[] blockMaxima = findBlockMaximaNxNInternal(data, maxx, maxy, n, border);
		int nMaxima = 0;
		int heightThreshold = (int) getHeightThreshold();
		int floatBackground = (int) background;

		boolean[] maximaFlag = null;
		if (isNeighbourCheck())
			maximaFlag = getFlagBuffer(data.length);

		FIND_MAXIMUM: for (int index : blockMaxima)
		{
			int v = data[index];

			if (v < heightThreshold)
				continue;

			int x = index % maxx;
			int y = index / maxx;

			// Compare the maxima to the surroundings. Ignore the block region already processed.
			//
			//.......(mi-n,mj-n)----------------------------------(mi+n,mj-n)
			//............|aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa|
			//............|aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa|
			//............|aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa|
			//............|aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa|
			//............|aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa|
			//............|bbbbbbbb.(i,j)-------------(i+n,j).ccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbbbbb|.....(mi,mj)......|ccccccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbbbbb|..................|ccccccccccccc|
			//............|bbbbbbbb.(i,j+n)-----------(i+n,j+n).ccccccc|
			//............|dddddddddddddddddddddddddddddddddddddddddddd|
			//............|dddddddddddddddddddddddddddddddddddddddddddd|
			//............|dddddddddddddddddddddddddddddddddddddddddddd|
			//......(mi-n,mj+n)----------------------------------(mi+n,mj+n)
			//
			// No check for over-running boundaries since this is the internal version
			//int steps = 0;
			int mi = x;
			int mj = y;
			int i = (n + 1) * ((mi - border) / (n + 1)) + border; // Blocks n+1 wide
			int j = (n + 1) * ((mj - border) / (n + 1)) + border; // Blocks n+1 wide
			// The block boundaries will have been truncated on the final block. Ensure this is swept
			int i_plus_n = FastMath.min(i + n + 1, maxx - border - 1);
			int j_plus_n = FastMath.min(j + n + 1, maxy - border - 1);
			int mi_minus_n = FastMath.max(mi - n, 0);
			int mi_plus_n = FastMath.min(mi + n, maxx - 1);
			int mj_minus_n = FastMath.max(mj - n, 0);
			int mj_plus_n = FastMath.min(mj + n, maxy - 1);

			//System.out.printf("Block [%d,%d] => [%d,%d]\n", x, y, i, j);

			// A
			for (int jj = mj_minus_n; jj < j; jj++)
			{
				int indexStart = jj * maxx + mi_minus_n;
				int indexEnd = jj * maxx + mi_plus_n;
				for (; indexStart <= indexEnd; indexStart++)
				{
					//steps++;
					if (v < data[indexStart])
						continue FIND_MAXIMUM;
				}
			}
			for (int jj = j; jj < j_plus_n; jj++)
			{
				// B
				{
					int indexStart = jj * maxx + mi_minus_n;
					int indexEnd = jj * maxx + i;
					for (; indexStart < indexEnd; indexStart++)
					{
						//steps++;
						if (v < data[indexStart])
							continue FIND_MAXIMUM;
					}
				}

				// C
				{
					int indexStart = jj * maxx + i_plus_n;
					int indexEnd = jj * maxx + mi_plus_n;
					for (; indexStart <= indexEnd; indexStart++)
					{
						//steps++;
						if (v < data[indexStart])
							continue FIND_MAXIMUM;
					}
				}
			}
			// D
			for (int jj = j_plus_n; jj <= mj_plus_n; jj++)
			{
				int indexStart = jj * maxx + mi_minus_n;
				int indexEnd = jj * maxx + mi_plus_n;
				for (; indexStart <= indexEnd; indexStart++)
				{
					//steps++;
					if (v < data[indexStart])
						continue FIND_MAXIMUM;
				}
			}
			//System.out.printf("%.2f @ %d,%d. Steps = %d, max = %b\n", v, x, y, steps, isMax);

			if (isNeighbourCheck())
			{
				// A
				for (int jj = mj_minus_n; jj < j; jj++)
				{
					int indexStart = jj * maxx + mi_minus_n;
					int indexEnd = jj * maxx + mi_plus_n;
					for (; indexStart <= indexEnd; indexStart++)
					{
						//steps++;
						if (maximaFlag[indexStart])
							continue FIND_MAXIMUM;
					}
				}
				for (int jj = j; jj < j_plus_n; jj++)
				{
					// B
					{
						int indexStart = jj * maxx + mi_minus_n;
						int indexEnd = jj * maxx + i;
						for (; indexStart < indexEnd; indexStart++)
						{
							//steps++;
							if (maximaFlag[indexStart])
								continue FIND_MAXIMUM;
						}
					}

					// C
					{
						int indexStart = jj * maxx + i_plus_n;
						int indexEnd = jj * maxx + mi_plus_n;
						for (; indexStart <= indexEnd; indexStart++)
						{
							//steps++;
							if (maximaFlag[indexStart])
								continue FIND_MAXIMUM;
						}
					}
				}
				// D
				for (int jj = j_plus_n; jj <= mj_plus_n; jj++)
				{
					int indexStart = jj * maxx + mi_minus_n;
					int indexEnd = jj * maxx + mi_plus_n;
					for (; indexStart <= indexEnd; indexStart++)
					{
						//steps++;
						if (maximaFlag[indexStart])
							continue FIND_MAXIMUM;
					}
				}
			}

			// Check the maximum width
			if (minimumWidth > 0)
			{
				int v_half = intHalfMaximum(floatBackground, v);
				int index2;
				// Scan right
				int x1 = x + 1;
				index2 = index + 1;
				while (x1 < maxx && data[index2] > v_half)
				{
					x1++;
					index2++;
				}
				// Scan left
				int x2 = x - 1;
				index2 = index - 1;
				while (x2 >= 0 && data[index2] > v_half)
				{
					x2--;
					index2--;
				}
				if (x1 - x2 < minimumWidth)
					continue;
				// Scan up
				int y1 = y + 1;
				index2 = index + maxx;
				while (y1 < maxy && data[index2] > v_half)
				{
					y1++;
					index2 += maxx;
				}
				// Scan down
				int y2 = y - 1;
				index2 = index - 1;
				while (y2 >= 0 && data[index2] > v_half)
				{
					y2--;
					index2 -= maxx;
				}
				if (y1 - y2 < minimumWidth)
					continue;
			}

			//System.out.printf("blockFind [%d,%d]\n", mi, mj);

			// Re-use storage space
			blockMaxima[nMaxima++] = index;
			if (isNeighbourCheck())
				maximaFlag[index] = true;
		} // end FIND_MAXIMA
		  //System.out.printf("---\n");

		return truncate(blockMaxima, nMaxima);
	}

	/**
	 * Search the data for the index of the maximum in each block of size n*n.
	 * <p>
	 * E.g. Max [ (i,i+n) x (i,j+n) ] for (i=0; i&lt;maxx; i+=n) x (j=0, j&lt;maxy; j+=n)
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @return The maxima indices
	 */
	public int[] findBlockMaximaNxN(int[] data, int maxx, int maxy, int n)
	{
		// Include the actual pixel in the block
		// This makes the block search inclusive: i=0; i<=n; i++
		n++;

		// The number of blocks in x and y
		int xblocks = getBlocks(maxx, n);
		int yblocks = getBlocks(maxy, n);

		// The final index in each dimension (where an incomplete block is found).
		// This equals maxx/maxy if the number of blocks fits exactly.
		int xfinal = n * (maxx / n);
		int yfinal = n * (maxy / n);

		int[] maxima = new int[xblocks * yblocks];

		int block = 0;
		for (int y = 0; y < maxy; y += n)
		{
			for (int x = 0; x < maxx; x += n)
			{
				// Find the sweep size in each direction
				int xsize = (x != xfinal) ? n : maxx - xfinal;
				int ysize = (y != yfinal) ? n : maxy - yfinal;

				int index = y * maxx + x;
				int maxIndex = index;
				int max = data[maxIndex];

				while (ysize-- > 0)
				{
					for (int x2 = xsize; x2-- > 0;)
					{
						if (max < data[index])
						{
							max = data[index];
							maxIndex = index;
						}
						index++;
					}
					index += maxx - xsize;
				}

				maxima[block++] = maxIndex;
			}
		}

		return maxima;
	}

	/**
	 * Search the data for the index of the maximum in each block of size n*n.
	 * An inner boundary of N is ignored as potential maxima.
	 * <p>
	 * E.g. Max [ (i,i+n) x (i,j+n) ] for (i=n; i&lt;maxx-n; i+=n) x (j=n, j&lt;maxy-n; j+=n)
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @param border
	 *            The internal border
	 * @return The maxima indices
	 */
	public int[] findBlockMaximaNxNInternal(int[] data, int maxx, int maxy, int n, int border)
	{
		// Include the actual pixel in the block
		// This makes the block search inclusive: i=0; i<=n; i++
		n++;

		// The number of blocks in x and y. Subtract the boundary blocks.
		int xblocks = getBlocks(maxx - border, n);
		int yblocks = getBlocks(maxy - border, n);

		if (xblocks < 1 || yblocks < 1)
			return new int[0];

		// The final index in each dimension (where an incomplete block is found).
		// This equals maxx/maxy if the number of blocks fits exactly.
		int xfinal = n * ((maxx - 2 * border) / n) + border;
		int yfinal = n * ((maxy - 2 * border) / n) + border;

		int[] maxima = new int[xblocks * yblocks];

		int block = 0;
		for (int y = border; y < maxy - border; y += n)
		{
			for (int x = border; x < maxx - border; x += n)
			{
				//System.out.printf("Block [%d,%d]\n", x, y);

				// Find the sweep size in each direction
				int xsize = (x != xfinal) ? n : maxx - xfinal - border;
				int ysize = (y != yfinal) ? n : maxy - yfinal - border;

				int index = y * maxx + x;
				int maxIndex = index;
				int max = data[maxIndex];

				while (ysize-- > 0)
				{
					for (int x2 = xsize; x2-- > 0;)
					{
						if (max < data[index])
						{
							max = data[index];
							maxIndex = index;
						}
						index++;
					}
					index += maxx - xsize;
				}

				maxima[block++] = maxIndex;
			}
		}

		return truncate(maxima, block);
	}

	/**
	 * Search the data for the index of the maximum in each block of size n*n.
	 * <p>
	 * E.g. Max [ (i,i+n) x (i,j+n) ] for (i=0; i&lt;maxx; i+=n) x (j=0, j&lt;maxy; j+=n)
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @return The maxima indices
	 */
	public int[] findBlockMaximaOptimised(int[] data, int maxx, int maxy, int n)
	{
		// =-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=
		// The following code DOES NOT increase speed. 
		// I believe the Java optimiser cannot unroll all the nested loops.
		// =-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=

		// Include the actual pixel in the block
		// This makes the block search inclusive: i=0; i<=n; i++
		n++;

		// The number of blocks in x and y
		int xblocks = getBlocks(maxx, n);
		int yblocks = getBlocks(maxy, n);

		// The final index in each dimension (where an incomplete block is found).
		// This equals maxx/maxy if the number of blocks fits exactly.
		int xfinal = n * (maxx / n);
		int yfinal = n * (maxy / n);

		int[] maxima = new int[xblocks * yblocks];

		// Sweep 4 regions:
		// x.............xfinal  
		// |.............|..maxx
		// |.............|..|  
		// aaaaaaaaaaaaaaabbb-- y
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb--yfinal
		// cccccccccccccccddd
		// cccccccccccccccddd--maxy

		int block = 0;
		for (int y = 0; y < yfinal; y += n)
		{
			// A
			for (int x = 0; x < xfinal; x += n)
			{
				// Find the sweep size in each direction
				int xsize = n;
				int ysize = n;

				int index = y * maxx + x;
				int maxIndex = index;
				int max = data[maxIndex];

				while (ysize-- > 0)
				{
					for (int x2 = xsize; x2-- > 0;)
					{
						if (max < data[index])
						{
							max = data[index];
							maxIndex = index;
						}
						index++;
					}
					index += maxx - xsize;
				}

				maxima[block++] = maxIndex;
			}

			// B
			if (xfinal != maxx)
			{
				// Find the sweep size in each direction
				int xsize = maxx - xfinal;
				int ysize = n;

				int index = y * maxx + xfinal;
				int maxIndex = index;
				int max = data[maxIndex];

				while (ysize-- > 0)
				{
					for (int x2 = xsize; x2-- > 0;)
					{
						if (max < data[index])
						{
							max = data[index];
							maxIndex = index;
						}
						index++;
					}
					index += maxx - xsize;
				}

				maxima[block++] = maxIndex;
			}
		}
		if (yfinal != maxy)
		{
			// C
			for (int x = 0; x < xfinal; x += n)
			{
				// Find the sweep size in each direction
				int xsize = n;
				int ysize = maxy - yfinal;

				int index = yfinal * maxx + x;
				int maxIndex = index;
				int max = data[maxIndex];

				while (ysize-- > 0)
				{
					for (int x2 = xsize; x2-- > 0;)
					{
						if (max < data[index])
						{
							max = data[index];
							maxIndex = index;
						}
						index++;
					}
					index += maxx - xsize;
				}

				maxima[block++] = maxIndex;
			}

			// D
			if (xfinal != maxx)
			{
				// Find the sweep size in each direction
				int xsize = maxx - xfinal;
				int ysize = maxy - yfinal;

				int index = yfinal * maxx + xfinal;
				int maxIndex = index;
				int max = data[maxIndex];

				while (ysize-- > 0)
				{
					for (int x2 = xsize; x2-- > 0;)
					{
						if (max < data[index])
						{
							max = data[index];
							maxIndex = index;
						}
						index++;
					}
					index += maxx - xsize;
				}

				maxima[block++] = maxIndex;
			}
		}

		return maxima;
	}

	/**
	 * Search the data for the index of the maximum in each block of size 2*2.
	 * <p>
	 * E.g. Max [ (i,i+n) x (i,j+n) ] for (i=0; i&lt;maxx; i+=n) x (j=0, j&lt;maxy; j+=n)
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @return The maxima indices
	 */
	public int[] findBlockMaxima2x2(int[] data, int maxx, int maxy)
	{
		// Optimised for 2x2 block
		//final int n = 2;

		// The number of blocks in x and y
		int xblocks = getBlocks(maxx, 2);
		int yblocks = getBlocks(maxy, 2);

		// The final index in each dimension (where an incomplete block is found).
		// This equals maxx/maxy if the number of blocks fits exactly.
		int xfinal = 2 * (maxx / 2);
		int yfinal = 2 * (maxy / 2);

		int[] maxima = new int[xblocks * yblocks];

		// Sweep 4 regions:
		// x.............xfinal  
		// |.............|..maxx
		// |.............|..|  
		// aaaaaaaaaaaaaaabbb-- y
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb
		// aaaaaaaaaaaaaaabbb--yfinal
		// cccccccccccccccddd
		// cccccccccccccccddd--maxy

		int block = 0;
		for (int y = 0; y < yfinal; y += 2)
		{
			// A
			int xindex = y * maxx;
			for (int x = 0; x < xfinal; x += 2)
			{
				// Compare 2x2 block
				int max1 = (data[xindex] > data[xindex + 1]) ? xindex : xindex + 1;
				int max2 = (data[xindex + maxx] > data[xindex + maxx + 1]) ? xindex + maxx : xindex + maxx + 1;

				maxima[block++] = (data[max1] > data[max2]) ? max1 : max2;

				xindex += 2;
			}
			// B
			if (xfinal != maxx)
			{
				// Compare 1x2 block
				int index = y * maxx + xfinal;
				maxima[block++] = (data[index] > data[index + maxx]) ? index : index + maxx;
			}
		}
		if (yfinal != maxy)
		{
			// C
			int xindex = yfinal * maxx;
			for (int x = 0; x < xfinal; x += 2)
			{
				// Compare 2x1 block
				maxima[block++] = (data[xindex] > data[xindex + 1]) ? xindex : xindex + 1;
				xindex += 2;
			}
			// D
			if (xfinal != maxx)
			{
				// Compare 1x1 block
				maxima[block++] = yfinal * maxx + xfinal;
			}
		}

		return maxima;
	}

	/**
	 * Compute the local-maxima within a 3x3 block
	 * <p>
	 * Any maxima below the configured fraction above background are ignored. Fraction = (Max-background)/Max within the
	 * 3x3 neighbourhood. Maxima below the minimum height (above background) or the minimum peak-width at half maximum
	 * in any dimension are ignored.
	 * <p>
	 * Note: The height thresholding step is ignored if the height threshold is zero.
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @return Indices to the local maxima (index = maxx * y + x)
	 */
	public int[] blockFind3x3(int[] data, int maxx, int maxy)
	{
		// The number of blocks in x and y
		int xblocks = getBlocks(maxx, 2);
		int yblocks = getBlocks(maxy, 2);

		// The final index in each dimension (where an incomplete block is found).
		// This equals maxx/maxy if the number of blocks fits exactly.
		int xfinal = 2 * (maxx / 2);
		int yfinal = 2 * (maxy / 2);

		// Expand the canvas
		int newx = maxx + (maxx - xfinal) + 2;
		int newy = maxy + (maxy - yfinal) + 2;
		data = expand(data, maxx, maxy, newx, newy);

		int[] maxima = new int[xblocks * yblocks];

		int heightThreshold = (int) getHeightThreshold();
		int floatBackground = (int) background;
		boolean validations = (heightThreshold > 0 || minimumWidth > 0 || isNeighbourCheck());

		boolean[] maximaFlag = null;
		if (isNeighbourCheck())
			maximaFlag = getFlagBuffer(data.length);

		// Compare 2x2 block
		// ....
		// .AB.
		// .CD.
		// ....

		// Create the scan arrays to search the remaining 5 locations around the 2x2 block maxima
		// ....  aaa.  .bbb  ....  .... 
		// .AB.  aA..  ..Bb  c...  ...d 
		// .CD.  a...  ...b  cC..  ..Dd 
		// ....  ....  ....  ccc.  .ddd 
		int[] a = new int[] { -newx - 1, -newx, -newx + 1, -1, +newx - 1 };
		int[] b = new int[] { -newx - 1, -newx, -newx + 1, +1, +newx + 1 };
		int[] c = new int[] { -newx - 1, -1, +newx - 1, +newx, +newx + 1 };
		int[] d = new int[] { -newx + 1, +1, +newx - 1, +newx, +newx + 1 };

		int nMaxima = 0;
		for (int y = 0; y < maxy; y += 2)
		{
			int xindex = (y + 1) * newx + 1;
			FIND_MAXIMUM: for (int x = 0; x < maxx; x += 2, xindex += 2)
			{
				int[] scan = a;
				int maxIndex = xindex;
				if (data[maxIndex] < data[xindex + 1])
				{
					scan = b;
					maxIndex = xindex + 1;
				}
				if (data[maxIndex] < data[xindex + newx])
				{
					scan = c;
					maxIndex = xindex + newx;
				}
				if (data[maxIndex] < data[xindex + newx + 1])
				{
					scan = d;
					maxIndex = xindex + newx + 1;
				}

				// Check the remaining region
				for (int offset : scan)
				{
					if (data[maxIndex] < data[maxIndex + offset])
						continue FIND_MAXIMUM;
				}

				if (validations)
				{
					if (data[maxIndex] < heightThreshold)
						continue;

					if (isNeighbourCheck())
					{
						for (int offset : scan)
						{
							if (maximaFlag[maxIndex + offset])
								continue FIND_MAXIMUM;
						}
					}

					// Check the maximum width
					if (minimumWidth > 0)
					{
						int x0 = maxIndex % maxx;
						int y0 = maxIndex / maxx;

						// Get the width at half maximum.
						int v_half = intHalfMaximum(floatBackground, data[maxIndex]);
						int index2;
						// Scan right
						int x1 = x0 + 1;
						index2 = maxIndex + 1;
						while (x1 < maxx && data[index2] > v_half)
						{
							x1++;
							index2++;
						}
						// Scan left
						int x2 = x0 - 1;
						index2 = maxIndex - 1;
						while (x2 >= 0 && data[index2] > v_half)
						{
							x2--;
							index2--;
						}
						if (x1 - x2 < minimumWidth)
							continue;
						// Scan up
						int y1 = y0 + 1;
						index2 = maxIndex + newx;
						while (y1 < maxy && data[index2] > v_half)
						{
							y1++;
							index2 += newx;
						}
						// Scan down
						int y2 = y0 - 1;
						index2 = maxIndex - newx;
						while (y2 >= 0 && data[index2] > v_half)
						{
							y2--;
							index2 -= newx;
						}
						if (y1 - y2 < minimumWidth)
							continue;
					}

					if (isNeighbourCheck())
						maximaFlag[maxIndex] = true;
				}

				// Remap the maxima
				int xx = maxIndex % newx;
				int yy = maxIndex / newx;

				//System.out.printf("blockFind3x3 [%d,%d]\n", xx-1, yy-1);
				maxima[nMaxima++] = (yy - 1) * maxx + xx - 1;
			} // end FIND_MAXIMA
		}

		return truncate(maxima, nMaxima);
	}

	/**
	 * Compute the local-maxima within a 3x3 block.
	 * An inner boundary of 1 is ignored as potential maxima on the top and left, and a boundary of 1 or 2 on the right
	 * or bottom (depending if the image is even/odd dimensions).
	 * <p>
	 * Any maxima below the configured fraction above background are ignored. Fraction = (Max-background)/Max within the
	 * 3x3 neighbourhood. Maxima below the minimum height (above background) or the minimum peak-width at half maximum
	 * in any dimension are ignored.
	 * <p>
	 * Note: The height thresholding step is ignored if the height threshold is zero.
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param border
	 *            The internal border
	 * @return Indices to the local maxima (index = maxx * y + x)
	 */
	public int[] blockFind3x3Internal(int[] data, int maxx, int maxy, int border)
	{
		if (border < 1)
			return blockFind3x3(data, maxx, maxy);

		// The number of blocks in x and y
		int xblocks = getBlocks(maxx - border, 2);
		int yblocks = getBlocks(maxy - border, 2);

		int[] maxima = new int[xblocks * yblocks];

		int heightThreshold = (int) getHeightThreshold();
		int floatBackground = (int) background;
		boolean validations = (heightThreshold > 0 || minimumWidth > 0 || isNeighbourCheck());

		boolean[] maximaFlag = null;
		if (isNeighbourCheck())
			maximaFlag = getFlagBuffer(data.length);

		// Compare 2x2 block
		// ....
		// .AB.
		// .CD.
		// ....

		// Create the scan arrays to search the remaining 5 locations around the 2x2 block maxima
		// ....  aaa.  .bbb  ....  .... 
		// .AB.  aA..  ..Bb  c...  ...d 
		// .CD.  a...  ...b  cC..  ..Dd 
		// ....  ....  ....  ccc.  .ddd 
		int[] a = new int[] { -maxx - 1, -maxx, -maxx + 1, -1, +maxx - 1 };
		int[] b = new int[] { -maxx - 1, -maxx, -maxx + 1, +1, +maxx + 1 };
		int[] c = new int[] { -maxx - 1, -1, +maxx - 1, +maxx, +maxx + 1 };
		int[] d = new int[] { -maxx + 1, +1, +maxx - 1, +maxx, +maxx + 1 };

		int nMaxima = 0;
		for (int y = border; y < maxy - border - 1; y += 2)
		{
			int xindex = y * maxx + border;
			FIND_MAXIMUM: for (int x = border; x < maxx - border - 1; x += 2, xindex += 2)
			{
				int[] scan = a;
				int maxIndex = xindex;
				if (data[maxIndex] < data[xindex + 1])
				{
					scan = b;
					maxIndex = xindex + 1;
				}
				if (data[maxIndex] < data[xindex + maxx])
				{
					scan = c;
					maxIndex = xindex + maxx;
				}
				if (data[maxIndex] < data[xindex + maxx + 1])
				{
					scan = d;
					maxIndex = xindex + maxx + 1;
				}

				// Check the remaining region
				for (int offset : scan)
				{
					if (data[maxIndex] < data[maxIndex + offset])
						continue FIND_MAXIMUM;
				}

				if (validations)
				{
					if (data[maxIndex] < heightThreshold)
						continue;

					if (isNeighbourCheck())
					{
						for (int offset : scan)
						{
							if (maximaFlag[maxIndex + offset])
								continue FIND_MAXIMUM;
						}
					}

					// Check the maximum width
					if (minimumWidth > 0)
					{
						int x0 = maxIndex % maxx;
						int y0 = maxIndex / maxx;

						// Get the width at half maximum.
						int v_half = intHalfMaximum(floatBackground, data[maxIndex]);
						int index2;
						// Scan right
						int x1 = x0 + 1;
						index2 = maxIndex + 1;
						while (x1 < maxx && data[index2] > v_half)
						{
							x1++;
							index2++;
						}
						// Scan left
						int x2 = x0 - 1;
						index2 = maxIndex - 1;
						while (x2 >= 0 && data[index2] > v_half)
						{
							x2--;
							index2--;
						}
						if (x1 - x2 < minimumWidth)
							continue;
						// Scan up
						int y1 = y0 + 1;
						index2 = maxIndex + maxx;
						while (y1 < maxy && data[index2] > v_half)
						{
							y1++;
							index2 += maxx;
						}
						// Scan down
						int y2 = y0 - 1;
						index2 = maxIndex - maxx;
						while (y2 >= 0 && data[index2] > v_half)
						{
							y2--;
							index2 -= maxx;
						}
						if (y1 - y2 < minimumWidth)
							continue;
					}

					if (isNeighbourCheck())
						maximaFlag[maxIndex] = true;
				}

				maxima[nMaxima++] = maxIndex;
			} // end FIND_MAXIMA
		}

		return truncate(maxima, nMaxima);
	}

	/**
	 * Search the data for the index of the maximum in each block of size n*n.
	 * <p>
	 * E.g. Max [ (i,i+n) x (i,j+n) ] for (i=0; i&lt;maxx; i+=n) x (j=0, j&lt;maxy; j+=n)
	 * 
	 * @param data
	 *            The input data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @return The maxima indices
	 */
	public int[] findBlockMaxima(int[] data, int maxx, int maxy, int n)
	{
		if (n == 1)
			return findBlockMaxima2x2(data, maxx, maxy);

		return findBlockMaximaNxN(data, maxx, maxy, n);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public NonMaximumSuppression clone()
	{
		try
		{
			NonMaximumSuppression o = (NonMaximumSuppression) super.clone();
			o.newDataFloat = null;
			o.newDataInt = null;
			o.maximaFlagBuffer = null;
			return o;
		}
		catch (CloneNotSupportedException e)
		{
			// Ignore
		}
		return null;
	}
}