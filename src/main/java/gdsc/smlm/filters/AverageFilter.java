package gdsc.smlm.filters;

import org.apache.commons.math3.util.FastMath;

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
 * Computes the block average for each point within the array.
 * <p>
 * block algorithm sweeps the entire (2n+1)*(2n+1) region around each pixel. Speed ~ Order(N*N).
 * <p>
 * stripedBlock algorithm uses two consecutive 1D passes using (2n+1) strips. Totals from the first pass are used in the
 * second pass resulting in a speed increase. Speed ~ Order(N). Results should match exactly the block algorithm.
 * <p>
 * rollingBlock algorithm uses two consecutive 1D passes using (2n+1) strips. Each pass is computed using a rolling
 * total thus each pixel sum can be computed using a single addition and subtraction of the end pixels of the strip. Due
 * to cumulative error of the rolling sum the results may differ from the other algorithms for large images (applies to
 * the the float version since integer arithmetic should be robust within Integer.MAX bounds). Speed ~ Order(1).
 * <p>
 * Note: Due to lack of small dimension checking the routines will fail if maxx or maxy are less than 2. All routines
 * are OK for 3x3 images and larger.
 */
public class AverageFilter implements Cloneable
{
	private float[] floatDataBuffer = null;
	private float[] floatRowBuffer = null;

	/**
	 * Compute the block average within a 2n+1 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 */
	public void rollingBlockAverageInternal(float[] data, final int maxx, final int maxy, final int n)
	{
		// Note: Speed tests show that this method is only marginally faster than rollingBlockAverageNxNInternal.
		// Sometimes it is slower. The intricacies of the java optimiser escape me.
		if (n == 1)
			rollingBlockAverage3x3Internal(data, maxx, maxy);
		else
			rollingBlockAverageNxNInternal(data, maxx, maxy, n);
	}

	/**
	 * Compute the block average within a 2n+1 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 */
	public void rollingBlockAverageNxNInternal(float[] data, final int maxx, final int maxy, final int n)
	{
		final int blockSize = 2 * n + 1;
		if (maxx < blockSize || maxy < blockSize)
			return;

		float[] newData = floatBuffer(data.length);

		final float divisor = (float) (1.0 / (blockSize * blockSize));

		// X-direction
		for (int y = 0; y < maxy; y++)
		{
			// Initialise the rolling sum
			float sum = 0;

			int endIndex = y * maxx;
			int x = 0;
			while (x < blockSize)
			{
				sum += data[endIndex];
				endIndex++;
				x++;
			}

			// Rolling sum over the X-direction
			int startIndex = y * maxx;
			int centreIndex = startIndex + n;

			newData[centreIndex] = sum;

			while (x < maxx)
			{
				centreIndex++;

				sum += data[endIndex] - data[startIndex];

				newData[centreIndex] = sum;

				x++;
				startIndex++;
				endIndex++;
			}
		}

		// Y-direction. 
		// Only sweep over the interior
		for (int x = n; x < maxx - n; x++)
		{
			// Initialise the rolling sum
			float sum = 0;

			int endIndex = x;
			int y = 0;
			while (y < blockSize)
			{
				sum += newData[endIndex];
				endIndex += maxx;
				y++;
			}

			// Rolling sum over the Y-direction
			int startIndex = x;
			int centreIndex = startIndex + n * maxx;

			data[centreIndex] = sum * divisor;

			while (y < maxy)
			{
				centreIndex += maxx;

				sum += newData[endIndex] - newData[startIndex];

				data[centreIndex] = sum * divisor;

				y++;
				startIndex += maxx;
				endIndex += maxx;
			}
		}
	}

	/**
	 * Compute the block average within a 3x3 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 */
	public void rollingBlockAverage3x3Internal(float[] data, final int maxx, final int maxy)
	{
		if (maxx < 3 || maxy < 3)
			return;

		float[] newData = floatBuffer(data.length);

		final float divisor = (float) (1.0 / 9);

		// X-direction
		for (int y = 0; y < maxy; y++)
		{
			// Initialise the rolling sum
			int startIndex = y * maxx;
			int centreIndex = startIndex + 1;
			int endIndex = centreIndex + 1;
			float sum = data[startIndex] + data[centreIndex] + data[endIndex];

			// Rolling sum over the X-direction
			newData[centreIndex++] = sum;

			for (int x = 0; x < maxx - 3; x++)
			{
				sum += data[++endIndex] - data[startIndex++];
				newData[centreIndex++] = sum;
			}
		}

		// Y-direction. 
		// Only sweep over the interior
		for (int x = 1; x < maxx - 1; x++)
		{
			// Initialise the rolling sum
			int startIndex = x;
			int centreIndex = startIndex + maxx;
			int endIndex = centreIndex + maxx;
			float sum = newData[startIndex] + newData[centreIndex] + newData[endIndex];

			// Rolling sum over the Y-direction
			data[centreIndex] = sum * divisor;

			for (int y = 0; y < maxy - 3; y++)
			{
				centreIndex += maxx;
				endIndex += maxx;
				sum += newData[endIndex] - newData[startIndex];
				data[centreIndex] = sum * divisor;
				startIndex += maxx;
			}
		}
	}

	/**
	 * Compute the block average within a 2n+1 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 */
	@Deprecated
	public void rollingBlockAverageNxNInternalTransposed(float[] data, final int maxx, final int maxy, final int n)
	{
		int blockSize = 2 * n + 1;
		if (maxx < blockSize || maxy < blockSize)
			return;

		float[] newData = floatBuffer(data.length);

		final float divisor = (float) (1.0 / ((2 * n + 1) * (2 * n + 1)));

		// NOTE: 
		// To increase speed when sweeping the arrays:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		for (int y = 0; y < maxy; y++)
		{
			// Initialise the rolling sum
			float sum = 0;

			int endIndex = y * maxx;
			int x = 0;
			while (x < blockSize)
			{
				sum += data[endIndex];
				endIndex++;
				x++;
			}

			// Rolling sum over the X-direction
			int startIndex = y * maxx;
			int newCentreIndex = y + n * maxy;

			newData[newCentreIndex] = sum;

			while (x < maxx)
			{
				newCentreIndex += maxy;

				sum += data[endIndex] - data[startIndex];

				newData[newCentreIndex] = sum;

				x++;
				startIndex++;
				endIndex++;
			}
		}

		// Y-direction. 
		// Only sweep over the interior
		for (int x = n; x < maxx - n; x++)
		{
			// Initialise the rolling sum
			float sum = 0;

			int endIndex = x * maxy;
			int y = 0;
			while (y < blockSize)
			{
				sum += newData[endIndex];
				endIndex++;
				y++;
			}

			// Rolling sum over the Y-direction
			int startIndex = x * maxy;
			int centreIndex = x + n * maxx;

			data[centreIndex] = sum * divisor;

			while (y < maxy)
			{
				centreIndex += maxx;

				sum += newData[endIndex] - newData[startIndex];

				data[centreIndex] = sum * divisor;

				y++;
				startIndex++;
				endIndex++;
			}
		}
	}

	/**
	 * Compute the block average within a 2n+1 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 */
	public void stripedBlockAverageInternal(float[] data, final int maxx, final int maxy, final int n)
	{
		if (n == 1)
			stripedBlockAverage3x3Internal(data, maxx, maxy);
		else
			stripedBlockAverageNxNInternal(data, maxx, maxy, n);
	}

	/**
	 * Compute the block average within a 2n+1 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 */
	public void stripedBlockAverageNxNInternal(float[] data, final int maxx, final int maxy, final int n)
	{
		final int blockSize = 2 * n + 1;
		if (maxx < blockSize || maxy < blockSize)
			return;

		float[] newData = floatBuffer(data.length);

		final float divisor = (float) (1.0 / (blockSize * blockSize));

		// NOTE: 
		// To increase speed when sweeping the arrays:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		for (int y = 0; y < maxy; y++)
		{
			int index = y * maxx;
			for (int x = 0; x <= maxx - blockSize; x++, index++)
			{
				float sum = 0;
				for (int x2 = 0; x2 < blockSize; x2++)
				{
					sum += data[index + x2];
				}
				newData[(x + n) * maxy + y] = sum;
			}
		}

		// Y-direction. 
		// Only sweep over the interior
		for (int x = n; x < maxx - n; x++)
		{
			int index = x * maxy;
			for (int y = 0; y <= maxy - blockSize; y++, index++)
			{
				float sum = 0;
				for (int y2 = 0; y2 < blockSize; y2++)
				{
					sum += newData[index + y2];
				}
				data[x + (y + n) * maxx] = sum * divisor;
			}
		}
	}

	/**
	 * Compute the block average within a 3x3 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 */
	public void stripedBlockAverage3x3Internal(float[] data, final int maxx, final int maxy)
	{
		if (maxx < 3 || maxy < 3)
			return;

		float[] newData = floatBuffer(data.length);

		final float divisor = (float) (1.0 / 9);

		// NOTE: 
		// To increase speed when sweeping the arrays:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		for (int y = 0; y < maxy; y++)
		{
			int index = y * maxx;
			int index2 = maxy + y;
			for (int x = 0; x <= maxx - 3; x++, index++)
			{
				float sum = data[index] + data[index + 1] + data[index + 2];
				newData[index2] = sum;
				index2 += maxy;
			}
		}

		// Y-direction. 
		// Only sweep over the interior
		for (int x = 1; x < maxx - 1; x++)
		{
			int index = x * maxy;
			int index2 = x + maxx;
			for (int y = 0; y <= maxy - 3; y++, index++)
			{
				float sum = newData[index] + newData[index + 1] + newData[index + 2];
				data[index2] = sum * divisor;
				index2 += maxx;
			}
		}
	}

	/**
	 * Compute the block average within a 2n+1 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 */
	public void blockAverageInternal(float[] data, final int maxx, final int maxy, final int n)
	{
		if (n == 1)
			blockAverage3x3Internal(data, maxx, maxy);
		else
			blockAverageNxNInternal(data, maxx, maxy, n);
	}

	/**
	 * Compute the block average within a 2w+1 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param w
	 *            The block size
	 */
	public void blockAverageInternal(float[] data, final int maxx, final int maxy, final float w)
	{
		if (w < 1)
			blockAverage3x3Internal(data, maxx, maxy, w);
		else
			blockAverageNxNInternal(data, maxx, maxy, w);
	}

	/**
	 * Compute the block average within a 2n+1 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 */
	public void blockAverageNxNInternal(float[] data, final int maxx, final int maxy, final int n)
	{
		// The size of the region 
		final int nX = (2 * n + 1);
		final int nY = (2 * n + 1);

		if (maxx < nX || maxy < nY)
			return;

		float[] newData = floatBuffer(data.length);

		int[] offset = new int[nX * nY - 1];
		for (int y = -n, d = 0; y <= n; y++)
			for (int x = -n; x <= n; x++)
				if (x != 0 || y != 0)
				{
					offset[d] = maxx * y + x;
					d++;
				}

		final float divisor = (float) (1.0 / ((2 * n + 1) * (2 * n + 1)));

		for (int y = n; y < maxy - n; y++)
		{
			int index = y * maxx + n;
			for (int x = n; x < maxx - n; x++, index++)
			{
				float sum = data[index];

				// Sweep neighbourhood - 
				// No check for boundaries as this should be an internal sweep.
				for (int offset_d : offset)
				{
					sum += data[index + offset_d];
				}

				newData[index] = sum * divisor;
			}
		}

		// Copy back
		for (int y = n; y < maxy - n; y++)
		{
			int index = y * maxx + n;
			for (int x = n; x < maxx - n; x++, index++)
			{
				data[index] = newData[index];
			}
		}
	}

	/**
	 * Compute the block average within a 2w+1 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Uses a normalised [[w*w, w, ..., w, w*w], [w, 1, ..., 1, w], [w*w, w, ..., w, w*w]] convolution kernel.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param w
	 *            The block size
	 */
	public void blockAverageNxNInternal(float[] data, final int maxx, final int maxy, final float w)
	{
		final int n = (int) w;
		final int n1 = (n == w) ? n : n + 1;

		if (n == n1)
		{
			// There is no edge
			blockAverage(data, maxx, maxy, n);
			return;
		}

		// The size of the region 
		final int nX = (2 * n1 + 1);
		final int nY = (2 * n1 + 1);

		if (maxx < nX || maxy < nY)
			return;

		final float[] newData = floatBuffer(data.length);

		// Boundary control
		final int xwidth = n1;
		final int ywidth = n1;
		final int xlimit = maxx - xwidth;
		final int ylimit = maxy - ywidth;

		// Inner block
		final int[] offset = new int[(2 * xwidth - 1) * (2 * ywidth - 1) - 1];
		final int[] xoffset = new int[offset.length];
		final int[] yoffset = new int[offset.length];
		for (int y = -ywidth + 1, d = 0; y < ywidth; y++)
			for (int x = -xwidth + 1; x < xwidth; x++)
				if (x != 0 || y != 0)
				{
					offset[d] = maxx * y + x;
					xoffset[d] = x;
					yoffset[d] = y;
					d++;
				}

		// Edges
		int j = 0;
		final int size = 2 * ((2 * xwidth - 1) + (2 * ywidth - 1));
		final int[] xoffset1 = new int[size];
		final int[] yoffset1 = new int[size];
		final int[] offset1 = new int[size];
		for (int y = -ywidth + 1; y < ywidth; y++)
		{
			yoffset1[j] = yoffset1[j + 1] = y;
			xoffset1[j] = -xwidth;
			xoffset1[j + 1] = xwidth;
			offset1[j] = maxx * y - xwidth;
			offset1[j + 1] = maxx * y + xwidth;
			j += 2;
		}
		for (int x = -xwidth + 1; x < xwidth; x++)
		{
			xoffset1[j] = xoffset1[j + 1] = x;
			yoffset1[j] = -ywidth;
			yoffset1[j + 1] = ywidth;
			offset1[j] = maxx * -ywidth + x;
			offset1[j + 1] = maxx * ywidth + x;
			j += 2;
		}

		// Corners
		final int[] xoffset2 = new int[] { -xwidth, -xwidth, xwidth, xwidth };
		final int[] yoffset2 = new int[] { -ywidth, ywidth, -ywidth, ywidth };
		final int[] offset2 = new int[xoffset2.length];
		for (int d = xoffset2.length; d-- > 0;)
			offset2[d] = maxx * yoffset2[d] + xoffset2[d];

		final float divisor = (float) (1.0 / ((2 * w + 1) * (2 * w + 1)));

		final float w1 = w - n;
		final float w2 = w1 * w1;
		for (int y = n1; y < ylimit; y++)
		{
			int index = y * maxx + n1;
			for (int x = n1; x < xlimit; x++, index++)
			{
				float sum = data[index];
				float sum1 = 0;
				float sum2 = 0;

				// Sweep neighbourhood
				// No check for boundaries as this should be an internal sweep.
				for (int d = offset.length; d-- > 0;)
					sum += data[index + offset[d]];
				for (int d = offset1.length; d-- > 0;)
					sum1 += data[index + offset1[d]];
				for (int d = offset2.length; d-- > 0;)
					sum2 += data[index + offset2[d]];
				
				newData[index] = (sum + sum1 * w1 + sum2 * w2) * divisor;
			}
		}

		// Copy back
		for (int y = n1; y < ylimit; y++)
		{
			int index = y * maxx + n1;
			for (int x = n1; x < xlimit; x++, index++)
			{
				data[index] = newData[index];
			}
		}
	}

	/**
	 * Compute the block average within a 3x3 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 */
	public void blockAverage3x3Internal(float[] data, final int maxx, final int maxy)
	{
		float[] newData = floatBuffer(data.length);

		final float divisor = (float) (1.0 / 9);

		for (int y = 1; y < maxy - 1; y++)
		{
			int index0 = (y - 1) * maxx + 1;
			int index1 = y * maxx + 1;
			int index2 = (y + 1) * maxx + 1;
			for (int x = 1; x < maxx - 1; x++)
			{
				float sum = data[index0 - 1] + data[index0] + data[index0 + 1] + data[index1 - 1] + data[index1] +
						data[index1 + 1] + data[index2 - 1] + data[index2] + data[index2 + 1];
				newData[index1] = sum * divisor;
				index0++;
				index1++;
				index2++;
			}
		}

		// Copy back
		for (int y = 1; y < maxy - 1; y++)
		{
			int index = y * maxx + 1;
			for (int x = 1; x < maxx - 1; x++, index++)
			{
				data[index] = newData[index];
			}
		}
	}

	/**
	 * Compute the weighted block average within a 3x3 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Uses a normalised [[w*w, w, w*w], [w, 1, w], [w*w, w, w*w]] convolution kernel.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param w
	 *            The weight
	 */
	public void blockAverage3x3Internal(float[] data, final int maxx, final int maxy, final float w)
	{
		float[] newData = floatBuffer(data.length);

		final float w2 = w * w;
		final float divisor = (float) (1.0 / (1 + 4 * w + 4 * w2));

		for (int y = 1; y < maxy - 1; y++)
		{
			int index0 = (y - 1) * maxx + 1;
			int index1 = y * maxx + 1;
			int index2 = (y + 1) * maxx + 1;
			for (int x = 1; x < maxx - 1; x++)
			{
				// Edges
				float sum1 = data[index0] + data[index1 - 1] + data[index1 + 1] + data[index2];
				// Corners
				float sum2 = data[index0 - 1] + data[index0 + 1] + data[index2 - 1] + data[index2 + 1];

				newData[index1] = (data[index1] + sum1 * w + sum2 * w2) * divisor;
				index0++;
				index1++;
				index2++;
			}
		}

		// Copy back
		for (int y = 1; y < maxy - 1; y++)
		{
			int index = y * maxx + 1;
			for (int x = 1; x < maxx - 1; x++, index++)
			{
				data[index] = newData[index];
			}
		}
	}

	/**
	 * Compute an approximate Gaussian convolution within a 3x3 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Uses a normalised [[1, 2, 1], [2, 4, 2], [1, 2, 1]] convolution kernel.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 */
	public void blockGaussian3x3Internal(float[] data, final int maxx, final int maxy)
	{
		float[] newData = floatBuffer(data.length);

		final float divisor = (float) (1.0 / 16);

		for (int y = 1; y < maxy - 1; y++)
		{
			int index0 = (y - 1) * maxx + 1;
			int index1 = y * maxx + 1;
			int index2 = (y + 1) * maxx + 1;
			for (int x = 1; x < maxx - 1; x++)
			{
				float sum = data[index0 - 1] + 2 * data[index0] + data[index0 + 1] + 2 * data[index1 - 1] + 4 *
						data[index1] + 2 * data[index1 + 1] + data[index2 - 1] + 2 * data[index2] + data[index2 + 1];
				newData[index1] = sum * divisor;
				index0++;
				index1++;
				index2++;
			}
		}

		// Copy back
		for (int y = 1; y < maxy - 1; y++)
		{
			int index = y * maxx + 1;
			for (int x = 1; x < maxx - 1; x++, index++)
			{
				data[index] = newData[index];
			}
		}
	}

	private float[] floatBuffer(int size)
	{
		if (floatDataBuffer == null || floatDataBuffer.length < size)
		{
			floatDataBuffer = new float[size];
		}
		return floatDataBuffer;
	}

	/**
	 * Compute the block average within a 2n+1 size block around each point.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 */
	public void rollingBlockAverage(float[] data, final int maxx, final int maxy, final int n)
	{
		if (n == 1)
			rollingBlockAverage3x3(data, maxx, maxy);
		else
			rollingBlockAverageNxN(data, maxx, maxy, n);
	}

	/**
	 * Compute the block average within a 2n+1 size block around each point.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 */
	public void rollingBlockAverageNxN(float[] data, final int maxx, final int maxy, final int n)
	{
		float[] newData = floatBuffer(data.length);

		final float divisor = (float) (1.0 / ((2 * n + 1) * (2 * n + 1)));

		// NOTE: 
		// To increase speed when sweeping the arrays and allow for reusing code:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		int width = maxx;
		int height = maxy;
		float[] inData = data;
		float[] outData = newData;

		float[] row = floatRowBuffer(width + 2 * n);
		for (int y = 0; y < height; y++)
		{
			extractRow(inData, y, width, n, row);

			// Initialise rolling sum
			float sum = (n + 1) * row[0];
			int endIndex = n + 1;
			for (int i = 0; i < n; i++)
			{
				sum += row[endIndex++];
			}

			int centreIndex = y;
			outData[centreIndex] = sum;

			// Rolling sum over the X-direction
			for (int x = 0; x < width - 1; x++)
			{
				sum += row[endIndex++] - row[x];
				centreIndex += height;
				outData[centreIndex] = sum;
			}
		}

		// Y-direction. 
		width = maxy;
		height = maxx;
		inData = newData;
		outData = data;

		row = floatRowBuffer(width + 2 * n);
		for (int y = 0; y < height; y++)
		{
			extractRow(inData, y, width, n, row);

			// Initialise rolling sum
			float sum = (n + 1) * row[0];
			int endIndex = n + 1;
			for (int i = 0; i < n; i++)
			{
				sum += row[endIndex++];
			}

			int centreIndex = y;
			outData[centreIndex] = sum * divisor;

			// Rolling sum over the X-direction
			for (int x = 0; x < width - 1; x++)
			{
				sum += row[endIndex++] - row[x];
				centreIndex += height;
				outData[centreIndex] = sum * divisor;
			}
		}
	}

	/**
	 * Compute the block average within a 3x3 size block around each point.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 */
	public void rollingBlockAverage3x3(float[] data, final int maxx, final int maxy)
	{
		float[] newData = floatBuffer(data.length);

		final float divisor = (float) (1.0 / 9);

		// NOTE: 
		// To increase speed when sweeping the arrays and allow for reusing code:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		int width = maxx;
		int height = maxy;
		float[] inData = data;
		float[] outData = newData;

		float[] row = floatRowBuffer(width + 2);
		for (int y = 0; y < height; y++)
		{
			extractRow1(inData, y, width, row);

			// Initialise rolling sum
			float sum = 2 * row[0] + row[2];
			int endIndex = 3;

			int centreIndex = y;
			outData[centreIndex] = sum;

			// Rolling sum over the X-direction
			for (int x = 0; x < width - 1; x++)
			{
				sum += row[endIndex++] - row[x];
				centreIndex += height;
				outData[centreIndex] = sum;
			}
		}

		// Y-direction. 
		width = maxy;
		height = maxx;
		inData = newData;
		outData = data;

		row = floatRowBuffer(width + 2);
		for (int y = 0; y < height; y++)
		{
			extractRow1(inData, y, width, row);

			// Initialise rolling sum
			float sum = 2 * row[0] + row[2];
			int endIndex = 3;

			int centreIndex = y;
			outData[centreIndex] = sum * divisor;

			// Rolling sum over the X-direction
			for (int x = 0; x < width - 1; x++)
			{
				sum += row[endIndex++] - row[x];
				centreIndex += height;
				outData[centreIndex] = sum * divisor;
			}
		}
	}

	/**
	 * Compute the block average within a 2n+1 size block around each point.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 */
	public void stripedBlockAverage(float[] data, final int maxx, final int maxy, final int n)
	{
		if (n == 1)
			stripedBlockAverage3x3(data, maxx, maxy);
		else
			stripedBlockAverageNxN(data, maxx, maxy, n);
	}

	/**
	 * Compute the block average within a 2n+1 size block around each point.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 */
	public void stripedBlockAverageNxN(float[] data, final int maxx, final int maxy, final int n)
	{
		int blockSize = 2 * n + 1;

		float[] newData = floatBuffer(data.length);

		final float divisor = (float) (1.0 / ((2 * n + 1) * (2 * n + 1)));

		// NOTE: 
		// To increase speed when sweeping the arrays and allow for reusing code:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		int width = maxx;
		int height = maxy;
		float[] inData = data;
		float[] outData = newData;

		float[] row = floatRowBuffer(width + 2 * n);
		for (int y = 0; y < height; y++)
		{
			extractRow(inData, y, width, n, row);

			int centreIndex = y;
			for (int x = 0; x < width; x++)
			{
				// Sum strips
				float sum = 0;

				for (int j = 0; j < blockSize; j++)
				{
					sum += row[x + j];
				}

				// Store result in transpose
				outData[centreIndex] = sum;
				centreIndex += height;
			}
		}

		// Y-direction. 
		width = maxy;
		height = maxx;
		inData = newData;
		outData = data;

		row = floatRowBuffer(width + 2 * n);
		for (int y = 0; y < height; y++)
		{
			// Extract row (pad ends)
			extractRow(inData, y, width, n, row);

			int centreIndex = y;
			for (int x = 0; x < width; x++)
			{
				// Sum strips
				float sum = 0;

				for (int j = 0; j < blockSize; j++)
				{
					sum += row[x + j];
				}

				// Store result in transpose
				outData[centreIndex] = sum * divisor;
				centreIndex += height;
			}
		}
	}

	/**
	 * Compute the block average within a 3x3 size block around each point.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 */
	public void stripedBlockAverage3x3(float[] data, final int maxx, final int maxy)
	{
		float[] newData = floatBuffer(data.length);

		final float divisor = (float) (1.0 / 9);

		// NOTE: 
		// To increase speed when sweeping the arrays and allow for reusing code:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		int width = maxx;
		int height = maxy;
		float[] inData = data;
		float[] outData = newData;

		float[] row = floatRowBuffer(width + 2);
		for (int y = 0; y < height; y++)
		{
			extractRow1(inData, y, width, row);

			int centreIndex = y;
			for (int x = 0; x < width; x++)
			{
				// Sum strips
				float sum = row[x] + row[x + 1] + row[x + 2];

				// Store result in transpose
				outData[centreIndex] = sum;
				centreIndex += height;
			}
		}

		// Y-direction. 
		width = maxy;
		height = maxx;
		inData = newData;
		outData = data;

		row = floatRowBuffer(width + 2);
		for (int y = 0; y < height; y++)
		{
			// Extract row (pad ends)
			extractRow1(inData, y, width, row);

			int centreIndex = y;
			for (int x = 0; x < width; x++)
			{
				// Sum strips
				float sum = row[x] + row[x + 1] + row[x + 2];

				// Store result in transpose
				outData[centreIndex] = sum * divisor;
				centreIndex += height;
			}
		}
	}

	private float[] floatRowBuffer(int size)
	{
		if (floatRowBuffer == null || floatRowBuffer.length < size)
		{
			floatRowBuffer = new float[size];
		}
		return floatRowBuffer;
	}

	private void extractRow(float[] inData, int y, int width, final int n, float[] row)
	{
		final int index = y * width;

		// Pad ends
		for (int i = 0; i < n; i++)
		{
			row[i] = inData[index];
			row[i + n + width] = inData[index + width - 1];
		}

		// Fill in data
		System.arraycopy(inData, index, row, n, width);
		//for (int x = 0, i = n; x < width; x++)
		//{
		//	row[i++] = inData[index++];
		//}
	}

	private void extractRow1(float[] inData, int y, int width, float[] row)
	{
		final int index = y * width;

		// Pad ends
		row[0] = inData[index];
		row[1 + width] = inData[index + width - 1];

		// Fill in data
		System.arraycopy(inData, index, row, 1, width);
		//for (int x = 0, i = 1; x < width; x++)
		//{
		//	row[i++] = inData[index++];
		//}
	}

	/**
	 * Compute the block average within a 2n+1 size block around each point.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 */
	public void blockAverage(float[] data, final int maxx, final int maxy, final int n)
	{
		if (n == 1)
			blockAverage3x3(data, maxx, maxy);
		else
			blockAverageNxN(data, maxx, maxy, n);
	}

	/**
	 * Compute the block average within a 2w+1 size block around each point.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param w
	 *            The block size
	 */
	public void blockAverage(float[] data, final int maxx, final int maxy, final float w)
	{
		if (w < 1)
			blockAverage3x3(data, maxx, maxy, w);
		else
			blockAverageNxN(data, maxx, maxy, w);
	}

	/**
	 * Compute the block average within a 2n+1 size block around each point.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 */
	public void blockAverageNxN(float[] data, final int maxx, final int maxy, final int n)
	{
		float[] newData = floatBuffer(data.length);

		// Boundary control
		final int xwidth = FastMath.min(n, maxx - 1);
		final int ywidth = FastMath.min(n, maxy - 1);
		final int xlimit = maxx - xwidth;
		final int ylimit = maxy - ywidth;

		int[] offset = new int[(2 * xwidth + 1) * (2 * ywidth + 1) - 1];
		int[] xoffset = new int[offset.length];
		int[] yoffset = new int[offset.length];
		for (int y = -ywidth, d = 0; y <= ywidth; y++)
			for (int x = -xwidth; x <= xwidth; x++)
				if (x != 0 || y != 0)
				{
					offset[d] = maxx * y + x;
					xoffset[d] = x;
					yoffset[d] = y;
					d++;
				}

		final float divisor = (float) (1.0 / ((2 * n + 1) * (2 * n + 1)));

		int index = 0;
		for (int y = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++, index++)
			{
				float sum = data[index];

				// Flag to indicate this pixels has a complete (2n+1) neighbourhood 
				boolean isInnerXY = (y >= ywidth && y < ylimit) && (x >= xwidth && x < xlimit);

				// Sweep neighbourhood
				if (isInnerXY)
				{
					for (int d = offset.length; d-- > 0;)
					{
						sum += data[index + offset[d]];
					}
				}
				else
				{
					for (int d = offset.length; d-- > 0;)
					{
						// Get the pixel with boundary checking
						int yy = y + yoffset[d];
						int xx = x + xoffset[d];
						if (xx <= 0)
							xx = 0;
						else if (xx >= maxx)
							xx = maxx - 1;
						if (yy <= 0)
							yy = 0;
						else if (yy >= maxy)
							yy = maxy - 1;
						sum += data[xx + yy * maxx];
					}
				}
				newData[index] = sum * divisor;
			}
		}

		// Copy back
		System.arraycopy(newData, 0, data, 0, data.length);
	}

	/**
	 * Compute the block average within a 2w+1 size block around each point.
	 * <p>
	 * Uses a normalised [[w*w, w, ..., w, w*w], [w, 1, ..., 1, w], [w*w, w, ..., w, w*w]] convolution kernel.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param w
	 *            The block size
	 */
	public void blockAverageNxN(float[] data, final int maxx, final int maxy, final float w)
	{
		final int n = (int) w;
		final int n1 = (n == w) ? n : n + 1;

		if (n == n1)
		{
			// There is no edge
			blockAverage(data, maxx, maxy, n);
			return;
		}

		final float[] newData = floatBuffer(data.length);

		// Boundary control
		final int xwidth = FastMath.min(n1, maxx - 1);
		final int ywidth = FastMath.min(n1, maxy - 1);
		final int xlimit = maxx - xwidth;
		final int ylimit = maxy - ywidth;

		// Inner block
		final int[] offset = new int[(2 * xwidth - 1) * (2 * ywidth - 1) - 1];
		final int[] xoffset = new int[offset.length];
		final int[] yoffset = new int[offset.length];
		for (int y = -ywidth + 1, d = 0; y < ywidth; y++)
			for (int x = -xwidth + 1; x < xwidth; x++)
				if (x != 0 || y != 0)
				{
					offset[d] = maxx * y + x;
					xoffset[d] = x;
					yoffset[d] = y;
					d++;
				}

		// Edges
		int j = 0;
		final int size = 2 * ((2 * xwidth - 1) + (2 * ywidth - 1));
		final int[] xoffset1 = new int[size];
		final int[] yoffset1 = new int[size];
		final int[] offset1 = new int[size];
		for (int y = -ywidth + 1; y < ywidth; y++)
		{
			yoffset1[j] = yoffset1[j + 1] = y;
			xoffset1[j] = -xwidth;
			xoffset1[j + 1] = xwidth;
			offset1[j] = maxx * y - xwidth;
			offset1[j + 1] = maxx * y + xwidth;
			j += 2;
		}
		for (int x = -xwidth + 1; x < xwidth; x++)
		{
			xoffset1[j] = xoffset1[j + 1] = x;
			yoffset1[j] = -ywidth;
			yoffset1[j + 1] = ywidth;
			offset1[j] = maxx * -ywidth + x;
			offset1[j + 1] = maxx * ywidth + x;
			j += 2;
		}

		// Corners
		final int[] xoffset2 = new int[] { -xwidth, -xwidth, xwidth, xwidth };
		final int[] yoffset2 = new int[] { -ywidth, ywidth, -ywidth, ywidth };
		final int[] offset2 = new int[xoffset2.length];
		for (int d = xoffset2.length; d-- > 0;)
			offset2[d] = maxx * yoffset2[d] + xoffset2[d];

		final float divisor = (float) (1.0 / ((2 * w + 1) * (2 * w + 1)));

		final float w1 = w - n;
		final float w2 = w1 * w1;
		int index = 0;
		for (int y = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++, index++)
			{
				float sum = data[index];
				float sum1 = 0;
				float sum2 = 0;

				// Flag to indicate this pixels has a complete (2n1+1) neighbourhood 
				boolean isInnerXY = (y >= ywidth && y < ylimit) && (x >= xwidth && x < xlimit);

				// Sweep neighbourhood
				if (isInnerXY)
				{
					for (int d = offset.length; d-- > 0;)
						sum += data[index + offset[d]];
					for (int d = offset1.length; d-- > 0;)
						sum1 += data[index + offset1[d]];
					for (int d = offset2.length; d-- > 0;)
						sum2 += data[index + offset2[d]];
				}
				else
				{
					// Get the pixel with boundary checking

					// Inner block
					for (int d = offset.length; d-- > 0;)
					{
						int yy = y + yoffset[d];
						int xx = x + xoffset[d];
						if (xx <= 0)
							xx = 0;
						else if (xx >= maxx)
							xx = maxx - 1;
						if (yy <= 0)
							yy = 0;
						else if (yy >= maxy)
							yy = maxy - 1;
						sum += data[xx + yy * maxx];
					}
					// Edges
					for (int d = offset1.length; d-- > 0;)
					{
						int yy = y + yoffset1[d];
						int xx = x + xoffset1[d];
						if (xx <= 0)
							xx = 0;
						else if (xx >= maxx)
							xx = maxx - 1;
						if (yy <= 0)
							yy = 0;
						else if (yy >= maxy)
							yy = maxy - 1;
						sum1 += data[xx + yy * maxx];
					}
					// Corners
					for (int d = offset2.length; d-- > 0;)
					{
						int yy = y + yoffset2[d];
						int xx = x + xoffset2[d];
						if (xx <= 0)
							xx = 0;
						else if (xx >= maxx)
							xx = maxx - 1;
						if (yy <= 0)
							yy = 0;
						else if (yy >= maxy)
							yy = maxy - 1;
						sum2 += data[xx + yy * maxx];
					}
				}
				newData[index] = (sum + sum1 * w1 + sum2 * w2) * divisor;
			}
		}

		// Copy back
		System.arraycopy(newData, 0, data, 0, data.length);
	}

	/**
	 * Compute the block average within a 3x3 size block around each point.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 */
	public void blockAverage3x3(float[] data, final int maxx, final int maxy)
	{
		float[] newData = floatBuffer(data.length);

		// Boundary control
		final int xwidth = 1;
		final int ywidth = 1;
		final int xlimit = maxx - xwidth;
		final int ylimit = maxy - ywidth;

		int[] offset = new int[8];
		int[] xoffset = new int[offset.length];
		int[] yoffset = new int[offset.length];
		for (int y = -ywidth, d = 0; y <= ywidth; y++)
			for (int x = -xwidth; x <= xwidth; x++)
				if (x != 0 || y != 0)
				{
					offset[d] = maxx * y + x;
					xoffset[d] = x;
					yoffset[d] = y;
					d++;
				}

		final float divisor = (float) (1.0 / 9);

		for (int y = 0; y < maxy; y++)
		{
			int index0 = (y - 1) * maxx;
			int index1 = y * maxx;
			int index2 = (y + 1) * maxx;
			for (int x = 0; x < maxx; x++)
			{
				// Flag to indicate this pixels has a complete (2n+1) neighbourhood 
				boolean isInnerXY = (y >= ywidth && y < ylimit) && (x >= xwidth && x < xlimit);

				// Sweep neighbourhood
				if (isInnerXY)
				{
					float sum = data[index0 - 1] + data[index0] + data[index0 + 1] + data[index1 - 1] + data[index1] +
							data[index1 + 1] + data[index2 - 1] + data[index2] + data[index2 + 1];
					newData[index1] = sum * divisor;
				}
				else
				{
					float sum = data[index1];

					for (int d = offset.length; d-- > 0;)
					{
						// Get the pixel with boundary checking
						int yy = y + yoffset[d];
						int xx = x + xoffset[d];
						if (xx <= 0)
							xx = 0;
						else if (xx >= maxx)
							xx = maxx - 1;
						if (yy <= 0)
							yy = 0;
						else if (yy >= maxy)
							yy = maxy - 1;
						sum += data[xx + yy * maxx];
					}
					newData[index1] = sum * divisor;
				}
				index0++;
				index1++;
				index2++;
			}
		}

		// Copy back
		System.arraycopy(newData, 0, data, 0, data.length);
	}

	/**
	 * Compute the weighted block average within a 3x3 size block around each point.
	 * <p>
	 * Uses a normalised [[w*w, w, w*w], [w, 1, w], [w*w, w, w*w]] convolution kernel.
	 * <p>
	 * Note: the input data is destructively modified.
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param w
	 *            The weight
	 */
	public void blockAverage3x3(float[] data, final int maxx, final int maxy, final float w)
	{
		float[] newData = floatBuffer(data.length);

		// Boundary control
		final int xlimit = maxx - 1;
		final int ylimit = maxy - 1;

		// Edges
		final int[] xoffset = new int[] { -1, 0, 0, 1 };
		final int[] yoffset = new int[] { 0, -1, 1, 0 };
		// Corners
		final int[] xoffset2 = new int[] { -1, -1, 1, 1 };
		final int[] yoffset2 = new int[] { -1, 1, -1, 1 };

		final float w2 = w * w;
		final float divisor = (float) (1.0 / (1 + 4 * w + 4 * w2));

		for (int y = 0; y < maxy; y++)
		{
			int index0 = (y - 1) * maxx;
			int index1 = y * maxx;
			int index2 = (y + 1) * maxx;
			for (int x = 0; x < maxx; x++)
			{
				// Flag to indicate this pixels has a complete (2n+1) neighbourhood 
				boolean isInnerXY = (y > 0 && y < ylimit) && (x > 0 && x < xlimit);

				// Sweep neighbourhood
				float sum1 = 0;
				float sum2 = 0;
				if (isInnerXY)
				{
					// Edges
					sum1 = data[index0] + data[index1 - 1] + data[index1 + 1] + data[index2];
					// Corners
					sum2 = data[index0 - 1] + data[index0 + 1] + data[index2 - 1] + data[index2 + 1];
				}
				else
				{
					// Get the pixel with boundary checking

					// Edges
					for (int d = xoffset.length; d-- > 0;)
					{
						int yy = y + yoffset[d];
						int xx = x + xoffset[d];
						if (xx <= 0)
							xx = 0;
						else if (xx >= maxx)
							xx = maxx - 1;
						if (yy <= 0)
							yy = 0;
						else if (yy >= maxy)
							yy = maxy - 1;
						sum1 += data[xx + yy * maxx];
					}
					// Corners
					for (int d = xoffset2.length; d-- > 0;)
					{
						int yy = y + yoffset2[d];
						int xx = x + xoffset2[d];
						if (xx <= 0)
							xx = 0;
						else if (xx >= maxx)
							xx = maxx - 1;
						if (yy <= 0)
							yy = 0;
						else if (yy >= maxy)
							yy = maxy - 1;
						sum2 += data[xx + yy * maxx];
					}
				}
				newData[index1] = (data[index1] + sum1 * w + sum2 * w2) * divisor;

				index0++;
				index1++;
				index2++;
			}
		}

		// Copy back
		System.arraycopy(newData, 0, data, 0, data.length);
	}

	/**
	 * Compute an approximate Gaussian convolution within a 3x3 size block around each point.
	 * <p>
	 * Uses a normalised [[1, 2, 1], [2, 4, 2], [1, 2, 1]] convolution kernel.
	 * <p>
	 * Note: the input data is destructively modified.
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 */
	public void blockGaussian3x3(float[] data, final int maxx, final int maxy)
	{
		float[] newData = floatBuffer(data.length);

		// Boundary control
		final int xwidth = 1;
		final int ywidth = 1;
		final int xlimit = maxx - xwidth;
		final int ylimit = maxy - ywidth;

		int[] offset = new int[(2 * xwidth + 1) * (2 * ywidth + 1) - 1];
		int[] xoffset = new int[offset.length];
		int[] yoffset = new int[offset.length];
		for (int y = -ywidth, d = 0; y <= ywidth; y++)
			for (int x = -xwidth; x <= xwidth; x++)
				if (x != 0 || y != 0)
				{
					offset[d] = maxx * y + x;
					xoffset[d] = x;
					yoffset[d] = y;
					d++;
				}

		float[] kernel = new float[] { 1f / 16, 2f / 16, 1f / 16, 2f / 16, /* 4f / 16, */2f / 16, 1f / 16, 2f / 16,
				1f / 16 };
		final float divisor = (float) (1.0 / 16.0);

		for (int y = 0; y < maxy; y++)
		{
			int index0 = (y - 1) * maxx;
			int index1 = y * maxx;
			int index2 = (y + 1) * maxx;
			for (int x = 0; x < maxx; x++)
			{
				// Flag to indicate this pixels has a complete (2n+1) neighbourhood 
				boolean isInnerXY = (y >= ywidth && y < ylimit) && (x >= xwidth && x < xlimit);

				// Sweep neighbourhood
				if (isInnerXY)
				{
					float sum = data[index0 - 1] + 2 * data[index0] + data[index0 + 1] + 2 * data[index1 - 1] + 4 *
							data[index1] + 2 * data[index1 + 1] + data[index2 - 1] + 2 * data[index2] +
							data[index2 + 1];
					newData[index1] = sum * divisor;
				}
				else
				{
					float sum = data[index1] * 0.25f; // = 4 / 16

					for (int d = offset.length; d-- > 0;)
					{
						// Get the pixel with boundary checking
						int yy = y + yoffset[d];
						int xx = x + xoffset[d];
						if (xx <= 0)
							xx = 0;
						else if (xx >= maxx)
							xx = maxx - 1;
						if (yy <= 0)
							yy = 0;
						else if (yy >= maxy)
							yy = maxy - 1;
						sum += data[xx + yy * maxx] * kernel[d];
					}
					newData[index1] = sum; // Kernel already normalised
				}
				index0++;
				index1++;
				index2++;
			}
		}

		// Copy back
		for (int index = data.length; index-- > 0;)
		{
			data[index] = newData[index];
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public Object clone()
	{
		try
		{
			AverageFilter o = (AverageFilter) super.clone();
			o.floatDataBuffer = null;
			o.floatRowBuffer = null;
			return o;
		}
		catch (CloneNotSupportedException e)
		{
			// Ignore
		}
		return null;
	}
}