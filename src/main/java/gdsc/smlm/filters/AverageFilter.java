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
	private int[] intDataBuffer = null;
	private int[] intRowBuffer = null;

	private float floatGetAverage(float sum, final float divisor)
	{
		return sum * divisor;
	}

	// When copying code changes must be made for intGetAverage
	//		int boxArea = (2 * n + 1) * (2 * n + 1);
	//      int halfBoxArea = (boxArea + 1) / 2;
	//      Since: average = (sum + halfBoxArea) / boxArea;
	//      //Initialise sum = halfBoxArea
	//      data[centreIndex] = intGetAverage(sum, boxArea);

	private int intGetAverage(int sum, int boxArea)
	{
		return sum / boxArea;
	}

	private int intGetAverage(float sum, final float divisor)
	{
		return (int) (sum * divisor);
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

		final float divisor = (float) (1.0 / ((2 * n + 1) * (2 * n + 1)));

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

			data[centreIndex] = floatGetAverage(sum, divisor);

			while (y < maxy)
			{
				centreIndex += maxx;

				sum += newData[endIndex] - newData[startIndex];

				data[centreIndex] = floatGetAverage(sum, divisor);

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
			data[centreIndex] = floatGetAverage(sum, divisor);

			for (int y = 0; y < maxy - 3; y++)
			{
				centreIndex += maxx;
				endIndex += maxx;
				sum += newData[endIndex] - newData[startIndex];
				data[centreIndex] = floatGetAverage(sum, divisor);
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

			data[centreIndex] = floatGetAverage(sum, divisor);

			while (y < maxy)
			{
				centreIndex += maxx;

				sum += newData[endIndex] - newData[startIndex];

				data[centreIndex] = floatGetAverage(sum, divisor);

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
				data[x + (y + n) * maxx] = floatGetAverage(sum, divisor);
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
				data[index2] = floatGetAverage(sum, divisor);
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
		float[] newData = floatBuffer(data.length);

		// Boundary control
		final int xwidth = FastMath.min(n, maxx - 1);
		final int ywidth = FastMath.min(n, maxy - 1);

		int[] offset = new int[(2 * xwidth + 1) * (2 * ywidth + 1) - 1];
		for (int y = -ywidth, d = 0; y <= ywidth; y++)
			for (int x = -xwidth; x <= xwidth; x++)
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

				newData[index] = floatGetAverage(sum, divisor);
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
				newData[index1] = floatGetAverage(sum, divisor);
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
	 * Uses a normalised [[w, w, w], [w, 1, w], [w, w, w]] convolution kernel.
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

		final float divisor = (float) (1.0 / (1 + 8 * w));

		for (int y = 1; y < maxy - 1; y++)
		{
			int index0 = (y - 1) * maxx + 1;
			int index1 = y * maxx + 1;
			int index2 = (y + 1) * maxx + 1;
			for (int x = 1; x < maxx - 1; x++)
			{
				float sum = data[index0 - 1] + data[index0] + data[index0 + 1] + data[index1 - 1] + data[index1 + 1] +
						data[index2 - 1] + data[index2] + data[index2 + 1];
				sum *= w;
				sum += data[index1];
				newData[index1] = floatGetAverage(sum, divisor);
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
				newData[index1] = floatGetAverage(sum, divisor);
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
			outData[centreIndex] = floatGetAverage(sum, divisor);

			// Rolling sum over the X-direction
			for (int x = 0; x < width - 1; x++)
			{
				sum += row[endIndex++] - row[x];
				centreIndex += height;
				outData[centreIndex] = floatGetAverage(sum, divisor);
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
			outData[centreIndex] = floatGetAverage(sum, divisor);

			// Rolling sum over the X-direction
			for (int x = 0; x < width - 1; x++)
			{
				sum += row[endIndex++] - row[x];
				centreIndex += height;
				outData[centreIndex] = floatGetAverage(sum, divisor);
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
				outData[centreIndex] = floatGetAverage(sum, divisor);
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
				outData[centreIndex] = floatGetAverage(sum, divisor);
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
		int index = y * width;

		// Pad ends
		for (int i = 0; i < n; i++)
		{
			row[i] = inData[index];
			row[i + n + width] = inData[index + width - 1];
		}

		// Fill in data
		for (int x = 0, i = n; x < width; x++)
		{
			row[i++] = inData[index++];
		}
	}

	private void extractRow1(float[] inData, int y, int width, float[] row)
	{
		int index = y * width;

		// Pad ends
		row[0] = inData[index];
		row[1 + width] = inData[index + width - 1];

		// Fill in data
		for (int x = 0, i = 1; x < width; x++)
		{
			row[i++] = inData[index++];
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
	public void blockAverage(float[] data, final int maxx, final int maxy, final int n)
	{
		if (n == 1)
			blockAverage3x3(data, maxx, maxy);
		else
			blockAverageNxN(data, maxx, maxy, n);
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
				newData[index] = floatGetAverage(sum, divisor);
			}
		}

		// Copy back
		for (index = data.length; index-- > 0;)
		{
			data[index] = newData[index];
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
					newData[index1] = floatGetAverage(sum, divisor);
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
					newData[index1] = floatGetAverage(sum, divisor);
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

	/**
	 * Compute the weighted block average within a 3x3 size block around each point.
	 * <p>
	 * Uses a normalised [[w, w, w], [w, 1, w], [w, w, w]] convolution kernel.
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

		final float divisor = (float) (1.0 / (1 + 8 * w));

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
				float sum = 0;
				if (isInnerXY)
				{
					sum = data[index0 - 1] + data[index0] + data[index0 + 1] + data[index1 - 1] + data[index1 + 1] +
							data[index2 - 1] + data[index2] + data[index2 + 1];
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
				sum *= w;
				sum += data[index1];
				newData[index1] = floatGetAverage(sum, divisor);

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
					newData[index1] = floatGetAverage(sum, divisor);
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

	// ----------------------------------------------------
	// XXX
	// NOTE:
	// The following code is copied directly from above. 
	// All 'float' have been replaced with 'int'. Changes must then be made for intGetAverage(...)
	// ----------------------------------------------------

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
	public void rollingBlockAverageInternal(int[] data, final int maxx, final int maxy, final int n)
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
	public void rollingBlockAverageNxNInternal(int[] data, final int maxx, final int maxy, final int n)
	{
		final int blockSize = 2 * n + 1;
		if (maxx < blockSize || maxy < blockSize)
			return;

		int[] newData = intBuffer(data.length);

		final int divisor = (2 * n + 1) * (2 * n + 1);

		// X-direction
		for (int y = 0; y < maxy; y++)
		{
			// Initialise the rolling sum
			int sum = 0;

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
			int sum = 0;

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

			data[centreIndex] = intGetAverage(sum, divisor);

			while (y < maxy)
			{
				centreIndex += maxx;

				sum += newData[endIndex] - newData[startIndex];

				data[centreIndex] = intGetAverage(sum, divisor);

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
	public void rollingBlockAverage3x3Internal(int[] data, final int maxx, final int maxy)
	{
		if (maxx < 3 || maxy < 3)
			return;

		int[] newData = intBuffer(data.length);

		final int divisor = 9;

		// X-direction
		for (int y = 0; y < maxy; y++)
		{
			// Initialise the rolling sum
			int startIndex = y * maxx;
			int centreIndex = startIndex + 1;
			int endIndex = centreIndex + 1;
			int sum = data[startIndex] + data[centreIndex] + data[endIndex];

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
			int sum = newData[startIndex] + newData[centreIndex] + newData[endIndex];

			// Rolling sum over the Y-direction
			data[centreIndex] = intGetAverage(sum, divisor);

			for (int y = 0; y < maxy - 3; y++)
			{
				centreIndex += maxx;
				endIndex += maxx;
				sum += newData[endIndex] - newData[startIndex];
				data[centreIndex] = intGetAverage(sum, divisor);
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
	public void rollingBlockAverageNxNInternalTransposed(int[] data, final int maxx, final int maxy, final int n)
	{
		int blockSize = 2 * n + 1;
		if (maxx < blockSize || maxy < blockSize)
			return;

		int[] newData = intBuffer(data.length);

		final int divisor = (2 * n + 1) * (2 * n + 1);

		// NOTE: 
		// To increase speed when sweeping the arrays:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		for (int y = 0; y < maxy; y++)
		{
			// Initialise the rolling sum
			int sum = 0;

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
			int sum = 0;

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

			data[centreIndex] = intGetAverage(sum, divisor);

			while (y < maxy)
			{
				centreIndex += maxx;

				sum += newData[endIndex] - newData[startIndex];

				data[centreIndex] = intGetAverage(sum, divisor);

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
	public void stripedBlockAverageInternal(int[] data, final int maxx, final int maxy, final int n)
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
	public void stripedBlockAverageNxNInternal(int[] data, final int maxx, final int maxy, final int n)
	{
		int blockSize = 2 * n + 1;
		if (maxx < blockSize || maxy < blockSize)
			return;

		int[] newData = intBuffer(data.length);

		final int divisor = (2 * n + 1) * (2 * n + 1);

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
				int sum = 0;
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
				int sum = 0;
				for (int y2 = 0; y2 < blockSize; y2++)
				{
					sum += newData[index + y2];
				}
				data[x + (y + n) * maxx] = intGetAverage(sum, divisor);
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
	public void stripedBlockAverage3x3Internal(int[] data, final int maxx, final int maxy)
	{
		if (maxx < 3 || maxy < 3)
			return;

		int[] newData = intBuffer(data.length);

		final int divisor = 9;

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
				int sum = data[index] + data[index + 1] + data[index + 2];
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
				int sum = newData[index] + newData[index + 1] + newData[index + 2];
				data[index2] = intGetAverage(sum, divisor);
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
	public void blockAverageInternal(int[] data, final int maxx, final int maxy, final int n)
	{
		if (n == 1)
			blockAverage3x3Internal(data, maxx, maxy);
		else
			blockAverageNxNInternal(data, maxx, maxy, n);
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
	public void blockAverageNxNInternal(int[] data, final int maxx, final int maxy, final int n)
	{
		int[] newData = intBuffer(data.length);

		// Boundary control
		final int xwidth = FastMath.min(n, maxx - 1);
		final int ywidth = FastMath.min(n, maxy - 1);

		int[] offset = new int[(2 * xwidth + 1) * (2 * ywidth + 1) - 1];
		for (int y = -ywidth, d = 0; y <= ywidth; y++)
			for (int x = -xwidth; x <= xwidth; x++)
				if (x != 0 || y != 0)
				{
					offset[d] = maxx * y + x;
					d++;
				}

		final int divisor = (2 * n + 1) * (2 * n + 1);

		for (int y = n; y < maxy - n; y++)
		{
			int index = y * maxx + n;
			for (int x = n; x < maxx - n; x++, index++)
			{
				int sum = data[index];

				// Sweep neighbourhood - 
				// No check for boundaries as this should be an internal sweep.
				for (int offset_d : offset)
				{
					sum += data[index + offset_d];
				}

				newData[index] = intGetAverage(sum, divisor);
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
	public void blockAverage3x3Internal(int[] data, final int maxx, final int maxy)
	{
		int[] newData = intBuffer(data.length);

		final int divisor = 9;

		for (int y = 1; y < maxy - 1; y++)
		{
			int index0 = (y - 1) * maxx + 1;
			int index1 = y * maxx + 1;
			int index2 = (y + 1) * maxx + 1;
			for (int x = 1; x < maxx - 1; x++)
			{
				int sum = data[index0 - 1] + data[index0] + data[index0 + 1] + data[index1 - 1] + data[index1] +
						data[index1 + 1] + data[index2 - 1] + data[index2] + data[index2 + 1];
				newData[index1] = intGetAverage(sum, divisor);
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
	 * Uses a normalised [[w, w, w], [w, 1, w], [w, w, w]] convolution kernel.
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
	public void blockAverage3x3Internal(int[] data, final int maxx, final int maxy, final int w)
	{
		int[] newData = intBuffer(data.length);

		final float divisor = (float) (1.0 / (1 + 8 * w));

		for (int y = 1; y < maxy - 1; y++)
		{
			int index0 = (y - 1) * maxx + 1;
			int index1 = y * maxx + 1;
			int index2 = (y + 1) * maxx + 1;
			for (int x = 1; x < maxx - 1; x++)
			{
				int sum = data[index0 - 1] + data[index0] + data[index0 + 1] + data[index1 - 1] + data[index1 + 1] +
						data[index2 - 1] + data[index2] + data[index2 + 1];
				float fsum = sum * w + data[index1];
				newData[index1] = intGetAverage(fsum, divisor);
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
	public void blockGaussian3x3Internal(int[] data, final int maxx, final int maxy)
	{
		int[] newData = intBuffer(data.length);

		final int divisor = 16;

		for (int y = 1; y < maxy - 1; y++)
		{
			int index0 = (y - 1) * maxx + 1;
			int index1 = y * maxx + 1;
			int index2 = (y + 1) * maxx + 1;
			for (int x = 1; x < maxx - 1; x++)
			{
				int sum = data[index0 - 1] + 2 * data[index0] + data[index0 + 1] + 2 * data[index1 - 1] + 4 *
						data[index1] + 2 * data[index1 + 1] + data[index2 - 1] + 2 * data[index2] + data[index2 + 1];
				newData[index1] = intGetAverage(sum, divisor);
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

	private int[] intBuffer(int size)
	{
		if (intDataBuffer == null || intDataBuffer.length < size)
		{
			intDataBuffer = new int[size];
		}
		return intDataBuffer;
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
	public void rollingBlockAverage(int[] data, final int maxx, final int maxy, final int n)
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
	public void rollingBlockAverageNxN(int[] data, final int maxx, final int maxy, final int n)
	{
		int[] newData = intBuffer(data.length);

		final int divisor = (2 * n + 1) * (2 * n + 1);

		// NOTE: 
		// To increase speed when sweeping the arrays and allow for reusing code:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		int width = maxx;
		int height = maxy;
		int[] inData = data;
		int[] outData = newData;

		int[] row = intRowBuffer(width + 2 * n);
		for (int y = 0; y < height; y++)
		{
			extractRow(inData, y, width, n, row);

			// Initialise rolling sum
			int sum = (n + 1) * row[0];
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

		row = intRowBuffer(width + 2 * n);
		for (int y = 0; y < height; y++)
		{
			extractRow(inData, y, width, n, row);

			// Initialise rolling sum
			int sum = (n + 1) * row[0];
			int endIndex = n + 1;
			for (int i = 0; i < n; i++)
			{
				sum += row[endIndex++];
			}

			int centreIndex = y;
			outData[centreIndex] = intGetAverage(sum, divisor);

			// Rolling sum over the X-direction
			for (int x = 0; x < width - 1; x++)
			{
				sum += row[endIndex++] - row[x];
				centreIndex += height;
				outData[centreIndex] = intGetAverage(sum, divisor);
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
	public void rollingBlockAverage3x3(int[] data, final int maxx, final int maxy)
	{
		int[] newData = intBuffer(data.length);

		final int divisor = 9;

		// NOTE: 
		// To increase speed when sweeping the arrays and allow for reusing code:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		int width = maxx;
		int height = maxy;
		int[] inData = data;
		int[] outData = newData;

		int[] row = intRowBuffer(width + 2);
		for (int y = 0; y < height; y++)
		{
			extractRow1(inData, y, width, row);

			// Initialise rolling sum
			int sum = 2 * row[0] + row[2];
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

		row = intRowBuffer(width + 2);
		for (int y = 0; y < height; y++)
		{
			extractRow1(inData, y, width, row);

			// Initialise rolling sum
			int sum = 2 * row[0] + row[2];
			int endIndex = 3;

			int centreIndex = y;
			outData[centreIndex] = intGetAverage(sum, divisor);

			// Rolling sum over the X-direction
			for (int x = 0; x < width - 1; x++)
			{
				sum += row[endIndex++] - row[x];
				centreIndex += height;
				outData[centreIndex] = intGetAverage(sum, divisor);
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
	public void stripedBlockAverage(int[] data, final int maxx, final int maxy, final int n)
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
	public void stripedBlockAverageNxN(int[] data, final int maxx, final int maxy, final int n)
	{
		int blockSize = 2 * n + 1;

		int[] newData = intBuffer(data.length);

		final int divisor = (2 * n + 1) * (2 * n + 1);

		// NOTE: 
		// To increase speed when sweeping the arrays and allow for reusing code:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		int width = maxx;
		int height = maxy;
		int[] inData = data;
		int[] outData = newData;

		int[] row = intRowBuffer(width + 2 * n);
		for (int y = 0; y < height; y++)
		{
			extractRow(inData, y, width, n, row);

			int centreIndex = y;
			for (int x = 0; x < width; x++)
			{
				// Sum strips
				int sum = 0;

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

		row = intRowBuffer(width + 2 * n);
		for (int y = 0; y < height; y++)
		{
			// Extract row (pad ends)
			extractRow(inData, y, width, n, row);

			int centreIndex = y;
			for (int x = 0; x < width; x++)
			{
				// Sum strips
				int sum = 0;

				for (int j = 0; j < blockSize; j++)
				{
					sum += row[x + j];
				}

				// Store result in transpose
				outData[centreIndex] = intGetAverage(sum, divisor);
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
	public void stripedBlockAverage3x3(int[] data, final int maxx, final int maxy)
	{
		int[] newData = intBuffer(data.length);

		final int divisor = 9;

		// NOTE: 
		// To increase speed when sweeping the arrays and allow for reusing code:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		int width = maxx;
		int height = maxy;
		int[] inData = data;
		int[] outData = newData;

		int[] row = intRowBuffer(width + 2);
		for (int y = 0; y < height; y++)
		{
			extractRow1(inData, y, width, row);

			int centreIndex = y;
			for (int x = 0; x < width; x++)
			{
				// Sum strips
				int sum = row[x] + row[x + 1] + row[x + 2];

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

		row = intRowBuffer(width + 2);
		for (int y = 0; y < height; y++)
		{
			// Extract row (pad ends)
			extractRow1(inData, y, width, row);

			int centreIndex = y;
			for (int x = 0; x < width; x++)
			{
				// Sum strips
				int sum = row[x] + row[x + 1] + row[x + 2];

				// Store result in transpose
				outData[centreIndex] = intGetAverage(sum, divisor);
				centreIndex += height;
			}
		}
	}

	private int[] intRowBuffer(int size)
	{
		if (intRowBuffer == null || intRowBuffer.length < size)
		{
			intRowBuffer = new int[size];
		}
		return intRowBuffer;
	}

	private void extractRow(int[] inData, int y, int width, final int n, int[] row)
	{
		int index = y * width;

		// Pad ends
		for (int i = 0; i < n; i++)
		{
			row[i] = inData[index];
			row[i + n + width] = inData[index + width - 1];
		}

		// Fill in data
		for (int x = 0, i = n; x < width; x++)
		{
			row[i++] = inData[index++];
		}
	}

	private void extractRow1(int[] inData, int y, int width, int[] row)
	{
		int index = y * width;

		// Pad ends
		row[0] = inData[index];
		row[1 + width] = inData[index + width - 1];

		// Fill in data
		for (int x = 0, i = 1; x < width; x++)
		{
			row[i++] = inData[index++];
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
	public void blockAverage(int[] data, final int maxx, final int maxy, final int n)
	{
		if (n == 1)
			blockAverage3x3(data, maxx, maxy);
		else
			blockAverageNxN(data, maxx, maxy, n);
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
	public void blockAverageNxN(int[] data, final int maxx, final int maxy, final int n)
	{
		int[] newData = intBuffer(data.length);

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

		final int divisor = (2 * n + 1) * (2 * n + 1);

		int index = 0;
		for (int y = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++, index++)
			{
				int sum = data[index];

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
				newData[index] = intGetAverage(sum, divisor);
			}
		}

		// Copy back
		for (index = data.length; index-- > 0;)
		{
			data[index] = newData[index];
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
	public void blockAverage3x3(int[] data, final int maxx, final int maxy)
	{
		int[] newData = intBuffer(data.length);

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

		final int divisor = 9;

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
					int sum = data[index0 - 1] + data[index0] + data[index0 + 1] + data[index1 - 1] + data[index1] +
							data[index1 + 1] + data[index2 - 1] + data[index2] + data[index2 + 1];
					newData[index1] = intGetAverage(sum, divisor);
				}
				else
				{
					int sum = data[index1];

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
					newData[index1] = intGetAverage(sum, divisor);
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

	/**
	 * Compute the weighted block average within a 3x3 size block around each point.
	 * <p>
	 * Uses a normalised [[w, w, w], [w, 1, w], [w, w, w]] convolution kernel.
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
	public void blockAverage3x3(int[] data, final int maxx, final int maxy, final float w)
	{
		int[] newData = intBuffer(data.length);

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

		final float divisor = (float) (1.0 / (1 + 8 * w));

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
				int sum = 0;
				if (isInnerXY)
				{
					sum = data[index0 - 1] + data[index0] + data[index0 + 1] + data[index1 - 1] + data[index1 + 1] +
							data[index2 - 1] + data[index2] + data[index2 + 1];
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
				float fsum = sum * w + data[index1];
				newData[index1] = intGetAverage(fsum, divisor);

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
	public void blockGaussian3x3(int[] data, final int maxx, final int maxy)
	{
		int[] newData = intBuffer(data.length);

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

		int[] kernel = new int[] { 1, 2, 1, 2, /* 4, */2, 1, 2, 1 };
		final int divisor = 16;

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
					int sum = data[index0 - 1] + 2 * data[index0] + data[index0 + 1] + 2 * data[index1 - 1] + 4 *
							data[index1] + 2 * data[index1 + 1] + data[index2 - 1] + 2 * data[index2] +
							data[index2 + 1];
					newData[index1] = intGetAverage(sum, divisor);
				}
				else
				{
					int sum = data[index1] * 4;

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
					newData[index1] = intGetAverage(sum, divisor);
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
			o.intDataBuffer = null;
			o.intRowBuffer = null;
			return o;
		}
		catch (CloneNotSupportedException e)
		{
			// Ignore
		}
		return null;
	}
}