package gdsc.smlm.filters;

import org.apache.commons.math3.util.FastMath;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Computes the block filter for each point within the array.
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
public abstract class BlockFilter extends BaseWeightedFilter
{
	private float[] floatDataBuffer = null;
	private float[] floatRowBuffer = null;

	private Normaliser normaliser;
	private Normaliser weightedNormaliser = null;
	private float weightedNormaliserN = 0;

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.BaseWeightedFilter#newWeights()
	 */
	@Override
	protected void newWeights()
	{
		weightedNormaliser = null;
	}

	/**
	 * Updates the weighted normaliser within a 2n+1 size block around each point.
	 *
	 * @param n
	 *            The block size
	 */
	private void updateWeightedNormaliser(final float n)
	{
		// Cache the normaliser
		if (weightedNormaliser == null || weightedNormaliserN != n)
		{
			weightedNormaliserN = n;
			weightedNormaliser = computeWeightedNormaliser(n);
		}
		normaliser = weightedNormaliser;
	}

	/**
	 * Computes the weighted normaliser within a 2n+1 size block around each point.
	 *
	 * @param n
	 *            The block size
	 * @return the weighted normaliser
	 */
	protected abstract Normaliser computeWeightedNormaliser(final float n);

	/**
	 * Computes the normaliser within a 2n+1 size block around each point.
	 *
	 * @param n
	 *            The block size
	 * @return the weighted normaliser
	 */
	protected abstract Normaliser computeNormaliser(final float n);

	/**
	 * Compute the filter within a 2n+1 size block around each point.
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
	public void rollingBlockFilterInternal(float[] data, final int maxx, final int maxy, final int n)
	{
		if (n == 1)
			rollingBlockFilter3x3Internal(data, maxx, maxy);
		else
			rollingBlockFilterNxNInternal(data, maxx, maxy, n);
	}

	/**
	 * Compute the filter within a 2n+1 size block around each point.
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
	void rollingBlockFilterNxNInternal(float[] data, final int maxx, final int maxy, final int n)
	{
		final int blockSize = 2 * n + 1;
		if (maxx < blockSize || maxy < blockSize)
			return;

		float[] newData = initialise(data, maxx, maxy, n);

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

			data[centreIndex] = normaliser.normalise(sum, centreIndex);

			while (y < maxy)
			{
				centreIndex += maxx;

				sum += newData[endIndex] - newData[startIndex];

				data[centreIndex] = normaliser.normalise(sum, centreIndex);

				y++;
				startIndex += maxx;
				endIndex += maxx;
			}
		}
	}

	/**
	 * Compute the filter within a 3x3 size block around each point.
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
	void rollingBlockFilter3x3Internal(float[] data, final int maxx, final int maxy)
	{
		if (maxx < 3 || maxy < 3)
			return;

		float[] newData = initialise(data, maxx, maxy, 1);

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
			data[centreIndex] = normaliser.normalise(sum, centreIndex);

			for (int y = 0; y < maxy - 3; y++)
			{
				centreIndex += maxx;
				endIndex += maxx;
				sum += newData[endIndex] - newData[startIndex];
				data[centreIndex] = normaliser.normalise(sum, centreIndex);
				startIndex += maxx;
			}
		}
	}

	/**
	 * Compute the filter within a 2n+1 size block around each point.
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
	public void stripedBlockFilterInternal(float[] data, final int maxx, final int maxy, final int n)
	{
		if (n == 1)
			stripedBlockFilter3x3Internal(data, maxx, maxy);
		else if (n == 2)
			stripedBlockFilter5x5Internal(data, maxx, maxy);
		else if (n == 3)
			stripedBlockFilter7x7Internal(data, maxx, maxy);
		else
			stripedBlockFilterNxNInternal(data, maxx, maxy, n);
	}

	/**
	 * Compute the filter within a 2w+1 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Uses a [[w*w, w, ..., w, w*w], [w, 1, ..., 1, w], [w*w, w, ..., w, w*w]] convolution kernel.
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
	public void stripedBlockFilterInternal(float[] data, final int maxx, final int maxy, final float w)
	{
		if (w <= 1)
			stripedBlockFilter3x3Internal(data, maxx, maxy, w);
		else if (w <= 2)
			stripedBlockFilter5x5Internal(data, maxx, maxy, w);
		else if (w <= 3)
			stripedBlockFilter7x7Internal(data, maxx, maxy, w);
		else
			stripedBlockFilterNxNInternal(data, maxx, maxy, w);
	}

	/**
	 * Compute the filter within a 2n+1 size block around each point.
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
	void stripedBlockFilterNxNInternal(float[] data, final int maxx, final int maxy, final int n)
	{
		final int blockSize = 2 * n + 1;
		if (maxx < blockSize || maxy < blockSize)
			return;

		float[] newData = initialise(data, maxx, maxy, n);

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
				int index2 = x + (y + n) * maxx;
				data[index2] = normaliser.normalise(sum, index2);
			}
		}
	}

	/**
	 * Compute the filter within a 2w+1 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Uses a [[w*w, w, ..., w, w*w], [w, 1, ..., 1, w], [w*w, w, ..., w, w*w]] convolution kernel.
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
	void stripedBlockFilterNxNInternal(float[] data, final int maxx, final int maxy, final float w)
	{
		final int n = (int) w;
		final int n1 = (n == w) ? n : n + 1;

		if (n == n1)
		{
			// There is no edge
			stripedBlockFilterInternal(data, maxx, maxy, n);
			return;
		}

		// The size of the region 
		final int nX = (2 * n1 + 1);
		final int nY = (2 * n1 + 1);

		if (maxx < nX || maxy < nY)
			return;

		final int blockSize = 2 * n1;

		float[] newData = initialise(data, maxx, maxy, w);

		final float w1 = w - n;

		// NOTE: 
		// To increase speed when sweeping the arrays:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		for (int y = 0; y < maxy; y++)
		{
			int index = y * maxx;
			for (int x = 0; x < maxx - blockSize; x++, index++)
			{
				float sum = data[index] * w1;
				for (int x2 = 1; x2 < blockSize; x2++)
				{
					sum += data[index + x2];
				}
				sum += data[index + blockSize] * w1;
				newData[(x + n1) * maxy + y] = sum;
			}
		}

		// Y-direction. 
		// Only sweep over the interior
		for (int x = n1; x < maxx - n1; x++)
		{
			int index = x * maxy;
			for (int y = 0; y < maxy - blockSize; y++, index++)
			{
				float sum = newData[index] * w1;
				for (int y2 = 1; y2 < blockSize; y2++)
				{
					sum += newData[index + y2];
				}
				sum += newData[index + blockSize] * w1;
				int index2 = x + (y + n1) * maxx;
				data[index2] = normaliser.normalise(sum, index2);
			}
		}
	}

	/**
	 * Compute the filter within a 3x3 size block around each point.
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
	void stripedBlockFilter3x3Internal(float[] data, final int maxx, final int maxy)
	{
		if (maxx < 3 || maxy < 3)
			return;

		float[] newData = initialise(data, maxx, maxy, 1);

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
				newData[index2] = data[index] + data[index + 1] + data[index + 2];
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
				data[index2] = normaliser.normalise(newData[index] + newData[index + 1] + newData[index + 2], index2);
				index2 += maxx;
			}
		}
	}

	/**
	 * Compute the filter within a 3x3 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Uses a [[w*w, w, w*w], [w, 1, w], [w*w, w, w*w]] convolution kernel.
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
	void stripedBlockFilter3x3Internal(float[] data, final int maxx, final int maxy, final float w)
	{
		if (maxx < 3 || maxy < 3)
			return;

		float[] newData = initialise(data, maxx, maxy, w);

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
				newData[index2] = w * (data[index] + data[index + 2]) + data[index + 1];
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
				data[index2] = normaliser.normalise(w * (newData[index] + newData[index + 2]) + newData[index + 1],
						index2);
				index2 += maxx;
			}
		}
	}

	/**
	 * Compute the filter within a 5x5 size block around each point.
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
	void stripedBlockFilter5x5Internal(float[] data, final int maxx, final int maxy)
	{
		if (maxx < 5 || maxy < 5)
			return;

		float[] newData = initialise(data, maxx, maxy, 2);

		// NOTE: 
		// To increase speed when sweeping the arrays:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		for (int y = 0; y < maxy; y++)
		{
			int index = y * maxx;
			int index2 = 2 * maxy + y;
			for (int x = 0; x <= maxx - 5; x++, index++)
			{
				float sum = data[index] + data[index + 1] + data[index + 2] + data[index + 3] + data[index + 4];
				newData[index2] = sum;
				index2 += maxy;
			}
		}

		// Y-direction. 
		// Only sweep over the interior
		for (int x = 2; x < maxx - 2; x++)
		{
			int index = x * maxy;
			int index2 = x + 2 * maxx;
			for (int y = 0; y <= maxy - 5; y++, index++)
			{
				float sum = newData[index] + newData[index + 1] + newData[index + 2] + newData[index + 3] +
						newData[index + 4];
				data[index2] = normaliser.normalise(sum, index2);
				index2 += maxx;
			}
		}
	}

	/**
	 * Compute the filter within a 5x5 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Uses a [[w*w, w, ..., w, w*w], [w, 1, ..., 1, w], [w*w, w, ..., w, w*w]] convolution kernel.
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
	 *            The weight (should be between 1 and 2)
	 */
	void stripedBlockFilter5x5Internal(float[] data, final int maxx, final int maxy, final float w)
	{
		if (maxx < 5 || maxy < 5)
			return;

		float[] newData = initialise(data, maxx, maxy, w);

		final float w1 = (w < 2) ? w - (int) w : 1;

		// NOTE: 
		// To increase speed when sweeping the arrays:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		for (int y = 0; y < maxy; y++)
		{
			int index = y * maxx;
			int index2 = 2 * maxy + y;
			for (int x = 0; x <= maxx - 5; x++, index++)
			{
				newData[index2] = w1 * (data[index] + data[index + 4]) + data[index + 1] + data[index + 2] +
						data[index + 3];
				index2 += maxy;
			}
		}

		// Y-direction. 
		// Only sweep over the interior
		for (int x = 2; x < maxx - 2; x++)
		{
			int index = x * maxy;
			int index2 = x + 2 * maxx;
			for (int y = 0; y <= maxy - 5; y++, index++)
			{
				data[index2] = normaliser.normalise(w1 * (newData[index] + newData[index + 4]) + newData[index + 1] +
						newData[index + 2] + newData[index + 3], index2);
				index2 += maxx;
			}
		}
	}

	/**
	 * Compute the filter within a 7x7 size block around each point.
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
	void stripedBlockFilter7x7Internal(float[] data, final int maxx, final int maxy)
	{
		if (maxx < 7 || maxy < 7)
			return;

		float[] newData = initialise(data, maxx, maxy, 3);

		// NOTE: 
		// To increase speed when sweeping the arrays:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		for (int y = 0; y < maxy; y++)
		{
			int index = y * maxx;
			int index2 = 3 * maxy + y;
			for (int x = 0; x <= maxx - 7; x++, index++)
			{
				float sum = data[index] + data[index + 1] + data[index + 2] + data[index + 3] + data[index + 4] +
						data[index + 5] + data[index + 6];
				newData[index2] = sum;
				index2 += maxy;
			}
		}

		// Y-direction. 
		// Only sweep over the interior
		for (int x = 3; x < maxx - 3; x++)
		{
			int index = x * maxy;
			int index2 = x + 3 * maxx;
			for (int y = 0; y <= maxy - 7; y++, index++)
			{
				float sum = newData[index] + newData[index + 1] + newData[index + 2] + newData[index + 3] +
						newData[index + 4] + newData[index + 5] + newData[index + 6];
				data[index2] = normaliser.normalise(sum, index2);
				index2 += maxx;
			}
		}
	}

	/**
	 * Compute the filter within a 7x7 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Uses a [[w*w, w, ..., w, w*w], [w, 1, ..., 1, w], [w*w, w, ..., w, w*w]] convolution kernel.
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
	 *            The weight (should be between 2 and 3)
	 */
	void stripedBlockFilter7x7Internal(float[] data, final int maxx, final int maxy, final float w)
	{
		if (maxx < 7 || maxy < 7)
			return;

		float[] newData = initialise(data, maxx, maxy, w);

		final float w1 = (w < 3) ? w - (int) w : 1;

		// NOTE: 
		// To increase speed when sweeping the arrays:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		for (int y = 0; y < maxy; y++)
		{
			int index = y * maxx;
			int index2 = 3 * maxy + y;
			for (int x = 0; x <= maxx - 7; x++, index++)
			{
				newData[index2] = w1 * (data[index] + data[index + 6]) + data[index + 1] + data[index + 2] +
						data[index + 3] + data[index + 4] + data[index + 5];
				index2 += maxy;
			}
		}

		// Y-direction. 
		// Only sweep over the interior
		for (int x = 3; x < maxx - 3; x++)
		{
			int index = x * maxy;
			int index2 = x + 3 * maxx;
			for (int y = 0; y <= maxy - 7; y++, index++)
			{
				data[index2] = normaliser.normalise(w1 * (newData[index] + newData[index + 6]) + newData[index + 1] +
						newData[index + 2] + newData[index + 3] + newData[index + 4] + newData[index + 5], index2);
				index2 += maxx;
			}
		}
	}

	/**
	 * Compute the filter within a 2n+1 size block around each point.
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
	public void blockFilterInternal(float[] data, final int maxx, final int maxy, final int n)
	{
		if (n == 1)
			blockFilter3x3Internal(data, maxx, maxy);
		else
			blockFilterNxNInternal(data, maxx, maxy, n);
	}

	/**
	 * Compute the filter within a 2w+1 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Uses a [[w*w, w, ..., w, w*w], [w, 1, ..., 1, w], [w*w, w, ..., w, w*w]] convolution kernel.
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
	public void blockFilterInternal(float[] data, final int maxx, final int maxy, final float w)
	{
		if (w < 1)
			blockFilter3x3Internal(data, maxx, maxy, w);
		else
			blockFilterNxNInternal(data, maxx, maxy, w);
	}

	/**
	 * Compute the filter within a 2n+1 size block around each point.
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
	void blockFilterNxNInternal(float[] data, final int maxx, final int maxy, final int n)
	{
		// The size of the region 
		final int nX = (2 * n + 1);
		final int nY = (2 * n + 1);

		if (maxx < nX || maxy < nY)
			return;

		float[] newData = initialise(data, maxx, maxy, n);

		int[] offset = new int[nX * nY - 1];
		for (int y = -n, d = 0; y <= n; y++)
			for (int x = -n; x <= n; x++)
				if (x != 0 || y != 0)
				{
					offset[d] = maxx * y + x;
					d++;
				}

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

				newData[index] = normaliser.normalise(sum, index);
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
	 * Compute the filter within a 2w+1 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Uses a [[w*w, w, ..., w, w*w], [w, 1, ..., 1, w], [w*w, w, ..., w, w*w]] convolution kernel.
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
	void blockFilterNxNInternal(float[] data, final int maxx, final int maxy, final float w)
	{
		final int n = (int) w;
		final int n1 = (n == w) ? n : n + 1;

		if (n == n1)
		{
			// There is no edge
			blockFilterInternal(data, maxx, maxy, n);
			return;
		}

		// The size of the region 
		final int nX = (2 * n1 + 1);
		final int nY = (2 * n1 + 1);

		if (maxx < nX || maxy < nY)
			return;

		final float[] newData = initialise(data, maxx, maxy, w);

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

				newData[index] = normaliser.normalise(sum + sum1 * w1 + sum2 * w2, index);
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
	 * Compute the filter within a 3x3 size block around each point.
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
	void blockFilter3x3Internal(float[] data, final int maxx, final int maxy)
	{
		float[] newData = initialise(data, maxx, maxy, 1);

		for (int y = 1; y < maxy - 1; y++)
		{
			int index0 = (y - 1) * maxx + 1;
			int index1 = y * maxx + 1;
			int index2 = (y + 1) * maxx + 1;
			for (int x = 1; x < maxx - 1; x++)
			{
				float sum = data[index0 - 1] + data[index0] + data[index0 + 1] + data[index1 - 1] + data[index1] +
						data[index1 + 1] + data[index2 - 1] + data[index2] + data[index2 + 1];
				newData[index1] = normaliser.normalise(sum, index1);
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
	 * Compute the weighted filter within a 3x3 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Uses a [[w*w, w, w*w], [w, 1, w], [w*w, w, w*w]] convolution kernel.
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
	void blockFilter3x3Internal(float[] data, final int maxx, final int maxy, final float w)
	{
		float[] newData = initialise(data, maxx, maxy, w);

		final float w2 = w * w;

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

				newData[index1] = normaliser.normalise(data[index1] + sum1 * w + sum2 * w2, index1);
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
	 * Initialise for filtering.
	 *
	 * @param data
	 *            the data
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 * @return the float buffer
	 */
	private float[] initialise(float[] data, final int maxx, final int maxy, final float n)
	{
		int size = data.length;
		if (hasWeights())
		{
			if (weights.length != size || this.weightWidth != maxx || this.weightHeight != maxy)
				throw new IllegalStateException("Weights are not the correct size");
			updateWeightedNormaliser(n);
			// Apply weights
			for (int i = 0; i < size; i++)
				data[i] *= weights[i];
		}
		else
		{
			normaliser = computeNormaliser(n);
		}
		return createFloatBuffer(size);
	}

	/**
	 * Creates the float buffer.
	 *
	 * @param size
	 *            the size
	 * @return the float buffer
	 */
	private float[] createFloatBuffer(int size)
	{
		if (floatDataBuffer == null || floatDataBuffer.length < size)
		{
			floatDataBuffer = new float[size];
		}
		return floatDataBuffer;
	}

	/**
	 * Compute the filter within a 2n+1 size block around each point.
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
	public void rollingBlockFilter(float[] data, final int maxx, final int maxy, final int n)
	{
		if (n == 1)
			rollingBlockFilter3x3(data, maxx, maxy);
		else
			rollingBlockFilterNxN(data, maxx, maxy, n);
	}

	/**
	 * Compute the filter within a 2n+1 size block around each point.
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
	void rollingBlockFilterNxN(float[] data, final int maxx, final int maxy, final int n)
	{
		float[] newData = initialise(data, maxx, maxy, n);

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
			outData[centreIndex] = normaliser.normalise(sum, centreIndex);

			// Rolling sum over the X-direction
			for (int x = 0; x < width - 1; x++)
			{
				sum += row[endIndex++] - row[x];
				centreIndex += height;
				outData[centreIndex] = normaliser.normalise(sum, centreIndex);
			}
		}
	}

	/**
	 * Compute the filter within a 3x3 size block around each point.
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
	void rollingBlockFilter3x3(float[] data, final int maxx, final int maxy)
	{
		float[] newData = initialise(data, maxx, maxy, 1);

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
			outData[centreIndex] = normaliser.normalise(sum, centreIndex);

			// Rolling sum over the X-direction
			for (int x = 0; x < width - 1; x++)
			{
				sum += row[endIndex++] - row[x];
				centreIndex += height;
				outData[centreIndex] = normaliser.normalise(sum, centreIndex);
			}
		}
	}

	/**
	 * Compute the filter within a 2n+1 size block around each point.
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
	public void stripedBlockFilter(float[] data, final int maxx, final int maxy, final int n)
	{
		if (n == 1)
			stripedBlockFilter3x3(data, maxx, maxy);
		else if (n == 2)
			stripedBlockFilter5x5(data, maxx, maxy, n);
		else if (n == 3)
			stripedBlockFilter7x7(data, maxx, maxy, n);
		else
			stripedBlockFilterNxN(data, maxx, maxy, n);
	}

	/**
	 * Compute the filter within a 2w+1 size block around each point.
	 * <p>
	 * Uses a [[w*w, w, ..., w, w*w], [w, 1, ..., 1, w], [w*w, w, ..., w, w*w]] convolution kernel.
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
	public void stripedBlockFilter(float[] data, final int maxx, final int maxy, final float w)
	{
		if (w <= 1)
			stripedBlockFilter3x3(data, maxx, maxy, w);
		else if (w <= 2)
			stripedBlockFilter5x5(data, maxx, maxy, w);
		else if (w <= 3)
			stripedBlockFilter7x7(data, maxx, maxy, w);
		else
			stripedBlockFilterNxN(data, maxx, maxy, w);
	}

	/**
	 * Compute the filter within a 2n+1 size block around each point.
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
	void stripedBlockFilterNxN(float[] data, final int maxx, final int maxy, final int n)
	{
		int blockSize = 2 * n + 1;

		float[] newData = initialise(data, maxx, maxy, n);

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
				outData[centreIndex] = normaliser.normalise(sum, centreIndex);
				centreIndex += height;
			}
		}
	}

	/**
	 * Compute the filter within a 2w+1 size block around each point.
	 * <p>
	 * Uses a [[w*w, w, ..., w, w*w], [w, 1, ..., 1, w], [w*w, w, ..., w, w*w]] convolution kernel.
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
	void stripedBlockFilterNxN(float[] data, final int maxx, final int maxy, final float w)
	{
		final int n = (int) w;
		final int n1 = (n == w) ? n : n + 1;

		if (n == n1)
		{
			// There is no edge
			stripedBlockFilter(data, maxx, maxy, n);
			return;
		}

		int blockSize = 2 * n1;

		float[] newData = initialise(data, maxx, maxy, w);

		final float w1 = w - n;

		// NOTE: 
		// To increase speed when sweeping the arrays and allow for reusing code:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		int width = maxx;
		int height = maxy;
		float[] inData = data;
		float[] outData = newData;

		float[] row = floatRowBuffer(width + 2 * n1);
		for (int y = 0; y < height; y++)
		{
			extractRow(inData, y, width, n1, row);

			int centreIndex = y;
			for (int x = 0; x < width; x++)
			{
				// Sum strips
				float sum = row[x] * w1;
				for (int j = 1; j < blockSize; j++)
				{
					sum += row[x + j];
				}
				sum += row[x + blockSize] * w1;

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

		row = floatRowBuffer(width + 2 * n1);
		for (int y = 0; y < height; y++)
		{
			// Extract row (pad ends)
			extractRow(inData, y, width, n1, row);

			int centreIndex = y;
			for (int x = 0; x < width; x++)
			{
				// Sum strips
				float sum = row[x] * w1;
				for (int j = 1; j < blockSize; j++)
				{
					sum += row[x + j];
				}
				sum += row[x + blockSize] * w1;

				// Store result in transpose
				outData[centreIndex] = normaliser.normalise(sum, centreIndex);
				centreIndex += height;
			}
		}
	}

	/**
	 * Compute the filter within a 3x3 size block around each point.
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
	void stripedBlockFilter3x3(float[] data, final int maxx, final int maxy)
	{
		float[] newData = initialise(data, maxx, maxy, 1);

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
				outData[centreIndex] = normaliser.normalise(sum, centreIndex);
				centreIndex += height;
			}
		}
	}

	/**
	 * Compute the filter within a 3x3 size block around each point.
	 * <p>
	 * Uses a [[w*w, w, w*w], [w, 1, w], [w*w, w, w*w]] convolution kernel.
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
	void stripedBlockFilter3x3(float[] data, final int maxx, final int maxy, final float w)
	{
		float[] newData = initialise(data, maxx, maxy, w);

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
				// Store result in transpose
				outData[centreIndex] = w * (row[x] + row[x + 2]) + row[x + 1];
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
				// Store result in transpose
				outData[centreIndex] = normaliser.normalise(w * (row[x] + row[x + 2]) + row[x + 1], centreIndex);
				centreIndex += height;
			}
		}
	}

	/**
	 * Compute the filter within a 5x5 size block around each point.
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
	void stripedBlockFilter5x5(float[] data, final int maxx, final int maxy)
	{
		float[] newData = initialise(data, maxx, maxy, 2);

		// NOTE: 
		// To increase speed when sweeping the arrays and allow for reusing code:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		int width = maxx;
		int height = maxy;
		float[] inData = data;
		float[] outData = newData;

		float[] row = floatRowBuffer(width + 4);
		for (int y = 0; y < height; y++)
		{
			extractRow2(inData, y, width, row);

			int centreIndex = y;
			for (int x = 0; x < width; x++)
			{
				// Sum strips
				// Store result in transpose
				outData[centreIndex] = row[x] + row[x + 1] + row[x + 2] + row[x + 3] + row[x + 4];
				centreIndex += height;
			}
		}

		// Y-direction. 
		width = maxy;
		height = maxx;
		inData = newData;
		outData = data;

		row = floatRowBuffer(width + 4);
		for (int y = 0; y < height; y++)
		{
			// Extract row (pad ends)
			extractRow2(inData, y, width, row);

			int centreIndex = y;
			for (int x = 0; x < width; x++)
			{
				// Sum strips
				// Store result in transpose
				outData[centreIndex] = normaliser.normalise(row[x] + row[x + 1] + row[x + 2] + row[x + 3] + row[x + 4],
						centreIndex);
				centreIndex += height;
			}
		}
	}

	/**
	 * Compute the filter within a 5x5 size block around each point.
	 * <p>
	 * Uses a [[w*w, w, ..., w, w*w], [w, 1, ..., 1, w], [w*w, w, ..., w, w*w]] convolution kernel.
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
	 *            The weight (should be between 1 and 2)
	 */
	void stripedBlockFilter5x5(float[] data, final int maxx, final int maxy, final float w)
	{
		float[] newData = initialise(data, maxx, maxy, w);

		final float w1 = (w < 2) ? w - (int) w : 1;

		// NOTE: 
		// To increase speed when sweeping the arrays and allow for reusing code:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		int width = maxx;
		int height = maxy;
		float[] inData = data;
		float[] outData = newData;

		float[] row = floatRowBuffer(width + 4);
		for (int y = 0; y < height; y++)
		{
			extractRow2(inData, y, width, row);

			int centreIndex = y;
			for (int x = 0; x < width; x++)
			{
				// Sum strips
				// Store result in transpose
				outData[centreIndex] = w1 * (row[x] + row[x + 4]) + row[x + 1] + row[x + 2] + row[x + 3];
				centreIndex += height;
			}
		}

		// Y-direction. 
		width = maxy;
		height = maxx;
		inData = newData;
		outData = data;

		row = floatRowBuffer(width + 4);
		for (int y = 0; y < height; y++)
		{
			// Extract row (pad ends)
			extractRow2(inData, y, width, row);

			int centreIndex = y;
			for (int x = 0; x < width; x++)
			{
				// Sum strips
				// Store result in transpose
				outData[centreIndex] = normaliser
						.normalise(w1 * (row[x] + row[x + 4]) + row[x + 1] + row[x + 2] + row[x + 3], centreIndex);
				centreIndex += height;
			}
		}
	}

	/**
	 * Compute the filter within a 7x7 size block around each point.
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
	void stripedBlockFilter7x7(float[] data, final int maxx, final int maxy)
	{
		float[] newData = initialise(data, maxx, maxy, 3);

		// NOTE: 
		// To increase speed when sweeping the arrays and allow for reusing code:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		int width = maxx;
		int height = maxy;
		float[] inData = data;
		float[] outData = newData;

		float[] row = floatRowBuffer(width + 6);
		for (int y = 0; y < height; y++)
		{
			extractRow3(inData, y, width, row);

			int centreIndex = y;
			for (int x = 0; x < width; x++)
			{
				// Sum strips
				// Store result in transpose
				outData[centreIndex] = row[x] + row[x + 1] + row[x + 2] + row[x + 3] + row[x + 4] + row[x + 5] +
						row[x + 6];
				centreIndex += height;
			}
		}

		// Y-direction. 
		width = maxy;
		height = maxx;
		inData = newData;
		outData = data;

		row = floatRowBuffer(width + 6);
		for (int y = 0; y < height; y++)
		{
			// Extract row (pad ends)
			extractRow3(inData, y, width, row);

			int centreIndex = y;
			for (int x = 0; x < width; x++)
			{
				// Sum strips
				// Store result in transpose
				outData[centreIndex] = normaliser.normalise(
						row[x] + row[x + 1] + row[x + 2] + row[x + 3] + row[x + 4] + row[x + 5] + row[x + 6],
						centreIndex);
				centreIndex += height;
			}
		}
	}

	/**
	 * Compute the filter within a 7x7 size block around each point.
	 * <p>
	 * Uses a [[w*w, w, ..., w, w*w], [w, 1, ..., 1, w], [w*w, w, ..., w, w*w]] convolution kernel.
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
	 *            The weight (should be between 2 and 3)
	 */
	void stripedBlockFilter7x7(float[] data, final int maxx, final int maxy, final float w)
	{
		float[] newData = initialise(data, maxx, maxy, w);

		final float w1 = (w < 3) ? w - (int) w : 1;

		// NOTE: 
		// To increase speed when sweeping the arrays and allow for reusing code:
		//   newData is XY ordinal => x * maxy + y
		//   data is YX ordinal    => y * maxx + x

		// X-direction
		int width = maxx;
		int height = maxy;
		float[] inData = data;
		float[] outData = newData;

		float[] row = floatRowBuffer(width + 6);
		for (int y = 0; y < height; y++)
		{
			extractRow3(inData, y, width, row);

			int centreIndex = y;
			for (int x = 0; x < width; x++)
			{
				// Sum strips
				// Store result in transpose
				outData[centreIndex] = w1 * (row[x] + row[x + 6]) + row[x + 1] + row[x + 2] + row[x + 3] + row[x + 4] +
						row[x + 5];
				centreIndex += height;
			}
		}

		// Y-direction. 
		width = maxy;
		height = maxx;
		inData = newData;
		outData = data;

		row = floatRowBuffer(width + 6);
		for (int y = 0; y < height; y++)
		{
			// Extract row (pad ends)
			extractRow3(inData, y, width, row);

			int centreIndex = y;
			for (int x = 0; x < width; x++)
			{
				// Sum strips
				// Store result in transpose
				outData[centreIndex] = normaliser.normalise(
						w1 * (row[x] + row[x + 6]) + row[x + 1] + row[x + 2] + row[x + 3] + row[x + 4] + row[x + 5],
						centreIndex);
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
	}

	private void extractRow1(float[] inData, int y, int width, float[] row)
	{
		final int index = y * width;

		// Pad ends
		row[0] = inData[index];
		row[1 + width] = inData[index + width - 1];

		// Fill in data
		System.arraycopy(inData, index, row, 1, width);
	}

	private void extractRow2(float[] inData, int y, int width, float[] row)
	{
		final int index = y * width;

		// Pad ends
		row[0] = row[1] = inData[index];
		row[2 + width] = row[3 + width] = inData[index + width - 1];

		// Fill in data
		System.arraycopy(inData, index, row, 2, width);
	}

	private void extractRow3(float[] inData, int y, int width, float[] row)
	{
		final int index = y * width;

		// Pad ends
		row[0] = row[1] = row[2] = inData[index];
		row[3 + width] = row[4 + width] = row[5 + width] = inData[index + width - 1];

		// Fill in data
		System.arraycopy(inData, index, row, 3, width);
	}

	/**
	 * Compute the filter within a 2n+1 size block around each point.
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
	public void blockFilter(float[] data, final int maxx, final int maxy, final int n)
	{
		if (n == 1)
			blockFilter3x3(data, maxx, maxy);
		else
			blockFilterNxN(data, maxx, maxy, n);
	}

	/**
	 * Compute the filter within a 2w+1 size block around each point.
	 * <p>
	 * Uses a [[w*w, w, ..., w, w*w], [w, 1, ..., 1, w], [w*w, w, ..., w, w*w]] convolution kernel.
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
	public void blockFilter(float[] data, final int maxx, final int maxy, final float w)
	{
		if (w < 1)
			blockFilter3x3(data, maxx, maxy, w);
		else
			blockFilterNxN(data, maxx, maxy, w);
	}

	/**
	 * Compute the filter within a 2n+1 size block around each point.
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
	void blockFilterNxN(float[] data, final int maxx, final int maxy, final int n)
	{
		float[] newData = initialise(data, maxx, maxy, n);

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
				newData[index] = normaliser.normalise(sum, index);
			}
		}

		// Copy back
		System.arraycopy(newData, 0, data, 0, data.length);
	}

	/**
	 * Compute the filter within a 2w+1 size block around each point.
	 * <p>
	 * Uses a [[w*w, w, ..., w, w*w], [w, 1, ..., 1, w], [w*w, w, ..., w, w*w]] convolution kernel.
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
	void blockFilterNxN(float[] data, final int maxx, final int maxy, final float w)
	{
		final int n = (int) w;
		final int n1 = (n == w) ? n : n + 1;

		if (n == n1)
		{
			// There is no edge
			blockFilter(data, maxx, maxy, n);
			return;
		}

		final float[] newData = initialise(data, maxx, maxy, w);

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
				newData[index] = normaliser.normalise(sum + sum1 * w1 + sum2 * w2, index);
			}
		}

		// Copy back
		System.arraycopy(newData, 0, data, 0, data.length);
	}

	/**
	 * Compute the filter within a 3x3 size block around each point.
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
	void blockFilter3x3(float[] data, final int maxx, final int maxy)
	{
		float[] newData = initialise(data, maxx, maxy, 1);

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
					newData[index1] = normaliser.normalise(sum, index1);
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
					newData[index1] = normaliser.normalise(sum, index1);
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
	 * Compute the weighted filter within a 3x3 size block around each point.
	 * <p>
	 * Uses a [[w*w, w, w*w], [w, 1, w], [w*w, w, w*w]] convolution kernel.
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
	void blockFilter3x3(float[] data, final int maxx, final int maxy, final float w)
	{
		float[] newData = initialise(data, maxx, maxy, w);

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
				newData[index1] = normaliser.normalise(data[index1] + sum1 * w + sum2 * w2, index1);

				index0++;
				index1++;
				index2++;
			}
		}

		// Copy back
		System.arraycopy(newData, 0, data, 0, data.length);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public BlockFilter clone()
	{
		BlockFilter o = (BlockFilter) super.clone();
		o.floatDataBuffer = null;
		o.floatRowBuffer = null;
		return o;
	}
}