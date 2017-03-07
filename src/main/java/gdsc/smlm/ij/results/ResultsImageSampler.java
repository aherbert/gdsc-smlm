package gdsc.smlm.ij.results;

import java.awt.Rectangle;
import java.util.Arrays;
import java.util.Comparator;

import gdsc.core.utils.Random;
import gdsc.core.utils.TurboList;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gnu.trove.list.array.TLongArrayList;
import gnu.trove.map.hash.TLongIntHashMap;
import gnu.trove.procedure.TLongIntProcedure;
import ij.ImageStack;
import ij.process.ImageProcessor;

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
 * Allows sampling sections from a source image for localisation results.
 */
public class ResultsImageSampler
{
	private static class LongComparator implements Comparator<long[]>
	{
		public int compare(long[] o1, long[] o2)
		{
			return Long.compare(o1[1], o2[1]);
		}
	}

	private static LongComparator lc = new LongComparator();

	private final MemoryPeakResults results;
	private final ImageStack stack;

	private final int minx, miny;
	private final int maxx, maxy;
	private final int maxx_maxy;

	/** The size for samples */
	public final int size;

	private long[] no;
	private long[][] data;
	private int lower, upper;
	private TurboList<long[]> list = new TurboList<long[]>();
	private Random r = new Random();

	/**
	 * Instantiates a new results image sampler.
	 *
	 * @param results
	 *            the results
	 * @param stack
	 *            the source stack for the results
	 * @param size
	 *            the size of the sample blocks
	 */
	public ResultsImageSampler(MemoryPeakResults results, ImageStack stack, int size)
	{
		this.results = results;
		this.stack = stack;
		this.size = size;

		Rectangle bounds = results.getBounds(true);
		// Round the image dimensions to the nearest block interval
		minx = size * (bounds.x / size);
		miny = size * (bounds.y / size);
		int ux = size * (int) Math.ceil((double) (bounds.x + bounds.width) / size);
		int uy = size * (int) Math.ceil((double) (bounds.y + bounds.height) / size);
		maxx = ux - minx;
		maxy = uy - miny;
		maxx_maxy = maxx * maxy;
	}

	/**
	 * Analyse the input data to allow sampling.
	 *
	 * @return true, if successful
	 */
	public boolean analyse()
	{
		no = null;
		if (results.isEmpty())
			return false;

		// Count the number of localisations in each region
		TLongIntHashMap map = new TLongIntHashMap(results.size());
		for (PeakResult p : results)
		{
			if (p.peak < 1)
				continue;
			long index = getIndex(p.getXPosition(), p.getYPosition(), p.peak - 1);
			map.adjustOrPutValue(index, 1, 1);
		}

		final long[][] data = new long[map.size()][];
		map.forEachEntry(new TLongIntProcedure()
		{
			int i = 0;

			public boolean execute(long a, int b)
			{
				data[i++] = new long[] { a, b };
				return true;
			}
		});
		this.data = data;

		// Split the image into frames with zero localisations, low density, high density

		// Sort
		Arrays.sort(data, lc);

		// For the low and high sample we just split in half. 
		lower = data.length / 2;
		upper = data.length - lower;

		// Do the empty blocks
		long max = 1 + data[data.length - 1][0];
		long empty = max - map.size();
		if (empty == 0)
		{
			// All indices are used 
			no = new long[0];
		}
		else
		{
			// Randomly sample indices that are not used.
			// Do this by picking blocks after those with localisations.
			TLongArrayList list = new TLongArrayList(data.length);
			long emptyCandidate = 1;
			for (int i = 0; i < data.length; i++)
			{
				long current = data[i][0];
				// If the current index is bigger than the candidate then it must be empty
				if (current > emptyCandidate)
					list.add(emptyCandidate);
				// Set the next candidate
				emptyCandidate = current + 1;
			}
			no = list.toArray();
		}

		return true;
	}

	/**
	 * Gets the index of the block from the floating point coordinates.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param t
	 *            the time frame
	 * @return the index
	 */
	private long getIndex(float x, float y, int t)
	{
		return (long) ((maxx_maxy) * t) + maxx * getY(y) + getX(x);
	}

	private int getX(float f)
	{
		return (int) ((f - minx) / size);
	}

	private int getY(float f)
	{
		return (int) ((f - miny) / size);
	}

	/**
	 * Convert the single index into x,y,z coords, Input array must be length >= 3.
	 * 
	 * @param index
	 * @param xyz
	 * @return The xyz array
	 */
	private int[] getXYZ(long index, int[] xyz)
	{
		// Frames start at 1
		xyz[2] = (int) (index / (maxx_maxy)) + 1;
		int mod = (int) (index % (maxx_maxy));
		xyz[1] = mod / maxx;
		xyz[0] = mod % maxx;
		// Convert back to real coords
		xyz[1] = (xyz[1] * size) + miny;
		xyz[0] = (xyz[0] * size) + minx;
		return xyz;
	}

	/**
	 * Gets the results used to construct this instance.
	 *
	 * @return the results
	 */
	public MemoryPeakResults getResults()
	{
		return results;
	}

	/**
	 * Checks if is valid (i.e. samples can be obtained).
	 *
	 * @return true, if is valid
	 */
	public boolean isValid()
	{
		return no != null;
	}

	/**
	 * Gets the number of empty samples.
	 *
	 * @return the number of empty samples
	 */
	public int getNumberOfEmptySamples()
	{
		return no.length;
	}

	/**
	 * Gets the number of low density samples.
	 *
	 * @return the number of low density samples
	 */
	public int getNumberOfLowDensitySamples()
	{
		return lower;
	}

	/**
	 * Gets the number of high density samples.
	 *
	 * @return the number of high density samples
	 */
	public int getNumberOfHighDensitySamples()
	{
		return upper;
	}

	/**
	 * Gets the sample.
	 *
	 * @param nNo
	 *            the number of samples with no localisations
	 * @param nLow
	 *            the number of samples with low localisations
	 * @param nHigh
	 *            the number of samples with high localisations
	 * @return the sample (could be empty if no samples were made)
	 */
	public ImageStack getSample(int nNo, int nLow, int nHigh)
	{
		ImageStack out = new ImageStack(size, size);
		if (!isValid())
			return out;

		list.clearf();

		// empty
		for (int i : r.sample(nNo, no.length))
			list.add(new long[] { no[i], 0 });
		// low
		for (int i : r.sample(nLow, lower))
			list.add(data[i]);
		// high
		for (int i : r.sample(nHigh, upper))
			list.add(data[i + lower]);

		if (list.isEmpty())
			return out;

		// Sort descending by number in the frame
		long[][] sample = list.toArray(new long[list.size()][]);
		Arrays.sort(sample, lc);

		int[] xyz = new int[3];
		Rectangle stackBounds = new Rectangle(stack.getWidth(), stack.getHeight());
		for (long[] s : sample)
		{
			getXYZ(s[0], xyz);

			// Construct the region to extract
			Rectangle target = new Rectangle(xyz[0], xyz[1], size, size);
			target = target.intersection(stackBounds);
			if (target.width == 0 || target.height == 0)
				continue;

			// Extract the frame
			int slice = xyz[2];
			ImageProcessor ip = stack.getProcessor(slice);

			// Cut out the desired pixels (some may be blank if the block overruns the source image)
			float[] pixels = new float[size * size];
			for (int y = 0; y < target.height; y++)
				for (int x = 0, i = y * size, index = target.y * ip.getWidth() +
						target.x; x < target.width; x++, i++, index++)
				{
					pixels[i] = ip.getf(index);
				}
			out.addSlice(String.format("Frame=%d (x=%d,y=%d) (n=%d)", slice, xyz[0], xyz[1], s[1]), pixels);
		}

		return out;
	}

	/**
	 * Set the random for use during sampling
	 * 
	 * @param r
	 *            the random to set (ignored if null)
	 */
	public void setRandom(Random r)
	{
		if (r != null)
			this.r = r;
	}
}
