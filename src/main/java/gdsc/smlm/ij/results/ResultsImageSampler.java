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
	private static class IndexComparator implements Comparator<long[]>
	{
		public int compare(long[] o1, long[] o2)
		{
			return Long.compare(o1[0], o2[0]);
		}
	}

	private static class CountComparator implements Comparator<long[]>
	{
		public int compare(long[] o1, long[] o2)
		{
			int result = Long.compare(o1[1], o2[1]);
			if (result == 0)
				// Use index if the same count
				return Long.compare(o1[0], o2[0]);
			return result;
		}
	}

	private static class ReverseCountComparator implements Comparator<long[]>
	{
		public int compare(long[] o1, long[] o2)
		{
			int result = Long.compare(o2[1], o1[1]);
			if (result == 0)
				// Use index if the same count
				return Long.compare(o1[0], o2[0]);
			return result;
		}
	}

	private static IndexComparator ic = new IndexComparator();
	private static CountComparator cc = new CountComparator();
	private static ReverseCountComparator rcc = new ReverseCountComparator();

	private final MemoryPeakResults results;
	private final ImageStack stack;

	private final int lx, ly;
	private final int xblocks, yblocks;
	private final int xy_blocks;

	/** The size for samples */
	public final int size;

	/**
	 * The max number of empty samples. Since they are empty it should not matter
	 * unless the noise characteristics change over the image duration. Set to 0 to sample throughout the lifetime of
	 * the localisation occurrences.
	 */
	public int maxNumberOfEmptySamples = 500;

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
		lx = size * (bounds.x / size);
		ly = size * (bounds.y / size);
		int ux = size * (int) Math.ceil((double) (bounds.x + bounds.width) / size);
		int uy = size * (int) Math.ceil((double) (bounds.y + bounds.height) / size);
		xblocks = (ux - lx) / size;
		yblocks = (uy - ly) / size;
		xy_blocks = xblocks * yblocks;
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
			// Avoid invalid slices
			if (p.peak < 1 || p.peak > stack.getSize())
				continue;
			long index = getIndex(p.getXPosition(), p.getYPosition(), p.peak);
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

		// Sort by index
		Arrays.sort(data, ic);

		// Do the empty blocks
		long total = stack.getSize() * xy_blocks;
		long empty = total - map.size();
		if (empty == 0)
		{
			// All indices are used 
			no = new long[0];
		}
		else
		{
			// Randomly sample indices that are not used.
			if (maxNumberOfEmptySamples > 0)
			{
				// Just enumerate the first N. Since they are empty it should not matter
				// unless the noise characteristics change over the image duration.
				long emptyCandidate = 0;
				long[] list = new long[maxNumberOfEmptySamples];
				int c = 0;
				OUTER: for (int i = 0; i < data.length; i++)
				{
					long current = data[i][0];
					// If the current index is bigger than the candidate then it must be empty
					while (current > emptyCandidate)
					{
						// Add all those that are empty
						list[c++] = emptyCandidate++;
						if (c == maxNumberOfEmptySamples)
							break OUTER;
					}
					// Set the next candidate
					emptyCandidate = current + 1;
				}
				no = Arrays.copyOf(list, c);
			}
			else
			{
				// Sample throughout the localisation time course
				TLongArrayList list = new TLongArrayList(data.length);
				if (empty < data.length)
				{
					// We can pick all the indices that are missing 
					long emptyCandidate = 0;
					for (int i = 0; i < data.length; i++)
					{
						long current = data[i][0];
						// If the current index is bigger than the candidate then it must be empty
						while (current > emptyCandidate)
						{
							// Add all those that are empty
							list.add(emptyCandidate++);
						}
						// Set the next candidate
						emptyCandidate = current + 1;
					}
				}
				else
				{
					// There are many empty blocks so just sample blocks 
					// after those with localisations.
					long emptyCandidate = 1;
					for (int i = 0; i < data.length; i++)
					{
						long current = data[i][0];
						// If the current index is bigger than the candidate then it must be empty
						if (current > emptyCandidate)
						{
							// Note: we only sample the next empty index after an index with data
							// This means the number of frames could be lower
							list.add(emptyCandidate);
						}
						// Set the next candidate
						emptyCandidate = current + 1;
					}
				}
				no = list.toArray();
			}
		}

		// For the low and high sample we just split in half.
		// The data must be sorted by count.
		Arrays.sort(data, cc);
		lower = data.length / 2;
		upper = data.length - lower;

		//// Debug indexing
		//int[] xyz = new int[3];
		//for (long l = 0; l < max; l++)
		//{
		//	getXYZ(l, xyz);
		//	if (xyz[0] < 0 || xyz[0] > stack.getWidth() || xyz[1] < 0 || xyz[1] > stack.getHeight() || xyz[2] < 1 ||
		//			xyz[2] > stack.getSize())
		//	{
		//		System.out.printf("Bad %d = %s\n", l, Arrays.toString(xyz));
		//	}
		//}

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
		// Make frames start at zero for the index
		return (long) ((xy_blocks) * (t - 1)) + xblocks * getY(y) + getX(x);
	}

	private int getX(float f)
	{
		return (int) ((f - lx) / size);
	}

	private int getY(float f)
	{
		return (int) ((f - ly) / size);
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
		xyz[2] = (int) (index / (xy_blocks));
		int mod = (int) (index % (xy_blocks));
		xyz[1] = mod / xblocks;
		xyz[0] = mod % xblocks;

		// Convert back to real coords
		xyz[2]++; // Frames start at 1
		xyz[1] = (xyz[1] * size) + ly;
		xyz[0] = (xyz[0] * size) + lx;

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
		Arrays.sort(sample, rcc);

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
				for (int x = 0, i = y * size, index = (y + target.y) * ip.getWidth() +
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
