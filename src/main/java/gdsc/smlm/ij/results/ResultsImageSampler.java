package gdsc.smlm.ij.results;

import java.awt.Rectangle;
import java.util.Arrays;
import java.util.Comparator;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

import gdsc.core.utils.Random;
import gdsc.core.utils.TurboList;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TLongArrayList;
import gnu.trove.map.hash.TLongObjectHashMap;
import gnu.trove.procedure.TLongObjectProcedure;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Overlay;
import ij.gui.PointRoi;
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
	private static class ResultsSample
	{
		final long index;
		final float[] data;

		ResultsSample(long index)
		{
			this.index = index;
			this.data = new float[0]; // makes things simple
		}

		ResultsSample(long index, float[] data)
		{
			this.index = index;
			this.data = data;
		}

		public int size()
		{
			return data.length / 2;
		}
	}

	private static class IndexComparator implements Comparator<ResultsSample>
	{
		public int compare(ResultsSample o1, ResultsSample o2)
		{
			return Long.compare(o1.index, o2.index);
		}
	}

	private static class CountComparator implements Comparator<ResultsSample>
	{
		public int compare(ResultsSample o1, ResultsSample o2)
		{
			int result = Integer.compare(o1.data.length, o2.data.length);
			if (result == 0)
				// Use index if the same count
				return Long.compare(o1.index, o2.index);
			return result;
		}
	}

	private static class ReverseCountComparator implements Comparator<ResultsSample>
	{
		public int compare(ResultsSample o1, ResultsSample o2)
		{
			int result = Integer.compare(o2.data.length, o1.data.length);
			if (result == 0)
				// Use index if the same count
				return Long.compare(o1.index, o2.index);
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
	private ResultsSample[] data;
	private int lower, upper;
	private TurboList<ResultsSample> list = new TurboList<ResultsSample>();
	private RandomGenerator r = new Well19937c();

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
		// Clear old samples
		data = null;
		no = null;

		if (results.isEmpty())
			return false;

		createResultSamples();

		// Split the image into frames with zero localisations, low density, high density

		createEmptySampleIndices();

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

	private void createEmptySampleIndices()
	{
		// Create the empty blocks
		long total = stack.getSize() * xy_blocks;
		long empty = total - data.length;
		if (empty == 0)
		{
			// All indices are used 
			no = new long[0];
		}
		else
		{
			// Sort by index
			Arrays.sort(data, ic);

			// Randomly sample indices that are not used.
			if (maxNumberOfEmptySamples > 0)
			{
				// Just enumerate the first N. Since they are empty it should not matter
				// unless the noise characteristics change over the image duration.
				long emptyCandidate = 0;
				long[] list = new long[(int) Math.min(empty, maxNumberOfEmptySamples)];
				int c = 0;
				OUTER: for (int i = 0; i < data.length; i++)
				{
					long current = data[i].index;
					// If the current index is bigger than the candidate then it must be empty
					while (current > emptyCandidate)
					{
						// Add all those that are empty
						list[c++] = emptyCandidate++;
						if (c == list.length)
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
						long current = data[i].index;
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
						long current = data[i].index;
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
	}

	/**
	 * Creates the result samples.
	 * Do this by storing the coordinates at the region index.
	 */
	private void createResultSamples()
	{
		TLongObjectHashMap<TFloatArrayList> map = new TLongObjectHashMap<TFloatArrayList>(results.size());
		TFloatArrayList next = new TFloatArrayList(2);
		for (PeakResult p : results)
		{
			// Avoid invalid slices
			if (p.getFrame() < 1 || p.getFrame() > stack.getSize())
				continue;
			long index = getIndex(p.getXPosition(), p.getYPosition(), p.getFrame());

			TFloatArrayList current = map.putIfAbsent(index, next);
			if (current == null)
			{
				// If the return value is null then this is a new insertion.
				// Set the current value as the one we just added and create the next insertion object.
				current = next;
				next = new TFloatArrayList(2);
			}
			current.add(p.getXPosition());
			current.add(p.getYPosition());
		}

		// Create an array of all the sample entries.
		// This is used to sample regions by density.
		data = new ResultsSample[map.size()];
		map.forEachEntry(new TLongObjectProcedure<TFloatArrayList>()
		{
			int i = 0;

			public boolean execute(long a, TFloatArrayList b)
			{
				data[i++] = new ResultsSample(a, b.toArray());
				b.clear(0); // Clear memory
				return true;
			}
		});
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
	 * Gets the sample image. The image is a stack of the samples with an overlay of the localisation positions.
	 *
	 * @param nNo
	 *            the number of samples with no localisations
	 * @param nLow
	 *            the number of samples with low localisations
	 * @param nHigh
	 *            the number of samples with high localisations
	 * @return the sample image (could be null if no samples were made)
	 */
	public ImagePlus getSample(int nNo, int nLow, int nHigh)
	{
		ImageStack out = new ImageStack(size, size);
		if (!isValid())
			return null;

		list.clearf();

		// empty
		for (int i : Random.sample(nNo, no.length, r))
			list.add(new ResultsSample(no[i]));
		// low
		for (int i : Random.sample(nLow, lower, r))
			list.add(data[i]);
		// high
		for (int i : Random.sample(nHigh, upper, r))
			list.add(data[i + lower]);

		if (list.isEmpty())
			return null;

		// Sort descending by number in the frame
		ResultsSample[] sample = list.toArray(new ResultsSample[list.size()]);
		Arrays.sort(sample, rcc);

		int[] xyz = new int[3];
		Rectangle stackBounds = new Rectangle(stack.getWidth(), stack.getHeight());
		Overlay overlay = new Overlay();
		float[] ox = new float[10], oy = new float[10];
		for (ResultsSample s : sample)
		{
			getXYZ(s.index, xyz);

			// Construct the region to extract
			Rectangle target = new Rectangle(xyz[0], xyz[1], size, size);
			target = target.intersection(stackBounds);
			if (target.width == 0 || target.height == 0)
				continue;

			// Extract the frame
			int slice = xyz[2];
			ImageProcessor ip = stack.getProcessor(slice);

			// Cut out the desired pixels (some may be blank if the block overruns the source image)
			ImageProcessor ip2 = ip.createProcessor(size, size);
			for (int y = 0; y < target.height; y++)
				for (int x = 0, i = y * size, index = (y + target.y) * ip.getWidth() +
						target.x; x < target.width; x++, i++, index++)
				{
					ip2.setf(i, ip.getf(index));
				}

			int size = s.size();
			if (size > 0)
			{
				// Create an ROI with the localisations
				for (int i = 0, j = 0; i < size; i++)
				{
					ox[i] = s.data[j++] - xyz[0];
					oy[i] = s.data[j++] - xyz[1];
				}
				PointRoi roi = new PointRoi(ox, oy, size);
				roi.setPosition(out.getSize() + 1);
				overlay.add(roi);
			}

			out.addSlice(String.format("Frame=%d (x=%d,y=%d) (n=%d)", slice, xyz[0], xyz[1], size), ip2.getPixels());
		}

		if (out.getSize() == 0)
			return null;

		ImagePlus imp = new ImagePlus("Sample", out);
		imp.setOverlay(overlay);
		return imp;
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
