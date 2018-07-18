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
package gdsc.smlm.ij.results;

import java.awt.Rectangle;
import java.util.Arrays;
import java.util.Comparator;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

import gdsc.core.ij.Utils;
import gdsc.core.utils.Random;
import gdsc.core.utils.TurboList;
import gdsc.smlm.data.config.CalibrationReader;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gnu.trove.list.array.TLongArrayList;
import gnu.trove.map.hash.TLongObjectHashMap;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.process.ImageProcessor;

/**
 * Allows sampling sections from a source image for localisation results.
 */
public class ResultsImageSampler
{
	private static class PeakResultList
	{
		int size = 0;
		PeakResult[] data;

		PeakResultList()
		{
			this(new PeakResult[1]);
		}

		PeakResultList(PeakResult[] data)
		{
			this.data = data;
		}

		void add(PeakResult p)
		{
			data[size++] = p;
			if (size == data.length)
			{
				final PeakResult[] data2 = new PeakResult[size * 2];
				System.arraycopy(data, 0, data2, 0, size);
				data = data2;
			}
		}

		PeakResult get(int i)
		{
			return data[i];
		}
	}

	private static class ResultsSample
	{
		long index;
		final PeakResultList list;

		static ResultsSample createEmpty(long index)
		{
			return new ResultsSample(index, new PeakResultList(null));
		}

		static ResultsSample create(long index)
		{
			return new ResultsSample(index, new PeakResultList());
		}

		ResultsSample(long index, PeakResultList data)
		{
			this.index = index;
			this.list = data;
		}

		public int size()
		{
			return list.size;
		}

		public void add(PeakResult p)
		{
			list.add(p);
		}
	}

	private static class IndexComparator implements Comparator<ResultsSample>
	{
		@Override
		public int compare(ResultsSample o1, ResultsSample o2)
		{
			return Long.compare(o1.index, o2.index);
		}
	}

	private static class CountComparator implements Comparator<ResultsSample>
	{
		@Override
		public int compare(ResultsSample o1, ResultsSample o2)
		{
			final int result = compareInt(o1.size(), o2.size());
			if (result == 0)
				// Use index if the same count
				return Long.compare(o1.index, o2.index);
			return result;
		}
	}

	private static class ReverseCountComparator implements Comparator<ResultsSample>
	{
		@Override
		public int compare(ResultsSample o1, ResultsSample o2)
		{
			final int result = compareInt(o2.size(), o1.size());
			if (result == 0)
				// Use index if the same count
				return Long.compare(o1.index, o2.index);
			return result;
		}
	}

	private static int compareInt(int x, int y)
	{
		return (x < y) ? -1 : ((x == y) ? 0 : 1);
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
	private final TurboList<ResultsSample> list = new TurboList<>();
	private RandomGenerator r = new Well19937c();

	/**
	 * Instantiates a new results image sampler.
	 *
	 * @param results
	 *            the results (must be in pixel units)
	 * @param stack
	 *            the source stack for the results
	 * @param size
	 *            the size of the sample blocks
	 */
	public ResultsImageSampler(MemoryPeakResults results, ImageStack stack, int size)
	{
		if (results.getDistanceUnit() != DistanceUnit.PIXEL)
			throw new IllegalArgumentException("Results must be in pixel units");

		this.results = results;
		this.stack = stack;
		this.size = size;

		final Rectangle bounds = results.getBounds(true);
		// Round the image dimensions to the nearest block interval
		lx = size * (bounds.x / size);
		ly = size * (bounds.y / size);
		final int ux = size * (int) Math.ceil((double) (bounds.x + bounds.width) / size);
		final int uy = size * (int) Math.ceil((double) (bounds.y + bounds.height) / size);
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
		final long total = stack.getSize() * xy_blocks;
		final long empty = total - data.length;
		if (empty == 0)
			// All indices are used
			no = new long[0];
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
				final long[] list = new long[(int) Math.min(empty, maxNumberOfEmptySamples)];
				int c = 0;
				OUTER: for (int i = 0; i < data.length; i++)
				{
					final long current = data[i].index;
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
				final TLongArrayList list = new TLongArrayList(data.length);
				if (empty < data.length)
				{
					// We can pick all the indices that are missing
					long emptyCandidate = 0;
					for (int i = 0; i < data.length; i++)
					{
						final long current = data[i].index;
						// If the current index is bigger than the candidate then it must be empty
						while (current > emptyCandidate)
							// Add all those that are empty
							list.add(emptyCandidate++);
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
						final long current = data[i].index;
						// If the current index is bigger than the candidate then it must be empty
						if (current > emptyCandidate)
							// Note: we only sample the next empty index after an index with data
							// This means the number of frames could be lower
							list.add(emptyCandidate);
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
		final TLongObjectHashMap<ResultsSample> map = new TLongObjectHashMap<>(results.size());
		ResultsSample next = ResultsSample.create(-1);
		for (final PeakResult p : results.toArray())
		{
			// Avoid invalid slices
			if (p.getFrame() < 1 || p.getFrame() > stack.getSize())
				continue;
			final long index = getIndex(p.getXPosition(), p.getYPosition(), p.getFrame());

			ResultsSample current = map.putIfAbsent(index, next);
			if (current == null)
			{
				// If the return value is null then this is a new insertion.
				// Set the current value as the one we just added and create the next insertion object.
				current = next;
				current.index = index;
				next = ResultsSample.create(-1);
			}
			current.add(p);
		}

		// Create an array of all the sample entries.
		// This is used to sample regions by density.
		data = new ResultsSample[map.size()];
		map.values(data);
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
	 *            the index
	 * @param xyz
	 *            the xyz
	 * @return The xyz array
	 */
	private int[] getXYZ(long index, int[] xyz)
	{
		xyz[2] = (int) (index / (xy_blocks));
		final int mod = (int) (index % (xy_blocks));
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
	 * Gets the sample image. The image is a stack of the samples with an overlay of the localisation positions. The
	 * info property is set with details of the localisations and the image is calibrated.
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
		final ImageStack out = new ImageStack(size, size);
		if (!isValid())
			return null;

		list.clearf();

		// empty
		for (final int i : Random.sample(nNo, no.length, r))
			list.add(ResultsSample.createEmpty(no[i]));
		// low
		for (final int i : Random.sample(nLow, lower, r))
			list.add(data[i]);
		// high
		for (final int i : Random.sample(nHigh, upper, r))
			list.add(data[i + lower]);

		if (list.isEmpty())
			return null;

		double nmPerPixel = 1;
		if (results.hasCalibration())
		{
			final CalibrationReader calibration = results.getCalibrationReader();
			if (calibration.hasNmPerPixel())
				nmPerPixel = calibration.getNmPerPixel();
		}

		// Sort descending by number in the frame
		final ResultsSample[] sample = list.toArray(new ResultsSample[list.size()]);
		Arrays.sort(sample, rcc);

		final int[] xyz = new int[3];
		final Rectangle stackBounds = new Rectangle(stack.getWidth(), stack.getHeight());
		final Overlay overlay = new Overlay();
		final float[] ox = new float[10], oy = new float[10];
		final StringBuilder sb = new StringBuilder();
		if (nmPerPixel == 1)
			sb.append("Sample X Y Z Signal\n");
		else
			sb.append("Sample X(nm) Y(nm) Z(nm) Signal\n");

		for (final ResultsSample s : sample)
		{
			getXYZ(s.index, xyz);

			// Construct the region to extract
			Rectangle target = new Rectangle(xyz[0], xyz[1], size, size);
			target = target.intersection(stackBounds);
			if (target.width == 0 || target.height == 0)
				continue;

			// Extract the frame
			final int slice = xyz[2];
			final ImageProcessor ip = stack.getProcessor(slice);

			// Cut out the desired pixels (some may be blank if the block overruns the source image)
			final ImageProcessor ip2 = ip.createProcessor(size, size);
			for (int y = 0; y < target.height; y++)
				for (int x = 0, i = y * size, index = (y + target.y) * ip.getWidth() +
						target.x; x < target.width; x++, i++, index++)
					ip2.setf(i, ip.getf(index));

			final int size = s.size();
			if (size > 0)
			{
				final int position = out.getSize() + 1;
				// Create an ROI with the localisations
				for (int i = 0; i < size; i++)
				{
					final PeakResult p = s.list.get(i);
					ox[i] = p.getXPosition() - xyz[0];
					oy[i] = p.getYPosition() - xyz[1];
					sb.append(position).append(' ');
					sb.append(Utils.rounded(ox[i] * nmPerPixel)).append(' ');
					sb.append(Utils.rounded(oy[i] * nmPerPixel)).append(' ');
					sb.append(Utils.rounded(p.getZPosition() * nmPerPixel)).append(' ');
					sb.append(Utils.rounded(p.getIntensity())).append('\n');
				}
				final PointRoi roi = new PointRoi(ox, oy, size);
				roi.setPosition(position);
				overlay.add(roi);
			}

			out.addSlice(String.format("Frame=%d @ %d,%d px (n=%d)", slice, xyz[0], xyz[1], size), ip2.getPixels());
		}

		if (out.getSize() == 0)
			return null;

		final ImagePlus imp = new ImagePlus("Sample", out);
		imp.setOverlay(overlay);
		// Note: Only the info property can be saved to a TIFF file
		imp.setProperty("Info", sb.toString());
		if (nmPerPixel != 1)
		{
			final ij.measure.Calibration cal = new ij.measure.Calibration();
			cal.setUnit("nm");
			cal.pixelHeight = cal.pixelWidth = nmPerPixel;
			imp.setCalibration(cal);
		}

		return imp;
	}

	/**
	 * Set the random generator for use during sampling
	 *
	 * @param r
	 *            the random generator to set (ignored if null)
	 */
	public void setRandom(RandomGenerator r)
	{
		if (r != null)
			this.r = r;
	}
}
