package gdsc.smlm.model;

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

import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

/**
 * Samples uniformly from the specified masks. All non-zero pixels are sampled. The centre of the
 * mask stack corresponds to XY=0. The z coordinate is randomly sampled from the slice depth offset
 * by the slice position in the stack. The distribution of Z is centred on zero.
 * <p>
 * X coordinates are returned in the interval -width/2 to width/2. These can be converted to different values using the
 * scale parameter. Likewise for the Y coordinates. E.g. a mask of 100x100 (range of -50:50) can be used to generate
 * coordinates in the range -100:100 using a scale of 2.
 * <p>
 * Sub pixel locations and z-depth are sampled from a uniform distribution. A Halton sequence is used by default but
 * this can be changed by setting a custom uniform distribution.
 */
public class MaskDistribution3D implements SpatialDistribution
{
	private RandomGenerator randomGenerator;
	private double sliceDepth;
	private double minDepth, depth;
	private MaskDistribution[] slices;
	private int size;
	private int[] cumulativeSize;
	private MaskDistribution projection = null;

	/**
	 * Create a distribution from the stack of mask images (packed in YX order)
	 * 
	 * @param masks
	 * @param width
	 *            The width of the mask in pixels
	 * @param height
	 *            the height of the mask in pixels
	 * @param sliceDepth
	 *            The depth of each slice
	 * @param scaleX
	 *            Used to scale the mask X-coordinate to a new value
	 * @param scaleY
	 *            Used to scale the mask Y-coordinate to a new value
	 */
	public MaskDistribution3D(List<int[]> masks, int width, int height, double sliceDepth, double scaleX, double scaleY)
	{
		this(masks, width, height, sliceDepth, scaleX, scaleY, null);
	}

	/**
	 * Create a distribution from the stack of mask images (packed in YX order)
	 * 
	 * @param masks
	 * @param width
	 *            The width of the mask in pixels
	 * @param height
	 *            the height of the mask in pixels
	 * @param sliceDepth
	 *            The depth of each slice
	 * @param scaleX
	 *            Used to scale the mask X-coordinate to a new value
	 * @param scaleY
	 *            Used to scale the mask Y-coordinate to a new value
	 * @param randomGenerator
	 *            Used to pick random pixels in the mask
	 */
	public MaskDistribution3D(List<int[]> masks, int width, int height, double sliceDepth, double scaleX,
			double scaleY, RandomGenerator randomGenerator)
	{
		this(masks, width, height, sliceDepth, scaleX, scaleY, randomGenerator, null);
	}

	/**
	 * Create a distribution from the stack of mask images (packed in YX order)
	 * 
	 * @param masks
	 * @param width
	 *            The width of the mask in pixels
	 * @param height
	 *            the height of the mask in pixels
	 * @param sliceDepth
	 *            The depth of each slice
	 * @param scaleX
	 *            Used to scale the mask X-coordinate to a new value
	 * @param scaleY
	 *            Used to scale the mask Y-coordinate to a new value
	 * @param randomGenerator
	 *            Used to pick random pixels in the mask
	 * @param uniformDistribution
	 *            Used for sub-pixel location and slice z-depth
	 */
	public MaskDistribution3D(List<int[]> masks, int width, int height, double sliceDepth, double scaleX,
			double scaleY, RandomGenerator randomGenerator, UniformDistribution uniformDistribution)
	{
		if (width < 1 || height < 1)
			throw new IllegalArgumentException("Dimensions must be above zero");
		if (sliceDepth <= 0)
			throw new IllegalArgumentException("Slice depth must be above zero");
		if (masks == null || masks.isEmpty())
			throw new IllegalArgumentException("Mask must not be null or empty");
		if (randomGenerator == null)
			randomGenerator = new Well19937c(System.currentTimeMillis() + System.identityHashCode(this));

		this.randomGenerator = randomGenerator;
		if (uniformDistribution == null)
			uniformDistribution = new UniformDistribution(null, new double[] { 1, 1, 1 }, randomGenerator.nextInt());

		// Create a mask for each distribution
		this.sliceDepth = sliceDepth;
		depth = sliceDepth * masks.size();
		minDepth = depth * -0.5;
		slices = new MaskDistribution[masks.size()];
		cumulativeSize = new int[masks.size()];
		int i = 0;
		final int length = masks.get(0).length;
		for (int[] mask : masks)
		{
			if (length != mask.length)
				throw new IllegalArgumentException("Masks must be the same size");
			slices[i] = new MaskDistribution(mask, width, height, sliceDepth, scaleX, scaleY, randomGenerator,
					uniformDistribution);
			size += slices[i].getSize();
			cumulativeSize[i] = size;
			i++;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialDistribution#next()
	 */
	public double[] next()
	{
		// Find a random pixel within the cumulative mask
		int n = randomGenerator.nextInt(size);
		for (int i = 0; i < cumulativeSize.length; i++)
		{
			// Find the correct mask
			if (n < cumulativeSize[i])
			{
				double[] xyz = slices[i].next();
				// Adjust the depth to the correct position.
				// Note the slice will return in the range -0.5*sliceDepth to 0.5*sliceDepth so
				// add half the slice depth to get the range 0-sliceDepth.
				xyz[2] += (i + 0.5) * sliceDepth + minDepth;
				return xyz;
			}
		}
		// This should not happen since n will be less than cumulativeSize[cumulativeSize.length-1]
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialDistribution#isWithin(double[])
	 */
	public boolean isWithin(double[] xyz)
	{
		MaskDistribution mask = getMask(xyz[2]);
		return (mask == null) ? false : mask.isWithinXY(xyz);
	}

	private MaskDistribution getMask(double z)
	{
		if (z < minDepth)
			return null;
		int slice = (int) ((z - minDepth) / sliceDepth);
		return (slice < slices.length) ? slices[slice] : null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialDistribution#isWithinXY(double[])
	 */
	public boolean isWithinXY(double[] xyz)
	{
		createProjection();
		return projection.isWithinXY(xyz);
	}

	private void createProjection()
	{
		// Create a projection of the masks
		if (projection == null)
		{
			int[] mask = Arrays.copyOf(slices[0].getMask(), slices[0].getMask().length);
			for (int i = 1; i < slices.length; i++)
			{
				int[] mask2 = slices[i].getMask();
				for (int j = 0; j < mask.length; j++)
				{
					if (mask2[j] != 0)
						mask[j] = 1;
				}
			}
			projection = new MaskDistribution(mask, slices[0].getWidth(), slices[0].getHeight(), sliceDepth,
					slices[0].getScaleX(), slices[0].getScaleY(), randomGenerator);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialDistribution#initialise(double[])
	 */
	public void initialise(double[] xyz)
	{
		MaskDistribution mask = getMask(xyz[2]);
		if (mask != null)
			mask.initialise(xyz);

		// Also initialise for isWithinXY()
		createProjection();
		projection.initialise(xyz);
	}

	/**
	 * The UniformDistribution to pick the sub pixel x,y coordinates and slice z-depth
	 * 
	 * @param uniformDistribution
	 *            the uniformDistribution to set
	 */
	public void setUniformDistribution(UniformDistribution uniformDistribution)
	{
		if (uniformDistribution == null)
			uniformDistribution = new UniformDistribution(null, new double[] { 1, 1, 1 }, randomGenerator.nextInt());
		for (MaskDistribution slice : slices)
			slice.setUniformDistribution(uniformDistribution);
	}
}
