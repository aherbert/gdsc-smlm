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

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

/**
 * Samples uniformly from the specified mask. All non-zero pixels are sampled. The centre of the
 * mask corresponds to XY=0. The z coordinate is randomly sampled from the image depth and centred on
 * zero.
 * <p>
 * X coordinates are returned in the interval -width/2 to width/2. These can be converted to different values using the
 * scale parameter. Likewise for the Y coordinates. E.g. a mask of 100x100 (range of -50:50) can be used to generate
 * coordinates in the range -100:100 using a scale of 2.
 * <p>
 * Sub pixel locations and z-depth are sampled from a uniform distribution. A Halton sequence is used by default but
 * this can be changed by setting a custom uniform distribution.
 */
public class MaskDistribution implements SpatialDistribution
{
	private RandomGenerator randomGenerator;
	private UniformDistribution uniformDistribution;
	private int[] mask;
	private int[] indices;
	private final int width, height, half_width, half_height;
	private final double min, depth;
	private int particle = 0;
	private final double scaleX, scaleY;

	/**
	 * Create a distribution from the mask image (packed in YX order)
	 * 
	 * @param mask
	 * @param width
	 *            The width of the mask in pixels
	 * @param height
	 *            the height of the mask in pixels
	 * @param depth
	 *            The mask depth
	 * @param scaleX
	 *            Used to scale the mask X-coordinate to a new value
	 * @param scaleY
	 *            Used to scale the mask Y-coordinate to a new value
	 */
	public MaskDistribution(byte[] mask, int width, int height, double depth, double scaleX, double scaleY)
	{
		this(mask, width, height, depth, scaleX, scaleY, null);
	}

	/**
	 * Create a distribution from the mask image (packed in YX order)
	 * 
	 * @param mask
	 * @param width
	 *            The width of the mask in pixels
	 * @param height
	 *            the height of the mask in pixels
	 * @param depth
	 *            The mask depth
	 * @param scaleX
	 *            Used to scale the mask X-coordinate to a new value
	 * @param scaleY
	 *            Used to scale the mask Y-coordinate to a new value
	 */
	public MaskDistribution(int[] mask, int width, int height, double depth, double scaleX, double scaleY)
	{
		this(mask, width, height, depth, scaleX, scaleY, null);
	}

	/**
	 * Create a distribution from the mask image (packed in YX order)
	 * 
	 * @param mask
	 * @param width
	 *            The width of the mask in pixels
	 * @param height
	 *            the height of the mask in pixels
	 * @param depth
	 *            The mask depth
	 * @param scaleX
	 *            Used to scale the mask X-coordinate to a new value
	 * @param scaleY
	 *            Used to scale the mask Y-coordinate to a new value
	 * @param randomGenerator
	 *            Used to pick random pixels in the mask
	 */
	public MaskDistribution(byte[] mask, int width, int height, double depth, double scaleX, double scaleY,
			RandomGenerator randomGenerator)
	{
		this(convert(mask), width, height, depth, scaleX, scaleY, randomGenerator);
	}

	private static int[] convert(byte[] mask)
	{
		if (mask == null)
			return null;
		int[] newMask = new int[mask.length];
		for (int i = 0; i < mask.length; i++)
		{
			if (mask[i] != 0)
				newMask[i] = 1;
		}
		return newMask;
	}

	/**
	 * Create a distribution from the mask image (packed in YX order)
	 * 
	 * @param mask
	 * @param width
	 *            The width of the mask in pixels
	 * @param height
	 *            the height of the mask in pixels
	 * @param depth
	 *            The mask depth
	 * @param scaleX
	 *            Used to scale the mask X-coordinate to a new value
	 * @param scaleY
	 *            Used to scale the mask Y-coordinate to a new value
	 * @param randomGenerator
	 *            Used to pick random pixels in the mask
	 */
	public MaskDistribution(int[] mask, int width, int height, double depth, double scaleX, double scaleY,
			RandomGenerator randomGenerator)
	{
		this(mask, width, height, depth, scaleX, scaleY, randomGenerator, null);
	}

	/**
	 * Create a distribution from the mask image (packed in YX order)
	 * 
	 * @param mask
	 * @param width
	 *            The width of the mask in pixels
	 * @param height
	 *            the height of the mask in pixels
	 * @param depth
	 *            The mask depth
	 * @param scaleX
	 *            Used to scale the mask X-coordinate to a new value
	 * @param scaleY
	 *            Used to scale the mask Y-coordinate to a new value
	 * @param randomGenerator
	 *            Used to pick random pixels in the mask
	 * @param uniformDistribution
	 *            Used for sub-pixel location and z-depth
	 */
	public MaskDistribution(int[] mask, int width, int height, double depth, double scaleX, double scaleY,
			RandomGenerator randomGenerator, UniformDistribution uniformDistribution)
	{
		if (width < 1 || height < 1)
			throw new IllegalArgumentException("Dimensions must be above zero");
		if (scaleX < 0 || scaleY < 0)
			throw new IllegalArgumentException("Scale must be above zero");
		if (mask == null || mask.length < width * height)
			throw new IllegalArgumentException("Mask must not be null and must at least (width * height) in size");
		if (randomGenerator == null)
			randomGenerator = new Well19937c(System.currentTimeMillis() + System.identityHashCode(this));

		this.randomGenerator = randomGenerator;
		setUniformDistribution(uniformDistribution);
		this.mask = mask;
		this.width = width;
		this.scaleX = scaleX;
		this.scaleY = scaleY;
		this.height = height;
		this.half_width = width / 2;
		this.half_height = height / 2;
		this.min = -depth / 2;
		this.depth = depth;

		final int size = width * height;
		int count = 0;
		for (int i = 0; i < size; i++)
		{
			if (mask[i] != 0)
				count++;
		}

		if (count == 0)
			throw new IllegalArgumentException("Mask must have non-zero pixels");

		indices = new int[count];
		count = 0;
		for (int i = 0; i < size; i++)
		{
			if (mask[i] != 0)
				indices[count++] = i;
		}

		// Fischer-Yates shuffle the indices to scramble the mask positions
		for (int i = indices.length; i-- > 1;)
		{
			final int j = randomGenerator.nextInt(i + 1);
			final int tmp = indices[i];
			indices[i] = indices[j];
			indices[j] = tmp;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialDistribution#next()
	 */
	public double[] next()
	{
		final int randomPosition = randomGenerator.nextInt(indices.length);
		// Ensure XY = 0 is the centre of the image
		final int x = indices[randomPosition] % width - half_width;
		final int y = indices[randomPosition] / width - half_height;
		final double[] d = uniformDistribution.nextUnit();
		d[0] = (x + d[0]) * scaleX;
		d[1] = (y + d[1]) * scaleY;
		d[2] = min + d[2] * depth;
		return d;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialDistribution#isWithin(double[])
	 */
	public boolean isWithin(double[] xyz)
	{
		if (!isWithinXY(xyz))
			return false;
		if (xyz[2] < min || xyz[2] > min + depth)
			return false;
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialDistribution#isWithinXY(double[])
	 */
	public boolean isWithinXY(double[] xyz)
	{
		// Ensure XY = 0 is the centre of the image
		int index = getIndex(xyz);
		if (index < 0 || index >= mask.length || mask[index] == 0)
			return false;
		// Check if the search was initialised in a particle 
		if (particle == 0)
		{
			// No starting particle so just accept the position
			return true;
		}
		// Must be in the same particle as the initial position
		return mask[index] == particle;
	}

	private int getIndex(double[] xyz)
	{
		int x = (int) (xyz[0] / scaleX) + half_width;
		int y = (int) (xyz[1] / scaleY) + half_height;
		if (x < 0 || x >= width || y < 0 || y >= height)
			return -1;
		int index = y * width + x;
		return index;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialDistribution#initialise(double[])
	 */
	public void initialise(double[] xyz)
	{
		findParticles();

		// Now store the particle that contains the position
		int index = getIndex(xyz);
		particle = (index < 0 || index >= mask.length) ? 0 : mask[index];
	}

	// Used for a particle search
	private static final int[] DIR_X_OFFSET = new int[] { 0, 1, 1, 1, 0, -1, -1, -1 };
	private static final int[] DIR_Y_OFFSET = new int[] { -1, -1, 0, 1, 1, 1, 0, -1 };
	private int maxx = 0, maxy;
	private int xlimit, ylimit;
	private int[] offset = null;

	/**
	 * Convert the mask to connected particles, each with a unique number.
	 * This allows the within function to restrict movement to the particle of origin
	 */
	private void findParticles()
	{
		// Check if already initialised
		if (maxx > 0)
			return;

		maxx = width;
		maxy = height;

		xlimit = maxx - 1;
		ylimit = maxy - 1;

		// Create the offset table (for single array 2D neighbour comparisons)
		offset = new int[DIR_X_OFFSET.length];
		for (int d = offset.length; d-- > 0;)
		{
			offset[d] = maxx * DIR_Y_OFFSET[d] + DIR_X_OFFSET[d];
		}

		int[] pList = new int[mask.length];

		// Store all the non-zero positions
		boolean[] binaryMask = new boolean[mask.length];
		for (int i = 0; i < mask.length; i++)
			binaryMask[i] = (mask[i] != 0);

		// Find particles 
		int particles = 0;
		for (int i = 0; i < binaryMask.length; i++)
		{
			if (binaryMask[i])
			{
				expandParticle(binaryMask, mask, pList, i, ++particles);

				// Debug: Show the particles
				//float[] pixels = new float[mask.length];
				//for (int j = 0; j < mask.length; j++)
				//	pixels[j] = (binaryMask[j]) ? 0 : mask[j];
				//Utils.display("Particle", new FloatProcessor(maxx, maxy, pixels));
			}
		}

		// Free memory
		offset = null;
	}

	/**
	 * Searches from the specified point to find all connected points and assigns them to given particle.
	 */
	private void expandParticle(boolean[] binaryMask, int[] mask, int[] pList, int index0, final int particle)
	{
		binaryMask[index0] = false; // mark as processed
		int listI = 0; // index of current search element in the list
		int listLen = 1; // number of elements in the list

		// we create a list of connected points and start the list at the particle start position
		pList[listI] = index0;

		do
		{
			int index1 = pList[listI];
			// Mark this position as part of the particle
			mask[index1] = particle;

			// Search the 8-connected neighbours 
			int x1 = index1 % maxx;
			int y1 = index1 / maxx;

			boolean isInnerXY = (x1 != 0 && x1 != xlimit) && (y1 != 0 && y1 != ylimit);

			if (isInnerXY)
			{
				for (int d = 8; d-- > 0;)
				{
					int index2 = index1 + offset[d];
					if (binaryMask[index2])
					{
						binaryMask[index2] = false; // mark as processed
						// Add this to the search
						pList[listLen++] = index2;
					}
				}
			}
			else
			{
				for (int d = 8; d-- > 0;)
				{
					if (isInside(x1, y1, d))
					{
						int index2 = index1 + offset[d];
						if (binaryMask[index2])
						{
							binaryMask[index2] = false; // mark as processed
							// Add this to the search
							pList[listLen++] = index2;
						}
					}
				}
			}

			listI++;

		} while (listI < listLen);
	}

	/**
	 * Returns whether the neighbour in a given direction is within the image. NOTE: it is assumed that the pixel x,y
	 * itself is within the image! Uses class variables xlimit, ylimit: (dimensions of the image)-1
	 * 
	 * @param x
	 *            x-coordinate of the pixel that has a neighbour in the given direction
	 * @param y
	 *            y-coordinate of the pixel that has a neighbour in the given direction
	 * @param direction
	 *            the direction from the pixel towards the neighbour
	 * @return true if the neighbour is within the image (provided that x, y is within)
	 */
	private boolean isInside(int x, int y, int direction)
	{
		switch (direction)
		{
			case 0:
				return (y > 0);
			case 1:
				return (y > 0 && x < xlimit);
			case 2:
				return (x < xlimit);
			case 3:
				return (y < ylimit && x < xlimit);
			case 4:
				return (y < ylimit);
			case 5:
				return (y < ylimit && x > 0);
			case 6:
				return (x > 0);
			case 7:
				return (y > 0 && x > 0);
		}
		return false;
	}

	/**
	 * @return The number of non-zero pixels in the mask
	 */
	public int getSize()
	{
		return indices.length;
	}

	/**
	 * @return The width of the mask
	 */
	public int getWidth()
	{
		return width;
	}

	/**
	 * @return The height of the mask
	 */
	public int getHeight()
	{
		return height;
	}

	/**
	 * @return The X-scale
	 */
	public double getScaleX()
	{
		return scaleX;
	}

	/**
	 * @return The Y-scale
	 */
	public double getScaleY()
	{
		return scaleY;
	}

	/**
	 * @return The mask
	 */
	protected int[] getMask()
	{
		return mask;
	}

	/**
	 * The UniformDistribution to pick the sub pixel x,y coordinates and z-depth
	 * 
	 * @param uniformDistribution
	 *            the uniformDistribution to set
	 */
	public void setUniformDistribution(UniformDistribution uniformDistribution)
	{
		if (uniformDistribution == null)
			uniformDistribution = new UniformDistribution(null, new double[] { 1, 1, 1 }, randomGenerator.nextInt());
		this.uniformDistribution = uniformDistribution;
	}
}
