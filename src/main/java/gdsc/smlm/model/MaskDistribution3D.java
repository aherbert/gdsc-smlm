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
package gdsc.smlm.model;


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
	private UniformDistribution uniformDistribution;
	private int[] mask;
	private int[] indices;
	private final int maxx, maxy, maxz, maxx_maxy;
	private final double halfWidth, halfHeight;
	private double minDepth, depth;
	private int particle = 0;
	private final double scaleX, scaleY;

	private double sliceDepth;
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
		setUniformDistribution(uniformDistribution);
		maxx = width;
		maxy = height;
		maxz = masks.size();
		this.scaleX = scaleX;
		this.scaleY = scaleY;
		halfWidth = width * 0.5;
		halfHeight = height * 0.5;
		this.sliceDepth = sliceDepth;
		depth = sliceDepth * maxz;
		minDepth = depth * -0.5;

		maxx_maxy = maxx * maxy;
		mask = new int[maxz * maxx_maxy];
		indices = new int[mask.length];

		int count = 0, index = 0;
		for (int[] mask : masks)
		{
			if (mask.length < maxx_maxy)
				throw new IllegalArgumentException("Masks must be the same size");
			for (int i = 0; i < maxx_maxy; i++)
			{
				this.mask[index] = mask[i];
				if (mask[i] != 0)
					indices[count++] = index;
				index++;
			}
		}

		if (count == 0)
			throw new IllegalArgumentException("Mask must have non-zero pixels");

		indices = Arrays.copyOf(indices, count);

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
		int[] xyz = new int[3];
		getXYZ(indices[randomPosition], xyz);
		final double[] d = uniformDistribution.nextUnit();

		// Ensure XY = 0 is the centre of the image by subtracting half the width/height
		d[0] = (xyz[0] + d[0] - halfWidth) * scaleX;
		d[1] = (xyz[1] + d[1] - halfHeight) * scaleY;
		d[2] = (xyz[2] + d[2]) * sliceDepth + minDepth;

		return d;
	}

	private int getIndex(double[] xyz)
	{
		int x = (int) (halfWidth + xyz[0] / scaleX);
		int y = (int) (halfHeight + xyz[1] / scaleY);
		int z = (int) ((xyz[2] - minDepth) / sliceDepth);
		if (x < 0 || x >= maxx || y < 0 || y >= maxy || z < 0 || z >= maxz)
			return -1;
		return getIndex(x, y, z);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialDistribution#isWithin(double[])
	 */
	public boolean isWithin(double[] xyz)
	{
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
			int[] mask2 = new int[maxx_maxy];
			for (int z = 0, index = 0; z < maxz; z++)
			{
				for (int j = 0; j < maxx_maxy; j++)
				{
					if (mask[index++] != 0)
						mask2[j] = 1;
				}
			}
			projection = new MaskDistribution(mask2, maxx, maxy, sliceDepth, scaleX, scaleY, randomGenerator);
		}
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
		final int index = getIndex(xyz);
		particle = (index < 0 || index >= mask.length) ? 0 : mask[index];

		// Also initialise for isWithinXY()
		createProjection();
		projection.initialise(xyz);
	}

	// Used for a particle search
	private final int[] DIR_X_OFFSET = new int[] { 0, 1, 1, 1, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, -1, -1, 0, 0, 1, 1, 1,
			0, -1, -1, -1, 0 };
	private final int[] DIR_Y_OFFSET = new int[] { -1, -1, 0, 1, 1, 1, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, 0, -1, -1, 0,
			1, 1, 1, 0, -1, 0 };
	private final int[] DIR_Z_OFFSET = new int[] { 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1,
			1, 1, 1, 1, 1, 1 };
	private int xlimit = -1, ylimit, zlimit;
	private int[] offset = null;

	/**
	 * Convert the mask to connected particles, each with a unique number.
	 * This allows the within function to restrict movement to the particle of origin
	 */
	private void findParticles()
	{
		// Check if already initialised
		if (xlimit != -1)
			return;

		xlimit = maxx - 1;
		ylimit = maxy - 1;
		zlimit = maxz - 1;

		// Create the offset table (for single array 3D neighbour comparisons)
		offset = new int[DIR_X_OFFSET.length];
		for (int d = offset.length; d-- > 0;)
		{
			offset[d] = getIndex(DIR_X_OFFSET[d], DIR_Y_OFFSET[d], DIR_Z_OFFSET[d]);
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
			}
		}

		//// Debug 
		//ImageStack stack = new ImageStack(maxx, maxy, maxz);
		//for (int i=0, index=0; i<maxz; i++)
		//{
		//	short[] pixels = new short[maxx_maxy];
		//	for (int j=0; j<maxx_maxy; j++)
		//		pixels[j] = (short) mask[index++];
		//	stack.setPixels(pixels, i+1);
		//}
		//Utils.display("objects", stack);

		// Free memory
		offset = null;
	}

	/**
	 * Return the single index associated with the x,y,z coordinates
	 * 
	 * @param x
	 * @param y
	 * @param z
	 * @return The index
	 */
	private int getIndex(int x, int y, int z)
	{
		return (maxx_maxy) * z + maxx * y + x;
	}

	/**
	 * Convert the single index into x,y,z coords, Input array must be length >= 3.
	 * 
	 * @param index
	 * @param xyz
	 * @return The xyz array
	 */
	private int[] getXYZ(int index, int[] xyz)
	{
		xyz[2] = index / (maxx_maxy);
		int mod = index % (maxx_maxy);
		xyz[1] = mod / maxx;
		xyz[0] = mod % maxx;
		return xyz;
	}

	/**
	 * Searches from the specified point to find all connected points and assigns them to given particle.
	 */
	private void expandParticle(boolean[] binaryMask, int[] mask, int[] pList, int index0, final int particle)
	{
		binaryMask[index0] = false; // mark as processed
		int listI = 0; // index of current search element in the list
		int listLen = 1; // number of elements in the list

		int[] xyz = new int[3];

		// we create a list of connected points and start the list at the particle start position
		pList[listI] = index0;

		do
		{
			int index1 = pList[listI];
			// Mark this position as part of the particle
			mask[index1] = particle;

			getXYZ(index1, xyz);

			final int x1 = xyz[0];
			final int y1 = xyz[1];
			final int z1 = xyz[2];

			final boolean isInnerXY = (y1 != 0 && y1 != ylimit) && (x1 != 0 && x1 != xlimit);
			final boolean isInnerXYZ = (zlimit == 0) ? isInnerXY : isInnerXY && (z1 != 0 && z1 != zlimit);

			// Search the neighbours 
			for (int d = 26; d-- > 0;)
			{
				if (isInnerXYZ || (isInnerXY && isWithinZ(z1, d)) || isWithinXYZ(x1, y1, z1, d))
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

			listI++;

		} while (listI < listLen);
	}

	/**
	 * returns whether the neighbour in a given direction is within the image. NOTE: it is assumed that the pixel x,y,z
	 * itself is within the image! Uses class variables xlimit, ylimit, zlimit: (dimensions of the image)-1
	 * 
	 * @param x
	 *            x-coordinate of the pixel that has a neighbour in the given direction
	 * @param y
	 *            y-coordinate of the pixel that has a neighbour in the given direction
	 * @param z
	 *            z-coordinate of the pixel that has a neighbour in the given direction
	 * @param direction
	 *            the direction from the pixel towards the neighbour
	 * @return true if the neighbour is within the image (provided that x, y, z is within)
	 */
	private boolean isWithinXYZ(int x, int y, int z, int direction)
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
			case 8:
				return (z > 0 && y > 0);
			case 9:
				return (z > 0 && y > 0 && x < xlimit);
			case 10:
				return (z > 0 && x < xlimit);
			case 11:
				return (z > 0 && y < ylimit && x < xlimit);
			case 12:
				return (z > 0 && y < ylimit);
			case 13:
				return (z > 0 && y < ylimit && x > 0);
			case 14:
				return (z > 0 && x > 0);
			case 15:
				return (z > 0 && y > 0 && x > 0);
			case 16:
				return (z > 0);
			case 17:
				return (z < zlimit && y > 0);
			case 18:
				return (z < zlimit && y > 0 && x < xlimit);
			case 19:
				return (z < zlimit && x < xlimit);
			case 20:
				return (z < zlimit && y < ylimit && x < xlimit);
			case 21:
				return (z < zlimit && y < ylimit);
			case 22:
				return (z < zlimit && y < ylimit && x > 0);
			case 23:
				return (z < zlimit && x > 0);
			case 24:
				return (z < zlimit && y > 0 && x > 0);
			case 25:
				return (z < zlimit);
		}
		return false;
	}

	/**
	 * returns whether the neighbour in a given direction is within the image. NOTE: it is assumed that the pixel z
	 * itself is within the image! Uses class variables zlimit: (dimensions of the image)-1
	 * 
	 * @param z
	 *            z-coordinate of the pixel that has a neighbour in the given direction
	 * @param direction
	 *            the direction from the pixel towards the neighbour
	 * @return true if the neighbour is within the image (provided that z is within)
	 */
	private boolean isWithinZ(int z, int direction)
	{
		// z = 0
		if (direction < 8)
			return true;
		// z = -1
		if (direction < 17)
			return (z > 0);
		// z = 1
		return z < zlimit;
	}

	/**
	 * @return The number of non-zero pixels in the mask
	 */
	public int getSize()
	{
		return indices.length;
	}

	/**
	 * @return The width of the mask in pixels
	 */
	public int getWidth()
	{
		return maxx;
	}

	/**
	 * @return The height of the mask in pixels
	 */
	public int getHeight()
	{
		return maxy;
	}

	/**
	 * @return The depth of the mask in pixels
	 */
	public int getDepth()
	{
		return maxz;
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
	 * @return The mask (packed in ZYX order)
	 */
	protected int[] getMask()
	{
		return mask;
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
		this.uniformDistribution = uniformDistribution;
	}
}
