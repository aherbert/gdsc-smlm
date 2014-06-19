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

import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * Populates an image with well spaced unary or binary localisations.
 * <p>
 * Creates a grid layout using cells of the specified size within the image area. Each cell can have one or two
 * localisations. The first localisation is placed within the central 50% of the cell. The second localisation (if
 * present) is placed randomly up to a maximum distance away. When all cells have been sampled then no more
 * localisations are generated.
 */
public class GridDistribution implements SpatialDistribution
{
	private RandomGenerator randomGenerator;
	private RandomDataGenerator dataGenerator;
	private int size, cellSize;
	private double pBinary;
	private double minBinaryDistance, maxBinaryDistance;
	private double min, depth;

	private int cell = -1;
	private int nCellsPerRow, nCells;
	private double[] previous = null;

	/**
	 * Create a distribution with the binary spots placed from 0 - distance  
	 * 
	 * @param size
	 * @param depth
	 * @param cellSize
	 * @param pBinary
	 * @param binaryDistance
	 */
	public GridDistribution(int size, double depth, int cellSize, double pBinary, double binaryDistance)
	{
		this(size, depth, cellSize, pBinary, 0, binaryDistance, null);
	}

	/**
	 * Create a distribution with the binary spots placed from min - max distance  
	 * 
	 * @param size
	 * @param depth
	 * @param cellSize
	 * @param pBinary
	 * @param minBinaryDistance
	 * @param maxBinaryDistance
	 */
	public GridDistribution(int size, double depth, int cellSize, double pBinary, double minBinaryDistance, double maxBinaryDistance)
	{
		this(size, depth, cellSize, pBinary, minBinaryDistance, maxBinaryDistance, null);
	}

	/**
	 * Create a distribution with the binary spots placed from min - max distance  
	 * 
	 * @param size
	 * @param depth
	 * @param cellSize
	 * @param pBinary
	 * @param minBinaryDistance
	 * @param maxBinaryDistance
	 * @param randomGenerator
	 */
	public GridDistribution(int size, double depth, int cellSize, double pBinary, double minBinaryDistance, double maxBinaryDistance,
			RandomGenerator randomGenerator)
	{
		if (size < 1)
			throw new IllegalArgumentException("Size must be above zero");
		if (size < cellSize)
			throw new IllegalArgumentException("Size must be >= cell size");
		if (pBinary < 0 || pBinary > 1)
			throw new IllegalArgumentException("Probability must be between 0 and 1");
		if (maxBinaryDistance < 0)
			throw new IllegalArgumentException("Max distance must be positive");
		if (minBinaryDistance > maxBinaryDistance)
			throw new IllegalArgumentException("Min distance must be below max distance");
		if (randomGenerator == null)
			randomGenerator = new JDKRandomGenerator();

		this.randomGenerator = randomGenerator;
		this.dataGenerator = new RandomDataGenerator(randomGenerator);
		this.size = size;
		this.min = -depth / 2;
		this.depth = depth;
		this.cellSize = cellSize;
		this.pBinary = pBinary;
		this.minBinaryDistance = minBinaryDistance;
		this.maxBinaryDistance = maxBinaryDistance;

		nCellsPerRow = size / cellSize;
		nCells = nCellsPerRow * nCellsPerRow;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialDistribution#next()
	 */
	public double[] next()
	{
		if (previous != null)
		{
			// See if a binary localisation should be created near the previous spot
			if (randomGenerator.nextDouble() < pBinary)
			{
				double[] xyz = Arrays.copyOf(previous, 3);

				// Create a random unit vector
				double x = dataGenerator.nextGaussian(0, 1);
				double y = dataGenerator.nextGaussian(0, 1);
				double z = dataGenerator.nextGaussian(0, 1);
				final double length = Math.sqrt(x * x + y * y + z * z);
				if (length != 0)
				{
					// Shift by a random distance
					final double distance = (maxBinaryDistance == minBinaryDistance) ? maxBinaryDistance : 
							dataGenerator.nextUniform(minBinaryDistance, maxBinaryDistance, true);
					final double d = distance / length;
					x *= d;
					y *= d;
					z *= d;
				}
				xyz[0] += x;
				xyz[1] += y;
				xyz[2] += z;
				previous = null;
				return xyz;
			}
		}
		previous = null;
		// See if any more localisations will fit in the grid
		if (++cell < nCells)
		{
			int cellx = cell % nCellsPerRow;
			int celly = cell / nCellsPerRow;

			previous = new double[3];
			// Ensure the centre of the distribution is [0,0,0]
			previous[0] = cellx * cellSize - size / 2 + cellSize * dataGenerator.nextUniform(0.25, 0.75);
			previous[1] = celly * cellSize - size / 2 + cellSize * dataGenerator.nextUniform(0.25, 0.75);
			previous[2] = min + randomGenerator.nextDouble() * depth;
		}
		return previous;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialDistribution#isWithin(double[])
	 */
	public boolean isWithin(double[] xyz)
	{
		for (int i = 0; i < xyz.length; i++)
			if (xyz[i] < 0 || xyz[i] > size)
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
		for (int i = 0; i < 2; i++)
			if (xyz[i] < 0 || xyz[i] > size)
				return false;
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialDistribution#initialise(double[])
	 */
	public void initialise(double[] xyz)
	{
		// Ignore		
	}
}
