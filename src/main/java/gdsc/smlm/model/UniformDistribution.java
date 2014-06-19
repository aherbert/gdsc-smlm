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

import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * Samples uniformly from the specified volume
 */
public class UniformDistribution implements SpatialDistribution
{
	private double[] min, max, range;
	private RandomGenerator randomGenerator;

	public UniformDistribution(double[] max)
	{
		this(null, max, null);
	}

	public UniformDistribution(double[] min, double[] max)
	{
		this(min, max, null);
	}

	public UniformDistribution(double[] max, RandomGenerator randomGenerator)
	{
		this(null, max, randomGenerator);
	}

	public UniformDistribution(double[] min, double[] max, RandomGenerator randomGenerator)
	{
		if (randomGenerator == null)
			randomGenerator = new JDKRandomGenerator();
		if (min == null)
			min = new double[0];
		if (max == null)
			max = new double[0];

		this.randomGenerator = randomGenerator;

		this.min = new double[3];
		for (int i = 0; i < min.length; i++)
			this.min[i] = min[i];
		this.max = new double[3];
		for (int i = 0; i < max.length; i++)
		{
			if (max[i] < this.min[i])
				throw new IllegalArgumentException(String.format("Max %f must be greater than min %f", max[i],
						this.min[i]));
			this.max[i] = max[i];
		}

		this.range = new double[3];
		for (int i = 0; i < this.max.length; i++)
			this.range[i] = this.max[i] - this.min[i];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialDistribution#next()
	 */
	public double[] next()
	{
		double[] d = new double[3];
		for (int i = 0; i < 3; i++)
			d[i] = min[i] + randomGenerator.nextDouble() * range[i];
		return d;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialDistribution#isWithin(double[])
	 */
	public boolean isWithin(double[] xyz)
	{
		for (int i = 0; i < xyz.length; i++)
			if (xyz[i] < min[i] || xyz[i] > max[i])
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
			if (xyz[i] < min[i] || xyz[i] > max[i])
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
