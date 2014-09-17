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

import org.apache.commons.math3.random.HaltonSequenceGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomVectorGenerator;
import org.apache.commons.math3.random.Well44497b;

/**
 * Samples uniformly from the specified volume.
 * <p>
 * Uses a Halton sequence generator by default. This can be overriden by passing in a vector sequence generator or by
 * setting a random generator (which should produce uniform equi-distributed numbers in the domain [0,1]).
 */
public class UniformDistribution implements SpatialDistribution
{
	/**
	 * Wrap a standard random generator to create a vector generator for 3 dimensions
	 */
	private class VectorGeneratorWrapper implements RandomVectorGenerator
	{
		private RandomGenerator randomGenerator;

		public VectorGeneratorWrapper(RandomGenerator randomGenerator)
		{
			if (randomGenerator == null)
				randomGenerator = new Well44497b(System.currentTimeMillis() + System.identityHashCode(this));
			this.randomGenerator = randomGenerator;
		}

		/**
		 * @return
		 */
		@Override
		public double[] nextVector()
		{
			return new double[] { randomGenerator.nextDouble(), randomGenerator.nextDouble(),
					randomGenerator.nextDouble() };
		}
	}

	private double[] min, max, range;
	private RandomVectorGenerator vectorGenerator;

	public UniformDistribution(double[] max)
	{
		init(min, max, null);
	}

	public UniformDistribution(double[] min, double[] max)
	{
		init(min, max, null);
	}

	/**
	 * @param min
	 * @param max
	 * @param seed Start at the i-th point in the Halton sequence
	 */
	public UniformDistribution(double[] min, double[] max, int seed)
	{
		HaltonSequenceGenerator randomVectorGenerator = new HaltonSequenceGenerator(3);
		randomVectorGenerator.skipTo(Math.abs(seed));
		init(min, max, randomVectorGenerator);
	}
	
	/**
	 * @param min
	 * @param max
	 * @param randomVectorGenerator
	 *            Must produce vectors with dimension 3 (or above)
	 */
	public UniformDistribution(double[] min, double[] max, RandomVectorGenerator randomVectorGenerator)
	{
		init(min, max, randomVectorGenerator);
	}

	private void init(double[] min, double[] max, RandomVectorGenerator randomVectorGenerator)
	{
		if (min == null)
			min = new double[0];
		if (max == null)
			max = new double[0];

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
		
		if (randomVectorGenerator == null)
			randomVectorGenerator = new HaltonSequenceGenerator(3);
		this.vectorGenerator = randomVectorGenerator;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialDistribution#next()
	 */
	public double[] next()
	{
		double[] d = vectorGenerator.nextVector();
		for (int i = 0; i < 3; i++)
			d[i] = min[i] + d[i] * range[i];
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

	/**
	 * Set the random generator. The generator should produce uniform equi-distributed numbers in the domain [0,1].
	 * 
	 * @param randomGenerator If null then a Well generator is used
	 */
	public void setRandomGenerator(RandomGenerator randomGenerator)
	{
		this.vectorGenerator = new VectorGeneratorWrapper(randomGenerator);
	}
}
