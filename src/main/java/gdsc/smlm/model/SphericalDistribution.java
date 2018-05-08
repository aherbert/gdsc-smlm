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
import org.apache.commons.math3.util.FastMath;

/**
 * Samples uniformly from the specified spherical volume
 */
public class SphericalDistribution implements SpatialDistribution
{
	private final double radius, r2, range;
	private RandomGenerator randomGenerator;
	private boolean useRejectionMethod = true;
	private final double[] origin = new double[3];

	public SphericalDistribution(double radius)
	{
		this(radius, null);
	}

	public SphericalDistribution(double radius, RandomGenerator randomGenerator)
	{
		if (randomGenerator == null)
			randomGenerator = new JDKRandomGenerator();
		if (radius < 0)
			throw new IllegalArgumentException("Radius must be positive: {0}");
		this.radius = radius;
		this.r2 = radius * radius;
		this.range = 2 * radius;
		this.randomGenerator = randomGenerator;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialDistribution#next()
	 */
	public double[] next()
	{
		double[] xyz = new double[3];
		if (radius > 0)
		{
			// See: http://math.stackexchange.com/questions/87230/
			// picking-random-points-in-the-volume-of-sphere-with-uniform-probability

			if (useRejectionMethod)
			{
				// -=-=-=-
				// Rejection method:
				// Sample from a cube and then check if within a sphere 
				// -=-=-=-
				double d2 = 0;
				do
				{
					for (int i = 0; i < 3; i++)
					{
						//xyz[i] = randomGenerator.nextDouble() * ((randomGenerator.nextBoolean()) ? -radius : radius);
						// Avoid extra call to the random generator
						xyz[i] = randomGenerator.nextDouble() * range - radius;
					}
					d2 = xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2];
				} while (d2 > r2);
			}
			else
			{
				// -=-=-=-
				// Transformation method:
				// Generate a random point on the surface of the sphere and then sample within.
				// -=-=-=-

				// Generate a random unit vector: X1, X2, X3 sampled with mean 0 and variance 1
				for (int i = 0; i < 3; i++)
				{
					xyz[i] = randomGenerator.nextGaussian();
				}

				// Calculate the distance: RsU^1/3 / length			
				final double d = (radius * FastMath.cbrt(randomGenerator.nextDouble())) /
						Math.sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
				for (int i = 0; i < 3; i++)
					xyz[i] *= d;
			}
		}
		return xyz;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialDistribution#isWithin(double[])
	 */
	public boolean isWithin(double[] xyz)
	{
		final double[] delta = { xyz[0] - origin[0], xyz[1] - origin[1], xyz[2] - origin[2] };
		return (delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]) < r2;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialDistribution#isWithinXY(double[])
	 */
	public boolean isWithinXY(double[] xyz)
	{
		final double[] delta = { xyz[0] - origin[0], xyz[1] - origin[1] };
		return (delta[0] * delta[0] + delta[1] * delta[1]) < r2;
	}

	/**
	 * @return If true then sample from the distribution using the rejection method. The alternative is a transformation
	 *         method.
	 */
	public boolean isUseRejectionMethod()
	{
		return useRejectionMethod;
	}

	/**
	 * @param useRejectionMethod
	 *            If true then sample from the distribution using the rejection method. The alternative is a
	 *            transformation
	 *            method.
	 */
	public void setUseRejectionMethod(boolean useRejectionMethod)
	{
		this.useRejectionMethod = useRejectionMethod;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialDistribution#initialise(double[])
	 */
	public void initialise(double[] xyz)
	{
		if (xyz != null && xyz.length > 2)
		{
			for (int i = 0; i < 3; i++)
				origin[i] = xyz[i];
		}
	}
}
