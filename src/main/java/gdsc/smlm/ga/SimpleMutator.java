package gdsc.smlm.ga;

import org.apache.commons.math3.random.RandomDataGenerator;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2015 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Mutates the sequence by selecting random positions and random shifts.
 */
public class SimpleMutator extends Randomiser implements Mutator
{
	final double fraction;

	private boolean override = false;
	private double[] stepSize, lower, upper;

	/**
	 * @param random
	 * @param fraction
	 *            The fraction of the sequence positions to mutate on average
	 */
	public SimpleMutator(RandomDataGenerator random, double fraction)
	{
		super(random);
		this.fraction = fraction;
	}

	/**
	 * Override the mutation parameters that are obtained from the Chromosome interface.
	 * The arrays must match the fixed size of the Chromosome sequences to be mutated.
	 * <p>
	 * All settings are overriden even if null arrays are passed for some arguments.
	 * 
	 * @param stepSize
	 *            The mutation step size
	 * @param lower
	 *            The lower limit for the sequence positions
	 * @param upper
	 *            The upper limit for the sequence positions
	 */
	public void overrideChromosomeSettings(double[] stepSize, double[] lower, double[] upper)
	{
		this.stepSize = stepSize;
		this.lower = lower;
		this.upper = upper;
		override = true; // (stepSize != null || lower != null || upper != null);
	}

	/**
	 * Mutates the chromosome to form a new sequence.
	 * <p>
	 * The number of positions are chosen from a Poisson distribution with an average using a fraction of the total
	 * positions. The positions are then chosen randomly. Note that the same position may be chosen multiple times. The
	 * random shifts for each mutation are taken from a Gaussian using the chromosome mutation step range as the
	 * standard deviation.
	 * 
	 * @see gdsc.smlm.ga.Mutator#mutate(gdsc.smlm.ga.Chromosome)
	 */
	@Override
	public Chromosome mutate(Chromosome chromosome)
	{
		final double[] sequence = chromosome.sequence().clone();

		final double mean = fraction * chromosome.length();
		if (mean > 0)
		{
			int count = (int) random.nextPoisson(mean);
			final double[] step, min, max;
			if (override)
			{
				step = stepSize;
				min = lower;
				max = upper;
			}
			else
			{
				step = chromosome.mutationStepRange();
				min = chromosome.lowerLimit();
				max = chromosome.upperLimit();
			}
			// Override individually
			//final double[] step = (stepSize == null) ? chromosome.mutationStepRange() : stepSize;
			//final double[] min = (lower == null) ? chromosome.lowerLimit() : lower;
			//final double[] max = (upper == null) ? chromosome.upperLimit() : upper;

			while (count-- > 0)
			{
				int i = random.nextInt(0, chromosome.length());
				sequence[i] = random.nextGaussian(sequence[i], step[i]);
				// Check limits
				if (min != null)
				{
					if (sequence[i] < min[i])
						sequence[i] = min[i];
				}
				if (max != null)
				{
					if (sequence[i] > max[i])
						sequence[i] = max[i];
				}
			}
		}

		return chromosome.newChromosome(sequence);
	}
}
