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
			final double[] step = chromosome.mutationStepRange();
			final double[] min = chromosome.lowerLimit();
			final double[] max = chromosome.upperLimit();

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

		return chromosome.create(sequence);
	}
}
