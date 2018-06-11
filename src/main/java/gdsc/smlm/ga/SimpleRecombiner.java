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
package gdsc.smlm.ga;

import java.util.Arrays;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;

import gdsc.core.utils.SimpleArrayUtils;


/**
 * Recombine sequence by selecting random positions for crossover.
 */
public class SimpleRecombiner<T extends Comparable<T>> extends Randomiser implements Recombiner<T>
{
	final double fraction;
	final double meanChildren;

	/**
	 * @param random
	 * @param fraction
	 *            The fraction of the sequence positions to recombine on average
	 * @param meanChildren
	 *            The mean number of additional children
	 */
	public SimpleRecombiner(RandomDataGenerator random, double fraction, double meanChildren)
	{
		super(random);
		this.fraction = fraction;
		this.meanChildren = meanChildren;
	}

	/**
	 * Crossover the chromosome to form a new sequences.
	 * <p>
	 * The number of children are chosen from a Poisson distribution with an average using the mean children. A minimum
	 * of 1 child is selected to ensure crossover.
	 * <p>
	 * The number of positions are chosen from a Poisson distribution with an average using a fraction of the total
	 * positions. A minimum of 1 crossover position is selected to ensure crossover.
	 * <p>
	 * The positions are then chosen randomly and the new chromosome generated.
	 * 
	 * @see gdsc.smlm.ga.Recombiner#cross(gdsc.smlm.ga.Chromosome<T>, gdsc.smlm.ga.Chromosome<T>)
	 */
	public Chromosome<T>[] cross(Chromosome<T> chromosome1, Chromosome<T> chromosome2)
	{
		int nChildren = 1;
		if (meanChildren > 0)
			nChildren = Math.max(1, (int) random.nextPoisson(meanChildren));

		@SuppressWarnings("unchecked")
		Chromosome<T>[] children = new Chromosome[nChildren];
		int count = 0;
		double[] s1 = chromosome1.sequence();
		double[] s2 = chromosome2.sequence();
		while (count < nChildren)
		{
			ChromosomePair<T> pair = recombine(chromosome1, chromosome2, s1, s2);
			children[count++] = pair.c1;
			if (count == nChildren)
				break;
			children[count++] = pair.c2;
		}

		return children;
	}

	private ChromosomePair<T> recombine(Chromosome<T> chromosome1, Chromosome<T> chromosome2, double[] s1, double[] s2)
	{
		int nCrossovers = 1;
		if (fraction > 0)
			nCrossovers = Math.max(1, (int) random.nextPoisson(fraction * chromosome1.length()));

		// Randomly select positions using a partial Fischer-Yates shuffle
		int[] positions = new int[s1.length];
		for (int i = 0; i < positions.length; i++)
			positions[i] = i;

		RandomGenerator ran = random.getRandomGenerator();
		for (int i = positions.length, n = nCrossovers; i-- > 1 && n-- > 0;)
		{
			int j = (int) (ran.nextInt(i + 1));
			int tmp = positions[i];
			positions[i] = positions[j];
			positions[j] = tmp;
		}

		// Reverse the array because the end is random
		SimpleArrayUtils.reverse(positions);
		positions = Arrays.copyOf(positions, nCrossovers);
		// Get the positions in order
		if (nCrossovers != 1)
			Arrays.sort(positions);

		int nextSwap = 0;

		// Create the children by copying the parent, swapping at each crossover position
		double[] n1 = new double[s1.length];
		double[] n2 = new double[n1.length];
		for (int i = 0; i < n1.length; i++)
		{
			if (positions[nextSwap] == i)
			{
				double[] tmp = s1;
				s1 = s2;
				s2 = tmp;
				nextSwap++;
				// Avoid index out of bounds
				if (nextSwap == nCrossovers)
					nextSwap--;
			}
			n1[i] = s1[i];
			n2[i] = s2[i];
		}

		// Create the new chromosome using the correct parent, i.e.
		// If the first swap position was at the start then reverse them.
		Chromosome<T> c1 = (positions[0] == 0) ? chromosome2 : chromosome1;
		Chromosome<T> c2 = (positions[0] == 0) ? chromosome1 : chromosome2;
		c1 = c1.newChromosome(n1);
		c2 = c2.newChromosome(n2);

		// Ensure the child order is random
		return (ran.nextDouble() < 0.5) ? new ChromosomePair<T>(c1, c2) : new ChromosomePair<T>(c2, c1);
	}
}
