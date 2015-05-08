package gdsc.smlm.ga;

import java.util.ArrayList;
import java.util.List;

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
 * Selects the individuals using a linear ramp from the highest rank to the lowest to set the probability of selection.
 * Selection of breeding couples is uniformly random as per the SimpleSelectionStrategy.
 */
public class RampedSelectionStrategy extends SimpleSelectionStrategy implements SelectionStrategy
{
	/**
	 * @param random
	 * @param fraction
	 *            The fraction of the individuals to select (set between 0 and 1)
	 */
	public RampedSelectionStrategy(RandomDataGenerator random, double fraction)
	{
		super(random, fraction);
	}

	/**
	 * Select the top individual and then the rest using a probability set by their rank. The resulting subset will be
	 * at least size 2 (unless the input is smaller).
	 * 
	 * @param individuals
	 * @return the subset
	 * @see gdsc.smlm.ga.SelectionStrategy#select(java.util.List)
	 */
	@Override
	public List<? extends Chromosome> select(List<? extends Chromosome> individuals)
	{
		if (individuals == null || individuals.size() < 3)
			return individuals;
		
		// Sort the list
		ArrayList<Chromosome> sorted = new ArrayList<Chromosome>();
		sorted.addAll(individuals);
		ChromosomeComparator.sort(sorted);
		int size = (int) Math.round(sorted.size() * fraction);
		if (size < 2)
			size = 2;

		// Create the output subset
		ArrayList<Chromosome> subset = new ArrayList<Chromosome>(size);
		
		// Get the number of point available for selection:
		// n in this case is (size-1) since we include the top ranking individual.
		int n = sorted.size() - 1;
		
		// Add the top individual
		subset.add(sorted.get(0));

		// Get the cumulative total of the rank: total = n(n+1)/2
		long cumulative = ((long) n * (n + 1l)) / 2l;

		// Build an array of rank weighting. The highest ranked starts at n.
		// The first index is ignored since this has been included already.
		int[] rank = new int[sorted.size()];
		for (int i = 1; i < rank.length; i++)
			rank[i] = n--;

		// Now pick chromosomes using the cumulative as the upper limit
		while (subset.size() < size)
		{
			// Used to check we pick something
			final long previous = cumulative;

			// Generate a random positive within the cumulative range - note the end points are inclusive
			long next = random.nextLong(1l, cumulative);
			// Find the random position
			long sum = 0;
			for (int i = 1; i < rank.length; i++)
			{
				sum += rank[i];
				if (next <= sum)
				{
					// Pick this point
					subset.add(sorted.get(i));
					// Update the cumulative then eliminate from future selection  
					cumulative -= rank[i];
					rank[i] = 0;
					break;
				}
			}

			// Check we chose something 
			if (previous == cumulative)
				throw new RuntimeException("Failed to select a candidate. Size = " + subset.size() + " / " + size);
		}

		return subset;
	}
}
