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
 * Selection of breeding couples is uniform from the top n (n is incremented from 1 each selection) crossed with a
 * sample from the entire population using a linear ramped weighting.
 */
public class RampedSelectionStrategy extends SimpleSelectionStrategy implements SelectionStrategy
{
	private List<? extends Chromosome> sorted = null;
	private int n;
	private long[] sum;
	private long upper;

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
	 * at least size 2 (unless the input is smaller or there are not enough valid individuals (fitness above zero)).
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

		ArrayList<Chromosome> sorted = new ArrayList<Chromosome>(individuals.size());
		// Add only those with a fitness score
		for (Chromosome c : individuals)
			if (c.getFitness() > 0)
				sorted.add(c);
		if (sorted.size() < 3)
			return sorted;

		int size = (int) Math.round(sorted.size() * fraction);
		if (size < 2)
			size = 2;
		
		if (tracker!= null)
			tracker.progress(0, size);
		
		// Sort the list
		ChromosomeComparator.sort(sorted);

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
			if (tracker!= null)
				tracker.progress(subset.size(), size);
			
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
		
		if (tracker!= null)
			tracker.progress(1);

		return subset;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.SelectionStrategy#initialiseBreeding(java.util.List)
	 */
	@Override
	public void initialiseBreeding(List<? extends Chromosome> individuals)
	{
		if (individuals != null && individuals.size() < 2)
		{
			individuals = null;
			return;
		}
		this.sorted = new ArrayList<Chromosome>(individuals);
		ChromosomeComparator.sort(this.sorted);
		n = 0;

		// Build a cumulative array of rank weighting. The highest ranked starts at n.
		sum = new long[sorted.size()];
		int rank = sorted.size();
		sum[0] = rank--;
		for (int i = 1; i < sum.length; i++)
			sum[i] = rank-- + sum[i - 1];

		// Bounds are inclusive so subtract 1
		upper = sum[sum.length - 1] - 1;
	}

	/**
	 * Select pairs randomly from the population. The first is selected from the top n individuals with n starting at 1
	 * and incrementing for each call. The second is selected from the entire population with the weighting equal to
	 * their ranking by fitness.
	 * 
	 * @see gdsc.smlm.ga.SelectionStrategy#next()
	 */
	@Override
	public ChromosomePair next()
	{
		if (sorted == null)
			return null;
		int first, second;
		if (sorted.size() == 2)
		{
			first = 0;
			second = 1;
		}
		else
		{
			// Pick from the first n
			if (n == 0)
			{
				// This is the first call so select the top individual
				first = 0;
			}
			else
			{
				// Restrict the upper limit to the population range
				first = random.nextInt(0, Math.min(n, sorted.size() - 1));
			}
			n++;

			// Select the second from the ramped cumulative
			second = nextSample();
			// Avoid crossover with the same parent
			while (second == first)
				second = nextSample();
		}
		//System.out.printf("Next [%d] %d x %d\n", n, first, second);
		return new ChromosomePair(sorted.get(first), sorted.get(second));
	}

	private int nextSample()
	{
		final long next = random.nextLong(0, upper);
		for (int i = 0; i < sum.length; i++)
			if (next < sum[i])
				return i;
		return 0; // This should not happen
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.SelectionStrategy#finishBreeding()
	 */
	@Override
	public void finishBreeding()
	{
		// Free memory
		sorted = null;
		sum = null;
	}
}
