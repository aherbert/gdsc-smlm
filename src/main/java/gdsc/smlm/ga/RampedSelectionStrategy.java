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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.random.RandomDataGenerator;

/**
 * Selects the individuals using a linear ramp from the highest rank to the lowest to set the probability of selection.
 * Selection of breeding couples is uniform from the top n (n is incremented from 1 each selection) crossed with a
 * sample from the entire population using a linear ramped weighting.
 */
public class RampedSelectionStrategy<T extends Comparable<T>> extends SimpleSelectionStrategy<T>
		implements SelectionStrategy<T>
{
	private List<? extends Chromosome<T>> sorted = null;
	private int n;
	private long[] sum;
	private long upper;

	/**
	 * Instantiates a new ramped selection strategy.
	 *
	 * @param random
	 *            the random
	 * @param fraction
	 *            The fraction of the individuals to select (set between 0 and 1)
	 * @param max
	 *            The maximum number of individuals to select
	 */
	public RampedSelectionStrategy(RandomDataGenerator random, double fraction, int max)
	{
		super(random, fraction, max);
	}

	/**
	 * Select the top individual and then the rest using a probability set by their rank. The resulting subset will be
	 * at least size 2 (unless the input is smaller or there are not enough valid individuals (fitness not null)).
	 *
	 * @param individuals
	 * @return the subset
	 * @see gdsc.smlm.ga.SelectionStrategy#select(java.util.List)
	 */
	@Override
	public List<? extends Chromosome<T>> select(List<? extends Chromosome<T>> individuals)
	{
		if (individuals == null || individuals.size() < 3)
			return individuals;

		ArrayList<Chromosome<T>> sorted = new ArrayList<>(individuals.size());
		// Add only those with a fitness score
		for (Chromosome<T> c : individuals)
			if (c.getFitness() != null)
				sorted.add(c);
		if (sorted.size() < 3)
			return sorted;

		// Get the fraction relative to the input list size
		//final int size = getSize(sorted.size());
		final int size = getSize(individuals.size());

		if (tracker != null)
			tracker.progress(0);

		// Sort the list
		ChromosomeComparator.sort(sorted);

		// Create the output subset
		ArrayList<Chromosome<T>> subset = new ArrayList<>(size);

		// Get the number of point available for selection:
		// n in this case is (size-1) since we include the top ranking individual.
		int n = sorted.size() - 1;

		// Add the top individual
		subset.add(sorted.get(0));

		// Get the cumulative total of the rank: total = n(n+1)/2
		long cumulative = (n * (n + 1l)) / 2l;

		// Build an array of rank weighting. The highest ranked starts at n.
		// The first index is ignored since this has been included already.
		int[] rank = new int[sorted.size()];
		for (int i = 1; i < rank.length; i++)
			rank[i] = n--;

		// Now pick chromosomes using the cumulative as the upper limit
		while (subset.size() < size)
		{
			if (tracker != null)
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

		if (tracker != null)
			tracker.progress(1);

		return subset;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.ga.SelectionStrategy#initialiseBreeding(java.util.List)
	 */
	@Override
	public void initialiseBreeding(List<? extends Chromosome<T>> individuals)
	{
		sorted = null;
		if (individuals != null && individuals.size() < 2)
		{
			return;
		}

		// This method may be called when fitness is unknown.
		// Only sort those with a fitness score.
		ArrayList<Chromosome<T>> list = null;
		int count = 0;
		for (Chromosome<T> c : individuals)
			if (c.getFitness() == null)
				count++;

		if (count != 0)
		{
			int toSort = individuals.size() - count;
			if (toSort == 0)
			{
				// No sort possible. Revert to the simple selection strategy
				super.initialiseBreeding(individuals);
				return;
			}
			else
			{
				// A mixed population, some of which we can sort and some we cannot.
				list = new ArrayList<>(toSort);
				ArrayList<Chromosome<T>> subset = new ArrayList<>(count);
				for (Chromosome<T> c : individuals)
				{
					if (c.getFitness() == null)
						subset.add(c);
					else
						list.add(c);
				}

				ChromosomeComparator.sort(list);

				// Create a ramped sum for all those we can sort
				sum = createSum(list.size());
				// Extend the sum linearly for those we cannot sort (i.e. they have the same selection chance)
				sum = extendSum(sum, subset.size());

				list.addAll(subset);
			}
		}
		else
		{
			list = new ArrayList<>(individuals);
			ChromosomeComparator.sort(list);

			// Build a cumulative array of rank weighting. The highest ranked starts at n.
			sum = createSum(list.size());
		}

		this.sorted = list;

		n = 0;

		// Bounds are inclusive so subtract 1
		upper = sum[sum.length - 1] - 1;
	}

	public static long[] createSum(int size)
	{
		long[] sum = new long[size];
		int rank = size;
		sum[0] = rank--;
		for (int i = 1; i < sum.length; i++)
			sum[i] = rank-- + sum[i - 1];
		return sum;
	}

	public static long[] extendSum(long[] sum, int size)
	{
		long[] sum2 = Arrays.copyOf(sum, sum.length + size);
		long s = sum[sum.length - 1];
		for (int i = sum.length; i < sum2.length; i++)
			sum2[i] = ++s;
		return sum2;
	}

	/**
	 * Select pairs randomly from the population. The first is selected from the top n individuals with n starting at 1
	 * and incrementing for each call. The second is selected from the entire population with the weighting equal to
	 * their ranking by fitness.
	 *
	 * @see gdsc.smlm.ga.SelectionStrategy#next()
	 */
	@Override
	public ChromosomePair<T> next()
	{
		if (sorted == null)
			return super.next();
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
		return new ChromosomePair<>(sorted.get(first), sorted.get(second));
	}

	private int nextSample()
	{
		final long next = random.nextLong(0, upper);
		// JUnit test shows that the simple scan through the array is faster for small arrays
		if (sum.length < 100)
			return search(sum, next);
		else
			return binarySearch(sum, next);
	}

	/**
	 * Find the index such that sum[index-1] <= key < sum[index]
	 *
	 * @param sum
	 * @param p
	 * @return the index (or -1)
	 */
	public static int search(final long[] sum, final long key)
	{
		for (int i = 0; i < sum.length; i++)
			if (key < sum[i])
				return i;
		return 0; // This should not happen
	}

	/**
	 * Find the index such that sum[index-1] <= key < sum[index]
	 *
	 * @param sum
	 * @param p
	 * @return the index (or -1)
	 */
	public static int binarySearch(final long[] sum, final long key)
	{
		if (key < sum[0])
			return 0;

		// Adapted from Arrays.binarySearch which
		// finds the actual key value or returns a negative insertion point:

		// If we find the actual key return the next index above.
		// If we do not find the actual key return the insertion point

		int low = 0;
		int high = sum.length - 1;

		while (low <= high)
		{
			int mid = (low + high) >>> 1;
			long midVal = sum[mid];

			if (midVal < key)
				low = mid + 1;
			else if (midVal > key)
				high = mid - 1;
			else
				// Changed from Arrays.binarySearch
				//return mid; // key found
				// If we find the actual key return the next index above.
				return mid + 1;
		}

		// Changed from Arrays.binarySearch
		//return -(low + 1);  // key not found.
		// If we do not find the actual key return the insertion point
		return low;
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
