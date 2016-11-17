package gdsc.smlm.ga;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.random.RandomDataGenerator;

import gdsc.core.logging.TrackProgress;

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
 * Selects the top individuals
 */
public class SimpleSelectionStrategy<T extends Comparable<T>> extends Randomiser implements SelectionStrategy<T>
{
	final double fraction;
	final int max;

	private List<? extends Chromosome<T>> individuals = null;
	TrackProgress tracker = null;

	/**
	 * Instantiates a new simple selection strategy.
	 *
	 * @param random
	 *            the random
	 * @param fraction
	 *            The fraction of the individuals to select (set between 0 and 1)
	 * @param max
	 *            The maximum number of individuals to select
	 */
	public SimpleSelectionStrategy(RandomDataGenerator random, double fraction, int max)
	{
		super(random);
		if (fraction > 1)
			fraction = 1;
		this.fraction = fraction;
		this.max = max;
	}

	/**
	 * Select the top individuals using the configured fraction. The resulting subset will be at least size 2 (unless
	 * the input is smaller or there are not enough valid individuals (fitness not null)).
	 * 
	 * @param individuals
	 * @return the subset
	 * @see gdsc.smlm.ga.SelectionStrategy#select(java.util.List)
	 */
	public List<? extends Chromosome<T>> select(List<? extends Chromosome<T>> individuals)
	{
		if (individuals == null || individuals.size() < 2)
			return individuals;
		@SuppressWarnings("unchecked")
		Chromosome<T>[] subset = new Chromosome[individuals.size()];
		int size = 0;
		// Add only those with a fitness score
		for (Chromosome<T> c : individuals)
			if (c.getFitness() != null)
			{
				subset[size++] = c;
			}
		if (size < 3)
			return Arrays.asList(Arrays.copyOf(subset, size));
		if (tracker != null)
			tracker.progress(0.5);
		ChromosomeComparator.sort(subset, 0, size);
		
		// Get the fraction relative to the input list size
		//size = getSize(size);
		size = Math.min(size, getSize(individuals.size()));
		if (tracker != null)
			tracker.progress(1);
		return Arrays.asList(Arrays.copyOf(subset, size));
	}

	/**
	 * Calculate the new size of the population after selection
	 * 
	 * @param size
	 *            The current size of the population before selection
	 * @return The new size of the population
	 */
	protected int getSize(int size)
	{
		// Get the size using the fraction
		size = (int) Math.round(size * fraction);
		// Check against the max number to select
		if (max > 2 && size > max)
			size = max;
		// Check the size is at least 2
		if (size < 2)
			size = 2;
		return size;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.SelectionStrategy#initialiseBreeding(java.util.List)
	 */
	public void initialiseBreeding(List<? extends Chromosome<T>> individuals)
	{
		if (individuals != null && individuals.size() < 2)
			individuals = null;
		this.individuals = individuals;
	}

	/**
	 * Select pairs randomly from the population
	 * 
	 * @see gdsc.smlm.ga.SelectionStrategy#next()
	 */
	public ChromosomePair<T> next()
	{
		if (individuals == null)
			return null;
		int first, second;
		if (individuals.size() == 2)
		{
			first = 0;
			second = 1;
		}
		else
		{
			// Bounds are inclusive so subtract 1
			final int upper = individuals.size() - 1;
			first = random.nextInt(0, upper);
			second = random.nextInt(0, upper);
			// Avoid crossover with the same parent
			while (second == first)
				second = random.nextInt(0, upper);
		}
		return new ChromosomePair<T>(individuals.get(first), individuals.get(second));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.SelectionStrategy#finishBreeding()
	 */
	public void finishBreeding()
	{
		// Free memory
		individuals = null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.SelectionStrategy#setTracker(gdsc.core.logging.TrackProgress)
	 */
	public void setTracker(TrackProgress tracker)
	{
		this.tracker = tracker;
	}
}
