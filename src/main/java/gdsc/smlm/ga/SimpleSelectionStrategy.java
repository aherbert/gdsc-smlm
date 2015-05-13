package gdsc.smlm.ga;

import gdsc.smlm.results.TrackProgress;

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
 * Selects the top individuals
 */
public class SimpleSelectionStrategy extends Randomiser implements SelectionStrategy
{
	final double fraction;

	private List<? extends Chromosome> individuals = null;
	TrackProgress tracker = null;

	/**
	 * @param random
	 * @param fraction
	 *            The fraction of the individuals to select (set between 0 and 1)
	 */
	public SimpleSelectionStrategy(RandomDataGenerator random, double fraction)
	{
		super(random);
		if (fraction > 1)
			fraction = 1;
		this.fraction = fraction;
	}

	/**
	 * Select the top individuals using the configured fraction. The resulting subset will be at least size 2 (unless
	 * the input is smaller or there are not enough valid individuals (fitness above zero)).
	 * 
	 * @param individuals
	 * @return the subset
	 * @see gdsc.smlm.ga.SelectionStrategy#select(java.util.List)
	 */
	@Override
	public List<? extends Chromosome> select(List<? extends Chromosome> individuals)
	{
		if (individuals == null || individuals.size() < 2)
			return individuals;
		ArrayList<Chromosome> subset = new ArrayList<Chromosome>();
		// Add only those with a fitness score
		for (Chromosome c : individuals)
			if (c.getFitness() > 0)
				subset.add(c);
		if (subset.size() < 3)
			return subset;
		if (tracker != null)
			tracker.progress(0.5);
		ChromosomeComparator.sort(subset);
		int size = (int) Math.round(subset.size() * fraction);
		if (size < 2)
			size = 2;
		if (tracker != null)
			tracker.progress(1);
		return subset.subList(0, size);
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
			individuals = null;
		this.individuals = individuals;
	}

	/**
	 * Select pairs randomly from the population
	 * 
	 * @see gdsc.smlm.ga.SelectionStrategy#next()
	 */
	@Override
	public ChromosomePair next()
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
		return new ChromosomePair(individuals.get(first), individuals.get(second));
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
		individuals = null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.SelectionStrategy#setTracker(gdsc.smlm.results.TrackProgress)
	 */
	public void setTracker(TrackProgress tracker)
	{
		this.tracker = tracker;
	}
}
