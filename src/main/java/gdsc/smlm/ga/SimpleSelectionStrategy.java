package gdsc.smlm.ga;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
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
	 * the input is smaller).
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
		subset.addAll(individuals);
		ChromosomeComparator.sort(subset);
		int size = (int) Math.round(subset.size() * fraction);
		if (size < 2)
			size = 2;
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
		if (individuals != null && individuals.size() == 1)
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
			first = random.nextInt(0, individuals.size());
			second = random.nextInt(0, individuals.size());
			while (second != first)
				second = random.nextInt(0, individuals.size());
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
}
