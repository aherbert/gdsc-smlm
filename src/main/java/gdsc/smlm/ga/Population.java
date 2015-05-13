package gdsc.smlm.ga;

import gdsc.smlm.results.TrackProgress;

import java.util.ArrayList;
import java.util.List;

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
 * Contains a population of individuals that may crossover and mutate to evolve.
 * <p>
 * For simplicity the individuals have one chromosome sequence.
 */
/*
 * An extension would be to have an Individual class that has many Chromosomes, each is allowed to crossover with its
 * matching pair and then segregation occurs to new individuals.
 */
public class Population
{
	private List<? extends Chromosome> individuals;
	private int populationSize = 500;
	private int failureLimit = 3;
	private int iteration = 0;
	// This introduces a dependency on another gdsc.smlm package
	private TrackProgress tracker = null;

	/**
	 * Create a population of individuals
	 * 
	 * @param individuals
	 *            The population of individuals
	 * @throws InvalidPopulationSize
	 *             if the population is less than 2
	 */
	public Population(List<? extends Chromosome> individuals)
	{
		if (individuals == null)
			throw new InvalidPopulationSize(0, 1);
		checkSize(individuals.size());
		this.individuals = individuals;
	}

	private void checkSize(int size)
	{
		if (size < 1)
			throw new InvalidPopulationSize(0, 1);
	}

	/**
	 * @return the individuals
	 */
	public List<? extends Chromosome> getIndividuals()
	{
		return individuals;
	}

	/**
	 * Evolve the population of individuals until convergence of the most fit individual in the population.
	 * <p>
	 * The population will grow until the desired population size by recombination of individual pairs chosen from the
	 * population by the selection strategy. Child sequences will be subject to mutation. The fitness of all the
	 * individuals in the new population is evaluated and convergence checked for the fittest individual. If the initial
	 * population is small (<2 or <Chromosome.length()) then mutation will be used to expand it before recombination.
	 * <p>
	 * The process of grow, evaluate, select is repeated until convergence.
	 * <p>
	 * Note: the subset of individuals selected for the next generation by the selection strategy will be unchanged
	 * (i.e. no mutation). This allows the fittest individuals to remain unchanged.
	 * 
	 * @param mutator
	 * @param recombiner
	 * @param checker
	 * @throws InvalidPopulationSize
	 *             if the population is less than 2 (this can occur after selection)
	 * @return The best individual
	 */
	public Chromosome evolve(Mutator mutator, Recombiner recombiner, FitnessFunction fitnessFunction,
			SelectionStrategy selectionStrategy, ConvergenceChecker checker)
	{
		// Find the best individual
		grow(selectionStrategy, mutator, recombiner);
		Chromosome current = evaluateFitness(fitnessFunction);
		Chromosome previous;

		boolean converged = false;
		while (!converged)
		{
			previous = current;
			// Select the best individuals and expand the population
			select(selectionStrategy);
			grow(selectionStrategy, mutator, recombiner);
			// Evaluate the fitness and check convergence
			current = evaluateFitness(fitnessFunction);
			converged = checker.converged(previous, current);
		}
		if (tracker != null)
			tracker.status("Converged [%d]", iteration);
		return current;
	}

	private void grow(SelectionStrategy selectionStrategy, Mutator mutator, Recombiner recombiner)
	{
		iteration++;
		if (tracker != null)
			tracker.status("Grow [%d]", iteration);

		if (individuals.size() >= populationSize)
			return;

		ArrayList<Chromosome> newIndividuals = new ArrayList<Chromosome>(populationSize - individuals.size());

		// Check for a minimum population size & mutate the individuals to achieve it.
		// This allows a seed population of 1 to evolve.
		final int minSize = Math.max(2, individuals.get(0).length());
		int target = minSize - individuals.size();
		if (target > 0)
		{
			if (tracker != null)
				tracker.progress(individuals.size(), populationSize);

			int next = 0;
			int fails = 0;
			while (newIndividuals.size() < target && fails < failureLimit)
			{
				Chromosome c = mutator.mutate(individuals.get(next++ % individuals.size()));
				if (c != null && !isDuplicate(newIndividuals, c))
				{
					newIndividuals.add(c);
					fails = 0;
					if (tracker != null)
						tracker.progress(newIndividuals.size() + individuals.size(), populationSize);
				}
				else
				{
					fails++;
				}
			}

			// Combine the lists
			newIndividuals.addAll(individuals);
			individuals = newIndividuals;

			if (individuals.size() < 2)
			{
				if (tracker != null)
					tracker.progress(1);
				return; // Failed to mutate anything to achieve a breeding population
			}
		}

		// Now breed the population
		selectionStrategy.initialiseBreeding(individuals);
		target = populationSize - individuals.size();
		int previousSize = -1;
		int fails = 0;
		while (newIndividuals.size() < target && fails < failureLimit)
		{
			previousSize = newIndividuals.size();

			// Select two individuals for recombination
			ChromosomePair pair = selectionStrategy.next();
			Chromosome[] children = recombiner.cross(pair.c1, pair.c2);
			if (children != null && children.length != 0)
			{
				// New children have been generated so mutate them
				for (int i = 0; i < children.length && newIndividuals.size() < target; i++)
				{
					Chromosome c = mutator.mutate(children[i]);
					if (c == null)
						continue;

					// Ignore duplicates
					if (isDuplicate(newIndividuals, c))
						continue;

					newIndividuals.add(c);
				}
			}

			if (previousSize == newIndividuals.size())
				fails++;
			else
			{
				fails = 0;
				if (tracker != null)
					tracker.progress(newIndividuals.size() + individuals.size(), populationSize);
			}
		}
		selectionStrategy.finishBreeding();

		// Combine the lists
		newIndividuals.addAll(individuals);
		individuals = newIndividuals;
	}

	/**
	 * Check for duplicates in the current and new populations
	 * 
	 * @param newIndividuals
	 *            The new population
	 * @param c
	 *            The chromosome
	 * @return true if a duplicate
	 */
	private boolean isDuplicate(ArrayList<? extends Chromosome> newIndividuals, Chromosome c)
	{
		final double[] s = c.sequence();
		for (Chromosome i : this.individuals)
			if (match(i, s))
				return true;
		for (Chromosome i : newIndividuals)
			if (match(i, s))
				return true;
		return false;
	}

	/**
	 * Check if a chromosome matches the sequence
	 * 
	 * @param c
	 *            The chromosome
	 * @param s
	 *            The sequence
	 * @return True if a match
	 */
	private boolean match(Chromosome c, double[] s)
	{
		final double[] s2 = c.sequence();
		for (int i = 0; i < s.length; i++)
			if (s[i] != s2[i])
				return false;
		return true;
	}

	/**
	 * Calculate the fitness of the population
	 * 
	 * @param fitnessFunction
	 * @return The fittest individual
	 */
	private Chromosome evaluateFitness(FitnessFunction fitnessFunction)
	{
		if (tracker != null)
			tracker.status("Score [%d]", iteration);

		Chromosome best = null;
		double max = Double.NEGATIVE_INFINITY;

		// Subset only those with no fitness score (the others must be unchanged)
		ArrayList<Chromosome> subset = new ArrayList<Chromosome>(individuals.size());
		long count = 0;
		for (Chromosome c : individuals)
		{
			final double f = c.getFitness();
			if (f == 0)
			{
				subset.add(c);
			}
			else
			{
				if (tracker != null)
					tracker.progress(++count, individuals.size());
				if (max < f)
				{
					max = f;
					best = c;
				}
			}
		}

		fitnessFunction.initialise(subset);
		for (Chromosome c : subset)
		{
			final double f = fitnessFunction.fitness(c);
			c.setFitness(f);
			if (max < f)
			{
				max = f;
				best = c;
			}
			if (tracker != null)
				tracker.progress(++count, individuals.size());
		}
		fitnessFunction.shutdown();

		return best;
	}

	/**
	 * Select a subset of the population
	 * 
	 * @param selection
	 *            The selection strategy
	 */
	private void select(SelectionStrategy selection)
	{
		if (tracker != null)
		{
			tracker.status("Select [%d]", iteration);
			//selection.setTracker(tracker);
		}
		individuals = selection.select(individuals);
		checkSize(individuals.size());
	}

	/**
	 * Get the population size limit to achieve when growing the population
	 * 
	 * @return the populationSize
	 */
	public int getPopulationSize()
	{
		return populationSize;
	}

	/**
	 * Set the population size limit to achieve when growing the population
	 * 
	 * @param populationSize
	 *            the population size to set
	 */
	public void setPopulationSize(int populationSize)
	{
		checkSize(populationSize);
		this.populationSize = populationSize;
	}

	/**
	 * Get the number of failed recombinations/mutations to allow before the stopping attempts to grow the population
	 * 
	 * @return the failure limit
	 */
	public int getFailureLimit()
	{
		return failureLimit;
	}

	/**
	 * Set the number of failed recombinations/mutations to allow before the stopping attempts to grow the population
	 * 
	 * @param failureLimit
	 *            the failure limit
	 */
	public void setFailureLimit(int failureLimit)
	{
		this.failureLimit = failureLimit;
	}

	/**
	 * @return the tracker
	 */
	public TrackProgress getTracker()
	{
		return tracker;
	}

	/**
	 * Set a tracker to allow the progress to be followed
	 * 
	 * @param tracker
	 *            the tracker to set
	 */
	public void setTracker(TrackProgress tracker)
	{
		this.tracker = tracker;
	}

	/**
	 * Get the iteration. The iteration is increased each time the population grows as part of the [grow, evaluate,
	 * select] cycle.
	 * 
	 * @return the iteration
	 */
	public int getIteration()
	{
		return iteration;
	}
}
