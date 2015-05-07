package gdsc.smlm.ga;

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
			throw new InvalidPopulationSize();
		checkSize(individuals.size());
		this.individuals = individuals;
	}

	private void checkSize(int size)
	{
		if (size < 2)
			throw new InvalidPopulationSize();
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
	 * individuals in the new population is evaluated and convergence checked for the fittest individual.
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
		return current;
	}

	private void grow(SelectionStrategy selectionStrategy, Mutator mutator, Recombiner recombiner)
	{
		if (individuals.size() >= populationSize)
			return;

		selectionStrategy.initialiseBreeding(individuals);
		int fails = 0;
		int target = populationSize - individuals.size();
		ArrayList<Chromosome> newIndividuals = new ArrayList<Chromosome>(target);
		while (newIndividuals.size() < target && fails < failureLimit)
		{
			// Select two individuals for recombination
			ChromosomePair pair = selectionStrategy.next();
			Chromosome[] children = recombiner.cross(pair.c1, pair.c2);
			if (children == null || children.length == 0)
			{
				fails++;
				continue;
			}

			// New children have been generated so mutate them
			for (int i = 0; i < children.length && newIndividuals.size() < target; i++)
			{
				Chromosome c = mutator.mutate(children[i]);
				if (c == null)
				{
					fails++;
					continue;
				}

				// Ignore duplicates
				if (isDuplicate(newIndividuals, c))
				{
					fails++;
					continue;
				}

				newIndividuals.add(c);
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
		for (Chromosome i : this.individuals)
			if (i.distance(c) == 0)
				return true;
		for (Chromosome i : newIndividuals)
			if (i.distance(c) == 0)
				return true;
		return false;
	}

	/**
	 * Calculate the fitness of the population
	 * 
	 * @param fitnessFunction
	 * @return The fittest individual
	 */
	private Chromosome evaluateFitness(FitnessFunction fitnessFunction)
	{
		fitnessFunction.initialise(individuals);
		Chromosome best = null;
		double max = Double.NEGATIVE_INFINITY;
		for (Chromosome c : individuals)
		{
			final double f = fitnessFunction.fitness(c);
			c.setFitness(f);
			if (max < f)
			{
				max = f;
				best = c;
			}
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
}
