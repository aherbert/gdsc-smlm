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
import java.util.List;

import gdsc.core.logging.TrackProgress;

/**
 * Contains a population of individuals that may crossover and mutate to evolve.
 * <p>
 * For simplicity the individuals have one chromosome sequence.
 *
 * @param <T>
 *            the generic type
 */
/*
 * An extension would be to have an Individual class that has many Chromosomes, each is allowed to crossover with its
 * matching pair and then segregation occurs to new individuals.
 */
public class Population<T extends Comparable<T>>
{
	private List<? extends Chromosome<T>> individuals;
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
	public Population(List<? extends Chromosome<T>> individuals)
	{
		if (individuals == null)
			throw new InvalidPopulationSize(0, 1);
		checkSize(individuals.size());
		this.individuals = individuals;
	}

	private static void checkSize(int size)
	{
		if (size < 1)
			throw new InvalidPopulationSize(0, 1);
	}

	/**
	 * @return the individuals
	 */
	public List<? extends Chromosome<T>> getIndividuals()
	{
		return individuals;
	}

	/**
	 * Evolve the population of individuals until convergence of the most fit individual in the population.
	 * <p>
	 * The population will grow until the desired population size by recombination of individual pairs chosen from the
	 * population by the selection strategy. Child sequences will be subject to mutation. The fitness of all the
	 * individuals in the new population is evaluated and convergence checked for the fittest individual. If the initial
	 * population is small (<2 or <Chromosome<T>.length()) then mutation will be used to expand it before recombination.
	 * <p>
	 * The process of grow, evaluate, select is repeated until convergence.
	 * <p>
	 * Note: the subset of individuals selected for the next generation by the selection strategy will be unchanged
	 * (i.e. no mutation). This allows the fittest individuals to remain unchanged.
	 *
	 * @param mutator
	 *            the mutator
	 * @param recombiner
	 *            the recombiner
	 * @param fitnessFunction
	 *            the fitness function
	 * @param selectionStrategy
	 *            the selection strategy
	 * @param checker
	 *            the checker
	 * @return The best individual
	 * @throws InvalidPopulationSize
	 *             if the population is less than 2 (this can occur after selection)
	 */
	public Chromosome<T> evolve(Mutator<T> mutator, Recombiner<T> recombiner, FitnessFunction<T> fitnessFunction,
			SelectionStrategy<T> selectionStrategy, ConvergenceChecker<T> checker)
	{
		// Reset the fitness
		for (Chromosome<T> c : individuals)
			c.setFitness(null);

		// Find the best individual
		grow(selectionStrategy, mutator, recombiner);
		Chromosome<T> current = evaluateFitness(fitnessFunction);
		Chromosome<T> previous;

		boolean converged = false;
		while (!converged)
		{
			previous = current;
			// Select the best individuals and expand the population
			if (!select(selectionStrategy))
			{
				current = null;
				break;
			}
			grow(selectionStrategy, mutator, recombiner);
			// Evaluate the fitness and check convergence
			current = evaluateFitness(fitnessFunction);
			converged = checker.converged(previous, current);
		}
		if (tracker != null)
			tracker.status("Converged [%d]", iteration);
		return current;
	}

	private void grow(SelectionStrategy<T> selectionStrategy, Mutator<T> mutator, Recombiner<T> recombiner)
	{
		iteration++;
		start("Grow");

		if (individuals.size() >= populationSize)
			return;

		ArrayList<Chromosome<T>> newIndividuals = new ArrayList<>(populationSize - individuals.size());

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
				Chromosome<T> c = mutator.mutate(individuals.get(next++ % individuals.size()));
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
				end();
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
			ChromosomePair<T> pair = selectionStrategy.next();
			Chromosome<T>[] children = recombiner.cross(pair.c1, pair.c2);
			if (children != null && children.length != 0)
			{
				// New children have been generated so mutate them
				for (int i = 0; i < children.length && newIndividuals.size() < target; i++)
				{
					Chromosome<T> c = mutator.mutate(children[i]);
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

		end();
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
	private boolean isDuplicate(ArrayList<? extends Chromosome<T>> newIndividuals, Chromosome<T> c)
	{
		final double[] s = c.sequence();
		for (Chromosome<T> i : this.individuals)
			if (match(i, s))
				return true;
		for (Chromosome<T> i : newIndividuals)
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
	private boolean match(Chromosome<T> c, double[] s)
	{
		final double[] s2 = c.sequence();
		for (int i = 0; i < s.length; i++)
			if (s[i] != s2[i])
				return false;
		return true;
	}

	/**
	 * Calculate the fitness of the population.
	 *
	 * @param fitnessFunction
	 *            the fitness function
	 * @return The fittest individual
	 */
	private Chromosome<T> evaluateFitness(FitnessFunction<T> fitnessFunction)
	{
		start("Score");

		Chromosome<T> best = null;
		T max = null;

		// Subset only those with no fitness score (the others must be unchanged)
		ArrayList<Chromosome<T>> subset = new ArrayList<>(individuals.size());
		long count = 0;
		for (Chromosome<T> c : individuals)
		{
			final T f = c.getFitness();
			if (f == null)
			{
				subset.add(c);
			}
			else
			{
				if (tracker != null)
					tracker.progress(++count, individuals.size());
				if (f.compareTo(max) < 0)
				{
					max = f;
					best = c;
				}
			}
		}

		fitnessFunction.initialise(subset);
		for (Chromosome<T> c : subset)
		{
			final T f = fitnessFunction.fitness(c);
			c.setFitness(f);
			if (f != null && f.compareTo(max) < 0)
			{
				max = f;
				best = c;
			}
			if (tracker != null)
				tracker.progress(++count, individuals.size());
		}
		fitnessFunction.shutdown();

		end();

		return best;
	}

	/**
	 * Select a subset of the population
	 *
	 * @param selection
	 *            The selection strategy
	 * @return True if a valid population was selected (size>=1)
	 */
	private boolean select(SelectionStrategy<T> selection)
	{
		start("Select");
		individuals = selection.select(individuals);
		end();
		return !individuals.isEmpty();
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

	private void start(String stage)
	{
		if (tracker != null)
		{
			tracker.status(stage + " [%d]", iteration);
			tracker.progress(0);
		}
	}

	private void end()
	{
		if (tracker != null)
			tracker.progress(1);
	}
}
