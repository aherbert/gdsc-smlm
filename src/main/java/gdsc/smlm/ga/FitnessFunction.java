package gdsc.smlm.ga;

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
 * Calculate the fitness of a chromosome
 */
public interface FitnessFunction
{
	/**
	 * Initialise the fitness function using a population of individuals. This can be used to pre-process the population
	 * before the {@link #fitness(Chromosome)} method is run on each individual.
	 * 
	 * @param individuals The population of individuals that will be assessed
	 */
	void initialise(List<? extends Chromosome> individuals);

	/**
	 * Calculate the fitness
	 * 
	 * @param chromosome
	 * @return the fitness
	 */
	double fitness(Chromosome chromosome);
	
	/**
	 * Shutdown the fitness function. This can be used to post-process the population
	 * after the {@link #fitness(Chromosome)} method is run on each individual.
	 */
	void shutdown();
}
