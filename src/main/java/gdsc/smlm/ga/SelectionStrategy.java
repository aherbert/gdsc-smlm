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
 * Defines a selection strategy of a population of individuals
 */
public interface SelectionStrategy
{
	/**
	 * Select a subset from the population using the fitness
	 * 
	 * @param individuals
	 * @return a selection of individuals
	 */
	List<? extends Chromosome> select(List<? extends Chromosome> individuals);

	/**
	 * Initialise the selection of pairs for breeding using the fitness
	 * 
	 * @param individuals
	 *            the population of individuals
	 */
	void initialiseBreeding(List<? extends Chromosome> individuals);

	/**
	 * Get the next pair of individuals for breeding. Must be called after {@link #initialiseBreeding(List)}.
	 * 
	 * @return The next pair
	 */
	ChromosomePair next();
	
	/**
	 * Finish selection of pairs for breeding
	 */
	void finishBreeding();
}
