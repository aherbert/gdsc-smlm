package gdsc.smlm.ga;

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
 * Defines recombination crossover of a chromosome pair
 */
public interface Recombiner<T extends Comparable<T>>
{
	/**
	 * Crossover the provided chromosomes to produce one or more new sequences
	 * 
	 * @param chromosome1
	 * @param chromosome2
	 * @return one or more new sequences
	 */
	Chromosome<T>[] cross(Chromosome<T> chromosome1, Chromosome<T> chromosome2);
}
