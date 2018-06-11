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

import java.util.List;


/**
 * Calculate the fitness of a chromosome
 */
public interface FitnessFunction<T extends Comparable<T>>
{
	/**
	 * Initialise the fitness function using a population of individuals. This can be used to pre-process the population
	 * before the {@link #fitness(Chromosome)} method is run on each individual.
	 * 
	 * @param individuals The population of individuals that will be assessed
	 */
	void initialise(List<? extends Chromosome<T>> individuals);

	/**
	 * Calculate the fitness
	 * 
	 * @param chromosome
	 * @return the fitness
	 */
	T fitness(Chromosome<T> chromosome);
	
	/**
	 * Shutdown the fitness function. This can be used to post-process the population
	 * after the {@link #fitness(Chromosome)} method is run on each individual.
	 */
	void shutdown();
}
