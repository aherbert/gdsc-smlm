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
 * Define the genetic sequence that can evolve
 */
public interface Chromosome
{
	/**
	 * Get the chromosome length
	 * 
	 * @return The chromosome length
	 */
	int length();

	/**
	 * Get the chromosome sequence
	 * 
	 * @return the chromosome sequence (must equal the length)
	 */
	double[] sequence();

	/**
	 * Create a new chromosome
	 * 
	 * @param sequence
	 *            the chromosome sequence (must equal the current length)
	 * @return A new chromosome with the given sequence
	 */
	Chromosome newChromosome(double[] sequence);

	/**
	 * Get the range for mutation at each position in the sequence. This defines how far each position in the sequence
	 * can mutate in a single step.
	 * 
	 * @return The range for mutation at each position in the sequence (must equal length)
	 */
	double[] mutationStepRange();

	/**
	 * Get the lower limit at each position in the sequence. It is valid to return negative infinity for any position or
	 * null for no limit.
	 * 
	 * @return The lower limit for each position in the sequence (must equal length)
	 */
	double[] lowerLimit();

	/**
	 * Get the upper limit at each position in the sequence. It is valid to return positive infinity for any position or
	 * null for no limit.
	 * 
	 * @return The upper limit for each position in the sequence (must equal length)
	 */
	double[] upperLimit();

	// Note: Default implementation of the getter/setter to store the double would require using Java 8. 

	/**
	 * Set the fitness
	 * 
	 * @param fitness
	 *            The fitness of the sequence
	 */
	void setFitness(double fitness);

	/**
	 * Get the fitness
	 * 
	 * @return The fitness of the sequence
	 */
	double getFitness();

	/**
	 * Calculate the distance to another chromosome
	 * @param other the other chromosome
	 * @return the distance (zero is a match)
	 */
	double distance(Chromosome other);

	/**
	 * Calculate if equal to another chromosome
	 * @param other the other chromosome
	 * @return true if the same
	 */
	boolean equals(Chromosome other);
}
