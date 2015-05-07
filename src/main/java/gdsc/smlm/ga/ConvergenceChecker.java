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
 * Defines convergence of a chromosome
 */
public interface ConvergenceChecker
{
	/**
	 * Check if the chromosome has converged
	 * 
	 * @param previous
	 * @param current
	 * @return true if the chromosome has converged
	 */
	boolean converged(Chromosome previous, Chromosome current);
}
