package gdsc.smlm.model;

import org.apache.commons.math3.random.RandomGenerator;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Specify methods for building RandomGenerators.
 */
public interface RandomGeneratorFactory
{
	/**
	 * Create a new random generator
	 * 
	 * @return a new random generator
	 */
	public RandomGenerator createRandomGenerator();
}
