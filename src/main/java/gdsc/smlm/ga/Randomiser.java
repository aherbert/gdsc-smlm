package gdsc.smlm.ga;

import org.apache.commons.math3.random.RandomDataGenerator;

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
 * Base class for data generation using randomness
 */
public abstract class Randomiser
{
	final RandomDataGenerator random;

	/**
	 * @param random
	 *            the random generator
	 * @throws IllegalArgumentException
	 *             if the random generator is null
	 */
	public Randomiser(RandomDataGenerator random)
	{
		if (random == null)
			throw new IllegalArgumentException("Random data generator cannot be null");
		this.random = random;
	}
}
