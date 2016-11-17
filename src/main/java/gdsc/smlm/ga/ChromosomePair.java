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
 * Define a pair of chromosomes
 */
public class ChromosomePair<T extends Comparable<T>>
{
	final Chromosome<T> c1, c2;

	/**
	 * Create a pair.
	 *
	 * @param c1
	 *            Chromosome 1
	 * @param c2
	 *            Chromosome 1
	 * @throws IllegalArgumentException
	 *             if either chromosome is null
	 */
	public ChromosomePair(Chromosome<T> c1, Chromosome<T> c2)
	{
		if (c1 == null || c2 == null)
			throw new IllegalArgumentException("Chromosomes must not be null");
		this.c1 = c1;
		this.c2 = c2;
	}
}
