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
 * Exception to throw if the population size is invalid (<2 individuals)
 */
public class InvalidPopulationSize extends RuntimeException
{
	private static final long serialVersionUID = 1L;

	public InvalidPopulationSize(int size, int minSize)
	{
		super("Population size (" + size + ") must be at least " + minSize);
	}
}
