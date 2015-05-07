package gdsc.smlm.ga;

import java.util.Comparator;

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
 * Mutates the sequence by selecting random positions and random shifts.
 */
public class ChromosomeComparator implements Comparator<Chromosome>
{
	/*
	 * (non-Javadoc)
	 * 
	 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
	 */
	@Override
	public int compare(Chromosome chromosome1, Chromosome chromosome2)
	{
		if (chromosome1.getFitness() > chromosome2.getFitness())
			return -1;
		if (chromosome1.getFitness() < chromosome2.getFitness())
			return 1;
		return 0;
	}
}
