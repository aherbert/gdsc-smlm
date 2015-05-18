package gdsc.smlm.ga;

import java.util.Collections;
import java.util.Comparator;
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
 * Sorts chromosome using the fitness, highest first.
 */
public class ChromosomeComparator implements Comparator<Chromosome>
{
	/*
	 * (non-Javadoc)
	 * 
	 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
	 */
	public int compare(Chromosome chromosome1, Chromosome chromosome2)
	{
		if (chromosome1.getFitness() > chromosome2.getFitness())
			return -1;
		if (chromosome1.getFitness() < chromosome2.getFitness())
			return 1;
		return 0;
	}

	/**
	 * Sort the list (highest first)
	 * 
	 * @param list
	 */
	public static void sort(List<? extends Chromosome> individuals)
	{
		Collections.sort(individuals, new ChromosomeComparator());
	}

	/**
	 * Sort the list (lowest first)
	 * 
	 * @param list
	 */
	public static void sortAscending(List<? extends Chromosome> individuals)
	{
		Collections.sort(individuals, new ChromosomeComparator());
		Collections.reverse(individuals);
	}
}
