package gdsc.smlm.ga;

import java.util.Arrays;
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
 * Sorts chromosome using the fitness, highest fitness first.
 */
public class ChromosomeComparator<T extends Comparable<T>> implements Comparator<Chromosome<T>>
{
	/*
	 * (non-Javadoc)
	 * 
	 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
	 */
	public int compare(Chromosome<T> chromosome1, Chromosome<T> chromosome2)
	{
		return chromosome1.getFitness().compareTo(chromosome2.getFitness());
	}

	/**
	 * Sort the list (highest fitness first)
	 * 
	 * @param list
	 */
	public static <T extends Comparable<T>> void sort(List<? extends Chromosome<T>> list)
	{
		Collections.sort(list, new ChromosomeComparator<T>());
	}

	/**
	 * Sort the list (highest fitness first)
	 * 
	 * @param list
	 */
	public static <T extends Comparable<T>> void sort(Chromosome<T>[] list)
	{
		sort(list, 0, list.length);
	}

	/**
	 * Sort the list (highest fitness first)
	 * 
	 * @param list
	 */
	public static <T extends Comparable<T>> void sort(Chromosome<T>[] list, int fromIndex, int toIndex)
	{
		Arrays.sort(list, fromIndex, toIndex, new ChromosomeComparator<T>());
	}
}
