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
public class ChromosomeComparator implements Comparator<Chromosome>
{
	private static class NullComparator implements Comparator<Chromosome>
	{
		public int compare(Chromosome o1, Chromosome o2)
		{
			return 0;
		}
	}

	final Comparator<Chromosome> comparator;

	/**
	 * Instantiates a new chromosome comparator.
	 */
	public ChromosomeComparator()
	{
		this(null);
	}

	/**
	 * Instantiates a new chromosome comparator.
	 *
	 * @param comparator
	 *            the comparator used when the fitness is equal
	 */
	public ChromosomeComparator(Comparator<Chromosome> comparator)
	{
		this.comparator = (comparator == null) ? new NullComparator() : comparator;
	}

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
		return comparator.compare(chromosome1, chromosome2);
	}

	/**
	 * Sort the list (highest fitness first)
	 * 
	 * @param list
	 */
	public static void sort(List<? extends Chromosome> list)
	{
		Collections.sort(list, new ChromosomeComparator());
	}

	/**
	 * Sort the list (highest fitness first), then using the comparator.
	 *
	 * @param list
	 *            the list
	 * @param comparator
	 *            the comparator used when the fitness is equal
	 */
	public static void sort(List<? extends Chromosome> list, Comparator<Chromosome> comparator)
	{
		Collections.sort(list, new ChromosomeComparator(comparator));
	}

	/**
	 * Sort the list (highest fitness first)
	 * 
	 * @param list
	 */
	public static void sort(Chromosome[] list)
	{
		sort(list, list.length);
	}

	/**
	 * Sort the list (highest fitness first)
	 * 
	 * @param list
	 */
	public static void sort(Chromosome[] list, int size)
	{
		Arrays.sort(list, 0, size, new ChromosomeComparator());
	}

	/**
	 * Sort the list (highest fitness first).
	 *
	 * @param list
	 *            the list
	 * @param size
	 *            the size
	 * @param comparator
	 *            the comparator used when the fitness is equal
	 */
	public static void sort(Chromosome[] list, int size, Comparator<Chromosome> comparator)
	{
		Arrays.sort(list, 0, size, new ChromosomeComparator(comparator));
	}
}
