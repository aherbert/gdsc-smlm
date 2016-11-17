package gdsc.smlm.ga;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Defines the fitness of a chromosome
 */
public class Fitness<T extends Comparable<T>> implements Comparable<Fitness<T>>
{
	final T t;
	final double score;

	/**
	 * Instantiates a new fitness.
	 *
	 * @param t
	 *            the t
	 * @param score
	 *            the score
	 */
	public Fitness(T t, double score)
	{
		this.t = t;
		this.score = score;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(Fitness<T> o)
	{
		if (t == null)
			return 1;
		if (o.t == null)
			return -1;
		return t.compareTo(o.t);
	}
}
