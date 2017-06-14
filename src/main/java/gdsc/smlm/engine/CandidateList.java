package gdsc.smlm.engine;

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

import java.util.Arrays;
import java.util.Comparator;

/**
 * Stores a list of candidates
 */
class CandidateList
{
	private static class CandidateComparator implements Comparator<Candidate>
	{
		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
		 */
		public int compare(Candidate o1, Candidate o2)
		{
			return o1.index - o2.index;
		}
	}

	final static CandidateComparator comp;
	static
	{
		comp = new CandidateComparator();
	}

	int size = 0;
	Candidate[] list = null;

	/**
	 * Instantiates a new candidate list.
	 */
	CandidateList()
	{
	}

	/**
	 * Instantiates a new candidate list.
	 *
	 * @param size
	 *            the size
	 * @param list
	 *            the list
	 */
	CandidateList(int size, Candidate[] list)
	{
		this.size = size;
		this.list = list;
	}

	/**
	 * Instantiates a new candidate list.
	 *
	 * @param list
	 *            the list
	 */
	CandidateList(Candidate[] list)
	{
		this.size = list.length;
		this.list = list;
	}

	/**
	 * Add a candidate
	 * 
	 * @param candidate
	 */
	public void add(Candidate candidate)
	{
		if (list == null)
			list = new Candidate[4];
		else if (list.length == size)
		{
			final Candidate[] list2 = new Candidate[size * 2];
			System.arraycopy(list, 0, list2, 0, size);
			list = list2;
		}
		list[size++] = candidate;
	}

	/**
	 * Sort in ascending order of Id
	 */
	public void sort()
	{
		if (size != 0)
			Arrays.sort(list, 0, size, comp);
	}

	/**
	 * Gets the size.
	 *
	 * @return the size
	 */
	public int getSize()
	{
		return size;
	}

	/**
	 * Gets the length of the list. This may be larger than {@link #getSize()}. It is used when the list of candidates
	 * is larger than the max candidate to process.
	 *
	 * @return the length of the list
	 */
	int getlength()
	{
		return (list != null) ? list.length : 0;
	}

	/**
	 * Gets the candidate
	 *
	 * @param index
	 *            the index
	 * @return the candidate
	 */
	public Candidate get(int index)
	{
		return list[index];
	}

	/**
	 * Copy this list.
	 *
	 * @return the new candidate list
	 */
	public CandidateList copy()
	{
		return new CandidateList(size, Arrays.copyOf(list, size));
	}
}