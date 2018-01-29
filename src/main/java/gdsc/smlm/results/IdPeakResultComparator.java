package gdsc.smlm.results;

import java.util.Comparator;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Compares the results by Id
 */
public class IdPeakResultComparator implements Comparator<PeakResult>
{
	/** An instance of the comparator */
	public static final IdPeakResultComparator INSTANCE = new IdPeakResultComparator();

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
	 */b
	public int compare(PeakResult o1, PeakResult o2)
	{
		int id1 = o1.getId();
		int id2 = o2.getId();
		if (id1 < id2)
			return -1;
		if (id1 == id2)
			// Sort by frame
			return compare(o1.getFrame(), o2.getFrame());
		return 1;
	}

	/**
	 * Compare the integers. Copied from Java 1.7.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @return the comparison
	 */
	private static int compare(int x, int y)
	{
		return (x < y) ? -1 : ((x == y) ? 0 : 1);
	}
}
