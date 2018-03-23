package gdsc.smlm.results;

import java.util.Comparator;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Compares the results by frame
 */
public class FramePeakResultComparator implements Comparator<PeakResult>
{
	/** An instance of the comparator */
	public static final FramePeakResultComparator INSTANCE = new FramePeakResultComparator();

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
	 */
	public int compare(PeakResult o1, PeakResult o2)
	{
		// Sort by frame
		return compare(o1.getFrame(), o2.getFrame());
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
