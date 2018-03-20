package gdsc.smlm.results;

import java.util.Comparator;

import org.apache.commons.math3.random.RandomGenerator;

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
 * Stores peak results with list access.
 */
public interface PeakResultStoreList extends PeakResultStore
{
	/**
	 * Gets the result.
	 *
	 * @param index
	 *            the index
	 * @return the peak result
	 */
	public PeakResult get(int index);

	/**
	 * Removes the result.
	 *
	 * @param index the index
	 * @return the peak result removed
     * @throws IndexOutOfBoundsException If the index is invalid
	 */
	public PeakResult remove(int index);
	
	/**
	 * Sort the results.
	 */
	public void sort();

	/**
	 * Sort the results.
	 *
	 * @param comparator
	 *            the comparator
	 */
	public void sort(Comparator<PeakResult> comparator);

	/**
	 * Shuffle the results.
	 *
	 * @param randomGenerator
	 *            the random generator
	 */
	public void shuffle(final RandomGenerator randomGenerator);

	/**
	 * Returns the index of the first occurrence of the specified result
	 * in this store, or -1 if this list does not contain the element.
	 * More formally, returns the lowest index <tt>i</tt> such that
	 * <tt>(result==null&nbsp;?&nbsp;get(i)==null&nbsp;:&nbsp;result.equals(get(i)))</tt>,
	 * or -1 if there is no such index.
	 *
	 * @param result
	 *            the result
	 * @return the index (or -1)
	 */
	public int indexOf(PeakResult result);

	/**
	 * Returns the index of the last occurrence of the specified result
	 * in this store, or -1 if this list does not contain the element.
	 * More formally, returns the highest index <tt>i</tt> such that
	 * <tt>(result==null&nbsp;?&nbsp;get(i)==null&nbsp;:&nbsp;result.equals(get(i)))</tt>,
	 * or -1 if there is no such index.
	 *
	 * @param result
	 *            the result
	 * @return the index (or -1)
	 */
	public int lastIndexOf(PeakResult result);
}
