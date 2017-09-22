package gdsc.smlm.results;

import java.util.Collection;
import java.util.Comparator;

import org.apache.commons.math3.random.RandomGenerator;

import gdsc.smlm.results.predicates.PeakResultPredicate;
import gdsc.smlm.results.procedures.PeakResultProcedure;

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
 * Stores peak results.
 */
public interface PeakResultStore
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
	 * Get the size.
	 *
	 * @return the size
	 */
	public int size();

	/**
	 * Add a result. Not synchronized.
	 *
	 * @param result
	 *            the result
	 */
	public void add(PeakResult result);

	/**
	 * Add all results.
	 *
	 * @param results
	 *            the results
	 */
	public void addCollection(Collection<PeakResult> results);

	/**
	 * Add all results.
	 *
	 * @param results
	 *            the results
	 */
	public void addArray(PeakResult[] results);

	/**
	 * Adds the results.
	 *
	 * @param results
	 *            the results
	 */
	public void addStore(PeakResultStore results);

	/**
	 * Clear the results.
	 */
	public void clear();

	/**
	 * Trims the capacity of this instance to be the current size. An application can use this operation to minimize
	 * the storage of an instance.
	 */
	public void trimToSize();

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
	 * Convert to an array. This is a new allocation of storage space.
	 *
	 * @return the peak result array
	 */
	public PeakResult[] toArray();

	/**
	 * Copy the results
	 */
	public PeakResultStore copy();

	/**
	 * Removes the result if it matches the filter. If objects are removed then the order of elements may change.
	 *
	 * @param filter
	 *            the filter
	 * @return true, if any were removed
	 */
	public boolean removeIf(PeakResultPredicate filter);

	/**
	 * Execute the procedure on each result in the store.
	 *
	 * @param procedure
	 *            the procedure
	 */
	public void forEach(PeakResultProcedure procedure);

	/**
	 * Get a subset of the results if they match the filter.
	 *
	 * @param filter
	 *            the filter
	 * @return the results
	 */
	public PeakResult[] subset(PeakResultPredicate filter);

	/**
	 * Shuffle the results.
	 *
	 * @param randomGenerator
	 *            the random generator
	 */
	public void shuffle(final RandomGenerator randomGenerator);
}
