package gdsc.smlm.results;

import java.util.Collection;
import java.util.Comparator;

import org.apache.commons.math3.random.RandomGenerator;

import gdsc.core.data.DataException;
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
 * Stores peak results and prevents modification.
 */
public class ImmutablePeakResultStore implements PeakResultStore
{
	private final PeakResultStore store;

	/**
	 * Instantiates a new immutable peak result store.
	 *
	 * @param store
	 *            the store
	 */
	public ImmutablePeakResultStore(PeakResultStore store)
	{
		if (store == null)
			throw new NullPointerException("Store must not be null");
		this.store = store;
	}

	public PeakResult get(int index)
	{
		return new ImmutablePeakResult(store.get(index));
	}

	public int size()
	{
		return store.size();
	}

	public void add(PeakResult result)
	{
		throw new DataException("This result store is immutable");
	}

	public void addCollection(Collection<PeakResult> results)
	{
		throw new DataException("This result store is immutable");
	}

	public void addArray(PeakResult[] results)
	{
		throw new DataException("This result store is immutable");
	}

	public void addStore(PeakResultStore results)
	{
		throw new DataException("This result store is immutable");
	}

	public void clear()
	{
		throw new DataException("This result store is immutable");
	}

	public void trimToSize()
	{
		store.trimToSize();
	}

	public void sort()
	{
		store.sort();
	}

	public void sort(Comparator<PeakResult> comparator)
	{
		store.sort(comparator);
	}

	public PeakResult[] toArray()
	{
		return makeImmutable(store.toArray());
	}

	/**
	 * Make the array a collection of immutable peak result objects.
	 *
	 * @param array
	 *            the array
	 * @return the array
	 */
	private static PeakResult[] makeImmutable(PeakResult[] array)
	{
		for (int i = 0; i < array.length; i++)
			array[i] = new ImmutablePeakResult(array[i]);
		return array;
	}

	public PeakResultStore copy()
	{
		return new ImmutablePeakResultStore(store.copy());
	}

	public PeakResultStore copy(boolean deepCopy)
	{
		return new ImmutablePeakResultStore(store.copy(deepCopy));
	}

	public boolean removeIf(PeakResultPredicate filter)
	{
		throw new DataException("This result store is immutable");
	}

	/**
	 * Used to wrap the results to make them immutable
	 */
	private static class ImmutablePeakResultProcedure implements PeakResultProcedure
	{
		PeakResultProcedure procedure;

		public ImmutablePeakResultProcedure(PeakResultProcedure procedure)
		{
			this.procedure = procedure;
		}

		public void execute(PeakResult peakResult)
		{
			procedure.execute(new ImmutablePeakResult(peakResult));
		}
	}

	public void forEach(PeakResultProcedure procedure)
	{
		store.forEach(new ImmutablePeakResultProcedure(procedure));
	}

	public PeakResult[] subset(PeakResultPredicate filter)
	{
		return makeImmutable(store.subset(filter));
	}

	public void shuffle(final RandomGenerator randomGenerator)
	{
		store.shuffle(randomGenerator);
	}
}