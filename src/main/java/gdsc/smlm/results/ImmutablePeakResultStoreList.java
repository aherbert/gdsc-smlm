package gdsc.smlm.results;

import java.util.Comparator;

import org.apache.commons.math3.random.RandomGenerator;

import gdsc.core.data.DataException;

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
 * Stores peak results and prevents modification.
 */
public class ImmutablePeakResultStoreList extends ImmutablePeakResultStore implements PeakResultStoreList
{
	private final PeakResultStoreList store;

	/**
	 * Instantiates a new immutable peak result store.
	 *
	 * @param store
	 *            the store
	 */
	public ImmutablePeakResultStoreList(PeakResultStoreList store)
	{
		super(store);
		this.store = store;
	}

	public PeakResult get(int index)
	{
		return new ImmutablePeakResult(store.get(index));
	}

	public PeakResult remove(int index)
	{
		throw new DataException("This result store is immutable");
	}

	public void sort()
	{
		store.sort();
	}

	public void sort(Comparator<PeakResult> comparator)
	{
		store.sort(comparator);
	}

	public PeakResultStoreList copy()
	{
		return new ImmutablePeakResultStoreList((PeakResultStoreList) store.copy());
	}

	public PeakResultStoreList copy(boolean deepCopy)
	{
		return new ImmutablePeakResultStoreList((PeakResultStoreList) store.copy(deepCopy));
	}

	public void shuffle(final RandomGenerator randomGenerator)
	{
		store.shuffle(randomGenerator);
	}

	public int indexOf(PeakResult result)
	{
		return store.indexOf(result);
	}

	public int lastIndexOf(PeakResult result)
	{
		return store.lastIndexOf(result);
	}
}
