package gdsc.smlm.results;

import java.util.Arrays;
import java.util.Collection;
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
 * Stores peak results using an array.
 */
public class ArrayPeakResultStore implements PeakResultStore
{
	/** The results. */
	private PeakResult[] results;

	/** The size. */
	private int size;

	/**
	 * Instantiates a new array list peak results store.
	 *
	 * @param capacity
	 *            the capacity
	 */
	public ArrayPeakResultStore(int capacity)
	{
		this.results = new PeakResult[Math.max(capacity, 0)];
		this.size = 0;
	}

	/**
	 * Instantiates a new array peak result store.
	 *
	 * @param store
	 *            the store to copy
	 */
	public ArrayPeakResultStore(ArrayPeakResultStore store)
	{
		this.results = store.toArray();
		this.size = store.size;
	}

	/**
	 * {@inheritDoc}
	 * <p>
	 * Note: This does not check against the current size so can return stale data.
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#get(int)
	 */
	public PeakResult get(int index)
	{
		return results[index];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#size()
	 */
	public int size()
	{
		return size;
	}

	/**
	 * Ensure that the specified number of elements can be added to the array.
	 * <p>
	 * This is not synchronized. However any class using the safeAdd() methods in different threads should be using the
	 * same synchronized method to add data thus this method will be within synchronized code.
	 * 
	 * @param length
	 */
	private void checkCapacity(int length)
	{
		final int minCapacity = size + length;
		final int oldCapacity = results.length;
		if (minCapacity > oldCapacity)
		{
			int newCapacity = (oldCapacity * 3) / 2 + 1;
			if (newCapacity < minCapacity)
				newCapacity = minCapacity;
			final PeakResult[] newResults = new PeakResult[newCapacity];
			System.arraycopy(results, 0, newResults, 0, size);
			results = newResults;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#add(gdsc.smlm.results.PeakResult)
	 */
	public void add(PeakResult result)
	{
		checkCapacity(1);
		results[size++] = result;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#addAll(java.util.Collection)
	 */
	public void addCollection(Collection<PeakResult> results)
	{
		addArray(results.toArray(new PeakResult[results.size()]));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#addAll(gdsc.smlm.results.PeakResult[])
	 */
	public void addArray(PeakResult[] results)
	{
		if (results == null)
			return;
		checkCapacity(results.length);
		System.arraycopy(results, 0, this.results, size, results.length);
		size += results.length;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#add(gdsc.smlm.results.PeakResultStore)
	 */
	public void addStore(PeakResultStore results)
	{
		if (results instanceof ArrayPeakResultStore)
		{
			addArray(((ArrayPeakResultStore) results).results);
		}
		else
		{
			addArray(results.toArray());
		}
	}

	/**
	 * {@inheritDoc}
	 * <p>
	 * Note: This does not remove the references to the underlying data or reallocate storage thus {@link #get(int)} can
	 * return stale data.
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#clear()
	 */
	public void clear()
	{
		size = 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#trimToSize()
	 */
	public void trimToSize()
	{
		if (size < results.length)
			results = toArray();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#sort()
	 */
	public void sort()
	{
		Arrays.sort(results, 0, size);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#sort(java.util.Comparator)
	 */
	public void sort(Comparator<PeakResult> comparator)
	{
		Arrays.sort(results, 0, size, comparator);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#toArray()
	 */
	public PeakResult[] toArray()
	{
		PeakResult[] array = new PeakResult[size];
		System.arraycopy(results, 0, array, 0, size);
		return array;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#copy()
	 */
	public PeakResultStore copy()
	{
		return new ArrayPeakResultStore(this);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#removeNullResults()
	 */
	public void removeNullResults()
	{
		int newSize = 0;
		for (int i = 0; i < size; i++)
		{
			if (results[i] != null)
				results[newSize++] = results[i];
		}
		size = newSize;
	}
}
