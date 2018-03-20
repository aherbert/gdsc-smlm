package gdsc.smlm.results;

import java.util.Arrays;
import java.util.BitSet;
import java.util.Collection;
import java.util.Comparator;
import java.util.Objects;

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
	 * Instantiates a new array peak result store.
	 *
	 * @param results
	 *            the results
	 */
	ArrayPeakResultStore(PeakResult[] results)
	{
		this.results = results;
		this.size = results.length;
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
	 * @see gdsc.smlm.results.PeakResultStore#copy(boolean)
	 */
	public PeakResultStore copy(boolean deepCopy)
	{
		if (deepCopy)
		{
			ArrayPeakResultStore copy = new ArrayPeakResultStore(size());
			for (int i = 0, size = size(); i < size; i++)
				copy.add(results[i].clone());
			return copy;
		}
		return copy();
	}

	/**
	 * {@inheritDoc}
	 * <p>
	 * Note: This does not remove the references to the underlying data or reallocate storage thus {@link #get(int)} can
	 * return stale data.
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#removeIf(gdsc.smlm.results.predicates.PeakResultPredicate)
	 */
	public boolean removeIf(PeakResultPredicate filter)
	{
		Objects.requireNonNull(filter);

		// Adapted from java.util.ArrayList (Java 1.8)

		// figure out which elements are to be removed
		// any exception thrown from the filter predicate at this stage
		// will leave the collection unmodified
		int removeCount = 0;
		final int size = this.size;
		final BitSet removeSet = new BitSet(size);
		for (int i = 0; i < size; i++)
		{
			if (filter.test(results[i]))
			{
				removeSet.set(i);
				removeCount++;
			}
		}

		// shift surviving elements left over the spaces left by removed elements
		final boolean anyToRemove = removeCount > 0;
		if (anyToRemove)
		{
			final int newSize = size - removeCount;
			for (int i = 0, j = 0; (i < size) && (j < newSize); i++, j++)
			{
				i = removeSet.nextClearBit(i);
				results[j] = results[i];
			}
			//for (int k = newSize; k < size; k++)
			//{
			//	results[k] = null; // Let gc do its work
			//}
			this.size = newSize;
		}

		return anyToRemove;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#forEach(gdsc.smlm.results.procedures.PeakResultProcedure)
	 */
	public void forEach(PeakResultProcedure procedure)
	{
		for (int i = 0; i < size; i++)
			procedure.execute(results[i]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#subset(gdsc.smlm.results.procedures.PeakResultPredicate)
	 */
	public PeakResult[] subset(PeakResultPredicate filter)
	{
		final ArrayPeakResultStore list = new ArrayPeakResultStore(10);
		for (int i = 0; i < size; i++)
			if (filter.test(results[i]))
				list.add(results[i]);
		return list.toArray();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#shuffle(org.apache.commons.math3.random.RandomGenerator)
	 */
	public void shuffle(RandomGenerator randomGenerator)
	{
		// Fisher-Yates shuffle
		for (int i = size; i-- > 1;)
		{
			int j = randomGenerator.nextInt(i + 1);
			PeakResult tmp = results[i];
			results[i] = results[j];
			results[j] = tmp;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#indexOf(gdsc.smlm.results.PeakResult)
	 */
	public int indexOf(PeakResult result)
	{
		if (result == null)
		{
			for (int i = 0; i < size; i++)
				if (results[i] == null)
					return i;
		}
		else
		{
			for (int i = 0; i < size; i++)
				if (result.equals(results[i]))
					return i;
		}
		return -1;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#lastIndexOf(gdsc.smlm.results.PeakResult)
	 */
	public int lastIndexOf(PeakResult result)
	{
		if (result == null)
		{
			for (int i = size; i-- > 0;)
				if (results[i] == null)
					return i;
		}
		else
		{
			for (int i = size; i-- > 0;)
				if (result.equals(results[i]))
					return i;
		}
		return -1;
	}
}
