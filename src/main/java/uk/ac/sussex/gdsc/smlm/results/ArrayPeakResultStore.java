/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package uk.ac.sussex.gdsc.smlm.results;

import java.util.Arrays;
import java.util.BitSet;
import java.util.Collection;
import java.util.Comparator;
import java.util.Objects;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.rng.UniformRandomProvider;

import uk.ac.sussex.gdsc.smlm.results.predicates.PeakResultPredicate;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

/**
 * Stores peak results using an array.
 */
public class ArrayPeakResultStore implements PeakResultStoreList
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
	 * @throws NullPointerException
	 *             if the store is null
	 */
	public ArrayPeakResultStore(ArrayPeakResultStore store) throws NullPointerException
	{
		this.results = store.toArray();
		this.size = store.size;
	}

	/**
	 * Instantiates a new array peak result store.
	 *
	 * @param results
	 *            the results
	 * @throws NullPointerException
	 *             if the results are null
	 */
	public ArrayPeakResultStore(PeakResult[] results) throws NullPointerException
	{
		this.results = results;
		this.size = results.length;
	}

	/**
	 * {@inheritDoc}
	 * <p>
	 * Note: This does not check against the current size so can return stale data.
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStoreList#get(int)
	 */
	@Override
	public PeakResult get(int index)
	{
		return results[index];
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#size()
	 */
	@Override
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
	 *            the length
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
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#add(uk.ac.sussex.gdsc.smlm.results.PeakResult)
	 */
	@Override
	public boolean add(PeakResult result)
	{
		checkCapacity(1);
		results[size++] = result;
		return true;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#addAll(java.util.Collection)
	 */
	@Override
	public boolean addCollection(Collection<PeakResult> results)
	{
		return addArray(results.toArray(new PeakResult[results.size()]));
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#addAll(uk.ac.sussex.gdsc.smlm.results.PeakResult[])
	 */
	@Override
	public boolean addArray(PeakResult[] results)
	{
		if (results == null)
			return false;
		return addArray(results, results.length);
	}

	private boolean addArray(PeakResult[] results, int length)
	{
		if (results == null || length == 0)
			return false;
		checkCapacity(length);
		System.arraycopy(results, 0, this.results, size, length);
		size += length;
		return true;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#add(uk.ac.sussex.gdsc.smlm.results.PeakResultStore)
	 */
	@Override
	public boolean addStore(PeakResultStore results)
	{
		if (results instanceof ArrayPeakResultStore)
		{
			final ArrayPeakResultStore store = (ArrayPeakResultStore) results;
			return addArray(store.results, store.size);
		}
		return addArray(results.toArray());
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStoreList#remove(int)
	 */
	@Override
	public PeakResult remove(int index)
	{
		rangeCheck(index);
		final PeakResult oldValue = results[index];
		fastRemove(index);
		return oldValue;
	}

	/**
	 * Checks if the given index is in range. If not, throws an appropriate
	 * runtime exception. This method does *not* check if the index is
	 * negative: It is always used immediately prior to an array access,
	 * which throws an ArrayIndexOutOfBoundsException if index is negative.
	 */
	private void rangeCheck(int index)
	{
		if (index >= size)
			throw new IndexOutOfBoundsException(outOfBoundsMsg(index));
	}

	/**
	 * A version of rangeCheck with lower bounds check
	 */
	private void rangeCheckWithLowerBounds(int index)
	{
		if (index > size || index < 0)
			throw new IndexOutOfBoundsException(outOfBoundsMsg(index));
	}

	/**
	 * Constructs an IndexOutOfBoundsException detail message.
	 */
	private String outOfBoundsMsg(int index)
	{
		return "Index: " + index + ", Size: " + size;
	}

	/*
	 * Private remove method that skips bounds checking and does not
	 * return the value removed.
	 */
	private void fastRemove(int index)
	{
		final int numMoved = size - index - 1;
		if (numMoved > 0)
			System.arraycopy(results, index + 1, results, index, numMoved);
		results[--size] = null; // Let gc do its work
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStoreList#remove(int, int)
	 */
	@Override
	public void remove(int fromIndex, int toIndex)
	{
		if (fromIndex > toIndex)
			throw new IllegalArgumentException("fromIndex must be <= toIndex");
		rangeCheckWithLowerBounds(fromIndex);
		rangeCheck(toIndex); // This is above fromIndex so ignore lower bounds check
		toIndex++; // Make exclusive
		final int numMoved = size - toIndex;
		if (numMoved > 0)
			System.arraycopy(results, toIndex, results, fromIndex, numMoved);
		// Let gc do its work
		while (fromIndex++ < toIndex)
			results[size--] = null;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#remove(uk.ac.sussex.gdsc.smlm.results.PeakResult)
	 */
	@Override
	public boolean remove(PeakResult result)
	{
		final int index = indexOf(result);
		if (index != -1)
		{
			fastRemove(index);
			return true;
		}
		return false;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#removeAll(java.util.Collection)
	 */
	@Override
	public boolean removeCollection(Collection<PeakResult> results)
	{
		return removeArray(results.toArray(new PeakResult[results.size()]));
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#removeAll(uk.ac.sussex.gdsc.smlm.results.PeakResult[])
	 */
	@Override
	public boolean removeArray(PeakResult[] results)
	{
		if (results == null || results.length == 0)
			return false;
		return batchRemove(results, false);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#remove(uk.ac.sussex.gdsc.smlm.results.PeakResultStore)
	 */
	@Override
	public boolean removeStore(PeakResultStore results)
	{
		if (results instanceof ArrayPeakResultStore)
			return removeArray(((ArrayPeakResultStore) results).results);
		return removeArray(results.toArray());
	}

	private boolean batchRemove(PeakResult[] results2, boolean complement)
	{
		// Adapted from java.utisl.ArrayList

		final ArrayPeakResultStore c = new ArrayPeakResultStore(results2);
		int r = 0, w = 0;
		boolean modified = false;
		try
		{
			for (; r < size; r++)
				if (c.contains(results[r]) == complement)
					results[w++] = results[r];
		}
		finally
		{
			// Preserve data even if c.contains() throws.
			if (r != size)
			{
				System.arraycopy(results, r, results, w, size - r);
				w += size - r;
			}
			if (w != size)
			{
				// clear to let GC do its work
				for (int i = w; i < size; i++)
					results[i] = null;
				size = w;
				modified = true;
			}
		}
		return modified;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#retainAll(java.util.Collection)
	 */
	@Override
	public boolean retainCollection(Collection<PeakResult> results)
	{
		return retainArray(results.toArray(new PeakResult[results.size()]));
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#retainArray(uk.ac.sussex.gdsc.smlm.results.PeakResult[])
	 */
	@Override
	public boolean retainArray(PeakResult[] results)
	{
		if (results == null || results.length == 0)
		{
			final boolean result = size != 0;
			clear();
			return result;
		}
		return batchRemove(results, true);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#retain(uk.ac.sussex.gdsc.smlm.results.PeakResultStore)
	 */
	@Override
	public boolean retainStore(PeakResultStore results)
	{
		if (results instanceof ArrayPeakResultStore)
			return retainArray(((ArrayPeakResultStore) results).results);
		return retainArray(results.toArray());
	}

	/**
	 * {@inheritDoc}
	 * <p>
	 * Note: This does not remove the references to the underlying data or reallocate storage thus {@link #get(int)} can
	 * return stale data.
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#clear()
	 */
	@Override
	public void clear()
	{
		size = 0;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#trimToSize()
	 */
	@Override
	public void trimToSize()
	{
		if (size < results.length)
			results = toArray();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStoreList#sort(java.util.Comparator)
	 */
	@Override
	public void sort(Comparator<PeakResult> comparator)
	{
		Arrays.sort(results, 0, size, comparator);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#toArray()
	 */
	@Override
	public PeakResult[] toArray()
	{
		final PeakResult[] array = new PeakResult[size];
		System.arraycopy(results, 0, array, 0, size);
		return array;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#copy()
	 */
	@Override
	public PeakResultStore copy()
	{
		return new ArrayPeakResultStore(this);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#copy(boolean)
	 */
	@Override
	public PeakResultStore copy(boolean deepCopy)
	{
		if (deepCopy)
		{
			final ArrayPeakResultStore copy = new ArrayPeakResultStore(size());
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
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#removeIf(uk.ac.sussex.gdsc.smlm.results.predicates.PeakResultPredicate)
	 */
	@Override
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
			if (filter.test(results[i]))
			{
				removeSet.set(i);
				removeCount++;
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
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#forEach(uk.ac.sussex.gdsc.smlm.results.procedures.
	 * PeakResultProcedure)
	 */
	@Override
	public void forEach(PeakResultProcedure procedure)
	{
		for (int i = 0; i < size; i++)
			procedure.execute(results[i]);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#subset(uk.ac.sussex.gdsc.smlm.results.procedures.
	 * PeakResultPredicate)
	 */
	@Override
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
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStoreList#shuffle(org.apache.commons.math3.random.RandomGenerator)
	 */
	@Override
	public void shuffle(RandomGenerator randomSource)
	{
		// Fisher-Yates shuffle
		for (int i = size; i-- > 1;)
		{
			final int j = randomSource.nextInt(i + 1);
			final PeakResult tmp = results[i];
			results[i] = results[j];
			results[j] = tmp;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStoreList#shuffle(org.apache.commons.rng.UniformRandomProvider)
	 */
	@Override
	public void shuffle(UniformRandomProvider randomSource)
	{
		// Fisher-Yates shuffle
		for (int i = size; i-- > 1;)
		{
			final int j = randomSource.nextInt(i + 1);
			final PeakResult tmp = results[i];
			results[i] = results[j];
			results[j] = tmp;
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStoreList#indexOf(uk.ac.sussex.gdsc.smlm.results.PeakResult)
	 */
	@Override
	public int indexOf(PeakResult result)
	{
		if (result == null)
		{
			for (int i = 0; i < size; i++)
				if (results[i] == null)
					return i;
		}
		else
			for (int i = 0; i < size; i++)
				if (result.equals(results[i]))
					return i;
		return -1;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStoreList#lastIndexOf(uk.ac.sussex.gdsc.smlm.results.PeakResult)
	 */
	@Override
	public int lastIndexOf(PeakResult result)
	{
		if (result == null)
		{
			for (int i = size; i-- > 0;)
				if (results[i] == null)
					return i;
		}
		else
			for (int i = size; i-- > 0;)
				if (result.equals(results[i]))
					return i;
		return -1;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#contains(uk.ac.sussex.gdsc.smlm.results.PeakResult)
	 */
	@Override
	public boolean contains(PeakResult result)
	{
		return indexOf(result) != -1;
	}

}
