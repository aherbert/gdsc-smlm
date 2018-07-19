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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;

import org.apache.commons.math3.random.RandomAdaptor;
import org.apache.commons.math3.random.RandomGenerator;

import uk.ac.sussex.gdsc.core.utils.TurboList;
import uk.ac.sussex.gdsc.core.utils.TurboList.SimplePredicate;
import uk.ac.sussex.gdsc.smlm.results.predicates.PeakResultPredicate;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

/**
 * Stores peak results using an ArrayList.
 */
public class ArrayListPeakResultStore implements PeakResultStoreList, PeakResultStoreCollection
{
	/** The results. */
	private ArrayList<PeakResult> results;

	/**
	 * Instantiates a new array list peak results store.
	 *
	 * @param capacity
	 *            the capacity
	 */
	public ArrayListPeakResultStore(int capacity)
	{
		this.results = new ArrayList<>(capacity);
	}

	/**
	 * Instantiates a new array list peak result store.
	 *
	 * @param store
	 *            the store to copy
	 */
	public ArrayListPeakResultStore(ArrayListPeakResultStore store)
	{
		this.results = new ArrayList<>(store.results);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStoreList#get(int)
	 */
	@Override
	public PeakResult get(int index)
	{
		return results.get(index);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#size()
	 */
	@Override
	public int size()
	{
		return results.size();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#add(uk.ac.sussex.gdsc.smlm.results.PeakResult)
	 */
	@Override
	public boolean add(PeakResult result)
	{
		return results.add(result);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#addCollection(java.util.Collection)
	 */
	@Override
	public boolean addCollection(Collection<PeakResult> results)
	{
		return this.results.addAll(results);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#addArray(uk.ac.sussex.gdsc.smlm.results.PeakResult[])
	 */
	@Override
	public boolean addArray(PeakResult[] results)
	{
		return this.results.addAll(Arrays.asList(results));
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#addStore(uk.ac.sussex.gdsc.smlm.results.PeakResultStore)
	 */
	@Override
	public boolean addStore(PeakResultStore results)
	{
		if (results instanceof PeakResultStoreCollection)
			return this.results.addAll(((PeakResultStoreCollection) results).getCollectionReference());
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
		return results.remove(index);
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
		for (int i = toIndex; i >= fromIndex; i--)
			results.remove(i);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#remove(uk.ac.sussex.gdsc.smlm.results.PeakResult)
	 */
	@Override
	public boolean remove(PeakResult result)
	{
		return results.remove(result);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#removeCollection(java.util.Collection)
	 */
	@Override
	public boolean removeCollection(Collection<PeakResult> results)
	{
		return this.results.removeAll(results);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#removeArray(uk.ac.sussex.gdsc.smlm.results.PeakResult[])
	 */
	@Override
	public boolean removeArray(PeakResult[] results)
	{
		return this.results.removeAll(Arrays.asList(results));
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#removeStore(uk.ac.sussex.gdsc.smlm.results.PeakResultStore)
	 */
	@Override
	public boolean removeStore(PeakResultStore results)
	{
		if (results instanceof PeakResultStoreCollection)
			return this.results.removeAll(((PeakResultStoreCollection) results).getCollectionReference());
		return removeArray(results.toArray());
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#retainCollection(java.util.Collection)
	 */
	@Override
	public boolean retainCollection(Collection<PeakResult> results)
	{
		return this.results.retainAll(results);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#retainArray(uk.ac.sussex.gdsc.smlm.results.PeakResult[])
	 */
	@Override
	public boolean retainArray(PeakResult[] results)
	{
		return this.results.retainAll(Arrays.asList(results));
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#retainStore(uk.ac.sussex.gdsc.smlm.results.PeakResultStore)
	 */
	@Override
	public boolean retainStore(PeakResultStore results)
	{
		if (results instanceof PeakResultStoreCollection)
			return this.results.retainAll(((PeakResultStoreCollection) results).getCollectionReference());
		return retainArray(results.toArray());
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#clear()
	 */
	@Override
	public void clear()
	{
		results.clear();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#trimToSize()
	 */
	@Override
	public void trimToSize()
	{
		results.trimToSize();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStoreList#sort(java.util.Comparator)
	 */
	@Override
	public void sort(Comparator<PeakResult> comparator)
	{
		Collections.sort(results, comparator);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#toArray()
	 */
	@Override
	public PeakResult[] toArray()
	{
		return results.toArray(new PeakResult[size()]);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#copy()
	 */
	@Override
	public PeakResultStore copy()
	{
		return new ArrayListPeakResultStore(this);
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
			final ArrayListPeakResultStore copy = new ArrayListPeakResultStore(size());
			for (int i = 0, size = size(); i < size; i++)
				copy.add(results.get(i).clone());
			return copy;
		}
		return copy();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#removeIf(uk.ac.sussex.gdsc.smlm.results.PeakResultPredicate)
	 */
	@Override
	public boolean removeIf(final PeakResultPredicate filter)
	{
		// Util we upgrade the Java version to 1.8 the ArrayList does not support
		// predicates so use a TurboList
		final TurboList<PeakResult> temp = new TurboList<>(this.results);
		if (temp.removeIf(new SimplePredicate<PeakResult>()
		{
			@Override
			public boolean test(PeakResult t)
			{
				return filter.test(t);
			}
		}))
		{
			this.results = new ArrayList<>(temp);
			return true;
		}
		return false;
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
		for (int i = 0, size = size(); i < size; i++)
			procedure.execute(results.get(i));
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
		for (int i = 0, size = size(); i < size; i++)
			if (filter.test(results.get(i)))
				list.add(results.get(i));
		return list.toArray();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStoreList#shuffle(org.apache.commons.math3.random.RandomGenerator)
	 */
	@Override
	public void shuffle(RandomGenerator randomGenerator)
	{
		Collections.shuffle(results, RandomAdaptor.createAdaptor(randomGenerator));
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStoreList#indexOf(uk.ac.sussex.gdsc.smlm.results.PeakResult)
	 */
	@Override
	public int indexOf(PeakResult result)
	{
		return results.indexOf(result);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStoreList#lastIndexOf(uk.ac.sussex.gdsc.smlm.results.PeakResult)
	 */
	@Override
	public int lastIndexOf(PeakResult result)
	{
		return results.lastIndexOf(result);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStore#contains(uk.ac.sussex.gdsc.smlm.results.PeakResult)
	 */
	@Override
	public boolean contains(PeakResult result)
	{
		return results.contains(result);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStoreCollection#getCollection()
	 */
	@Override
	@SuppressWarnings("unchecked")
	public Collection<PeakResult> getCollection()
	{
		return (Collection<PeakResult>) results.clone();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.PeakResultStoreCollection#getCollectionReference()
	 */
	@Override
	public Collection<PeakResult> getCollectionReference()
	{
		return results;
	}
}
