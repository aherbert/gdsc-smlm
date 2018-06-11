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
package gdsc.smlm.results;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;

import org.apache.commons.math3.random.RandomAdaptor;
import org.apache.commons.math3.random.RandomGenerator;

import gdsc.core.utils.TurboList;
import gdsc.core.utils.TurboList.SimplePredicate;
import gdsc.smlm.results.predicates.PeakResultPredicate;
import gdsc.smlm.results.procedures.PeakResultProcedure;


/**
 * Stores peak results using a TurboList. This is similar to an ArrayList but does not have concurrency checking.
 */
public class TurboListPeakResultStore implements PeakResultStoreList, PeakResultStoreCollection
{
	/** The results. */
	private TurboList<PeakResult> results;

	/**
	 * Instantiates a new array list peak results store.
	 *
	 * @param capacity
	 *            the capacity
	 */
	public TurboListPeakResultStore(int capacity)
	{
		this.results = new TurboList<PeakResult>(capacity);
	}

	/**
	 * Instantiates a new array list peak result store.
	 *
	 * @param store
	 *            the store to copy
	 */
	public TurboListPeakResultStore(TurboListPeakResultStore store)
	{
		this.results = new TurboList<PeakResult>(store.results);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStoreList#get(int)
	 */
	public PeakResult get(int index)
	{
		return results.get(index);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#size()
	 */
	public int size()
	{
		return results.size();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#add(gdsc.smlm.results.PeakResult)
	 */
	public boolean add(PeakResult result)
	{
		return results.add(result);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#addCollection(java.util.Collection)
	 */
	public boolean addCollection(Collection<PeakResult> results)
	{
		return this.results.addAll(results);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#addArray(gdsc.smlm.results.PeakResult[])
	 */
	public boolean addArray(PeakResult[] results)
	{
		return this.results.addAll(Arrays.asList(results));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#addStore(gdsc.smlm.results.PeakResultStore)
	 */
	public boolean addStore(PeakResultStore results)
	{
		if (results instanceof PeakResultStoreCollection)
		{
			return this.results.addAll(((PeakResultStoreCollection) results).getCollectionReference());
		}
		else
		{
			return addArray(results.toArray());
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStoreList#remove(int)
	 */
	public PeakResult remove(int index)
	{
		return results.remove(index);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStoreList#remove(int, int)
	 */
	public void remove(int fromIndex, int toIndex)
	{
		if (fromIndex > toIndex)
		{
			throw new IllegalArgumentException("fromIndex must be <= toIndex");
		}
		for (int i = toIndex; i >= fromIndex; i--)
		{
			results.remove(i);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#remove(gdsc.smlm.results.PeakResult)
	 */
	public boolean remove(PeakResult result)
	{
		return results.remove(result);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#removeCollection(java.util.Collection)
	 */
	public boolean removeCollection(Collection<PeakResult> results)
	{
		return this.results.removeAll(results);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#removeArray(gdsc.smlm.results.PeakResult[])
	 */
	public boolean removeArray(PeakResult[] results)
	{
		return this.results.removeAll(Arrays.asList(results));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#removeStore(gdsc.smlm.results.PeakResultStore)
	 */
	public boolean removeStore(PeakResultStore results)
	{
		if (results instanceof PeakResultStoreCollection)
		{
			return this.results.removeAll(((PeakResultStoreCollection) results).getCollectionReference());
		}
		else
		{
			return removeArray(results.toArray());
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#retainCollection(java.util.Collection)
	 */
	public boolean retainCollection(Collection<PeakResult> results)
	{
		return this.results.retainAll(results);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#retainArray(gdsc.smlm.results.PeakResult[])
	 */
	public boolean retainArray(PeakResult[] results)
	{
		return this.results.retainAll(Arrays.asList(results));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#retainStore(gdsc.smlm.results.PeakResultStore)
	 */
	public boolean retainStore(PeakResultStore results)
	{
		if (results instanceof PeakResultStoreCollection)
		{
			return this.results.retainAll(((PeakResultStoreCollection) results).getCollectionReference());
		}
		else
		{
			return retainArray(results.toArray());
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#clear()
	 */
	public void clear()
	{
		results.clear();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#trimToSize()
	 */
	public void trimToSize()
	{
		results.trimToSize();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStoreList#sort()
	 */
	public void sort()
	{
		Collections.sort(results);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStoreList#sort(java.util.Comparator)
	 */
	public void sort(Comparator<PeakResult> comparator)
	{
		Collections.sort(results, comparator);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#toArray()
	 */
	public PeakResult[] toArray()
	{
		return results.toArray(new PeakResult[size()]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#copy()
	 */
	public PeakResultStore copy()
	{
		return new TurboListPeakResultStore(this);
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
			TurboListPeakResultStore copy = new TurboListPeakResultStore(size());
			for (int i = 0, size = size(); i < size; i++)
				copy.add(results.getf(i).clone());
			return copy;
		}
		return copy();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#removeIf(gdsc.smlm.results.PeakResultPredicate)
	 */
	public boolean removeIf(final PeakResultPredicate filter)
	{
		// Delegate to the list implementation
		return this.results.removeIf(new SimplePredicate<PeakResult>()
		{
			public boolean test(PeakResult t)
			{
				return filter.test(t);
			}
		});
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#forEach(gdsc.smlm.results.procedures.PeakResultProcedure)
	 */
	public void forEach(PeakResultProcedure procedure)
	{
		for (int i = 0, size = size(); i < size; i++)
			procedure.execute(results.getf(i));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#subset(gdsc.smlm.results.procedures.PeakResultPredicate)
	 */
	public PeakResult[] subset(PeakResultPredicate filter)
	{
		final ArrayPeakResultStore list = new ArrayPeakResultStore(10);
		for (int i = 0, size = size(); i < size; i++)
			if (filter.test(results.getf(i)))
				list.add(results.getf(i));
		return list.toArray();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStoreList#shuffle(org.apache.commons.math3.random.RandomGenerator)
	 */
	public void shuffle(RandomGenerator randomGenerator)
	{
		Collections.shuffle(results, RandomAdaptor.createAdaptor(randomGenerator));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStoreList#indexOf(gdsc.smlm.results.PeakResult)
	 */
	public int indexOf(PeakResult result)
	{
		return results.indexOf(result);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStoreList#lastIndexOf(gdsc.smlm.results.PeakResult)
	 */
	public int lastIndexOf(PeakResult result)
	{
		return results.lastIndexOf(result);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStore#contains(gdsc.smlm.results.PeakResult)
	 */
	public boolean contains(PeakResult result)
	{
		return results.contains(result);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStoreCollection#getCollection()
	 */
	@SuppressWarnings("unchecked")
	public Collection<PeakResult> getCollection()
	{
		return (Collection<PeakResult>) results.clone();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultStoreCollection#getCollectionReference()
	 */
	public Collection<PeakResult> getCollectionReference()
	{
		return results;
	}
}
