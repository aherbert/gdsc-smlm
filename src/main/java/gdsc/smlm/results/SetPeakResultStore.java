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
import java.util.HashSet;

import gdsc.smlm.results.predicates.PeakResultPredicate;
import gdsc.smlm.results.procedures.PeakResultProcedure;

/**
 * Stores peak results using a set. This is similar to an HashSet but does not have concurrency checking.
 */
public class SetPeakResultStore implements PeakResultStore, PeakResultStoreCollection
{
	/** The results. */
	private HashSet<PeakResult> results;

	/**
	 * Instantiates a new set peak results store.
	 *
	 * @param capacity
	 *            the capacity
	 */
	public SetPeakResultStore(int capacity)
	{
		this.results = new HashSet<>(capacity);
	}

	/**
	 * Instantiates a new array list peak result store.
	 *
	 * @param store
	 *            the store to copy
	 */
	public SetPeakResultStore(SetPeakResultStore store)
	{
		this.results = new HashSet<>(store.results);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.PeakResultStore#size()
	 */
	@Override
	public int size()
	{
		return results.size();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.PeakResultStore#add(gdsc.smlm.results.PeakResult)
	 */
	@Override
	public boolean add(PeakResult result)
	{
		return results.add(result);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.PeakResultStore#addCollection(java.util.Collection)
	 */
	@Override
	public boolean addCollection(Collection<PeakResult> results)
	{
		return this.results.addAll(results);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.PeakResultStore#addArray(gdsc.smlm.results.PeakResult[])
	 */
	@Override
	public boolean addArray(PeakResult[] results)
	{
		return this.results.addAll(Arrays.asList(results));
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.PeakResultStore#addStore(gdsc.smlm.results.PeakResultStore)
	 */
	@Override
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
	 * @see gdsc.smlm.results.PeakResultStore#remove(gdsc.smlm.results.PeakResult)
	 */
	@Override
	public boolean remove(PeakResult result)
	{
		return results.remove(result);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.PeakResultStore#removeCollection(java.util.Collection)
	 */
	@Override
	public boolean removeCollection(Collection<PeakResult> results)
	{
		return this.results.removeAll(results);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.PeakResultStore#removeArray(gdsc.smlm.results.PeakResult[])
	 */
	@Override
	public boolean removeArray(PeakResult[] results)
	{
		return this.results.removeAll(Arrays.asList(results));
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.PeakResultStore#removeStore(gdsc.smlm.results.PeakResultStore)
	 */
	@Override
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
	@Override
	public boolean retainCollection(Collection<PeakResult> results)
	{
		return this.results.retainAll(results);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.PeakResultStore#retainArray(gdsc.smlm.results.PeakResult[])
	 */
	@Override
	public boolean retainArray(PeakResult[] results)
	{
		return this.results.retainAll(Arrays.asList(results));
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.PeakResultStore#retainStore(gdsc.smlm.results.PeakResultStore)
	 */
	@Override
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
	@Override
	public void clear()
	{
		results.clear();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.PeakResultStore#trimToSize()
	 */
	@Override
	public void trimToSize()
	{
		//results.trimToSize();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.PeakResultStore#toArray()
	 */
	@Override
	public PeakResult[] toArray()
	{
		return results.toArray(new PeakResult[size()]);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.PeakResultStore#copy()
	 */
	@Override
	public PeakResultStore copy()
	{
		return new SetPeakResultStore(this);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.PeakResultStore#copy(boolean)
	 */
	@Override
	public PeakResultStore copy(boolean deepCopy)
	{
		if (deepCopy)
		{
			SetPeakResultStore copy = new SetPeakResultStore(size());
			for (PeakResult r : results)
				copy.add(r.clone());
			return copy;
		}
		return copy();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.PeakResultStore#removeIf(gdsc.smlm.results.PeakResultPredicate)
	 */
	@Override
	public boolean removeIf(final PeakResultPredicate filter)
	{
		// Delegate to the list implementation
		final ArrayPeakResultStore list = new ArrayPeakResultStore(10);
		for (PeakResult r : results)
			if (filter.test(r))
				list.add(r);
		return this.results.removeAll(Arrays.asList(list.toArray()));
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.PeakResultStore#forEach(gdsc.smlm.results.procedures.PeakResultProcedure)
	 */
	@Override
	public void forEach(PeakResultProcedure procedure)
	{
		for (PeakResult r : results)
			procedure.execute(r);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.PeakResultStore#subset(gdsc.smlm.results.procedures.PeakResultPredicate)
	 */
	@Override
	public PeakResult[] subset(PeakResultPredicate filter)
	{
		final ArrayPeakResultStore list = new ArrayPeakResultStore(10);
		for (PeakResult r : results)
			if (filter.test(r))
				list.add(r);
		return list.toArray();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.PeakResultStore#contains(gdsc.smlm.results.PeakResult)
	 */
	@Override
	public boolean contains(PeakResult result)
	{
		return results.contains(result);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.PeakResultStoreCollection#getCollection()
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
	 * @see gdsc.smlm.results.PeakResultStoreCollection#getCollectionReference()
	 */
	@Override
	public Collection<PeakResult> getCollectionReference()
	{
		return results;
	}
}
