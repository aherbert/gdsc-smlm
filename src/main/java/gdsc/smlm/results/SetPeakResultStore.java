package gdsc.smlm.results;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;

import gdsc.smlm.results.predicates.PeakResultPredicate;
import gdsc.smlm.results.procedures.PeakResultProcedure;

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
		this.results = new HashSet<PeakResult>(capacity);
	}

	/**
	 * Instantiates a new array list peak result store.
	 *
	 * @param store
	 *            the store to copy
	 */
	public SetPeakResultStore(SetPeakResultStore store)
	{
		this.results = new HashSet<PeakResult>(store.results);
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
		//results.trimToSize();
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
		return new SetPeakResultStore(this);
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
