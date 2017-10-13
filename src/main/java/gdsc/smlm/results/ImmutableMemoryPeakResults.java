package gdsc.smlm.results;

import java.util.Collection;

import gdsc.core.data.DataException;
import gdsc.smlm.results.predicates.PeakResultPredicate;

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
 * Wraps peak results in memory and prevents modification of the results size.
 * <p>
 * Any method that modifies the size of the results set will throw a data exception.
 */
public class ImmutableMemoryPeakResults extends MemoryPeakResults
{
	/**
	 * Instantiates a new immutable memory peak results with the original results store.
	 *
	 * @param results
	 *            the results
	 */
	public ImmutableMemoryPeakResults(MemoryPeakResults results)
	{
		super(results.results);
	}

	/**
	 * Instantiates a new immutable memory peak results with an optional copy of the results store.
	 *
	 * @param results
	 *            the results
	 * @param copy
	 *            the copy
	 */
	public ImmutableMemoryPeakResults(MemoryPeakResults results, boolean copy)
	{
		super((copy) ? results.results.copy() : results.results);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.MemoryPeakResults#add(gdsc.smlm.results.PeakResult)
	 */
	public void add(PeakResult result)
	{
		throw new DataException("This results set is immutable");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.MemoryPeakResults#addAll(java.util.Collection)
	 */
	public void addAll(Collection<PeakResult> results)
	{
		throw new DataException("This results set is immutable");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.MemoryPeakResults#addAll(gdsc.smlm.results.PeakResult[])
	 */
	public void addAll(PeakResult[] results)
	{
		throw new DataException("This results set is immutable");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.MemoryPeakResults#addAll(gdsc.smlm.results.PeakResultStore)
	 */
	public void addAll(PeakResultStore results)
	{
		throw new DataException("This results set is immutable");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.MemoryPeakResults#add(gdsc.smlm.results.MemoryPeakResults)
	 */
	public void add(MemoryPeakResults results)
	{
		throw new DataException("This results set is immutable");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.MemoryPeakResults#removeNullResults()
	 */
	public void removeNullResults()
	{
		throw new DataException("This results set is immutable");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.MemoryPeakResults#removeIf(gdsc.smlm.results.predicates.PeakResultPredicate)
	 */
	public boolean removeIf(PeakResultPredicate filter)
	{
		throw new DataException("This results set is immutable");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.MemoryPeakResults#begin()
	 */
	public void begin()
	{
		throw new DataException("This results set is immutable");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.MemoryPeakResults#add(int, int, int, float, double, float, float[], float[])
	 */
	public void add(int peak, int origX, int origY, float origValue, double chiSquared, float noise, float[] params,
			float[] paramsStdDev)
	{
		throw new DataException("This results set is immutable");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.MemoryPeakResults#end()
	 */
	public void end()
	{
		throw new DataException("This results set is immutable");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.MemoryPeakResults#isActive()
	 */
	public boolean isActive()
	{
		throw new DataException("This results set is immutable");
	}
}
