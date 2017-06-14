package gdsc.smlm.results;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Wrapper class to output to multiple results destinations
 */
public class PeakResultsList extends AbstractPeakResults implements PeakResults
{
	private AtomicInteger size = new AtomicInteger(0);
	private List<PeakResults> results = new LinkedList<PeakResults>();

	/**
	 * Add a result format to the output. If a PeakResultsList is passed then it will be
	 * separated into the child PeakResults instances. This will break the size() function
	 * of any input PeakResultsList since only the children will remain within this list.
	 * <p>
	 * Sets the settings (source and configuration) of the child to the same as this list
	 * 
	 * @param peakResults
	 */
	public void addOutput(PeakResults peakResults)
	{
		if (peakResults instanceof PeakResultsList)
		{
			for (PeakResults r : ((PeakResultsList) peakResults).results)
				addOutput(r);
		}
		else
		{
			peakResults.copySettings(this);
			results.add(peakResults);
		}
	}

	/**
	 * @return The number of outputs contained in the list
	 */
	public int numberOfOutputs()
	{
		return results.size();
	}

	/**
	 * @return The outputs
	 */
	public PeakResults[] toArray()
	{
		return results.toArray(new PeakResults[results.size()]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#begin()
	 */
	public void begin()
	{
		size = new AtomicInteger(0);
		for (PeakResults peakResults : results)
			peakResults.begin();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#add(int, int, int, float, double, float, float[], float[])
	 */
	public void add(int peak, int origX, int origY, float origValue, double error, float noise, float[] params,
			float[] paramsStdDev)
	{
		size.incrementAndGet();
		for (PeakResults peakResults : results)
			peakResults.add(peak, origX, origY, origValue, error, noise, params, paramsStdDev);
	}

	public void add(PeakResult result)
	{
		size.incrementAndGet();
		for (PeakResults peakResults : results)
			peakResults.add(result);
	}

	public void addAll(Collection<PeakResult> results)
	{
		size.addAndGet(results.size());
		for (PeakResults peakResults : this.results)
			peakResults.addAll(results);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#size()
	 */
	public int size()
	{
		return size.intValue();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#end()
	 */
	public void end()
	{
		for (PeakResults peakResults : results)
			peakResults.end();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#isActive()
	 */
	public boolean isActive()
	{
		for (PeakResults peakResults : this.results)
			if (peakResults.isActive())
				return true;
		return false;
	}

	/**
	 * Checks all the results in the list. If any are not thread safe then they are wrapped with a
	 * SynchronizedPeakResults container.
	 *
	 * @return the thread safe list
	 */
	public PeakResultsList getThreadSafeList()
	{
		PeakResultsList newList = new PeakResultsList();
		for (PeakResults peakResults : this.results)
		{
			if (!(peakResults instanceof ThreadSafePeakResults))
			{
				peakResults = new SynchronizedPeakResults(peakResults);
			}
			newList.addOutput(peakResults);
		}
		return newList;
	}
}
