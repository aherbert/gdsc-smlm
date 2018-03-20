package gdsc.smlm.results.event;

import java.util.Arrays;
import java.util.LinkedList;

import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResultStore;
import gdsc.smlm.results.SetPeakResultStore;

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
 * Stores peak results and allows event propagation to listeners of the model.
 */
public class PeakResultModel
{
	final PeakResultStore results;
	final LinkedList<PeakResultModelListener> listeners;

	/**
	 * Instantiates a new peak result model.
	 */
	public PeakResultModel()
	{
		this(null);
	}

	/**
	 * Instantiates a new peak result model using the store.
	 * <p>
	 * The default store implementation is a set to avoid duplicate entries in the model.
	 *
	 * @param results
	 *            the results
	 */
	public PeakResultModel(PeakResultStore results)
	{
		if (results == null)
			results = new SetPeakResultStore(10);
		this.results = results;
		listeners = new LinkedList<PeakResultModelListener>();
	}

	/**
	 * Adds the peak result model listener.
	 *
	 * @param listener
	 *            the listener
	 * @return true, if successful
	 */
	public void addPeakResultModelListener(PeakResultModelListener listener)
	{
		if (listener == null || listeners.contains(listener))
			return;
		listeners.add(listener);
	}

	/**
	 * Adds the peak result model listener.
	 *
	 * @param listener
	 *            the listener
	 */
	public void removePeakResultModelListener(PeakResultModelListener listener)
	{
		if (listener == null)
			return;
		listeners.remove(listener);
	}

	/**
	 * Gets the peak result model listeners.
	 *
	 * @return the peak result model listeners
	 */
	public PeakResultModelListener[] getPeakResultModelListeners()
	{
		return listeners.toArray(new PeakResultModelListener[listeners.size()]);
	}

	/**
	 * Convert the model to an array.
	 *
	 * @return the peak result array
	 */
	public PeakResult[] toArray()
	{
		return results.toArray();
	}

	/**
	 * Adds the results to the model.
	 *
	 * @param source
	 *            the source of the event
	 * @param peakResults
	 *            the peak results
	 */
	public void add(Object source, PeakResult... peakResults)
	{
		if (peakResults.length == 0)
			return;
		// Reuse space
		int size = 0;
		for (int i = 0; i < peakResults.length; i++)
		{
			if (results.add(peakResults[i]))
				peakResults[size++] = peakResults[i];
		}
		fireAdded(source, peakResults, size);
	}

	/**
	 * Removes the results from the model.
	 *
	 * @param source
	 *            the source of the event
	 * @param peakResults
	 *            the peak results
	 */
	public void remove(Object source, PeakResult... peakResults)
	{
		if (peakResults.length == 0)
			return;
		// Reuse space
		int size = 0;
		for (int i = 0; i < peakResults.length; i++)
		{
			if (results.remove(peakResults[i]))
				peakResults[size++] = peakResults[i];
		}
		fireRemoved(source, peakResults, size);
	}
	
	/**
	 * Selects the results from the model.
	 *
	 * @param source
	 *            the source of the event
	 * @param peakResults
	 *            the peak results
	 */
	public void select(Object source, PeakResult... peakResults)
	{
		if (peakResults.length == 0)
			return;
		// Reuse space
		int size = 0;
		for (int i = 0; i < peakResults.length; i++)
		{
			if (results.contains(peakResults[i]))
				peakResults[size++] = peakResults[i];
		}
		fireSelected(source, peakResults, size);
	}	
	/**
	 * Checks for listeners.
	 *
	 * @return true, if successful
	 */
	public boolean hasListeners()
	{
		return !listeners.isEmpty();
	}
	
	private void fireAdded(Object source, PeakResult[] peakResults, int size)
	{
		if (size == 0 || listeners.isEmpty())
			return;
		if (size != peakResults.length)
			peakResults = Arrays.copyOf(peakResults, size);
		PeakResultModelEvent e = new PeakResultModelEvent(source, peakResults);
		for (PeakResultModelListener l : listeners)
			l.added(e);
	}

	private void fireRemoved(Object source, PeakResult[] peakResults, int size)
	{
		if (size == 0 || listeners.isEmpty())
			return;
		if (size != peakResults.length)
			peakResults = Arrays.copyOf(peakResults, size);
		PeakResultModelEvent e = new PeakResultModelEvent(source, peakResults);
		for (PeakResultModelListener l : listeners)
			l.removed(e);
	}

	private void fireSelected(Object source, PeakResult[] peakResults, int size)
	{
		// Don't check size is zero as we can select nothing
		if (listeners.isEmpty())
			return;
		if (size != peakResults.length)
			peakResults = Arrays.copyOf(peakResults, size);
		PeakResultModelEvent e = new PeakResultModelEvent(source, peakResults);
		for (PeakResultModelListener l : listeners)
			l.selected(e);
	}
}
