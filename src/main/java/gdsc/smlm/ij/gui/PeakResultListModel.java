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
package gdsc.smlm.ij.gui;

import java.util.Arrays;

import javax.swing.AbstractListModel;

import gdsc.smlm.results.ArrayPeakResultStore;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResultStoreList;

/**
 * Stores peak results and allows event propagation to listeners of the model.
 */
public class PeakResultListModel extends AbstractListModel<PeakResult>
{
	private static final long serialVersionUID = -7095827869962490626L;

	final PeakResultStoreList delegate;
	private boolean checkForDuplicates = false;

	/**
	 * Instantiates a new peak result model.
	 */
	public PeakResultListModel()
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
	public PeakResultListModel(PeakResultStoreList results)
	{
		if (results == null)
			results = new ArrayPeakResultStore(10);
		this.delegate = results;
	}

	/**
	 * Convert the model to an array.
	 *
	 * @return the peak result array
	 */
	public PeakResult[] toArray()
	{
		return delegate.toArray();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see javax.swing.ListModel#getSize()
	 */
	@Override
	public int getSize()
	{
		return delegate.size();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see javax.swing.ListModel#getElementAt(int)
	 */
	@Override
	public PeakResult getElementAt(int index)
	{
		return delegate.get(index);
	}

	/**
	 * Adds the results. Duplicates can be avoided using the check for duplicates property.
	 *
	 * @param source
	 *            the source
	 * @param peakResults
	 *            the peak results
	 * @see #isCheckDuplicates()
	 */
	public void add(Object source, PeakResult... peakResults)
	{
		if (peakResults.length == 0)
			return;
		int index0 = delegate.size();
		if (checkForDuplicates)
		{
			int size = 0;
			for (int i = 0; i < peakResults.length; i++)
			{
				if (!delegate.contains(peakResults[i]))
					peakResults[size++] = peakResults[i];
			}
			if (size == 0)
				return;
			if (size != peakResults.length)
				peakResults = Arrays.copyOf(peakResults, size);
		}
		delegate.addArray(peakResults);
		int index1 = delegate.size() - 1;
		fireIntervalAdded(source, index0, index1);
	}

	/**
	 * Removes the result
	 *
	 * @param source
	 *            the source
	 * @param peakResult
	 *            the peak result
	 */
	public void remove(Object source, PeakResult peakResult)
	{
		remove(delegate.indexOf(peakResult));
	}

	/**
	 * Removes the result.
	 *
	 * @param source
	 *            the source
	 * @param index
	 *            the index
	 */
	public void remove(Object source, int index)
	{
		if (index < 0 || index >= delegate.size())
			return;
		delegate.remove(index);
		fireIntervalRemoved(source, index, index);
	}

	/**
	 * Removes the results. This is a convenience method that calls remove for each result.
	 *
	 * @param source
	 *            the source
	 * @param peakResults
	 *            the peak results
	 */
	public void remove(Object source, PeakResult... peakResults)
	{
		if (peakResults.length == 0)
			return;
		for (int i = 0; i < peakResults.length; i++)
		{
			remove(peakResults[i]);
		}
	}

	/**
	 * Clear the results.
	 */
	public void clear()
	{
		int index1 = delegate.size() - 1;
		if (index1 >= 0)
		{
			delegate.clear();
			fireIntervalRemoved(this, 0, index1);
		}
	}

	/**
	 * Deletes the components at the specified range of indexes.
	 * The removal is inclusive, so specifying a range of (1,5)
	 * removes the component at index 1 and the component at index 5,
	 * as well as all components in between.
	 * <p>
	 * Throws an <code>ArrayIndexOutOfBoundsException</code>
	 * if the index was invalid.
	 *
	 * @param fromIndex
	 *            the index of the lower end of the range
	 * @param toIndex
	 *            the index of the upper end of the range
	 */
	public void removeRange(int fromIndex, int toIndex)
	{
		delegate.remove(fromIndex, toIndex);
		fireIntervalRemoved(this, fromIndex, toIndex);
	}

	/**
	 * If true then all results will be checked against the current contents before
	 * allowing an addition.
	 *
	 * @return true, if checking duplicates
	 */
	public boolean isCheckDuplicates()
	{
		return checkForDuplicates;
	}

	/**
	 * Sets the check duplicates flag. If true then all results will be checked against the current contents before
	 * allowing an addition.
	 *
	 * @param checkForDuplicates
	 *            the new check duplicates flag
	 */
	public void setCheckDuplicates(boolean checkForDuplicates)
	{
		this.checkForDuplicates = checkForDuplicates;
	}
}
