package gdsc.smlm.results.event;

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
 * The listener interface for receiving peakResultModel events.
 * The class that is interested in processing a peakResultModel
 * event implements this interface, and the object created
 * with that class is registered with a component using the
 * component's <code>addPeakResultModelListener<code> method. When
 * the peakResultModel event occurs, that object's appropriate
 * method is invoked.
 *
 * @see PeakResultModelEvent
 */
public interface PeakResultModelListener
{
	/**
	 * Invoked when results have been added to the model
	 *
	 * @param e
	 *            the event
	 */
	public void added(PeakResultModelEvent e);

	/**
	 * Invoked when results have been removed to the model
	 *
	 * @param e
	 *            the event
	 */
	public void removed(PeakResultModelEvent e);

	/**
	 * Invoked when results have been selected in the model
	 *
	 * @param e
	 *            the event
	 */
	public void selected(PeakResultModelEvent e);
}