package gdsc.smlm.results.event;

import java.util.EventObject;

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

import gdsc.smlm.results.PeakResult;

/**
 * The class containing information about the state changes in the PeakResultModel.
 */
public class PeakResultModelEvent extends EventObject
{
	private static final long serialVersionUID = -3191865242109194826L;
	
	private final PeakResult[] peakResults;

	/**
	 * Instantiates a new peak result model event.
	 *
	 * @param source the source
	 * @param peakResults the peak results
	 */
	public PeakResultModelEvent(Object source, PeakResult... peakResults)
	{
		super(source);
		this.peakResults = peakResults;
	}

	/**
	 * Gets the peak results.
	 *
	 * @return the peak results
	 */
	public PeakResult[] getPeakResults()
	{
		return peakResults;
	}
}