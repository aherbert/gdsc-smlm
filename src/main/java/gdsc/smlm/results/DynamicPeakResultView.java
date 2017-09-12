package gdsc.smlm.results;

import gdsc.smlm.results.predicates.FramePeakResultPredicate;
import gdsc.smlm.results.predicates.IdPeakResultPredicate;

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
 * Provides a dynamic view of the results. Changes to the underlying results store will be reflected in the view.
 */
public class DynamicPeakResultView implements PeakResultView
{
	private final PeakResultStore store;

	/**
	 * Instantiates a new cached peak result view.
	 *
	 * @param store
	 *            the store
	 */
	public DynamicPeakResultView(PeakResultStore store)
	{
		this.store = store;
	}

	/**
	 * Gets the results by frame.
	 *
	 * @param frame
	 *            the frame
	 * @return the results
	 */
	public PeakResult[] getResultsByFrame(int frame)
	{
		return store.subset(new FramePeakResultPredicate(frame));
	}

	/**
	 * Gets the results by id.
	 *
	 * @param id
	 *            the id
	 * @return the results
	 */
	public PeakResult[] getResultsById(int id)
	{
		return store.subset(new IdPeakResultPredicate(id));
	}
}
