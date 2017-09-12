package gdsc.smlm.results;

import gdsc.smlm.results.predicates.FramePeakResultPredicate;
import gdsc.smlm.results.predicates.IdPeakResultPredicate;
import gdsc.smlm.results.predicates.PeakResultPredicate;
import gnu.trove.map.hash.TIntObjectHashMap;

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
 * Provides a cache of the view of the results
 */
public class CachedPeakResultView implements PeakResultView
{
	private final PeakResultStore store;

	private TIntObjectHashMap<PeakResult[]> frameMap;
	private TIntObjectHashMap<PeakResult[]> idMap;

	/**
	 * Instantiates a new cached peak result view.
	 *
	 * @param store
	 *            the store
	 */
	public CachedPeakResultView(PeakResultStore store)
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
		if (frameMap == null)
		{
			frameMap = new TIntObjectHashMap<PeakResult[]>();
			return findResults(frameMap, frame, new FramePeakResultPredicate(frame));
		}
		else
		{
			PeakResult[] results = frameMap.get(frame);
			if (results == null)
				results = findResults(frameMap, frame, new FramePeakResultPredicate(frame));
			return results;
		}
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
		if (idMap == null)
		{
			idMap = new TIntObjectHashMap<PeakResult[]>();
			return findResults(idMap, id, new IdPeakResultPredicate(id));
		}
		else
		{
			PeakResult[] results = idMap.get(id);
			if (results == null)
				results = findResults(idMap, id, new IdPeakResultPredicate(id));
			return results;
		}
	}

	private PeakResult[] findResults(TIntObjectHashMap<PeakResult[]> map, int key, PeakResultPredicate filter)
	{
		PeakResult[] results = store.subset(filter);
		map.put(key, results);
		return results;
	}

	/**
	 * Clear the cache
	 */
	public void clear()
	{
		if (frameMap != null)
			frameMap.clear();
		if (idMap != null)
			idMap.clear();
	}
}
