package gdsc.smlm.results.filter;

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

import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.Trace;
import gdsc.smlm.results.TraceManager;
import gdsc.smlm.results.TraceManager.TraceMode;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

/**
 * Filter results using a precision threshold. Any results below the lower
 * precision limit are included. Any results above the upper precision limit are
 * excluded. Any results between the limits are included only if they can be
 * traced through time, optionally via other candidates, to a valid result.
 */
public abstract class HysteresisFilter extends Filter
{
	@XStreamAsAttribute
	double searchDistance = 2;
	@XStreamOmitField
	Set<PeakResult> ok;

	protected enum PeakStatus
	{
		OK, CANDIDATE, REJECT
	}

	public HysteresisFilter(double searchDistance)
	{
		this.searchDistance = searchDistance;
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		ok = new HashSet<PeakResult>();

		// Create a set of candidates and valid peaks
		MemoryPeakResults traceResults = new MemoryPeakResults();

		// Initialise peaks to check
		LinkedList<PeakResult> candidates = new LinkedList<PeakResult>();
		for (PeakResult result : peakResults.getResults())
		{
			switch (getStatus(result))
			{
				case OK:
					ok.add(result);
					traceResults.add(result);
					break;
				case CANDIDATE:
					candidates.add(result);
					traceResults.add(result);
					break;
				default:
					break;
			}
		}

		if (candidates.isEmpty())
			return;

		// Find average precision of the candidates and use it for the search
		// distance
		SummaryStatistics stats = new SummaryStatistics();
		for (PeakResult peakResult : candidates)
		{
			stats.addValue(peakResult.getPrecision(peakResults.getNmPerPixel(), peakResults.getGain()));
		}
		double distanceThreshold = stats.getMean() * searchDistance / peakResults.getNmPerPixel();

		// Trace through candidates
		TraceManager tm = new TraceManager(traceResults);
		tm.setTraceMode(TraceMode.LATEST_FORERUNNER);
		tm.traceMolecules(distanceThreshold, 1);
		Trace[] traces = tm.getTraces();

		for (Trace trace : traces)
		{
			if (trace.size() > 1)
			{
				// Check if the trace touches a valid point
				boolean isOk = false;
				for (PeakResult result : trace.getPoints())
				{
					if (ok.contains(result))
					{
						isOk = true;
						break;
					}
					ok.add(result);
				}
				// Add the entire trace to the OK points
				if (isOk)
				{
					for (PeakResult result : trace.getPoints())
					{
						ok.add(result);
					}
				}
			}
		}
	}

	protected abstract PeakStatus getStatus(PeakResult result);

	/**
	 * @throws NullPointerException
	 *             if not first initialised with a call to {@link #setup(MemoryPeakResults)}
	 * @see gdsc.smlm.results.filter.Filter#accept(gdsc.smlm.results.PeakResult)
	 */
	@Override
	public boolean accept(PeakResult peak) throws NullPointerException
	{
		return ok.contains(peak);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getDescription()
	 */
	@Override
	public String getDescription()
	{
		return "Any results between the limits (candidates) are included only if they can be traced "
				+ "through time, potentially via other candidates, to a valid result. The distance used for "
				+ "tracing is the search distance multiplied by the average precision of the candidates.";
	}
}