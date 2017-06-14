package gdsc.smlm.engine.filter;

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

import gdsc.smlm.fitting.FitResult;
import gdsc.smlm.results.PeakResult;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Filter the results using the distance to a set of coordinates. Any fitted position within the distance to the
 * target coordinates is accepted.
 * @deprecated Filtering of the results is no longer supported
 */
public class DistanceResultFilter extends ResultFilter
{
	public DistanceResultFilter(List<float[]> filter, float d, int nMaxima)
	{
		super(filter, d, nMaxima);
		filteredFitResults = new FitResult[nMaxima];
		filteredIndices = new int[nMaxima];
		peakResults = new ArrayList<PeakResult>(nMaxima);
	}

	/* (non-Javadoc)
	 * @see gdsc.smlm.engine.filter.ResultFilter#filter(gdsc.smlm.fitting.FitResult, int, gdsc.smlm.results.PeakResult[])
	 */
	@Override
	public void filter(FitResult fitResult, int maxIndex, PeakResult... results)
	{
		boolean found = false;
		for (PeakResult r : results)
		{
			if (r == null)
				continue;
			for (float[] coord : filter)
			{
				final float dx = r.getXPosition() - coord[0];
				final float dy = r.getYPosition() - coord[1];
				if (dx * dx + dy * dy < d2)
				{
					found = true;
					peakResults.add(r);
					break;
				}
			}
		}
		if (found)
		{
			// Add the result and the fitted index to the filtered results
			filteredFitResults[filteredCount] = fitResult;
			filteredIndices[filteredCount] = maxIndex;
			filteredCount++;
		}
	}

	/* (non-Javadoc)
	 * @see gdsc.smlm.engine.filter.ResultFilter#filter(gdsc.smlm.fitting.FitResult, int, float, float)
	 */
	@Override
	public void filter(FitResult fitResult, int maxIndex, float x, float y)
	{
		boolean found = false;
		for (float[] coord : filter)
		{
			final float dx = x - coord[0];
			final float dy = y - coord[1];
			if (dx * dx + dy * dy < d2)
			{
				found = true;
				break;
			}
		}
		if (found)
		{
			// Add the result and the fitted index to the filtered results
			filteredFitResults[filteredCount] = fitResult;
			filteredIndices[filteredCount] = maxIndex;
			filteredCount++;
		}
	}

	/* (non-Javadoc)
	 * @see gdsc.smlm.engine.filter.ResultFilter#finalise()
	 */
	@Override
	public void finalise()
	{
		filteredFitResults = Arrays.copyOf(filteredFitResults, filteredCount);
		filteredIndices = Arrays.copyOf(filteredIndices, filteredCount);
	}
}