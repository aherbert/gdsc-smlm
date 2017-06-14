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
import java.util.TreeSet;

/**
 * Filter the results using the distance to a set of coordinates. Positions must be within the distance threshold.
 * Fitted peaks are selected first and in the event of multiple results the peak with the strongest signal is
 * selected. Otherwise failed starting positions are selected and in the event of multiple results the closest
 * position will be chosen.
 * @deprecated Filtering of the results is no longer supported
 */
public class OptimumDistanceResultFilter extends ResultFilter
{
	private FitResult[] bestFitResults;
	private int[] bestIndices;
	private float[] bestD2;
	private float[] bestSignal;
	private PeakResult[] bestPeakResults;

	public OptimumDistanceResultFilter(List<float[]> filter, float d, int nMaxima)
	{
		super(filter, d, nMaxima);
		bestFitResults = new FitResult[filter.size()];
		bestIndices = new int[filter.size()];
		bestD2 = new float[filter.size()];
		Arrays.fill(bestD2, d2);
		bestSignal = new float[filter.size()];
		bestPeakResults = new PeakResult[filter.size()];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.engine.filter.ResultFilter#filter(gdsc.smlm.fitting.FitResult, int,
	 * gdsc.smlm.results.PeakResult[])
	 */
	@Override
	public void filter(FitResult fitResult, int maxIndex, PeakResult... results)
	{
		for (PeakResult r : results)
		{
			if (r == null)
				continue;
			for (int i = 0; i < filter.size(); i++)
			{
				float[] coord = filter.get(i);
				final float dx = r.getXPosition() - coord[0];
				final float dy = r.getYPosition() - coord[1];
				// Only check if within the distance threshold 
				if (dx * dx + dy * dy < d2)
				{
					// Then filter by signal strength
					float s = r.getSignal();
					if (s < bestSignal[i])
						continue;
					bestFitResults[i] = fitResult;
					bestIndices[i] = maxIndex;
					bestSignal[i] = s;
					bestPeakResults[i] = r;
				}
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.engine.filter.ResultFilter#filter(gdsc.smlm.fitting.FitResult, int, float, float)
	 */
	@Override
	public void filter(FitResult fitResult, int maxIndex, float x, float y)
	{
		for (int i = 0; i < filter.size(); i++)
		{
			// Skip if there is a peak result for this target coordinate
			if (bestPeakResults[i] != null)
				continue;
			float[] coord = filter.get(i);
			final float dx = x - coord[0];
			final float dy = y - coord[1];
			final float dd = dx * dx + dy * dy;
			// Check if this starting position is the closest
			if (dd < bestD2[i])
			{
				bestFitResults[i] = fitResult;
				bestIndices[i] = maxIndex;
				bestD2[i] = dd;
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.engine.filter.ResultFilter#finalise()
	 */
	@Override
	public void finalise()
	{
		// Note that there could be the same result allocated to two target positions
		// so find the unique results
		int[] uniqueIndices = new int[bestIndices.length];
		int unique = 0;
		for (int i = 0; i < bestIndices.length; i++)
		{
			if (bestFitResults[i] == null)
				continue;
			boolean found = false;
			for (int j = unique; j-- > 0;)
			{
				if (bestIndices[uniqueIndices[j]] == bestIndices[i])
				{
					found = true;
					break;
				}
			}
			if (!found)
				uniqueIndices[unique++] = i;
		}

		// The fit results and the indices must match so preserve the same order
		filteredCount = unique;
		filteredFitResults = new FitResult[unique];
		filteredIndices = new int[unique];
		for (int i = 0; i < unique; i++)
		{
			filteredFitResults[i] = bestFitResults[uniqueIndices[i]];
			filteredIndices[i] = bestIndices[uniqueIndices[i]];
		}

		// The peak results can be in any order so use a set to find the unique results
		if (unique > 0)
		{
			TreeSet<PeakResult> set = new TreeSet<PeakResult>();
			for (PeakResult r : bestPeakResults)
			{
				if (r != null)
					set.add(r);
			}

			peakResults = new ArrayList<PeakResult>(set.size());
			peakResults.addAll(set);
		}
		else
		{
			peakResults = new ArrayList<PeakResult>();
		}
	}
}