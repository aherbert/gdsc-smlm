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
package gdsc.smlm.engine.filter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import gdsc.smlm.fitting.FitResult;
import gdsc.smlm.results.PeakResult;

/**
 * Filter the results using the distance to a set of coordinates. Any fitted position within the distance to the
 * target coordinates is accepted.
 *
 * @deprecated Filtering of the results is no longer supported
 */
@Deprecated
public class DistanceResultFilter extends ResultFilter
{
	/**
	 * Instantiates a new distance result filter.
	 *
	 * @param filter
	 *            the filter
	 * @param d
	 *            the d
	 * @param nMaxima
	 *            the n maxima
	 */
	public DistanceResultFilter(List<float[]> filter, float d, int nMaxima)
	{
		super(filter, d, nMaxima);
		filteredFitResults = new FitResult[nMaxima];
		filteredIndices = new int[nMaxima];
		peakResults = new ArrayList<>(nMaxima);
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

	/*
	 * (non-Javadoc)
	 *
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

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.engine.filter.ResultFilter#finalise()
	 */
	@Override
	public void finalise()
	{
		filteredFitResults = Arrays.copyOf(filteredFitResults, filteredCount);
		filteredIndices = Arrays.copyOf(filteredIndices, filteredCount);
	}
}
