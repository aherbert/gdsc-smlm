package gdsc.smlm.results.filter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;

import gdsc.core.match.FractionClassificationResult;
import gdsc.core.match.FractionalAssignment;
import gdsc.core.match.RankedScoreCalculator;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Support direct filtering of FilterPeakResult objects.
 * <p>
 * Filter a multi-path set of peak results into accepted/rejected.
 */
public abstract class MultiPathFilter extends Filter
{
	/**
	 * Called before the accept method is called for PreprocessedPeakResult
	 * <p>
	 * This should be called once to initialise the filter before processing a batch of results.
	 * 
	 * @see #accept(PreprocessedPeakResult)
	 * @see #accept(MultiPathPeakResult)
	 */
	public void setup()
	{
	}

	/**
	 * Filter the peak result.
	 * 
	 * @param peak
	 *            The peak result
	 * @return true if the peak should be accepted, otherwise false to reject.
	 */
	public abstract boolean accept(PreprocessedPeakResult peak);

	/**
	 * Filter a multi-path set of peak results into a set that are accepted.
	 * 
	 * @param multiPathResult
	 * @return The peak results that are accepted (or null)
	 */
	final public PreprocessedPeakResult[] accept(MultiPathPeakResult multiPathResult)
	{
		// TODO - implement logic for filtering the multiPathResult

		return null;
	}

	/**
	 * Filter a set of multi-path results into a set of results
	 * 
	 * @param results
	 * @return the filtered results
	 */
	final public PreprocessedPeakResult[] filter(MultiPathPeakResult[] results)
	{
		setup();
		ArrayList<PreprocessedPeakResult> list = new ArrayList<PreprocessedPeakResult>(results.length);
		for (MultiPathPeakResult multiPathResult : results)
		{
			PreprocessedPeakResult[] result = accept(multiPathResult);
			if (result != null)
				list.addAll(Arrays.asList(result));
		}
		return list.toArray(new PreprocessedPeakResult[list.size()]);
	}

	/**
	 * Score a set of multi-path results. Filter each multi-path result. Any output results that are
	 * ClassifiedPeakResult are assumed to be positives and their assignments used to score the results per frame.
	 * <p>
	 * The number of consecutive rejections are counted per frame. When the configured number of failures is reached all
	 * remaining results for the frame are rejected. This assumes the results are ordered by the frame.
	 * <p>
	 * Note: The fractional scores are totalled as well as the integer tp/fp scores. These are returned in the positives
	 * and negatives fields of the result.
	 * 
	 * @param resultsList
	 *            a list of results to analyse
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @param n
	 *            The number of actual results
	 * @return the score
	 */
	public FractionClassificationResult fractionScore(MultiPathPeakResult[] results, final int failures, final int n)
	{
		final double[] score = new double[4];

		setup();
		final ArrayList<FractionalAssignment> assignments = new ArrayList<FractionalAssignment>();

		int frame = -1;
		int failCount = 0;
		int nPredicted = 0;
		for (MultiPathPeakResult multiPathResult : results)
		{
			// Reset fail count for new frames
			if (frame != multiPathResult.frame)
			{
				score(assignments, score, nPredicted);
				frame = multiPathResult.frame;
				failCount = 0;
				nPredicted = 0;
			}

			if (failCount <= failures)
			{
				// Assess the result if we are below the fail limit
				PreprocessedPeakResult[] result = accept(multiPathResult);
				final int size = nPredicted;
				if (result != null)
				{
					// For all the results that were returned, check if any are classified results
					// and store the classifications
					for (int i = 0; i < result.length; i++)
					{
						if (result[i] instanceof ClassifiedPeakResult)
						{
							nPredicted++;
							final FractionalAssignment[] a = ((ClassifiedPeakResult) result[i]).getAssignments();
							if (a != null && a.length > 0)
							{
								//list.addAll(Arrays.asList(a));
								assignments.addAll(new DummyCollection(a));
							}
						}
					}
				}
				if (size != nPredicted)
				{
					// More results were accepted
					failCount = 0;
					continue;
				}
			}

			// Nothing was accepted, increment fail count
			failCount++;
		}

		// Score final frame
		score(assignments, score, nPredicted);

		// Note: We are using the integer positives and negatives fields to actually store integer TP and FP
		return new FractionClassificationResult(score[0], score[1], 0, n - score[0], (int) score[2], (int) score[3]);
	}

	private void score(ArrayList<FractionalAssignment> assignments, double[] score, int nPredicted)
	{
		if (assignments.isEmpty())
			return;
		final RankedScoreCalculator calc = new RankedScoreCalculator(
				assignments.toArray(new FractionalAssignment[assignments.size()]));
		final double[] result = calc.score(nPredicted, false);
		score[0] += result[0];
		score[1] += result[1];
		score[2] += result[2];
		score[3] += result[3];
		assignments.clear();
	}

	/**
	 * Create a dummy collection that implements toArray() without cloning for the addAll() method in ArrayList
	 */
	private class DummyCollection implements Collection<FractionalAssignment>
	{
		final FractionalAssignment[] a;

		DummyCollection(FractionalAssignment[] a)
		{
			this.a = a;
		}

		public int size()
		{
			return a.length;
		}

		public boolean isEmpty()
		{
			return false;
		}

		public boolean contains(Object o)
		{
			return false;
		}

		public Iterator<FractionalAssignment> iterator()
		{
			return null;
		}

		public Object[] toArray()
		{
			return a;
		}

		public <T> T[] toArray(T[] a)
		{
			return null;
		}

		public boolean add(FractionalAssignment e)
		{
			return false;
		}

		public boolean remove(Object o)
		{
			return false;
		}

		public boolean containsAll(Collection<?> c)
		{
			return false;
		}

		public boolean addAll(Collection<? extends FractionalAssignment> c)
		{
			return false;
		}

		public boolean removeAll(Collection<?> c)
		{
			return false;
		}

		public boolean retainAll(Collection<?> c)
		{
			return false;
		}

		public void clear()
		{
		}
	}

	// TODO - Create a subset of multi-path results for scoring

	// TODO - Score a subset of multi-path results. Filter each multi-path result. Any output results that are ClassifiedFilterResults are accumulated. The accumulated set are scored once filtering has been performed.
}