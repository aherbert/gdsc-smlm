package gdsc.smlm.results.filter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;

import gdsc.core.match.FractionClassificationResult;
import gdsc.core.match.FractionalAssignment;
import gdsc.core.match.RankedScoreCalculator;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.results.filter.MultiPathFitResult.FitResult;

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
 * Filter a multi-path set of peak results into accepted/rejected.
 */
public class MultiPathFilter
{
	/**
	 * The direct filter to apply to the results
	 */
	@XStreamAsAttribute
	final DirectFilter filter;

	/**
	 * The residuals threshold to consider the residuals Quadrant Analysis (QA) score of a single for doublet fitting.
	 * The score should range from 0 to 1. A score of 1 will ignore doublet fitting
	 */
	@XStreamAsAttribute
	final double residualsThreshold;

	/**
	 * Create a new MultiPathFilter
	 * 
	 * @param filter
	 *            the direct filter for filtering the results
	 * @param residualsThreshold
	 *            The residuals threshold to consider a single fit for doublet fitting
	 */
	public MultiPathFilter(DirectFilter filter, double residualsThreshold)
	{
		this.filter = filter;
		this.residualsThreshold = residualsThreshold;
	}

	/**
	 * Called before the accept method is called for PreprocessedPeakResult. This calls the setup() method in the
	 * DirectFilter.
	 * <p>
	 * This should be called once to initialise the filter before processing a batch of results.
	 * 
	 * @see #accept(PreprocessedPeakResult)
	 */
	public void setup()
	{
		filter.setup();		
	}

	/**
	 * Called before the accept method is called for PreprocessedPeakResult. The flags can control the type of filtering
	 * requested. Filters are asked to respect the flags defined in this class. This calls the setup(int) method in the
	 * DirectFilter.
	 * <p>
	 * This should be called once to initialise the filter before processing a batch of results.
	 * 
	 * @param flags
	 *            Flags used to control the filter
	 * @see #accept(PreprocessedPeakResult)
	 */
	public void setup(final int flags)
	{
		filter.setup(flags);
	}

	/**
	 * Filter the peak result. This calls the accept() method in the DirectFilter.
	 * 
	 * @param peak
	 *            The peak result
	 * @return true if the peak should be accepted, otherwise false to reject.
	 */
	public boolean accept(final PreprocessedPeakResult peak)
	{
		return filter.accept(peak);
	}

	/**
	 * Filter a multi-path set of peak results into a set that are accepted.
	 * <p>
	 * Any existing or new results must pass the {@link #accept(PreprocessedPeakResult)} method. Any other
	 * results are
	 * assumed to be candidates that were fitted but will not be validated.
	 * 
	 * @param multiPathResult
	 * @return The new peak results that are accepted (or null)
	 */
	final public PreprocessedPeakResult[] accept(final MultiPathFitResult multiPathResult)
	{
		PreprocessedPeakResult[] results = null;

		// Filter multi-fit
		if ((results = acceptAll(multiPathResult.multiFitResult)) != null)
			return results;

		// Filter single-fit
		if ((results = acceptAll(multiPathResult.singleFitResult)) == null)
		{
			// The fit was not accepted. However it may have been rejected for being too wide
			// and is suitable for a doublet fit.

			// Check there is a result for the single spot
			if (multiPathResult.singleFitResult.status != FitStatus.OK)
				return null;

			// Check if the residuals score is below the configured threshold
			if (multiPathResult.singleQAScore < residualsThreshold)
				return null;

			// Get the single spot
			final PreprocessedPeakResult singleResult = extractFirstNew(multiPathResult.singleFitResult.results);
			if (singleResult == null)
				return null;

			// Check the width is reasonable given the size of the fitted region.
			//@formatter:off
			if (
					singleResult.getXSDFactor() < 1 || // Not a wide spot
					singleResult.getXSD() > multiPathResult.width * 0.5f || // width covers more than the region
					singleResult.getYSDFactor() < 1 || // Not a wide spot
					singleResult.getYSD() > multiPathResult.height * 0.5f // width covers more than the region
				)
				return null;
			//@formatter:on

			// We must validate the spot without width filtering
			setup(DirectFilter.NO_WIDTH);

			try
			{
				if (!accept(singleResult))
					// This is still a bad single result, without width filtering
					return null;
			}
			finally
			{
				// reset
				setup();
			}
		}
		else
		{
			// The single fit is OK.
			// If not eligible for a doublet fit then return.
			if (multiPathResult.singleQAScore < residualsThreshold)
				return results;
		}

		// We reached here with a single fit that is eligible for doublet fitting

		// Filter doublet fit
		final PreprocessedPeakResult[] doubletResults = acceptAny(multiPathResult.doubletFitResult);
		if (doubletResults != null)
			return doubletResults;

		return results;
	}

	/**
	 * Check all new and all existing results are valid. Returns the new results
	 * 
	 * @param fitResult
	 *            the results
	 * @return The new results that pass the filter
	 */
	private PreprocessedPeakResult[] acceptAll(final FitResult fitResult)
	{
		if (fitResult == null || fitResult.results == null)
			return null;
		final PreprocessedPeakResult[] results = fitResult.results;

		// All new and existing results should be valid
		int count = 0;
		for (int i = 0; i < results.length; i++)
		{
			if (results[i].isNewResult())
			{
				// All new results must pass
				if (!accept(results[i]))
					return null;
				count++;
			}
			else if (results[i].isExistingResult())
			{
				// All existing results must pass
				if (!accept(results[i]))
					return null;
			}
		}

		if (count == 0)
			return null;

		// Return the new results
		final PreprocessedPeakResult[] filtered = new PreprocessedPeakResult[count];
		count = 0;
		for (int i = 0; i < results.length; i++)
		{
			if (results[i].isNewResult())
				filtered[count++] = results[i];
		}
		return filtered;
	}

	/**
	 * Check any new and all existing results are valid. Returns the new results
	 * 
	 * @param fitResult
	 *            the results
	 * @return The new results that pass the filter
	 */
	private PreprocessedPeakResult[] acceptAny(final FitResult fitResult)
	{
		if (fitResult == null || fitResult.results == null)
			return null;
		final PreprocessedPeakResult[] results = fitResult.results;

		// Any new and all existing results should be valid
		int count = 0;
		final int[] ok = new int[results.length];
		for (int i = 0; i < results.length; i++)
		{
			if (results[i].isNewResult())
			{
				// Any new result that pass are OK
				if (accept(results[i]))
					ok[count++] = i;
			}
			else if (results[i].isExistingResult())
			{
				// All existing results must pass
				if (!accept(results[i]))
					return null;
			}
		}

		if (count == 0)
			return null;

		// Return the new results
		final PreprocessedPeakResult[] filtered = new PreprocessedPeakResult[count];
		for (int i = 0; i < count; i++)
		{
			filtered[i] = results[ok[i]];
		}
		return filtered;
	}

	private PreprocessedPeakResult extractFirstNew(PreprocessedPeakResult[] results)
	{
		if (results == null)
			return null;
		for (int i = 0; i < results.length; i++)
			if (results[i].isNewResult())
				return results[i];
		return null;
	}

	/**
	 * Filter a set of multi-path results into a set of results
	 * 
	 * @param results
	 * @return the filtered results
	 */
	final public PreprocessedPeakResult[] filter(final MultiPathFitResult[] results)
	{
		setup();
		ArrayList<PreprocessedPeakResult> list = new ArrayList<PreprocessedPeakResult>(results.length);
		for (MultiPathFitResult multiPathResult : results)
		{
			PreprocessedPeakResult[] result = accept(multiPathResult);
			if (result != null)
				list.addAll(Arrays.asList(result));
		}
		return list.toArray(new PreprocessedPeakResult[list.size()]);
	}

	/**
	 * Score a set of multi-path results. Filter each multi-path result. Any output results that are
	 * new results are assumed to be positives and their assignments used to score the results per frame.
	 * <p>
	 * The number of consecutive rejections are counted per frame. When the configured number of failures is reached all
	 * remaining results for the frame are rejected. This assumes the results are ordered by the frame.
	 * <p>
	 * Note: The fractional scores are totalled as well as the integer tp/fp scores. These are returned in the positives
	 * and negatives fields of the result.
	 * 
	 * @param results
	 *            a set of results to analyse
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @param n
	 *            The number of actual results
	 * @return the score
	 */
	public FractionClassificationResult fractionScore(final MultiPathFitResult[] results, final int failures,
			final int n)
	{
		return fractionScore(results, failures, n, true);
	}

	/**
	 * Create a subset of multi-path results, i.e. all those that pass the filter.
	 * <p>
	 * The number of consecutive rejections are counted per frame. When the configured number of failures is reached all
	 * remaining results for the frame are rejected. This assumes the results are ordered by the frame.
	 * <p>
	 * The number of failures before each result is stored in the failCount property of the MultiPathPeakResult.
	 * 
	 * @param results
	 *            a set of results to analyse
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @return the filtered results
	 */
	public MultiPathFitResult[] filterSubset(final MultiPathFitResult[] results, final int failures)
	{
		final MultiPathFitResult[] newResults = new MultiPathFitResult[results.length];
		int size = 0;

		setup();
		int frame = -1;
		int failCount = 0;
		for (MultiPathFitResult multiPathResult : results)
		{
			// Reset fail count for new frames
			if (frame != multiPathResult.frame)
			{
				frame = multiPathResult.frame;
				failCount = 0;
			}

			if (failCount <= failures)
			{
				// Assess the result if we are below the fail limit
				final PreprocessedPeakResult[] result = accept(multiPathResult);
				if (result != null)
				{
					// This has valid result so add to the output subset 
					newResults[size++] = multiPathResult;
					// Store the number of failures before this result
					multiPathResult.failCount = failCount;
					// Reset fail count
					failCount = 0;
				}
				else
				{
					// This was rejected, increment fail count
					failCount++;
				}
			}
		}

		return Arrays.copyOf(newResults, size);
	}

	/**
	 * Score a subset of multi-path results. The subset should be created with
	 * {@link #filterSubset(MultiPathFitResult[], int)}.
	 * <p>
	 * Filter each multi-path result. Any output results that are new results are assumed to be positives and
	 * their assignments used to score the results per frame.
	 * <p>
	 * The number of consecutive rejections are counted per frame. When the configured number of failures is reached all
	 * remaining results for the frame are rejected. This assumes the results are ordered by the frame.
	 * <p>
	 * Note: The fractional scores are totalled as well as the integer tp/fp scores. These are returned in the positives
	 * and negatives fields of the result.
	 * 
	 * @param results
	 *            a set of results to analyse
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @param n
	 *            The number of actual results
	 * @return the score
	 */
	public FractionClassificationResult fractionScoreSubset(final MultiPathFitResult[] results, final int failures,
			final int n)
	{
		return fractionScore(results, failures, n, true);
	}

	/**
	 * Score a set of multi-path results. If the set is a subset then the fail count will be accumulated using the
	 * failCount property of the MultiPathPeakResult.
	 * <p>
	 * Filter each multi-path result. Any output results that are new results are assumed to be positives and
	 * their assignments used to score the results per frame.
	 * <p>
	 * The number of consecutive rejections are counted per frame. When the configured number of failures is reached all
	 * remaining results for the frame are rejected. This assumes the results are ordered by the frame.
	 * <p>
	 * Note: The fractional scores are totalled as well as the integer tp/fp scores. These are returned in the positives
	 * and negatives fields of the result.
	 * 
	 * @param results
	 *            a set of results to analyse
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @param n
	 *            The number of actual results
	 * @param subset
	 *            True if a subset
	 * @return the score
	 */
	private FractionClassificationResult fractionScore(final MultiPathFitResult[] results, final int failures,
			final int n, final boolean subset)
	{
		final double[] score = new double[4];
		final ArrayList<FractionalAssignment> assignments = new ArrayList<FractionalAssignment>();

		setup();
		int frame = -1;
		int failCount = 0;
		int nPredicted = 0;
		for (MultiPathFitResult multiPathResult : results)
		{
			// Reset fail count for new frames
			if (frame != multiPathResult.frame)
			{
				score(assignments, score, nPredicted);
				frame = multiPathResult.frame;
				failCount = 0;
				nPredicted = 0;
			}

			// Include the number of failures before this result from the larger set
			if (subset)
				failCount += multiPathResult.failCount;

			if (failCount <= failures)
			{
				// Assess the result if we are below the fail limit
				final PreprocessedPeakResult[] result = accept(multiPathResult);
				final int size = nPredicted;
				if (result != null)
				{
					// For all the results that were returned, check if any are classified results
					// and store the classifications
					for (int i = 0; i < result.length; i++)
					{
						if (result[i].isNewResult())
						{
							final FractionalAssignment[] a = result[i].getAssignments(nPredicted++);
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
					// More results were accepted so reset the fail count
					failCount = 0;
				}
				else
				{
					// Nothing was accepted, increment fail count
					failCount++;
				}
			}
		}

		// Score final frame
		score(assignments, score, nPredicted);

		// Note: We are using the integer positives and negatives fields to actually store integer TP and FP
		return new FractionClassificationResult(score[0], score[1], 0, n - score[0], (int) score[2], (int) score[3]);
	}

	/**
	 * Score the assignments (TP/FP) and then clear the list
	 * 
	 * @param assignments
	 *            The assignments
	 * @param score
	 *            Scores array to accumilate TP/FP scores
	 * @param nPredicted
	 *            The number of predictions
	 */
	private void score(final ArrayList<FractionalAssignment> assignments, final double[] score, final int nPredicted)
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

		DummyCollection(final FractionalAssignment[] a)
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
}