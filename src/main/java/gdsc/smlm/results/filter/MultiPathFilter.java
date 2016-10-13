package gdsc.smlm.results.filter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;

import gdsc.core.match.FractionClassificationResult;
import gdsc.core.match.FractionalAssignment;
import gdsc.core.match.RankedScoreCalculator;
import gdsc.core.utils.NotImplementedException;
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
	 * Stores the results that were accepted when filtering a multi-path result. Also stores the fit result that was
	 * used to select the results.
	 */
	public class SelectedResult
	{
		final public PreprocessedPeakResult[] results;
		final public FitResult fitResult;

		public SelectedResult(PreprocessedPeakResult[] results, FitResult fitResult)
		{
			this.results = results;
			this.fitResult = fitResult;
		}
	}

	/**
	 * Allows storage of results that have been selected during multi-path filtering.
	 */
	public interface SelectedResultStore
	{
		/**
		 * Add a selected result to the store
		 * 
		 * @param selectedResult
		 */
		void add(SelectedResult selectedResult);

		/**
		 * Checks if is fit. Any candidate that has already been fit will not be stored.
		 *
		 * @param candidateId
		 *            the candidate id
		 * @return True if the candidate has been fit
		 */
		boolean isFit(int candidateId);

		/**
		 * Checks if is valid.
		 * <p>
		 * Return true if this candidate should definitely be filtered.
		 *
		 * @param candidateId
		 *            the candidate id
		 * @return true, if is valid
		 */
		boolean isValid(int candidateId);

		/**
		 * A result that passed the primary filter
		 * 
		 * @param result
		 */
		void pass(PreprocessedPeakResult result);

		/**
		 * A result that passed the minimal filter
		 * 
		 * @param result
		 */
		void passMin(PreprocessedPeakResult result);
	}

	/**
	 * Allow tracking of candidates that have been fit
	 */
	private class SimpleSelectedResultStore implements SelectedResultStore
	{
		boolean[] isFit;
		boolean[] isValid;

		SimpleSelectedResultStore()
		{
			isFit = new boolean[0];
			isValid = new boolean[0];
		}

		SimpleSelectedResultStore(int totalCandidates)
		{
			isFit = new boolean[totalCandidates];
			isValid = new boolean[totalCandidates];
		}

		public void add(SelectedResult selectedResult)
		{
		}

		public boolean isFit(int candidateId)
		{
			return isFit[candidateId];
		}

		public boolean isValid(int candidateId)
		{
			return isValid[candidateId];
		}

		public void pass(PreprocessedPeakResult result)
		{
			if (result.isNewResult())
				isFit[result.getCandidateId()] = true;
			// This an existing result or candidate. Mark as valid so candidates will be processed
			isValid[result.getCandidateId()] = true;
		}

		public void passMin(PreprocessedPeakResult result)
		{
			// Passing the minimal filter does not mean it is valid. This would be used to store
			// a fit estimate during processing for this candidate.
		}

		public void resize(int totalCandidates)
		{
			if (isFit.length < totalCandidates)
			{
				isFit = new boolean[totalCandidates];
				isValid = new boolean[totalCandidates];
			}
			else
			{
				Arrays.fill(isFit, 0, totalCandidates, false);
				Arrays.fill(isValid, 0, totalCandidates, false);
			}
		}
	}

	/**
	 * Used to return default behaviour for acceptAny/acceptAll
	 */
	private static class NullSelectedResultStore implements SelectedResultStore
	{
		public void add(SelectedResult selectedResult)
		{

		}

		public boolean isFit(int candidateId)
		{
			// Make sure non-candidate fits are ignored.
			return true;
		}

		public boolean isValid(int candidateId)
		{
			return false;
		}

		public void pass(PreprocessedPeakResult result)
		{

		}

		public void passMin(PreprocessedPeakResult result)
		{

		}
	}

	private static NullSelectedResultStore nullStore = new NullSelectedResultStore();

	/**
	 * The direct filter to apply to the results
	 */
	final IDirectFilter filter;

	/**
	 * The minimal direct filter to apply to the results.
	 * <p>
	 * This is applied if the result fails the primary filter. It is used to indicate that the result achieves a minimum
	 * set of criteria.
	 */
	final IDirectFilter minFilter;

	/**
	 * The residuals threshold to consider the residuals Quadrant Analysis (QA) score of a single for doublet fitting.
	 * The score should range from 0 to 1. A score equal or above 1 will ignore doublet fitting.
	 */
	@XStreamAsAttribute
	final public double residualsThreshold;

	/**
	 * Create a new MultiPathFilter.
	 *
	 * @param filter
	 *            the direct filter for filtering the results
	 * @param residualsThreshold
	 *            The residuals threshold to consider a single fit for doublet fitting
	 */
	public MultiPathFilter(IDirectFilter filter, double residualsThreshold)
	{
		this(filter, null, residualsThreshold);
	}

	/**
	 * Create a new MultiPathFilter.
	 *
	 * @param filter
	 *            the direct filter for filtering the results
	 * @param minFilter
	 *            the minimal direct filter for filtering the results
	 * @param residualsThreshold
	 *            The residuals threshold to consider a single fit for doublet fitting
	 */
	public MultiPathFilter(IDirectFilter filter, IDirectFilter minFilter, double residualsThreshold)
	{
		this.filter = filter;
		this.minFilter = minFilter;
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
		if (minFilter != null)
			minFilter.setup();
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
		if (minFilter != null)
			minFilter.setup(flags);
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
	 * Filter the peak result. This calls the accept() method in the minimal DirectFilter.
	 * 
	 * @param peak
	 *            The peak result
	 * @return true if the peak should be accepted, otherwise false to reject.
	 */
	private boolean minAccept(final PreprocessedPeakResult peak)
	{
		return minFilter.accept(peak);
	}

	/**
	 * Filter a multi-path set of peak results into a set that are accepted.
	 * <p>
	 * Any existing or new results must pass the {@link #accept(PreprocessedPeakResult)} method. Any other
	 * results are assumed to be candidates that were fitted but will not be validated unless required.
	 * <p>
	 * Note that new results may not be for the candidate identified by the MultiPathFitResult. This can
	 * happen when multi-fitting has fit another candidate that previously did not have a result. The
	 * SelectedResultStore is used to determine if that result has been fit already. If not it is added
	 * to the output list.
	 * <p>
	 * The SelectedResultStore will be passed any result that passes the configured filters.
	 *
	 * @param multiPathResult
	 *            the multi path result
	 * @param validateCandidates
	 *            Set to true to validate the candidates
	 * @param store
	 *            the store
	 * @return The new peak results that are accepted (and any valid candidates if found); or null
	 */
	final public PreprocessedPeakResult[] accept(final MultiPathFitResult multiPathResult, boolean validateCandidates,
			SelectedResultStore store)
	{
		PreprocessedPeakResult[] results = null;
		final int candidateId = multiPathResult.candidateId;

		// Ensure we don't have to check the store in acceptALl/acceptAny
		if (store == null)
			store = nullStore;

		// Filter multi-fit
		if ((results = acceptAll(candidateId, multiPathResult.getMultiFitResult(), validateCandidates, store)) != null)
			return results;

		// Filter single-fit
		if ((results = acceptAll(candidateId, multiPathResult.getSingleFitResult(), false, null)) == null)
		{
			// The fit was not accepted. However it may have been rejected for being too wide
			// and is suitable for a doublet fit.
			
			
			// TODO - What if the single fit drifts to another spot?
			

			// Check there is a result for the single spot
			if (multiPathResult.getSingleFitResult().status != 0)
				return null;

			// Check if the residuals score is below the configured threshold
			if (residualsThreshold >= 1 || multiPathResult.getSingleQAScore() < residualsThreshold)
				return null;

			// Get the single spot
			final PreprocessedPeakResult singleResult = extractFirstNew(multiPathResult.getSingleFitResult().results);
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
			if (residualsThreshold >= 1 || multiPathResult.getSingleQAScore() < residualsThreshold)
				return results;
		}

		// We reached here with a single fit that is eligible for doublet fitting

		// Filter doublet fit
		final PreprocessedPeakResult[] doubletResults = acceptAny(candidateId, multiPathResult.getDoubletFitResult(),
				validateCandidates, store);
		if (doubletResults != null)
			return doubletResults;

		return results;
	}

	/**
	 * Filter a multi-path set of peak results into a set that are accepted.
	 * <p>
	 * Any existing or new results must pass the {@link #accept(PreprocessedPeakResult)} method. Any other
	 * results are assumed to be candidates that were fitted but will not be validated unless required.
	 * <p>
	 * Note that new results may not be for the candidate identified by the MultiPathFitResult. This can
	 * happen when multi-fitting has fit another candidate that previously did not have a result. The
	 * SelectedResultStore is used to determine if that result has been fit already. If not it is added
	 * to the output list.
	 * <p>
	 * The method returns the the same results as {@link #accept(MultiPathFitResult, boolean)} but includes the
	 * FitResult that the data originated from.
	 * <p>
	 * The SelectedResultStore will be passed any result that passes the configured filters. It will not be passed the
	 * returned SelectedResult as the results will be duplicates of those passed to the store individually. They
	 * may also contain validated candidates. The returned results must thus be filtered for new results (e.g. not
	 * existing or candidate results).
	 *
	 * @param multiPathResult
	 *            the multi path result
	 * @param validateCandidates
	 *            Set to true to validate the candidates
	 * @param store
	 *            the store
	 * @return The results that are accepted; or null
	 */
	final public SelectedResult select(final MultiPathFitResult multiPathResult, boolean validateCandidates,
			SelectedResultStore store)
	{
		PreprocessedPeakResult[] results = null;
		final int candidateId = multiPathResult.candidateId;

		// Ensure we don't have to check the store in acceptALl/acceptAny
		if (store == null)
			store = nullStore;

		// Filter multi-fit
		if ((results = acceptAll(candidateId, multiPathResult.getMultiFitResult(), validateCandidates, store)) != null)
			return new SelectedResult(results, multiPathResult.getMultiFitResult());

		// Filter single-fit
		if ((results = acceptAll(candidateId, multiPathResult.getSingleFitResult(), false, store)) == null)
		{
			// The fit was not accepted. However it may have been rejected for being too wide
			// and is suitable for a doublet fit.

			// Check there is a result for the single spot
			if (multiPathResult.getSingleFitResult().status != 0)
				return null;

			// Check if the residuals score is below the configured threshold
			if (residualsThreshold >= 1 || multiPathResult.getSingleQAScore() < residualsThreshold)
				return null;

			// Get the single spot
			final PreprocessedPeakResult singleResult = extractFirstNew(multiPathResult.getSingleFitResult().results);
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
			if (residualsThreshold >= 1 || multiPathResult.getSingleQAScore() < residualsThreshold)
				return new SelectedResult(results, multiPathResult.getSingleFitResult());
		}

		// We reached here with a single fit that is eligible for doublet fitting

		// Filter doublet fit
		final PreprocessedPeakResult[] doubletResults = acceptAny(candidateId, multiPathResult.getDoubletFitResult(),
				validateCandidates, store);
		if (doubletResults != null)
			return new SelectedResult(doubletResults, multiPathResult.getDoubletFitResult());

		return new SelectedResult(results, multiPathResult.getSingleFitResult());
	}

	/**
	 * Select a set of peak results.
	 * <p>
	 * The number of consecutive rejections are counted. When the configured number of failures is reached all
	 * remaining results are rejected.
	 * <p>
	 * A selected result will be stored for each MultiPathFitResult that is assessed, even if the fitting failed. In
	 * this case the list of accepted results will be null.
	 * <p>
	 * The SelectedResultStore can be used to track results that pass validation. If this is null then the default
	 * behaviour is to track fitted candidates that pass validation. These will be processed even if the fail count has
	 * been reached.
	 *
	 * @param multiPathResults
	 *            the multi path results
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @param setup
	 *            Set to true to run the {@link #setup()} method
	 * @param store
	 *            the store (can be used to track results that pass validation)
	 * @return the results
	 */
	public void select(final IMultiPathFitResults multiPathResults, final int failures, boolean setup,
			SelectedResultStore store)
	{
		if (setup)
			setup();

		if (store == null)
			store = new SimpleSelectedResultStore(multiPathResults.getTotalCandidates());

		// TODO - this could be made iterative. Any pass through the data may store estimates 
		// using the SelectedResultStore and used to determine if 

		int failCount = 0;
		final int total = multiPathResults.getNumberOfResults();
		//while (multiPathResults.begin())
		//{
		for (int c = 0; c < total; c++)
		{
			final MultiPathFitResult multiPathResult = multiPathResults.getResult(c);

			if (multiPathResult == null)
				// Ignore this but do not count it as a failure
				continue;

			if (failCount <= failures || store.isValid(multiPathResult.candidateId))
			{
				// Assess the result if we are below the fail limit or have an estimate
				final SelectedResult result = select(multiPathResult, true, store);
				if (result != null)
				{
					final int[] ok = new int[result.results.length];
					int count = 0;
					for (int i = 0; i < ok.length; i++)
					{
						if (result.results[i].isNewResult())
							ok[count++] = i;
					}

					if (count != 0)
					{
						// This has valid results so add to the output subset only those that are new
						if (count == ok.length)
						{
							store.add(result);
						}
						else
						{
							final PreprocessedPeakResult[] filtered = new PreprocessedPeakResult[count];
							for (int i = 0; i < count; i++)
							{
								filtered[i] = result.results[ok[i]];
							}
							store.add(new SelectedResult(filtered, result.fitResult));
						}

						// More results were accepted so reset the fail count
						failCount = 0;
					}
					else
					{
						store.add(new SelectedResult(null, result.fitResult));

						// Nothing was accepted, increment fail count
						failCount++;
					}
				}
				else
				{
					// This failed. Just return the single result
					store.add(new SelectedResult(null, multiPathResult.getSingleFitResult()));

					// This was rejected, increment fail count
					failCount++;
				}
			}
		}
		//	multiPathResults.end();
		//}
	}

	/**
	 * Check all new and all existing results are valid. Returns the new results.
	 * <p>
	 * New results and validated candidates that fail the primary filter can be filtered using the minimal filter and
	 * sent to the store. The store can be used to determine if a fit for a different candidate has been performed
	 * already.
	 * 
	 * @param candidateId
	 *
	 * @param fitResult
	 *            the results
	 * @param validateCandidates
	 *            Set to true to validate the candidates
	 * @param store
	 *            the store
	 * @return The new results that pass the filter
	 */
	private PreprocessedPeakResult[] acceptAll(int candidateId, final FitResult fitResult, boolean validateCandidates,
			SelectedResultStore store)
	{
		if (fitResult == null || fitResult.results == null)
			return null;
		final PreprocessedPeakResult[] results = fitResult.results;

		// All new and existing results should be valid
		int count = 0;
		final int[] ok = new int[results.length];

		// Support for testing using the minimal filter.
		// Note: We do not check the store is not null. This is private method 
		// and we send in a null store if necessary.
		final boolean minimalFilter = minFilter != null;

		for (int i = 0; i < results.length; i++)
		{
			if (results[i].isNewResult())
			{
				if (results[i].getCandidateId() != candidateId)
				{
					// This is new result for a different candidate.
					// If a fit has already been accepted (or we don't know)
					// then it should be ignored.
					if (store.isFit(results[i].getCandidateId()))
						continue;
				}

				// All new results must pass
				if (accept(results[i]))
				{
					ok[count++] = i;
				}
				else
				{
					if (minimalFilter)
					{
						if (minAccept(results[i]))
							store.passMin(results[i]);
					}
					else
						return null;
				}
			}
			else if (results[i].isExistingResult())
			{
				// All existing results must pass
				if (!accept(results[i]))
				{
					return null;
				}
			}
			else if (validateCandidates)
			{
				// Optionally candidates must pass
				if (accept(results[i]))
				{
					ok[count++] = i;
				}
				else
				{
					if (minimalFilter)
					{
						if (minAccept(results[i]))
							store.passMin(results[i]);
					}
				}
			}
		}

		if (count == 0)
			return null;

		// Return the new results
		final PreprocessedPeakResult[] filtered = new PreprocessedPeakResult[count];
		for (int i = 0; i < count; i++)
		{
			filtered[i] = results[ok[i]];
			store.pass(filtered[i]);
		}

		return filtered;
	}

	/**
	 * Check any new and all existing results are valid. Returns the new results
	 * <p>
	 * New results and validated candidates that fail the primary filter can be filtered using the minimal filter and
	 * sent to the store. The store can be used to determine if a fit for a different candidate has been performed
	 * already.
	 * 
	 * @param candidateId
	 *
	 * @param fitResult
	 *            the results
	 * @param validateCandidates
	 *            Set to true to validate the candidates
	 * @param store
	 *            the store
	 * @return The new results that pass the filter
	 */
	private PreprocessedPeakResult[] acceptAny(int candidateId, final FitResult fitResult, boolean validateCandidates,
			SelectedResultStore store)
	{
		if (fitResult == null || fitResult.results == null)
			return null;
		final PreprocessedPeakResult[] results = fitResult.results;

		// Any new and all existing results should be valid
		int count = 0;
		final int[] ok = new int[results.length];

		// Support for testing using the minimal filter 
		// Note: We do not check the store is not null. This is private method 
		// and we send in a null store if necessary.
		final boolean minimalFilter = minFilter != null;

		for (int i = 0; i < results.length; i++)
		{
			if (results[i].isNewResult())
			{
				if (results[i].getCandidateId() != candidateId)
				{
					// This is new result for a different candidate.
					// If a fit has already been accepted (or we don't know)
					// then it should be ignored.
					if (store.isFit(results[i].getCandidateId()))
						continue;
				}

				// Any new result that pass are OK
				if (accept(results[i]))
				{
					ok[count++] = i;
				}
				else
				{
					if (minimalFilter)
					{
						if (minAccept(results[i]))
							store.passMin(results[i]);
					}
				}
			}
			else if (results[i].isExistingResult())
			{
				// All existing results must pass
				if (!accept(results[i]))
					return null;
			}
			else if (validateCandidates)
			{
				// Optionally candidates must pass
				if (accept(results[i]))
				{
					ok[count++] = i;
				}
				else
				{
					if (minimalFilter)
					{
						if (minAccept(results[i]))
							store.passMin(results[i]);
					}
				}
			}
		}

		if (count == 0)
			return null;

		// Return the new results
		final PreprocessedPeakResult[] filtered = new PreprocessedPeakResult[count];
		for (int i = 0; i < count; i++)
		{
			filtered[i] = results[ok[i]];
			store.pass(filtered[i]);
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
		final ArrayList<PreprocessedPeakResult> list = new ArrayList<PreprocessedPeakResult>(results.length);
		for (int i = 0; i < results.length; i++)
		{
			final PreprocessedPeakResult[] result = accept(results[i], false, null);
			if (result != null)
				list.addAll(Arrays.asList(result));
		}
		return list.toArray(new PreprocessedPeakResult[list.size()]);
	}

	/**
	 * Create a subset of multi-path results, i.e. all those that pass the filter.
	 * <p>
	 * The number of consecutive rejections are counted per frame. When the configured number of failures is reached all
	 * remaining results for the frame are rejected.
	 * <p>
	 * The number of failures before each result is stored in the failCount property of the MultiPathPeakResult. If the
	 * subset flag is set to true the current value of the failCount property is used to accumulate the fail count (as
	 * the current set is already a subset).
	 * 
	 * @param results
	 *            a set of results to analyse
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @param subset
	 *            True if a subset
	 * @return the filtered results
	 */
	public MultiPathFitResults[] filterSubset(final MultiPathFitResults[] results, final int failures, boolean subset)
	{
		final MultiPathFitResults[] newResults = new MultiPathFitResults[results.length];
		int size = 0;

		setup();
		for (int i = 0; i < results.length; i++)
		{
			final MultiPathFitResult[] newMultiPathResults = filter(results[i], failures, false, subset);
			if (newMultiPathResults != null)
				newResults[size++] = new MultiPathFitResults(results[i].frame, newMultiPathResults,
						results[i].totalCandidates);
		}

		return Arrays.copyOf(newResults, size);
	}

	/**
	 * Create a subset of multi-path results, i.e. all those that pass the filter.
	 * <p>
	 * The number of consecutive rejections are counted. When the configured number of failures is reached all
	 * remaining results for the frame are rejected.
	 * <p>
	 * The number of failures before each result is stored in the failCount property of the MultiPathPeakResult. If the
	 * subset flag is set to true the current value of the failCount property is used to accumulate the fail count (as
	 * the current set is already a subset).
	 * 
	 * @param results
	 *            a set of results to analyse
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @param subset
	 *            True if a subset
	 * @return the filtered results
	 */
	public MultiPathFitResult[] filter(final IMultiPathFitResults multiPathResults, final int failures, boolean subset)
	{
		return filter(multiPathResults, failures, true, subset);
	}

	/**
	 * Create a subset of multi-path results, i.e. all those that pass the filter.
	 * <p>
	 * The number of consecutive rejections are counted. When the configured number of failures is reached all
	 * remaining results for the frame are rejected.
	 * <p>
	 * The number of failures before each result is stored in the failCount property of the MultiPathPeakResult.
	 *
	 * @param multiPathResults
	 *            the multi path results
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @param setup
	 *            Set to true to run the {@link #setup()} method
	 * @param subset
	 *            True if a subset
	 * @return the filtered results
	 */
	private MultiPathFitResult[] filter(final IMultiPathFitResults multiPathResults, final int failures, boolean setup,
			boolean subset)
	{
		if (setup)
			setup();

		int failCount = 0;
		int size = 0;
		final MultiPathFitResult[] newMultiPathResults = new MultiPathFitResult[multiPathResults.getNumberOfResults()];
		final SimpleSelectedResultStore store = new SimpleSelectedResultStore(multiPathResults.getTotalCandidates());
		for (int c = 0; c < newMultiPathResults.length; c++)
		{
			final MultiPathFitResult multiPathResult = multiPathResults.getResult(c);

			// Include the number of failures before this result from the larger set
			if (subset)
				failCount += multiPathResult.failCount;

			if (failCount <= failures || store.isValid(multiPathResult.candidateId))
			{
				// Assess the result if we are below the fail limit or have an estimate
				final PreprocessedPeakResult[] result = accept(multiPathResult, true, store);
				if (result != null)
				{
					boolean isNew = false;
					for (int i = 0; i < result.length; i++)
					{
						if (result[i].isNewResult())
							isNew = true;
					}

					// Note: Even if the actual result failed, the candidate may have passed and so 
					// the entire multi-path result should be retained.

					// This has valid results so add to the output subset 
					newMultiPathResults[size++] = multiPathResult;
					// Store the number of failures before this result
					multiPathResult.failCount = failCount;

					if (isNew)
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
				else
				{
					// This was rejected, increment fail count
					failCount++;
				}
			}
		}

		if (size != 0)
			return Arrays.copyOf(newMultiPathResults, size);

		return null;
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
	public FractionClassificationResult fractionScore(final MultiPathFitResults[] results, final int failures,
			final int n)
	{
		return fractionScore(results, failures, n, false, null);
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
	public FractionClassificationResult fractionScoreSubset(final MultiPathFitResults[] results, final int failures,
			final int n)
	{
		return fractionScore(results, failures, n, true, null);
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
	 * @param assignments
	 *            the assignments
	 * @return the score
	 */
	public FractionClassificationResult fractionScore(final MultiPathFitResults[] results, final int failures,
			final int n, List<FractionalAssignment[]> assignments)
	{
		return fractionScore(results, failures, n, false, assignments);
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
	 * @param assignments
	 *            the assignments
	 * @return the score
	 */
	public FractionClassificationResult fractionScoreSubset(final MultiPathFitResults[] results, final int failures,
			final int n, List<FractionalAssignment[]> assignments)
	{
		return fractionScore(results, failures, n, true, assignments);
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
	 * @param assignments
	 *            the assignments
	 * @return the score
	 */
	private FractionClassificationResult fractionScore(final MultiPathFitResults[] results, final int failures,
			final int n, final boolean subset, List<FractionalAssignment[]> allAssignments)
	{
		final double[] score = new double[4];
		final ArrayList<FractionalAssignment> assignments = new ArrayList<FractionalAssignment>();

		final SimpleSelectedResultStore store = new SimpleSelectedResultStore();

		setup();
		for (MultiPathFitResults multiPathResults : results)
		{
			// Reset fail count for new frames
			int failCount = 0;
			int nPredicted = 0;
			store.resize(multiPathResults.totalCandidates);
			for (int c = 0; c < multiPathResults.multiPathFitResults.length; c++)
			{
				final MultiPathFitResult multiPathResult = multiPathResults.multiPathFitResults[c];

				// Include the number of failures before this result from the larger set
				if (subset)
					failCount += multiPathResult.failCount;

				if (failCount <= failures || store.isValid(multiPathResult.candidateId))
				{
					// Assess the result if we are below the fail limit or have an estimate
					final PreprocessedPeakResult[] result = accept(multiPathResult, true, store);
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

			final FractionalAssignment[] tmp = score(assignments, score, nPredicted);
			if (allAssignments != null)
				allAssignments.add(tmp);
		}

		// Note: We are using the integer positives and negatives fields to actually store integer TP and FP
		return new FractionClassificationResult(score[0], score[1], 0, n - score[0], (int) score[2], (int) score[3]);
	}

	/**
	 * Score the assignments (TP/FP) and then clear the list.
	 *
	 * @param assignments
	 *            The assignments
	 * @param score
	 *            Scores array to accumulate TP/FP scores
	 * @param nPredicted
	 *            The number of predictions
	 * @return the fractional assignments
	 */
	private FractionalAssignment[] score(final ArrayList<FractionalAssignment> assignments, final double[] score,
			final int nPredicted)
	{
		if (assignments.isEmpty())
			return null;
		final FractionalAssignment[] tmp = new FractionalAssignment[assignments.size()];
		final RankedScoreCalculator calc = new RankedScoreCalculator(assignments.toArray(tmp));
		final double[] result = calc.score(nPredicted, false);
		score[0] += result[0];
		score[1] += result[1];
		score[2] += result[2];
		score[3] += result[3];
		assignments.clear();
		return tmp;
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
			return size() == 0;
		}

		public boolean contains(Object o)
		{
			throw new NotImplementedException();
		}

		public Iterator<FractionalAssignment> iterator()
		{
			throw new NotImplementedException();
		}

		public Object[] toArray()
		{
			return a;
		}

		public <T> T[] toArray(T[] a)
		{
			return a;
		}

		public boolean add(FractionalAssignment e)
		{
			throw new NotImplementedException();
		}

		public boolean remove(Object o)
		{
			throw new NotImplementedException();
		}

		public boolean containsAll(Collection<?> c)
		{
			throw new NotImplementedException();
		}

		public boolean addAll(Collection<? extends FractionalAssignment> c)
		{
			throw new NotImplementedException();
		}

		public boolean removeAll(Collection<?> c)
		{
			throw new NotImplementedException();
		}

		public boolean retainAll(Collection<?> c)
		{
			throw new NotImplementedException();
		}

		public void clear()
		{
			throw new NotImplementedException();
		}
	}

	/**
	 * @return An XML representation of this object
	 */
	public String toXML()
	{
		return XStreamWrapper.toXML(this);
	}

	/**
	 * Create the filter from the XML representation
	 * 
	 * @param xml
	 * @return the filter
	 */
	public static MultiPathFilter fromXML(String xml)
	{
		try
		{
			return (MultiPathFilter) XStreamWrapper.fromXML(xml);
		}
		catch (ClassCastException ex)
		{
			//ex.printStackTrace();
		}
		return null;
	}
}