package gdsc.smlm.results.filter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

import gdsc.core.match.FractionClassificationResult;
import gdsc.core.match.FractionalAssignment;
import gdsc.core.match.RankedScoreCalculator;
import gdsc.core.utils.NotImplementedException;
import gdsc.smlm.results.filter.MultiPathFitResult.FitResult;

// TODO: Auto-generated Javadoc
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
public class MultiPathFilter implements Cloneable
{
	/**
	 * Stores the results that were accepted when filtering a multi-path result. Also stores the fit result that was
	 * used to select the results.
	 */
	public class SelectedResult
	{

		/** The results. */
		final public PreprocessedPeakResult[] results;

		/** The fit result. */
		final public FitResult fitResult;

		/**
		 * Instantiates a new selected result.
		 *
		 * @param results
		 *            the results
		 * @param fitResult
		 *            the fit result
		 */
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
		 * Add a selected result to the store.
		 *
		 * @param selectedResult
		 *            the selected result
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
		 * A result that passed the primary filter.
		 *
		 * @param result
		 *            the result
		 */
		void pass(PreprocessedPeakResult result);

		/**
		 * A result that passed the minimal filter.
		 *
		 * @param result
		 *            the result
		 */
		void passMin(PreprocessedPeakResult result);
	}

	/**
	 * Allow tracking of candidates that have been fit.
	 */
	private class SimpleSelectedResultStore implements SelectedResultStore
	{
		/** The is fit. */
		boolean[] isFit;

		/** The is valid. */
		boolean[] isValid;

		/**
		 * Instantiates a new simple selected result store.
		 */
		SimpleSelectedResultStore()
		{
			isFit = new boolean[0];
			isValid = new boolean[0];
		}

		/**
		 * Instantiates a new simple selected result store.
		 *
		 * @param totalCandidates
		 *            the total candidates
		 */
		SimpleSelectedResultStore(int totalCandidates)
		{
			isFit = new boolean[totalCandidates];
			isValid = new boolean[totalCandidates];
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see
		 * gdsc.smlm.results.filter.MultiPathFilter.SelectedResultStore#add(gdsc.smlm.results.filter.MultiPathFilter.
		 * SelectedResult)
		 */
		public void add(SelectedResult selectedResult)
		{
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see gdsc.smlm.results.filter.MultiPathFilter.SelectedResultStore#isFit(int)
		 */
		public boolean isFit(int candidateId)
		{
			return isFit[candidateId];
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see gdsc.smlm.results.filter.MultiPathFilter.SelectedResultStore#isValid(int)
		 */
		public boolean isValid(int candidateId)
		{
			return isValid[candidateId];
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see gdsc.smlm.results.filter.MultiPathFilter.SelectedResultStore#pass(gdsc.smlm.results.filter.
		 * PreprocessedPeakResult)
		 */
		public void pass(PreprocessedPeakResult result)
		{
			// This an existing result or candidate. Mark as valid so candidates will be processed
			isValid[result.getCandidateId()] = true;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see gdsc.smlm.results.filter.MultiPathFilter.SelectedResultStore#passMin(gdsc.smlm.results.filter.
		 * PreprocessedPeakResult)
		 */
		public void passMin(PreprocessedPeakResult result)
		{
			// Passing the minimal filter does not mean it is valid. This would be used to store
			// a fit estimate during processing for this candidate.
		}

		/**
		 * Resize.
		 *
		 * @param totalCandidates
		 *            the total candidates
		 */
		public void resize(int totalCandidates)
		{
			if (isFit.length < totalCandidates)
			{
				isFit = new boolean[totalCandidates];
				isValid = new boolean[totalCandidates];
			}
			else
			{
				for (int i = 0; i < totalCandidates; i++)
				{
					isFit[i] = false;
					isValid[i] = false;
				}
			}
		}
	}

	/**
	 * Used to return default behaviour for acceptAny/acceptAll.
	 */
	private static class NullSelectedResultStore implements SelectedResultStore
	{

		/*
		 * (non-Javadoc)
		 * 
		 * @see
		 * gdsc.smlm.results.filter.MultiPathFilter.SelectedResultStore#add(gdsc.smlm.results.filter.MultiPathFilter.
		 * SelectedResult)
		 */
		public void add(SelectedResult selectedResult)
		{

		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see gdsc.smlm.results.filter.MultiPathFilter.SelectedResultStore#isFit(int)
		 */
		public boolean isFit(int candidateId)
		{
			// Make sure non-candidate fits are ignored.
			return true;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see gdsc.smlm.results.filter.MultiPathFilter.SelectedResultStore#isValid(int)
		 */
		public boolean isValid(int candidateId)
		{
			return false;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see gdsc.smlm.results.filter.MultiPathFilter.SelectedResultStore#pass(gdsc.smlm.results.filter.
		 * PreprocessedPeakResult)
		 */
		public void pass(PreprocessedPeakResult result)
		{

		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see gdsc.smlm.results.filter.MultiPathFilter.SelectedResultStore#passMin(gdsc.smlm.results.filter.
		 * PreprocessedPeakResult)
		 */
		public void passMin(PreprocessedPeakResult result)
		{

		}
	}

	/** The null selected result store. */
	private static NullSelectedResultStore nullSelectedResultStore = new NullSelectedResultStore();

	/**
	 * Allows signalling of results that have been selected during multi-path filter scoring.
	 */
	public interface FractionScoreStore
	{

		/**
		 * Add the unique Id of a result that was selected
		 *
		 * @param uniqueId
		 *            the unique id
		 */
		void add(int uniqueId);
	}

	/**
	 * Used to return default behaviour
	 */
	private static class NullFractionScoreStore implements FractionScoreStore
	{
		public void add(int uniqueId)
		{
		}
	}

	/** The null fraction result store. */
	private static NullFractionScoreStore nullFractionScoreStore = new NullFractionScoreStore();

	/** The direct filter to apply to the results. */
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
	 * Return a deep copy of this object with a copy of the configured filters.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	public MultiPathFilter clone()
	{
		return new MultiPathFilter(copy(filter), copy(minFilter), residualsThreshold);
	}

	private IDirectFilter copy(IDirectFilter f)
	{
		return (f == null) ? null : f.copy();
	}

	/**
	 * Gets the filter.
	 *
	 * @return the filter
	 */
	public IDirectFilter getFilter()
	{
		return copy(filter);
	}

	/**
	 * Gets the minimal filter.
	 *
	 * @return the minimal filter
	 */
	public IDirectFilter getMinimalFilter()
	{
		return copy(minFilter);
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
		final int candidateId = multiPathResult.candidateId;

		// Ensure we don't have to check the store in acceptAll/acceptAny
		if (store == null)
			store = nullSelectedResultStore;

		// The aim is to obtain a new result for the current candidate Id. 
		// acceptAll/acceptAny will return all new results, even if they do not match the candidate.
		// So we check the candidate Id and return when we have a new result for the candidate. 
		// We accept the doublet fit over the single fit if we are performing doublet fitting.
		// If nothing matches then pick the result with the most new results, or use the default 
		// order we processed the fits. 

		boolean doDoublet = false;

		// Filter multi-fit

		// Accept all and then check if we can perform a doublet fit
		//		final PreprocessedPeakResult[] multiResults = acceptAll(candidateId, multiPathResult.getMultiFitResult(),
		//				validateCandidates, store);
		//		if (multiResults == null)
		//		{
		//			// The fit was not accepted. However it may have been rejected for being too wide
		//			// and is suitable for a doublet fit.
		//			doDoublet = isSuitableForDoubletFit(multiPathResult, multiPathResult.getMultiFitResult(), false);
		//		}
		//		else
		//		{
		//			doDoublet = (residualsThreshold < 1 && multiPathResult.getMultiQAScore() > residualsThreshold);
		//		}

		// Accept any and then check if we can perform a doublet fit
		final PreprocessedPeakResult[] multiResults = acceptAny(candidateId, multiPathResult.getMultiFitResult(),
				validateCandidates, store);
		doDoublet = isSuitableForDoubletFit(multiPathResult, multiPathResult.getMultiFitResult(), false);

		final PreprocessedPeakResult[] multiDoubletResults;
		if (doDoublet)
		{
			multiDoubletResults = acceptAny(candidateId, multiPathResult.getMultiDoubletFitResult(), validateCandidates,
					store);
			if (multiDoubletResults != null)
			{
				// Check we have a new result for the candidate
				if (contains(multiDoubletResults, candidateId))
					return multiDoubletResults;
			}
		}
		else
		{
			multiDoubletResults = null;
		}

		// Check if the multi result is to the correct candidate.
		if (multiResults != null && contains(multiResults, candidateId))
			return multiResults;

		// We reached here with:
		// a multi fit that failed or matched a different candidate
		// a doublet multi fit that failed or matched a different candidate

		doDoublet = false;

		// Filter single-fit
		final PreprocessedPeakResult[] singleResults = acceptAll(candidateId, multiPathResult.getSingleFitResult(),
				validateCandidates, store);
		if (singleResults == null)
		{
			// The fit was not accepted. However it may have been rejected for being too wide
			// and is suitable for a doublet fit.
			doDoublet = isSuitableForDoubletFit(multiPathResult, multiPathResult.getSingleFitResult(), true);
		}
		else
		{
			// The single fit is OK.
			doDoublet = (residualsThreshold < 1 && multiPathResult.getSingleQAScore() > residualsThreshold);
		}

		// We reached here with:
		// a multi fit that failed or matched a different candidate
		// a doublet multi fit that failed or matched a different candidate
		// a single fit that is eligible for doublet fitting, it may be null (if it passed without width filtering)

		final PreprocessedPeakResult[] singleDoubletResults;
		if (doDoublet)
		{
			singleDoubletResults = acceptAny(candidateId, multiPathResult.getDoubletFitResult(), validateCandidates,
					store);
			if (singleDoubletResults != null)
			{
				// Check we have a new result for the candidate
				if (contains(singleDoubletResults, candidateId))
					return singleDoubletResults;
			}
		}
		else
		{
			singleDoubletResults = null;
		}

		// Check if the single result is to the correct candidate.
		if (singleResults != null && contains(singleResults, candidateId))
			return singleResults;

		// We reached here with:
		// a multi fit that failed or matched a different candidate
		// a multi doublet fit that failed or matched a different candidate
		// a single fit that failed or matched a different candidate
		// a doublet fit that failed or matched a different candidate
		return rank(multiResults, multiDoubletResults, singleResults, singleDoubletResults);
	}

	/**
	 * Allows results to be ranked
	 */
	private class ResultRank implements Comparable<ResultRank>
	{
		/** The results. */
		final PreprocessedPeakResult[] results;

		/** The default rank (when the count of new result is the same). */
		final int rank;

		/** The count of new results. */
		final int count;

		/**
		 * Instantiates a new result rank.
		 *
		 * @param results
		 *            the results
		 * @param rank
		 *            the rank
		 */
		public ResultRank(PreprocessedPeakResult[] results, int rank)
		{
			this.results = results;
			this.rank = rank;
			if (results == null)
			{
				// Negative so null results are ranked below not-null results with no new results
				count = -1;
			}
			else
			{
				count = countNewResult(results);
			}
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		public int compareTo(ResultRank o)
		{
			final int result = o.count - count;
			if (result != 0)
				return result;
			return rank - o.rank;
		}
	}

	/**
	 * Rank the results. It is assumed that each result is either null or has results that do not match the current
	 * candidate Id. In this case we will return the result that has the highest number of new results. In the event of
	 * a tie we order as multi, multi-doublet, single, doublet. Results are discounted if null.
	 *
	 * @param multiResults
	 *            the multi results
	 * @param multiDoubletResults
	 *            the multi doublet results
	 * @param singleResults
	 *            the single results
	 * @param singleDoubletResults
	 *            the doublet results
	 * @return the preprocessed peak result[]
	 */
	private PreprocessedPeakResult[] rank(PreprocessedPeakResult[] multiResults,
			PreprocessedPeakResult[] multiDoubletResults, PreprocessedPeakResult[] singleResults,
			PreprocessedPeakResult[] singleDoubletResults)
	{
		if (multiResults == null && multiDoubletResults == null && singleResults == null &&
				singleDoubletResults == null)
			return null;
		final ResultRank[] rank = new ResultRank[4];
		rank[0] = new ResultRank(multiResults, 1);
		rank[1] = new ResultRank(multiDoubletResults, 2);
		rank[2] = new ResultRank(singleResults, 3);
		rank[3] = new ResultRank(singleDoubletResults, 4);
		Arrays.sort(rank);
		return rank[0].results;
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
		final int candidateId = multiPathResult.candidateId;

		// Ensure we don't have to check the store in acceptAll/acceptAny
		if (store == null)
			store = nullSelectedResultStore;

		// The aim is to obtain a new result for the current candidate Id. 
		// acceptAll/acceptAny will return all new results, even if they do not match the candidate.
		// So we check the candidate Id and return when we have a new result for the candidate. 
		// We accept the doublet fit over the single fit if we are performing doublet fitting.
		// If nothing matches then pick the result with the most new results, or use the default 
		// order we processed the fits. 

		boolean doDoublet = false;

		// Filter multi-fit
		// Accept all and then check if we can perform a doublet fit
		//		final PreprocessedPeakResult[] multiResults = acceptAll(candidateId, multiPathResult.getMultiFitResult(),
		//				validateCandidates, store);
		//		if (multiResults == null)
		//		{
		//			// The fit was not accepted. However it may have been rejected for being too wide
		//			// and is suitable for a doublet fit.
		//			doDoublet = isSuitableForDoubletFit(multiPathResult, multiPathResult.getMultiFitResult(), false);
		//		}
		//		else
		//		{
		//			doDoublet = (residualsThreshold < 1 && multiPathResult.getMultiQAScore() > residualsThreshold);
		//		}

		// Accept any and then check if we can perform a doublet fit
		final PreprocessedPeakResult[] multiResults = acceptAny(candidateId, multiPathResult.getMultiFitResult(),
				validateCandidates, store);
		doDoublet = isSuitableForDoubletFit(multiPathResult, multiPathResult.getMultiFitResult(), false);

		final PreprocessedPeakResult[] multiDoubletResults;
		if (doDoublet)
		{
			multiDoubletResults = acceptAny(candidateId, multiPathResult.getMultiDoubletFitResult(), validateCandidates,
					store);
			if (multiDoubletResults != null)
			{
				// Check we have a new result for the candidate
				if (contains(multiDoubletResults, candidateId))
					return new SelectedResult(multiDoubletResults, multiPathResult.getMultiDoubletFitResult());
			}
		}
		else
		{
			multiDoubletResults = null;
		}

		// Check if the multi result is to the correct candidate.
		if (multiResults != null && contains(multiResults, candidateId))
			return new SelectedResult(multiResults, multiPathResult.getMultiFitResult());

		// We reached here with:
		// a multi fit that failed or matched a different candidate
		// a doublet multi fit that failed or matched a different candidate

		doDoublet = false;

		// Filter single-fit
		final PreprocessedPeakResult[] singleResults = acceptAll(candidateId, multiPathResult.getSingleFitResult(),
				validateCandidates, store);
		if (singleResults == null)
		{
			// The fit was not accepted. However it may have been rejected for being too wide
			// and is suitable for a doublet fit.
			doDoublet = isSuitableForDoubletFit(multiPathResult, multiPathResult.getSingleFitResult(), true);
		}
		else
		{
			// The single fit is OK.
			doDoublet = (residualsThreshold < 1 && multiPathResult.getSingleQAScore() > residualsThreshold);
		}

		// We reached here with:
		// a multi fit that failed or matched a different candidate
		// a doublet multi fit that failed or matched a different candidate
		// a single fit that is eligible for doublet fitting, it may be null (if it passed without width filtering)

		final PreprocessedPeakResult[] singleDoubletResults;
		if (doDoublet)
		{
			singleDoubletResults = acceptAny(candidateId, multiPathResult.getDoubletFitResult(), validateCandidates,
					store);
			if (singleDoubletResults != null)
			{
				// Check we have a new result for the candidate
				if (contains(singleDoubletResults, candidateId))
					return new SelectedResult(singleDoubletResults, multiPathResult.getDoubletFitResult());
			}
		}
		else
		{
			singleDoubletResults = null;
		}

		// Check if the single result is to the correct candidate.
		if (singleResults != null && contains(singleResults, candidateId))
			return new SelectedResult(singleResults, multiPathResult.getSingleFitResult());

		// We reached here with:
		// a multi fit that failed or matched a different candidate
		// a multi doublet fit that failed or matched a different candidate
		// a single fit that failed or matched a different candidate
		// a doublet fit that failed or matched a different candidate
		final PreprocessedPeakResult[] result = rank(multiResults, multiDoubletResults, singleResults,
				singleDoubletResults);
		if (result == null)
			return null;
		//@formatter:off
		if (result == multiResults)	       return new SelectedResult(multiResults,        multiPathResult.getMultiFitResult());
		if (result == multiDoubletResults) return new SelectedResult(multiDoubletResults, multiPathResult.getMultiDoubletFitResult());
		if (result == singleResults)	   return new SelectedResult(singleResults,       multiPathResult.getSingleFitResult());
		                                   return new SelectedResult(singleDoubletResults,multiPathResult.getDoubletFitResult());
		//@formatter:on
	}

	private boolean isSuitableForDoubletFit(MultiPathFitResult multiPathResult, FitResult fitResult, boolean singleQA)
	{
		// Check there is a fit result
		if (fitResult == null || fitResult.status != 0 || fitResult.results == null)
			return false;

		// Check if the residuals score is below the configured threshold
		if (residualsThreshold >= 1)
			return false;

		// Check the other results are OK. Candidates are allowed to fail. New and existing results must pass.
		for (int i = 1; i < validationResults.length; i++)
			if ((fitResult.results[i].isNewResult() || fitResult.results[i].isExistingResult()) &&
					validationResults[i] != 0)
				return false;

		if (validationResults[0] == 0)
		{
			// The peak was valid so check the residuals
			return ((singleQA) ? multiPathResult.getSingleQAScore()
					: multiPathResult.getMultiQAScore()) > residualsThreshold;
		}

		// Check if it failed due to width
		if (!DirectFilter.anySet(validationResults[0], DirectFilter.V_X_SD_FACTOR | DirectFilter.V_X_SD_FACTOR))
			return false;

		// Get the first spot
		final PreprocessedPeakResult firstResult = fitResult.results[0];

		// Check the width is reasonable given the size of the fitted region.
		//@formatter:off
		if (	firstResult.getXSDFactor() < 1 || // Not a wide spot
				firstResult.getXSD() > multiPathResult.width || // width covers more than the region
				firstResult.getYSDFactor() < 1 || // Not a wide spot
				firstResult.getYSD() > multiPathResult.height // width covers more than the region
			)
			return false;
		//@formatter:on

		// Check the quadrant analysis on the fit residuals
		if (((singleQA) ? multiPathResult.getSingleQAScore() : multiPathResult.getMultiQAScore()) < residualsThreshold)
			return false;

		// We must validate the spot without width filtering. Do not change the min filter.
		filter.setup(DirectFilter.NO_WIDTH);

		try
		{
			if (!filter.accept(firstResult))
				// This is still a bad single result, without width filtering
				return false;
		}
		finally
		{
			// reset
			filter.setup();
		}

		return true;
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

			final boolean evaluateFit = failCount <= failures;
			if (evaluateFit || store.isValid(multiPathResult.candidateId))
			{
				// Assess the result if we are below the fail limit or have an estimate
				final SelectedResult result = select(multiPathResult, true, store);
				int size = 0;
				if (result != null)
				{
					final int[] ok = new int[result.results.length];
					for (int i = 0; i < ok.length; i++)
					{
						if (result.results[i].isNewResult())
							ok[size++] = i;
					}

					if (size != 0)
					{
						// This has valid results so add to the output subset only those that are new
						if (size == ok.length)
						{
							store.add(result);
						}
						else
						{
							final PreprocessedPeakResult[] filtered = new PreprocessedPeakResult[size];
							for (int i = 0; i < size; i++)
							{
								filtered[i] = result.results[ok[i]];
							}
							store.add(new SelectedResult(filtered, result.fitResult));
						}
					}
				}
				if (size == 0)
				{
					// This failed. Just return the single result
					store.add(new SelectedResult(null, multiPathResult.getSingleFitResult()));
				}
				if (evaluateFit && size != 0)
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
				// This failed. Just return the single result
				store.add(new SelectedResult(null, multiPathResult.getSingleFitResult()));

				// This was rejected, increment fail count
				failCount++;
			}
		}
		//	multiPathResults.end();
		//}
	}

	@XStreamOmitField
	private int[] validationResults;
	@XStreamOmitField
	private boolean failExisting;
	@XStreamOmitField
	private boolean failNew;

	/**
	 * Check all new and all existing results are valid. Returns the new results.
	 * <p>
	 * New results and validated candidates that fail the primary filter can be filtered using the minimal filter and
	 * sent to the store. The store can be used to determine if a fit for a different candidate has been performed
	 * already.
	 *
	 * @param candidateId
	 *            the candidate id
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

		// Validate the results
		validationResults = new int[results.length];
		for (int i = 0; i < results.length; i++)
		{
			validationResults[i] = filter.validate(results[i]);
		}

		// All new and existing results should be valid
		int count = 0;
		final int[] ok = new int[results.length];

		// Support for testing using the minimal filter.
		// Note: We do not check the store is not null. This is private method 
		// and we send in a null store if necessary.
		final boolean minimalFilter = minFilter != null;

		failExisting = false;
		failNew = false;

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

				if (validationResults[i] == 0)
				{
					ok[count++] = i;
				}
				else
				{
					failNew = true;
					if (minimalFilter)
					{
						if (minAccept(results[i]))
							store.passMin(results[i]);
					}
				}
			}
			else if (results[i].isExistingResult())
			{
				if (validationResults[i] != 0)
					failExisting = true;
			}
			else if (validateCandidates)
			{
				// Optionally candidates must pass
				if (validationResults[i] == 0)
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

		// All new results must pass
		// All existing results must pass
		if (count == 0 || failNew || failExisting)
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
	 *            the candidate id
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

		// Validate the results
		validationResults = new int[results.length];
		for (int i = 0; i < results.length; i++)
		{
			validationResults[i] = filter.validate(results[i]);
		}

		// Any new and all existing results should be valid
		int count = 0;
		final int[] ok = new int[results.length];

		// Support for testing using the minimal filter 
		// Note: We do not check the store is not null. This is private method 
		// and we send in a null store if necessary.
		final boolean minimalFilter = minFilter != null;

		failExisting = false;
		failNew = false;

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
				if (validationResults[i] == 0)
				{
					ok[count++] = i;
				}
				else
				{
					failNew = true;
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
				if (validationResults[i] != 0)
					failExisting = true;
			}
			else if (validateCandidates)
			{
				// Optionally candidates must pass
				if (validationResults[i] == 0)
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

		// All existing results must pass
		if (count == 0 || failExisting)
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
	 * Checks if there is the given candidate in the results.
	 *
	 * @param results
	 *            the results
	 * @param candidateId
	 *            the candidate id
	 * @return true, if there is the given candidate
	 */
	private boolean contains(final PreprocessedPeakResult[] results, final int candidateId)
	{
		for (int i = 0; i < results.length; i++)
			if (results[i].getCandidateId() == candidateId)
				return true;
		return false;
	}

	/**
	 * Counts the number of new results in the results.
	 *
	 * @param results
	 *            the results
	 * @return The count
	 */
	private int countNewResult(final PreprocessedPeakResult[] results)
	{
		int c = 0;
		if (results != null)
		{
			for (int i = 0; i < results.length; i++)
				if (results[i].isNewResult())
					c++;
		}
		return c;
	}

	/**
	 * Check if the results contain a new result
	 *
	 * @param results
	 *            the results
	 * @return True if a new result
	 */
	private boolean isNewResult(final PreprocessedPeakResult[] results)
	{
		if (results != null)
		{
			for (int i = 0; i < results.length; i++)
				if (results[i].isNewResult())
					return true;
		}
		return false;
	}

	/**
	 * Filter a set of multi-path results into a set of results.
	 *
	 * @param results
	 *            the results
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @param subset
	 *            True if a subset (the candidate Id will be used to determine the number of failed fits before the
	 *            current candidate)
	 * @return the filtered results
	 */
	final public PreprocessedPeakResult[] filter(final MultiPathFitResults[] results, final int failures,
			boolean subset)
	{
		setup();
		final SimpleSelectedResultStore store = new SimpleSelectedResultStore();

		final ArrayList<PreprocessedPeakResult> list = new ArrayList<PreprocessedPeakResult>(results.length);
		for (int i = 0; i < results.length; i++)
		{
			final MultiPathFitResults multiPathResults = results[i];

			int failCount = 0;
			int lastId = -1;
			int size = multiPathResults.getNumberOfResults();
			store.resize(multiPathResults.getTotalCandidates());

			for (int c = 0; c < size; c++)
			{
				final MultiPathFitResult multiPathResult = multiPathResults.getResult(c);

				// Include the number of failures before this result from the larger set
				if (subset)
				{
					failCount += (multiPathResult.candidateId - (lastId + 1));
					lastId = multiPathResult.candidateId;
				}

				final boolean evaluateFit = failCount <= failures;
				if (evaluateFit || store.isValid(multiPathResult.candidateId))
				{
					// Evaluate the result. 
					// This allows storing more estimates in the store even if we are past the failures limit.
					final PreprocessedPeakResult[] result = accept(multiPathResult, false, store);

					if (result != null)
						list.addAll(Arrays.asList(result));

					if (evaluateFit && isNewResult(result))
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
		return list.toArray(new PreprocessedPeakResult[list.size()]);
	}

	/**
	 * Create a subset of multi-path results, i.e. all those that pass the filter.
	 * <p>
	 * The number of consecutive rejections are counted per frame. When the configured number of failures is reached all
	 * remaining results for the frame are rejected.
	 * <p>
	 * If the subset flag is set to true the candidate Id will be used to determine the number of failed fits before the
	 * current candidate, assuming candidates start at zero and increment.
	 * 
	 * @param results
	 *            a set of results to analyse
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @param subset
	 *            True if a subset (the candidate Id will be used to determine the number of failed fits before the
	 *            current candidate)
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
						results[i].totalCandidates, results[i].nActual);
		}

		return Arrays.copyOf(newResults, size);
	}

	/**
	 * Create a subset of multi-path results, i.e. all those that pass the filter.
	 * <p>
	 * The number of consecutive rejections are counted. When the configured number of failures is reached all
	 * remaining results for the frame are rejected.
	 * <p>
	 * If the subset flag is set to true the candidate Id will be used to determine the number of failed fits before the
	 * current candidate, assuming candidates start at zero and increment.
	 *
	 * @param multiPathResults
	 *            the multi path results
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @param subset
	 *            True if a subset (the candidate Id will be used to determine the number of failed fits before the
	 *            current candidate)
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
	 * If the subset flag is set to true the candidate Id will be used to determine the number of failed fits before the
	 * current candidate, assuming candidates start at zero and increment.
	 *
	 * @param multiPathResults
	 *            the multi path results
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @param setup
	 *            Set to true to run the {@link #setup()} method
	 * @param subset
	 *            True if a subset (the candidate Id will be used to determine the number of failed fits before the
	 *            current candidate)
	 * @return the filtered results
	 */
	private MultiPathFitResult[] filter(final IMultiPathFitResults multiPathResults, final int failures, boolean setup,
			boolean subset)
	{
		if (setup)
			setup();

		int failCount = 0;
		int lastId = -1;
		int size = 0;
		final MultiPathFitResult[] newMultiPathResults = new MultiPathFitResult[multiPathResults.getNumberOfResults()];
		final SimpleSelectedResultStore store = new SimpleSelectedResultStore(multiPathResults.getTotalCandidates());
		//		if (multiPathResults.getFrame() == 12)
		//			System.out.println("Debug");
		for (int c = 0; c < newMultiPathResults.length; c++)
		{
			final MultiPathFitResult multiPathResult = multiPathResults.getResult(c);

			// Include the number of failures before this result from the larger set
			if (subset)
			{
				failCount += (multiPathResult.candidateId - (lastId + 1));
				lastId = multiPathResult.candidateId;
			}

			final boolean evaluateFit = failCount <= failures;
			if (evaluateFit || store.isValid(multiPathResult.candidateId))
			{
				// Evaluate the result. 
				// This allows storing more estimates in the store even if we are past the failures limit.
				final PreprocessedPeakResult[] result = accept(multiPathResult, false, store);

				// Note: Even if the actual result failed, the candidate may have passed and so 
				// the entire multi-path result should be retained.

				// Also note that depending on the filter, different results can be selected and pushed through
				// the store to set them valid. So we must push everything through the store to ensure nothing 
				// is removed that could be used.
				checkIsValid(multiPathResult.getSingleFitResult(), store);
				checkIsValid(multiPathResult.getMultiFitResult(), store);
				checkIsValid(multiPathResult.getDoubletFitResult(), store);
				checkIsValid(multiPathResult.getMultiDoubletFitResult(), store);

				// This has valid results so add to the output subset 
				newMultiPathResults[size++] = multiPathResult;

				if (evaluateFit && isNewResult(result))
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

		if (size != 0)
			return Arrays.copyOf(newMultiPathResults, size);

		return null;
	}

	private void checkIsValid(FitResult fitResult, SimpleSelectedResultStore store)
	{
		if (fitResult == null || fitResult.results == null)
			return;
		final PreprocessedPeakResult[] results = fitResult.results;

		for (int i = 0; i < results.length; i++)
		{
			if (!store.isValid[results[i].getCandidateId()])
			{
				if (accept(results[i]))
				{
					store.isValid[results[i].getCandidateId()] = true;
				}
			}
		}
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
		return fractionScore(results, failures, n, false, null, null);
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
		return fractionScore(results, failures, n, true, null, null);
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
	 * @param scoreStore
	 *            the score store
	 * @return the score
	 */
	public FractionClassificationResult fractionScore(final MultiPathFitResults[] results, final int failures,
			final int n, List<FractionalAssignment[]> assignments, FractionScoreStore scoreStore)
	{
		return fractionScore(results, failures, n, false, assignments, scoreStore);
	}

	/**
	 * Score a subset of multi-path results. The subset can be created with
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
	 * @param scoreStore
	 *            the score store
	 * @return the score
	 */
	public FractionClassificationResult fractionScoreSubset(final MultiPathFitResults[] results, final int failures,
			final int n, List<FractionalAssignment[]> assignments, FractionScoreStore scoreStore)
	{
		return fractionScore(results, failures, n, true, assignments, scoreStore);
	}

	String debugFilename;

	/**
	 * Sets the debug file for scoring.
	 *
	 * @param filename
	 *            the new debug file
	 */
	public void setDebugFile(String filename)
	{
		debugFilename = filename;
	}

	/**
	 * Score a set of multi-path results.
	 * <p>
	 * If the subset flag is set to true the candidate Id will be used to determine the number of failed fits before the
	 * current candidate, assuming candidates start at zero and increment.
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
	 *            True if a subset (the candidate Id will be used to determine the number of failed fits before the
	 *            current candidate)
	 * @param allAssignments
	 *            the assignments
	 * @return the score
	 */
	private FractionClassificationResult fractionScore(final MultiPathFitResults[] results, final int failures,
			final int n, final boolean subset, List<FractionalAssignment[]> allAssignments,
			FractionScoreStore scoreStore)
	{
		final double[] score = new double[4];
		final ArrayList<FractionalAssignment> assignments = new ArrayList<FractionalAssignment>();

		final SimpleSelectedResultStore store = new SimpleSelectedResultStore();
		if (scoreStore == null)
			scoreStore = nullFractionScoreStore;
		final boolean save = allAssignments != null;

		//		// Debugging the results that are scored
		//		java.io.OutputStreamWriter out = null;
		//		if (debugFilename != null)
		//		{
		//			try
		//			{
		//				out = new java.io.OutputStreamWriter(new java.io.FileOutputStream(debugFilename), "UTF-8");
		//			}
		//			catch (Exception e)
		//			{
		//			}
		//		}

		setup();
		for (int k = 0; k < results.length; k++)
		{
			final MultiPathFitResults multiPathResults = results[k];

			// Reset fail count for new frames
			int failCount = 0;
			int lastId = -1;
			int nPredicted = 0;
			store.resize(multiPathResults.totalCandidates);
			for (int c = 0; c < multiPathResults.multiPathFitResults.length; c++)
			{
				final MultiPathFitResult multiPathResult = multiPathResults.multiPathFitResults[c];

				// Include the number of failures before this result from the larger set
				if (subset)
				{
					failCount += (multiPathResult.candidateId - (lastId + 1));
					lastId = multiPathResult.candidateId;
				}

				final boolean evaluateFit = failCount <= failures;
				if (evaluateFit || store.isValid(multiPathResult.candidateId))
				{
					//					if (out != null)
					//					{
					//						try
					//						{
					//							out.write(String.format("[%d] %d : %d %b %b\n", multiPathResults.frame,
					//									multiPathResult.candidateId, failCount, store.isValid(multiPathResult.candidateId),
					//									isNewResult(accept(multiPathResult, true, null))));
					//						}
					//						catch (Exception e)
					//						{
					//							try
					//							{
					//								out.close();
					//							}
					//							catch (Exception ee)
					//							{
					//							}
					//							finally
					//							{
					//								out = null;
					//							}
					//						}
					//					}

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
								// This is a new fitted result
								scoreStore.add(result[i].getUniqueId());
								store.isFit[result[i].getCandidateId()] = true;

								final FractionalAssignment[] a = result[i].getAssignments(nPredicted++);
								if (a != null && a.length > 0)
								{
									//list.addAll(Arrays.asList(a));
									assignments.addAll(new DummyCollection(a));
								}
							}
						}
					}
					if (evaluateFit && size != nPredicted)
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

			final FractionalAssignment[] tmp = score(assignments, score, nPredicted, save, multiPathResults.nActual);
			if (save)
				allAssignments.add(tmp);

			//			if (out != null)
			//			{
			//				try
			//				{
			//					out.write(String.format("[%d] %s\n", multiPathResults.frame, Arrays.toString(score)));
			//				}
			//				catch (Exception e)
			//				{
			//					try
			//					{
			//						out.close();
			//					}
			//					catch (Exception ee)
			//					{
			//					}
			//					finally
			//					{
			//						out = null;
			//					}
			//				}
			//			}
		}

		//		if (out != null)
		//		{
		//			try
		//			{
		//				out.close();
		//			}
		//			catch (Exception ee)
		//			{
		//			}
		//		}

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
	 * @param save
	 *            Set to true to save the scored assignments
	 * @param nActual
	 *            The number of actual results in the frame
	 * @return the fractional assignments
	 */
	private FractionalAssignment[] score(final ArrayList<FractionalAssignment> assignments, final double[] score,
			final int nPredicted, boolean save, int nActual)
	{
		if (assignments.isEmpty())
			return null;
		final FractionalAssignment[] tmp = assignments.toArray(new FractionalAssignment[assignments.size()]);
		final RankedScoreCalculator calc = new RankedScoreCalculator(tmp, nActual, nPredicted);
		final double[] result = calc.score(nPredicted, false, save);
		score[0] += result[0];
		score[1] += result[1];
		score[2] += result[2];
		score[3] += result[3];
		assignments.clear();
		return calc.getScoredAssignments();
	}

	/**
	 * Create a dummy collection that implements toArray() without cloning for the addAll() method in ArrayList.
	 */
	private class DummyCollection implements Collection<FractionalAssignment>
	{

		/** The a. */
		final FractionalAssignment[] a;

		/**
		 * Instantiates a new dummy collection.
		 *
		 * @param a
		 *            the a
		 */
		DummyCollection(final FractionalAssignment[] a)
		{
			this.a = a;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Collection#size()
		 */
		public int size()
		{
			return a.length;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Collection#isEmpty()
		 */
		public boolean isEmpty()
		{
			return size() == 0;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Collection#contains(java.lang.Object)
		 */
		public boolean contains(Object o)
		{
			throw new NotImplementedException();
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Collection#iterator()
		 */
		public Iterator<FractionalAssignment> iterator()
		{
			throw new NotImplementedException();
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Collection#toArray()
		 */
		public Object[] toArray()
		{
			return a;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Collection#toArray(java.lang.Object[])
		 */
		public <T> T[] toArray(T[] a)
		{
			return a;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Collection#add(java.lang.Object)
		 */
		public boolean add(FractionalAssignment e)
		{
			throw new NotImplementedException();
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Collection#remove(java.lang.Object)
		 */
		public boolean remove(Object o)
		{
			throw new NotImplementedException();
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Collection#containsAll(java.util.Collection)
		 */
		public boolean containsAll(Collection<?> c)
		{
			throw new NotImplementedException();
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Collection#addAll(java.util.Collection)
		 */
		public boolean addAll(Collection<? extends FractionalAssignment> c)
		{
			throw new NotImplementedException();
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Collection#removeAll(java.util.Collection)
		 */
		public boolean removeAll(Collection<?> c)
		{
			throw new NotImplementedException();
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Collection#retainAll(java.util.Collection)
		 */
		public boolean retainAll(Collection<?> c)
		{
			throw new NotImplementedException();
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Collection#clear()
		 */
		public void clear()
		{
			throw new NotImplementedException();
		}
	}

	/**
	 * To XML.
	 *
	 * @return An XML representation of this object
	 */
	public String toXML()
	{
		return XStreamWrapper.toXML(this);
	}

	/**
	 * Create the filter from the XML representation.
	 *
	 * @param xml
	 *            the xml
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