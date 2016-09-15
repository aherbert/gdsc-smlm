package gdsc.smlm.results.filter;

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
 * Specifies a the result of fitting a frame using different fitting methods.
 * <p>
 * The multi-path results can be evaluated by the MultiPathFilter to determine which result from the different paths
 * should be accepted.
 * <p>
 * This class is used for benchmarking the fitting path options in the PeakFit algorithm.
 */
public class MultiPathFitResults implements IMultiPathFitResults
{
	/**
	 * The frame containing the results
	 */
	final public int frame;

	/**
	 * The multi-path results
	 */
	final public MultiPathFitResult[] multiPathFitResults;

	/**
	 * The total number of candidates. This may be greater than the size of the {@link #multiPathFitResults} array if
	 * this is a subset of the results, i.e. has been prefiltered.
	 */
	final public int totalCandidates;

	private boolean[] estimate = null;

	public MultiPathFitResults(int frame, MultiPathFitResult[] multiPathFitResults)
	{
		this(frame, multiPathFitResults, (multiPathFitResults == null) ? 0 : multiPathFitResults.length);
	}

	public MultiPathFitResults(int frame, MultiPathFitResult[] multiPathFitResults, int totalCandidates)
	{
		this.frame = frame;
		this.multiPathFitResults = multiPathFitResults;
		this.totalCandidates = totalCandidates;
	}

	public int getFrame()
	{
		return frame;
	}

	public int getNumberOfResults()
	{
		return multiPathFitResults.length;
	}

	public MultiPathFitResult getResult(int index)
	{
		return multiPathFitResults[index];
	}

	public int getTotalCandidates()
	{
		return totalCandidates;
	}

	public boolean isValid(int candidateId)
	{
		if (estimate == null)
			return false;
		return estimate[candidateId];
	}

	public void setValid(PreprocessedPeakResult preprocessedPeakResult)
	{
		if (estimate == null)
			estimate = new boolean[totalCandidates];
		estimate[preprocessedPeakResult.getCandidateId()] = true;
	}
}
