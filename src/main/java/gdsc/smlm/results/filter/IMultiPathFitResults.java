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
 * Specifies the result of fitting a frame using different fitting methods.
 * <p>
 * The multi-path results can be evaluated by the MultiPathFilter to determine which result from the different paths
 * should be accepted.
 */
public interface IMultiPathFitResults
{
	/**
	 * @return The frame containing the results
	 */
	int getFrame();

	
	/**
	 * @return The number of results
	 */
	int getNumberOfResults();
	
	
	/**
	 * Gets the result.
	 *
	 * @param index the index
	 * @return the result
	 */
	MultiPathFitResult getResult(int index);

	
	/**
	 * The total number of candidates. This may be greater than the size of the {@link #getNumberOfResults()} if
	 * this is a subset of the results, i.e. has been prefiltered.
	 *
	 * @return the total candidates
	 */
	int getTotalCandidates();


	/**
	 * Checks if is valid.
	 * <p>
	 * Return true if this candidate should definitely be filtered.
	 *
	 * @param candidateId the candidate id
	 * @return true, if is valid
	 */
	boolean isValid(int candidateId);


	/**
	 * Sets result as valid.
	 * <p>
	 * Should be called when the result passed validation and should be marked as valid.
	 *
	 * @param preprocessedPeakResult the result
	 */
	void setValid(PreprocessedPeakResult preprocessedPeakResult);
}
