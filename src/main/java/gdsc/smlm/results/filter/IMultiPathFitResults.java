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
	 * Get the number of results. The {@link #getResult(int)} method should support being called with any index up to
	 * the number of results (exclusive).
	 * 
	 * @return The number of results
	 */
	int getNumberOfResults();

	/**
	 * Gets the result.
	 *
	 * @param index
	 *            the index
	 * @return the result
	 */
	MultiPathFitResult getResult(int index);
	
	/**
	 * Called when the results that would be returned by {@link #getResult(int)} are no longer required
	 *
	 * @param index the index
	 */
	void complete(int index);

	/**
	 * The total number of candidates. This may be greater than the size of the {@link #getNumberOfResults()} if
	 * this is a subset of the results, i.e. has been prefiltered.
	 *
	 * @return the total candidates
	 */
	int getTotalCandidates();

	// Possible support for iteration
	//	/**
	//	 * Begin. Called before a pass through the results using {@link #getResult(int)}.
	//	 *
	//	 * @return true, if a pass through the results is possible
	//	 */
	//	boolean begin();
	//
	//	/**
	//	 * Called after a pass through the results. Returns a boolean indicating a repeat is possible.
	//	 * 
	//	 * @return true, if another pass through the results is possible
	//	 */
	//	boolean end();
}
