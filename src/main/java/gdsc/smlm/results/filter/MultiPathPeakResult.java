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
 * Specifies a the result of fitting a position using different fitting methods.
 * <p>
 * The multi-path result can be evaluated by the MultiPathFilter to determine which result from the different paths
 * should be accepted.
 * <p>
 * This class is used for benchmarking the fitting path options in the PeakFit algorithm.
 */
public class MultiPathPeakResult
{
	/**
	 * The number of failed results before this result
	 */
	public int failCount;
	
	/**
	 * The frame containing the result
	 */
	public int frame;
	
	// TODO - add multi-fit result, single fit result, doublet fit result
}
