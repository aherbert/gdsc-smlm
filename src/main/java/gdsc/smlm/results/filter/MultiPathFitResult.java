package gdsc.smlm.results.filter;

import gdsc.smlm.fitting.FitStatus;

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
public class MultiPathFitResult
{
	/**
	 * The number of failed results before this result
	 */
	public int failCount;

	/**
	 * The frame containing the result
	 */
	public int frame;

	/**
	 * The width of the fit region
	 */
	public int width;

	/**
	 * The height of the fit region
	 */
	public int height;

	/**
	 * The score from residuals analysis on the residuals of the single fit. This can be used to choose if the doublet
	 * fit should be considered.
	 */
	public double singleQAScore = 2;

	/**
	 * The results from the multi-fit. It is expected that one result will be true for isNewResult() and zero or more
	 * could be true for isExistingResult().
	 */
	public PreprocessedPeakResult[] multiFitResult;

	/**
	 * Fitting status of the multi-fit
	 */
	public FitStatus multiFitResultStatus;

	/**
	 * The results from the single-fit. It is expected that this should be one result that is true for isNewResult().
	 */
	public PreprocessedPeakResult[] singleFitResult;

	/**
	 * Fitting status of the single-fit
	 */
	public FitStatus singleFitResultStatus;
	
	/**
	 * The results from the doublet-fit. It is expected that this should be one or two results that are true for
	 * isNewResult().
	 */
	public PreprocessedPeakResult[] doubletFitResult;
	
	/**
	 * Fitting status of the doublet-fit
	 */
	public FitStatus doubletFitResultStatus;
}
