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
public class MultiPathFitResult implements Cloneable
{
	public static class FitResult
	{
		/**
		 * Fitting status of the fit. Zero for OK.
		 */
		final public int status;

		/**
		 * The results from the fit. It is expected that one or more results will be true for isNewResult() and zero or
		 * more could be true for isExistingResult().
		 */
		public PreprocessedPeakResult[] results;

		/**
		 * Allows storing any data associated with the fit result
		 */
		final public Object data;

		public FitResult(int status)
		{
			this(status, null);
		}

		public FitResult(int status, Object data)
		{
			this.status = status;
			this.data = data;
		}

		public int getStatus()
		{
			return status;
		}

		public PreprocessedPeakResult[] getResults()
		{
			return results;
		}

		public Object getData()
		{
			return data;
		}
	}

	/**
	 * The original value for the number of failed results before this result
	 */
	public int originalFailCount;

	/**
	 * The number of failed results before this result
	 */
	public int failCount;

	/**
	 * Reset the fail count to the original value.
	 * <p>
	 * Note that the MultiPathFilter will update the fail count when creating a subset. Multiple rounds of subsetting
	 * will continually update the fail count. This method can be used to reset the fail count after the subset has been
	 * processed.b
	 */
	public void resetFailCount()
	{
		this.failCount = this.originalFailCount;
	}

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
	 * Return the candidate Id of this result (i.e. the candidate used to identify this position for fitting)
	 */
	public int candidateId;

	/**
	 * The score from residuals analysis on the residuals of the single fit. This can be used to choose if the doublet
	 * fit should be considered.
	 */
	private double singleQAScore = -1;

	/**
	 * The results from the multi-fit. It is expected that one result will be true for isNewResult() and zero or more
	 * could be true for isExistingResult().
	 */
	private FitResult multiFitResult;

	/**
	 * The results from the single-fit. It is expected that this should be one result that is true for isNewResult().
	 */
	private FitResult singleFitResult;

	/**
	 * The results from the doublet-fit. It is expected that this should be one or two results that are true for
	 * isNewResult().
	 */
	private FitResult doubletFitResult;

	@Override
	public MultiPathFitResult clone()
	{
		try
		{
			return (MultiPathFitResult) super.clone();
		}
		catch (CloneNotSupportedException e)
		{
			return null;
		}
	}

	/**
	 * Copy the class level field values into a new object. Ignores the fail count fields.
	 * <p>
	 * To copy sub-class fields use {@link #clone()}.
	 * 
	 * @return A copy
	 */
	public MultiPathFitResult copy()
	{
		final MultiPathFitResult r = new MultiPathFitResult();
		r.candidateId = candidateId;
		r.frame = frame;
		r.width = width;
		r.height = height;
		r.candidateId = candidateId;
		r.singleQAScore = singleQAScore;
		r.multiFitResult = multiFitResult;
		r.singleFitResult = singleFitResult;
		r.doubletFitResult = doubletFitResult;
		return r;
	}

	/**
	 * @return the multiFitResult
	 */
	public FitResult getMultiFitResult()
	{
		return multiFitResult;
	}

	/**
	 * @param multiFitResult
	 *            the multiFitResult to set
	 */
	protected void setMultiFitResult(FitResult multiFitResult)
	{
		this.multiFitResult = multiFitResult;
	}

	/**
	 * @return the singleFitResult
	 */
	public FitResult getSingleFitResult()
	{
		return singleFitResult;
	}

	/**
	 * @param singleFitResult
	 *            the singleFitResult to set
	 */
	protected void setSingleFitResult(FitResult singleFitResult)
	{
		this.singleFitResult = singleFitResult;
	}

	/**
	 * @param residualsThreshold 
	 * @return the singleQAScore
	 */
	public double getSingleQAScore()
	{
		return singleQAScore;
	}

	/**
	 * @param singleQAScore
	 *            the singleQAScore to set
	 */
	protected void setSingleQAScore(double singleQAScore)
	{
		this.singleQAScore = singleQAScore;
	}

	/**
	 * @return the doubletFitResult
	 */
	public FitResult getDoubletFitResult()
	{
		return doubletFitResult;
	}

	/**
	 * @param doubletFitResult
	 *            the doubletFitResult to set
	 */
	protected void setDoubletFitResult(FitResult doubletFitResult)
	{
		this.doubletFitResult = doubletFitResult;
	}
}
