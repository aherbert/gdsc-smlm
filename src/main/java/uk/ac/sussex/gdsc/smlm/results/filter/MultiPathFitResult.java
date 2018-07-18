/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package uk.ac.sussex.gdsc.smlm.results.filter;

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
	/**
	 * The fit result.
	 */
	public static class FitResult implements Cloneable
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

		/**
		 * Instantiates a new fit result.
		 *
		 * @param status
		 *            the status
		 */
		public FitResult(int status)
		{
			this(status, null);
		}

		/**
		 * Instantiates a new fit result.
		 *
		 * @param status
		 *            the status
		 * @param data
		 *            the data
		 */
		public FitResult(int status, Object data)
		{
			this.status = status;
			this.data = data;
		}

		/**
		 * Gets the status.
		 *
		 * @return the status
		 */
		public int getStatus()
		{
			return status;
		}

		/**
		 * Gets the results.
		 *
		 * @return the results
		 */
		public PreprocessedPeakResult[] getResults()
		{
			return results;
		}

		/**
		 * Gets the data.
		 *
		 * @return the data
		 */
		public Object getData()
		{
			return data;
		}

		@Override
		public FitResult clone()
		{
			try
			{
				return (FitResult) super.clone();
			}
			catch (final CloneNotSupportedException e)
			{
				return null;
			}
		}
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
	 * The score from residuals analysis on the residuals of the multi fit. This can be used to choose if the doublet
	 * fit should be considered.
	 */
	private double multiQAScore = -1;

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
	 * The results from the doublet-fit on the multi-fit residuals. It is expected that this should be one or two
	 * results that are true for isNewResult().
	 */
	private FitResult multiDoubletFitResult;

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
		catch (final CloneNotSupportedException e)
		{
			return null;
		}
	}

	/**
	 * Copy the class level field values into a new object. Ignores the fail count fields.
	 * <p>
	 * To copy sub-class fields use {@link #clone()}.
	 *
	 * @param deep
	 *            Set to true to do a clone of the FitResult objects. Their array objects will not be copied.
	 * @return A copy
	 */
	public MultiPathFitResult copy(boolean deep)
	{
		final MultiPathFitResult r = new MultiPathFitResult();
		r.candidateId = candidateId;
		r.frame = frame;
		r.width = width;
		r.height = height;
		r.candidateId = candidateId;
		r.multiQAScore = multiQAScore;
		r.singleQAScore = singleQAScore;
		if (deep)
		{
			r.multiFitResult = clone(multiFitResult);
			r.multiDoubletFitResult = clone(multiDoubletFitResult);
			r.singleFitResult = clone(singleFitResult);
			r.doubletFitResult = clone(doubletFitResult);
		}
		else
		{
			r.multiFitResult = multiFitResult;
			r.multiDoubletFitResult = multiDoubletFitResult;
			r.singleFitResult = singleFitResult;
			r.doubletFitResult = doubletFitResult;
		}
		return r;
	}

	private static FitResult clone(FitResult f)
	{
		return (f == null) ? null : f.clone();
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
	 * @return the multiQAScore
	 */
	public double getMultiQAScore()
	{
		return multiQAScore;
	}

	/**
	 * @param multiQAScore
	 *            the multiQAScore to set
	 */
	protected void setMultiQAScore(double multiQAScore)
	{
		this.multiQAScore = multiQAScore;
	}

	/**
	 * @return the multiDoubletFitResult
	 */
	public FitResult getMultiDoubletFitResult()
	{
		return multiDoubletFitResult;
	}

	/**
	 * @param multiDoubletFitResult
	 *            the multiDoubletFitResult to set
	 */
	protected void setMultiDoubletFitResult(FitResult multiDoubletFitResult)
	{
		this.multiDoubletFitResult = multiDoubletFitResult;
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
