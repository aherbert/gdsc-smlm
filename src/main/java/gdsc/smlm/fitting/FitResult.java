package gdsc.smlm.fitting;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Contains the fitting result
 */
public class FitResult
{
	/**
	 * Provides a builder to allow simple adjustments to the immutable FitResult fields to create a new FitResult.
	 */
	public class Builder
	{
		private FitStatus status;
		private int degreesOfFreedom;
		private double error;
		private double[] initialParameters;
		private double[] parameters;
		private double[] parameterDevs;
		private int nPeaks;
		private int nFittedParameters;
		private Object data;
		private int iterations;
		private int evaluations;

		public Builder(FitStatus status, int degreesOfFreedom, double error, double[] initialParameters,
				double[] parameters, double[] parameterDevs, int nPeaks, int nFittedParameters, Object data,
				int iterations, int evaluations)
		{
			this.status = status;
			this.degreesOfFreedom = degreesOfFreedom;
			this.error = error;
			this.initialParameters = initialParameters;
			this.parameters = parameters;
			this.parameterDevs = parameterDevs;
			this.nPeaks = nPeaks;
			this.nFittedParameters = nFittedParameters;
			this.data = data;
			this.iterations = iterations;
			this.evaluations = evaluations;
		}

		public FitStatus getStatus()
		{
			return status;
		}

		public Builder setStatus(FitStatus status)
		{
			this.status = status;
			return this;
		}

		public int getDegreesOfFreedom()
		{
			return degreesOfFreedom;
		}

		public Builder setDegreesOfFreedom(int degreesOfFreedom)
		{
			this.degreesOfFreedom = degreesOfFreedom;
			return this;
		}

		public double getError()
		{
			return error;
		}

		public Builder setError(double error)
		{
			this.error = error;
			return this;
		}

		public double[] getInitialParameters()
		{
			return initialParameters;
		}

		public Builder setInitialParameters(double[] initialParameters)
		{
			this.initialParameters = initialParameters;
			return this;
		}

		public double[] getParameters()
		{
			return parameters;
		}

		public Builder setParameters(double[] parameters)
		{
			this.parameters = parameters;
			return this;
		}

		public double[] getParameterDeviations()
		{
			return parameterDevs;
		}

		public Builder setParameterDeviations(double[] parameterDevs)
		{
			this.parameterDevs = parameterDevs;
			return this;
		}

		public int getnPeaks()
		{
			return nPeaks;
		}

		public Builder setnPeaks(int nPeaks)
		{
			this.nPeaks = nPeaks;
			return this;
		}

		public int getnFittedParameters()
		{
			return nFittedParameters;
		}

		public Builder setnFittedParameters(int nFittedParameters)
		{
			this.nFittedParameters = nFittedParameters;
			return this;
		}

		public Object getData()
		{
			return data;
		}

		public Builder setData(Object data)
		{
			this.data = data;
			return this;
		}

		public int getIterations()
		{
			return iterations;
		}

		public Builder setIterations(int iterations)
		{
			this.iterations = iterations;
			return this;
		}

		public int getEvaluations()
		{
			return evaluations;
		}

		public Builder setEvaluations(int evaluations)
		{
			this.evaluations = evaluations;
			return this;
		}

		public FitResult build()
		{
			return new FitResult(status, degreesOfFreedom, error, initialParameters, initialParameters, parameterDevs,
					nPeaks, nFittedParameters, data, iterations, evaluations);
		}
	}

	private FitStatus status;
	private final int degreesOfFreedom;
	private double error;
	private final double[] initialParameters;
	private final double[] parameters;
	private final double[] parameterDevs;
	private final int nPeaks;
	private final int nFittedParameters;
	private Object data;
	private final int iterations, evaluations;

	/**
	 * Constructor
	 * 
	 * @param status
	 * @param degreesOfFreedom
	 * @param error
	 * @param initialParameters
	 * @param parameters
	 * @param parameterDevs
	 * @param nPeaks
	 * @param nFittedParameters
	 * @param data
	 */
	public FitResult(FitStatus status, int degreesOfFreedom, double error, double[] initialParameters,
			double[] parameters, double[] parameterDevs, int nPeaks, int nFittedParameters, Object data, int iterations,
			int evaluations)
	{
		this.status = status;
		this.degreesOfFreedom = degreesOfFreedom;
		this.error = error;
		this.initialParameters = initialParameters;
		this.parameters = parameters;
		this.parameterDevs = parameterDevs;
		this.nPeaks = nPeaks;
		this.nFittedParameters = nFittedParameters;
		this.data = data;
		this.iterations = iterations;
		this.evaluations = evaluations;
	}

	/**
	 * Constructor
	 * 
	 * @param result
	 */
	public FitResult(FitStatus result)
	{
		this.status = result;
		this.degreesOfFreedom = 0;
		this.error = 0;
		this.initialParameters = null;
		this.parameters = null;
		this.parameterDevs = null;
		this.nPeaks = 0;
		this.nFittedParameters = 0;
		this.data = null;
		this.iterations = 0;
		this.evaluations = 0;
	}

	/**
	 * @return the fit status
	 */
	public FitStatus getStatus()
	{
		return status;
	}

	/**
	 * @return the degrees of freedom for the fit
	 */
	public int getDegreesOfFreedom()
	{
		return degreesOfFreedom;
	}

	/**
	 * @return the error value for the fit (e.g. mean-squared error or the reduced Chi-squared measure)
	 */
	public double getError()
	{
		return error;
	}

	/**
	 * Sets the error. This can be used to update the error using a different metric.
	 *
	 * @param error
	 *            the new error
	 */
	public void setError(double error)
	{
		this.error = error;
	}

	/**
	 * Gets the initial parameters.
	 *
	 * @return the initial parameters
	 */
	public double[] getInitialParameters()
	{
		return initialParameters;
	}

	/**
	 * Get the fitted parameters
	 * 
	 * @return the fitted parameters
	 */
	public double[] getParameters()
	{
		return parameters;
	}

	/**
	 * Get the deviations (variances) for the fitted parameters
	 * 
	 * @return the parameter deviations
	 */
	public double[] getParameterDeviations()
	{
		return parameterDevs;
	}

	/**
	 * @return the number of fitted peaks
	 */
	public int getNumberOfPeaks()
	{
		return nPeaks;
	}

	/**
	 * @return the number of fitted parameters (e.g. fixed width fitting will reduce the number of fitted parameters)
	 */
	public int getNumberOfFittedParameters()
	{
		return nFittedParameters;
	}

	/**
	 * Returns an object containing data about the fit status. This is used to pass additional information depending on
	 * the fitting status.
	 * 
	 * @return the data
	 */
	public Object getStatusData()
	{
		return data;
	}

	/**
	 * @return the iterations
	 */
	public int getIterations()
	{
		return iterations;
	}

	/**
	 * @return the evaluations
	 */
	public int getEvaluations()
	{
		return evaluations;
	}

	/**
	 * Sets the status.
	 *
	 * @param fitStatus
	 *            the fit status
	 * @param data
	 *            the data
	 */
	public void setStatus(FitStatus fitStatus, Object data)
	{
		this.status = fitStatus;
		this.data = data;
	}

	/**
	 * Create a builder to allow creation of a new result. The builder is initialised using the current field values.
	 *
	 * @return the builder
	 */
	public Builder toBuilder()
	{
		return new Builder(status, degreesOfFreedom, error, initialParameters, parameters, parameterDevs, nPeaks,
				nFittedParameters, data, iterations, evaluations);
	}
}