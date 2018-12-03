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
package uk.ac.sussex.gdsc.smlm.fitting;

/**
 * Contains the fitting result.
 */
public class FitResult {
  /**
   * Provides a builder to allow simple adjustments to the immutable FitResult fields to create a
   * new FitResult.
   */
  public class Builder {
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

    /**
     * Instantiates a new builder.
     *
     * @param status the status
     * @param degreesOfFreedom the degrees of freedom
     * @param error the error
     * @param initialParameters the initial parameters
     * @param parameters the parameters
     * @param parameterDevs the parameter deviations
     * @param nPeaks the number of peaks
     * @param nFittedParameters the number of fitted parameters
     * @param data the data
     * @param iterations the iterations
     * @param evaluations the evaluations
     */
    public Builder(FitStatus status, int degreesOfFreedom, double error, double[] initialParameters,
        double[] parameters, double[] parameterDevs, int nPeaks, int nFittedParameters, Object data,
        int iterations, int evaluations) {
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
     * Gets the status.
     *
     * @return the status
     */
    public FitStatus getStatus() {
      return status;
    }

    /**
     * Sets the status.
     *
     * @param status the status
     * @return the builder
     */
    public Builder setStatus(FitStatus status) {
      this.status = status;
      return this;
    }

    /**
     * Gets the degrees of freedom.
     *
     * @return the degrees of freedom
     */
    public int getDegreesOfFreedom() {
      return degreesOfFreedom;
    }

    /**
     * Sets the degrees of freedom.
     *
     * @param degreesOfFreedom the degrees of freedom
     * @return the builder
     */
    public Builder setDegreesOfFreedom(int degreesOfFreedom) {
      this.degreesOfFreedom = degreesOfFreedom;
      return this;
    }

    /**
     * Gets the error.
     *
     * @return the error
     */
    public double getError() {
      return error;
    }

    /**
     * Sets the error.
     *
     * @param error the error
     * @return the builder
     */
    public Builder setError(double error) {
      this.error = error;
      return this;
    }

    /**
     * Gets the initial parameters.
     *
     * @return the initial parameters
     */
    public double[] getInitialParameters() {
      return initialParameters;
    }

    /**
     * Sets the initial parameters.
     *
     * @param initialParameters the initial parameters
     * @return the builder
     */
    public Builder setInitialParameters(double[] initialParameters) {
      this.initialParameters = initialParameters;
      return this;
    }

    /**
     * Gets the parameters.
     *
     * @return the parameters
     */
    public double[] getParameters() {
      return parameters;
    }

    /**
     * Sets the parameters.
     *
     * @param parameters the parameters
     * @return the builder
     */
    public Builder setParameters(double[] parameters) {
      this.parameters = parameters;
      return this;
    }

    /**
     * Gets the parameter deviations.
     *
     * @return the parameter deviations
     */
    public double[] getParameterDeviations() {
      return parameterDevs;
    }

    /**
     * Sets the parameter deviations.
     *
     * @param parameterDevs the parameter deviations
     * @return the builder
     */
    public Builder setParameterDeviations(double[] parameterDevs) {
      this.parameterDevs = parameterDevs;
      return this;
    }

    /**
     * Gets the number of peaks.
     *
     * @return the number of peaks
     */
    public int getNumberOfPeaks() {
      return nPeaks;
    }

    /**
     * Set the number of peaks.
     *
     * @param nPeaks the number of peaks
     * @return the builder
     */
    public Builder setNumberOfPeaks(int nPeaks) {
      this.nPeaks = nPeaks;
      return this;
    }

    /**
     * Gets the number of fitted parameters.
     *
     * @return the number of fitted parameters
     */
    public int getNumberOfFittedParameters() {
      return nFittedParameters;
    }

    /**
     * Set the number of fitted parameters.
     *
     * @param nFittedParameters the number of fitted parameters
     * @return the builder
     */
    public Builder setNumberOfFittedParameters(int nFittedParameters) {
      this.nFittedParameters = nFittedParameters;
      return this;
    }

    /**
     * Gets the data.
     *
     * @return the data
     */
    public Object getData() {
      return data;
    }

    /**
     * Sets the data.
     *
     * @param data the data
     * @return the builder
     */
    public Builder setData(Object data) {
      this.data = data;
      return this;
    }

    /**
     * Gets the iterations.
     *
     * @return the iterations
     */
    public int getIterations() {
      return iterations;
    }

    /**
     * Sets the iterations.
     *
     * @param iterations the iterations
     * @return the builder
     */
    public Builder setIterations(int iterations) {
      this.iterations = iterations;
      return this;
    }

    /**
     * Gets the evaluations.
     *
     * @return the evaluations
     */
    public int getEvaluations() {
      return evaluations;
    }

    /**
     * Sets the evaluations.
     *
     * @param evaluations the evaluations
     * @return the builder
     */
    public Builder setEvaluations(int evaluations) {
      this.evaluations = evaluations;
      return this;
    }

    /**
     * Builds the fit result.
     *
     * @return the fit result
     */
    public FitResult build() {
      return new FitResult(status, degreesOfFreedom, error, initialParameters, initialParameters,
          parameterDevs, nPeaks, nFittedParameters, data, iterations, evaluations);
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
   * Constructor.
   *
   * @param status the status
   * @param degreesOfFreedom the degrees of freedom
   * @param error the error
   * @param initialParameters the initial parameters
   * @param parameters the parameters
   * @param parameterDevs the parameter deviations
   * @param nPeaks the number of peaks
   * @param nFittedParameters the number of fitted parameters
   * @param data the data
   * @param iterations the iterations
   * @param evaluations the evaluations
   */
  public FitResult(FitStatus status, int degreesOfFreedom, double error, double[] initialParameters,
      double[] parameters, double[] parameterDevs, int nPeaks, int nFittedParameters, Object data,
      int iterations, int evaluations) {
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
   * Constructor.
   *
   * @param result the result
   */
  public FitResult(FitStatus result) {
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
   * Gets the status.
   *
   * @return the fit status
   */
  public FitStatus getStatus() {
    return status;
  }

  /**
   * Gets the degrees of freedom.
   *
   * @return the degrees of freedom for the fit
   */
  public int getDegreesOfFreedom() {
    return degreesOfFreedom;
  }

  /**
   * Gets the error.
   *
   * @return the error value for the fit (e.g. mean-squared error or the reduced Chi-squared
   *         measure)
   */
  public double getError() {
    return error;
  }

  /**
   * Sets the error. This can be used to update the error using a different metric.
   *
   * @param error the new error
   */
  public void setError(double error) {
    this.error = error;
  }

  /**
   * Gets the initial parameters.
   *
   * @return the initial parameters
   */
  public double[] getInitialParameters() {
    return initialParameters;
  }

  /**
   * Get the fitted parameters.
   *
   * @return the fitted parameters
   */
  public double[] getParameters() {
    return parameters;
  }

  /**
   * Get the deviations (variances) for the fitted parameters.
   *
   * @return the parameter deviations
   */
  public double[] getParameterDeviations() {
    return parameterDevs;
  }

  /**
   * Gets the number of peaks.
   *
   * @return the number of fitted peaks
   */
  public int getNumberOfPeaks() {
    return nPeaks;
  }

  /**
   * Gets the number of fitted parameters.
   *
   * @return the number of fitted parameters (e.g. fixed width fitting will reduce the number of
   *         fitted parameters)
   */
  public int getNumberOfFittedParameters() {
    return nFittedParameters;
  }

  /**
   * Returns an object containing data about the fit status. This is used to pass additional
   * information depending on the fitting status.
   *
   * @return the data
   */
  public Object getStatusData() {
    return data;
  }

  /**
   * Gets the iterations.
   *
   * @return the iterations
   */
  public int getIterations() {
    return iterations;
  }

  /**
   * Gets the evaluations.
   *
   * @return the evaluations
   */
  public int getEvaluations() {
    return evaluations;
  }

  /**
   * Sets the status.
   *
   * @param fitStatus the fit status
   * @param data the data
   */
  public void setStatus(FitStatus fitStatus, Object data) {
    this.status = fitStatus;
    this.data = data;
  }

  /**
   * Create a builder to allow creation of a new result. The builder is initialised using the
   * current field values.
   *
   * @return the builder
   */
  public Builder toBuilder() {
    return new Builder(status, degreesOfFreedom, error, initialParameters, parameters,
        parameterDevs, nPeaks, nFittedParameters, data, iterations, evaluations);
  }
}
