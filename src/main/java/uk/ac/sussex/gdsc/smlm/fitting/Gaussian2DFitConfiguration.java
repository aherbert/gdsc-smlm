/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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
 * Specifies the fitting configuration and validation for Gaussian 2D fitting.
 */
public interface Gaussian2DFitConfiguration {
  /**
   * Creates the appropriate stopping criteria and Gaussian function for the configuration.
   *
   * @param npeaks The number of peaks to fit
   * @param maxx The width of the XY data
   * @param maxy The height of the XY data
   */
  void initialise(int npeaks, int maxx, int maxy);

  /**
   * True if the fit residuals should be computed.
   *
   * @return true if the fit residuals should be computed
   */
  boolean isComputeResiduals();

  /**
   * Gets the function solver.
   *
   * @return the function solver
   */
  FunctionSolver getFunctionSolver();

  /**
   * True if the fit parameter deviations should be computed.
   *
   * @return true if the fit parameter deviations should be computed
   */
  boolean isComputeDeviations();

  /**
   * Checks if is background fitting.
   *
   * @return true, if is background fitting
   */
  boolean isBackgroundFitting();

  /**
   * Checks if is XSD fitting.
   *
   * @return true, if is XSD fitting
   */
  boolean isXSdFitting();

  /**
   * Checks if is YSD fitting.
   *
   * @return true, if is YSD fitting
   */
  boolean isYSdFitting();

  /**
   * Checks if is angle fitting.
   *
   * @return true, if is angle fitting
   */
  boolean isAngleFitting();

  /**
   * Checks if is z fitting.
   *
   * <p>If z fitting then it is assumed that the function is not fitting X SD or Y SD and the widths
   * are entirely defined by the z position. However the configuration should still return a value
   * for the initial X SD and Y SD, for example the width at z=0. This can be used to estimate the
   * centre of the Gaussian using its approximate size, or convert amplitude estimates to signal
   * (total volume of the Gaussian).
   *
   * @return true, if is z fitting
   */
  boolean isZFitting();

  /**
   * Gets the initial guess for the XSD parameter.
   *
   * @return the initial XSD
   */
  double getInitialXSd();

  /**
   * Gets the initial guess for the YSD parameter.
   *
   * @return the initial YSD
   */
  double getInitialYSd();

  /**
   * Gets the initial guess for the angle parameter.
   *
   * @return the initial angle
   */
  double getInitialAngle();

  /**
   * Gets the minimum width factor. This is used to limit the bounds of width fitting.
   *
   * @return the minimum width factor
   */
  double getMinWidthFactor();

  /**
   * Gets the maximum width factor. This is used to limit the bounds of width fitting.
   *
   * @return the maximum width factor
   */
  double getMaxWidthFactor();

  /**
   * Checks if fit validation should be performed using validateFit.
   *
   * @return true, if fit validation should be performed
   */
  boolean isFitValidation();

  /**
   * Check peaks to see if the fit was sensible. This is called after fitting so that the parameters
   * can be checked.
   *
   * @param peakCount The number of peaks
   * @param initialParams The initial peak parameters
   * @param params The fitted peak parameters
   * @param paramDevs the fitted peak parameters variances (can be null)
   * @return the fit status
   */
  FitStatus validateFit(int peakCount, double[] initialParams, double[] params, double[] paramDevs);

  /**
   * Gets the validation data. This can be set within the validateFit function.
   *
   * @return Data associated with the validation result
   */
  Object getValidationData();
}
