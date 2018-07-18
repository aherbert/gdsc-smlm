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
 * Specifies the fitting configuration for Gaussian 2D fitting.
 */
public interface Gaussian2DFitConfiguration
{
	/**
	 * Creates the appropriate stopping criteria and Gaussian function for the configuration.
	 *
	 * @param npeaks
	 *            The number of peaks to fit
	 * @param maxx
	 *            The height of the XY data
	 * @param maxy
	 *            the maxy
	 * @param params
	 *            The Gaussian parameters
	 */
	public void initialise(int npeaks, int maxx, int maxy, double[] params);

	/**
	 * True if the fit residuals should be computed
	 *
	 * @return true if the fit residuals should be computed
	 */
	public boolean isComputeResiduals();

	/**
	 * Gets the function solver.
	 *
	 * @return the function solver
	 */
	public FunctionSolver getFunctionSolver();

	/**
	 * True if the fit parameter deviations should be computed
	 *
	 * @return true if the fit parameter deviations should be computed
	 */
	public boolean isComputeDeviations();

	/**
	 * Checks if is background fitting.
	 *
	 * @return true, if is background fitting
	 */
	public boolean isBackgroundFitting();

	/**
	 * Checks if is XSD fitting.
	 *
	 * @return true, if is XSD fitting
	 */
	public boolean isXSDFitting();

	/**
	 * Checks if is YSD fitting.
	 *
	 * @return true, if is YSD fitting
	 */
	public boolean isYSDFitting();

	/**
	 * Checks if is angle fitting.
	 *
	 * @return true, if is angle fitting
	 */
	public boolean isAngleFitting();

	/**
	 * Checks if is z fitting.
	 * <p>
	 * If z fitting then it is assumed that the function is not fitting X SD or Y SD and the widths are entirely defined
	 * by the z position. However the configuration should still return a value for the initial X SD and Y SD, for
	 * example the width at z=0. This can be used to estimate the centre of the Gaussian using its approximate size, or
	 * convert amplitude estimates to signal (total volume of the Gaussian).
	 *
	 * @return true, if is z fitting
	 */
	public boolean isZFitting();

	/**
	 * Gets the initial guess for the XSD parameter.
	 *
	 * @return the initial XSD
	 */
	public double getInitialXSD();

	/**
	 * Gets the initial guess for the YSD parameter.
	 *
	 * @return the initial YSD
	 */
	public double getInitialYSD();

	/**
	 * Gets the initial guess for the angle parameter.
	 *
	 * @return the initial angle
	 */
	public double getInitialAngle();

	/**
	 * Gets the minimum width factor. This is used to limit the bounds of width fitting.
	 *
	 * @return the minimum width factor
	 */
	public double getMinWidthFactor();

	/**
	 * Gets the maximum width factor. This is used to limit the bounds of width fitting.
	 *
	 * @return the maximum width factor
	 */
	public double getMaxWidthFactor();

	/**
	 * Checks if fit validation should be performed using validateFit.
	 *
	 * @return true, if fit validation should be performed
	 */
	public boolean isFitValidation();

	/**
	 * Check peaks to see if the fit was sensible. This is called after fitting so that the parameters can be checked.
	 *
	 * @param nPeaks
	 *            The number of peaks
	 * @param initialParams
	 *            The initial peak parameters
	 * @param params
	 *            The fitted peak parameters
	 * @param paramDevs
	 *            the fitted peak parameters variances (can be null)
	 * @return the fit status
	 */
	public FitStatus validateFit(int nPeaks, double[] initialParams, double[] params, double[] paramDevs);

	/**
	 * Gets the validation data. This can be set within the validateFit function.
	 *
	 * @return Data associated with the validation result
	 */
	public Object getValidationData();
}
