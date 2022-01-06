/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear;

import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.function.GradientFunction;
import uk.ac.sussex.gdsc.smlm.function.NonLinearFunction;

/**
 * Uses Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a set of data
 * points (x, y).
 *
 * <p>Support bounded parameters using a hard-stop limit.
 *
 * <p>Support parameter clamping to prevent large parameter shifts. Optionally update the clamping
 * when the search direction changes.
 */
public class BoundedNonLinearFit extends NonLinearFit {
  private ParameterBounds bounds;

  /**
   * Default constructor.
   *
   * @param func The function to fit
   * @param bounds the bounds
   */
  public BoundedNonLinearFit(NonLinearFunction func, ParameterBounds bounds) {
    this(func, null, bounds);
  }

  /**
   * Default constructor.
   *
   * @param func The function to fit
   * @param sc The stopping criteria
   * @param bounds the bounds
   */
  public BoundedNonLinearFit(NonLinearFunction func, StoppingCriteria sc, ParameterBounds bounds) {
    super(func, sc);
    setBounds(bounds);
  }

  /**
   * Default constructor.
   *
   * @param func The function to fit
   * @param sc The stopping criteria
   * @param significantDigits Validate the Levenberg-Marquardt fit solution to the specified number
   *        of significant digits
   * @param maxAbsoluteError Validate the Levenberg-Marquardt fit solution using the specified
   *        maximum absolute error
   * @param bounds the bounds
   */
  public BoundedNonLinearFit(NonLinearFunction func, StoppingCriteria sc, int significantDigits,
      double maxAbsoluteError, ParameterBounds bounds) {
    super(func, sc, significantDigits, maxAbsoluteError);
    setBounds(bounds);
  }

  @Override
  protected boolean solve(double[] a, final int np) {
    if (super.solve(a, np)) {
      return true;
    }

    // If using a bounded LVM is there a chance that the gradient against the bounds will
    // be very large and effect the linear decomposition of the matrix?
    // If decomposition fails try again but set the bounded params to zero (these are
    // ignored by the solver), thus skipping these params for this iteration.

    if (bounds.atBounds(a)) {
      // System.out.printf("Failed when point was at the bounds\n");

      // This functionality has been removed since the solver no longer
      // extracts rows/columns from the matrix if the gradient vector is zero.
      // TODO - Determine if this support is needed.
    }

    return false;
  }

  @Override
  protected void updateFitParameters(double[] a, int[] gradientIndices, int np, double[] da,
      double[] ap) {
    bounds.applyBounds(a, da, ap);
  }

  @Override
  protected void accepted(double[] a, double[] ap, int np) {
    bounds.accepted(a, ap);
    super.accepted(a, ap, np);
  }

  @Override
  public FitStatus computeFit(double[] y, double[] fx, double[] a, double[] parametersVariance) {
    bounds.initialise();
    return super.computeFit(y, fx, a, parametersVariance);
  }

  @Override
  public boolean isBounded() {
    return true;
  }

  @Override
  public boolean isConstrained() {
    return false;
  }

  /**
   * {@inheritDoc}
   *
   * @throws IllegalArgumentException If the lower bound is above the upper bound
   */
  @Override
  public void setBounds(double[] lowerB, double[] upperB) {
    bounds.setBounds(lowerB, upperB);
  }

  /**
   * Sets the bounds. This can be used for dynamic clamping of the parameter updates.
   *
   * @param bounds the new bounds
   */
  public void setBounds(ParameterBounds bounds) {
    if (bounds == null) {
      bounds = new ParameterBounds(function);
    }
    this.bounds = bounds;
  }

  /**
   * Warning: If the function is changed then the clamp values may require updating. However setting
   * a new function does not set the clamp values to null to allow caching when the clamp values are
   * unchanged, e.g. evaluation of a different function in the same parameter space.
   *
   * <p>Setting a new function removes the current bounds.
   *
   * @param function the new gradient function
   */
  @Override
  public void setGradientFunction(GradientFunction function) {
    super.setGradientFunction(function);
    if (bounds != null) {
      bounds.setGradientFunction(function);
    }
  }
}
