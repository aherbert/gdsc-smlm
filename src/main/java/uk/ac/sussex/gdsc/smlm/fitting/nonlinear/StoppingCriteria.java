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

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear;

import java.util.logging.Logger;

/**
 * Defines the stopping criteria for the
 * {@link uk.ac.sussex.gdsc.smlm.fitting.nonlinear.NonLinearFit } class.
 *
 * <pre>
 * {@code
 * StoppingCritera sc; // Passed in
 *
 * double oldError = // ... Initialised;
 * sc.initialise(a);
 * while (sc.areNotSatisfied())
 * {
 *     // do fitting
 *     newError = fit(...);
 *
 *     sc.evaluate(oldError, newError, a);
 *
 *     // Set best error and update parameters a
 *     if (newError < oldError)
 *     {
 *         oldError = newError;
 *         ...
 *     }
 * }
 *
 * if (sc.areAchieved())
 *     // ... Do something
 * }
 * </pre>
 */
public abstract class StoppingCriteria {
  /** The logger. */
  protected Logger log;
  private int iteration;
  /** Set to true if iteration should continue (no reason to stop). */
  protected boolean notSatisfied;
  /** Set to true if the stopping criteria are achieved. */
  protected boolean areAchieved;
  /** The best parameters. */
  protected double[] bestA;
  private int minimumIterations;
  private int maximumIterations = 20;

  /**
   * Called at the start of the fit. It should be used to reset all iteration counters.
   *
   * @param a Set of m coefficients for the fit
   */
  public void initialise(double[] a) {
    iteration = 0;
    notSatisfied = true;
    areAchieved = false;
    copyCoefficients(a);
  }

  /**
   * Perform a deep copy of the fit coefficients array and store in a variable for later comparison.
   *
   * @param a the a
   */
  protected void copyCoefficients(double[] a) {
    this.bestA = a.clone();
  }

  /**
   * Called after each iteration of the fit. If the error value has reduced then the class should
   * decide if the fit is good enough (set areAchived to true) and can be stopped (set notSatisfied
   * to false).
   *
   * <p>If the fitting process appears to be failing then notSatisfied can be set to false to stop
   * iterating.
   *
   * <p>Note that if the error value is higher the NonLinearFit class will not update the fit
   * coefficients.
   *
   * @param oldError Previous error value for the fit
   * @param newError New error value for the fit
   * @param a Set of m coefficients for the current best fit (matches whichever error is lowest)
   */
  public abstract void evaluate(double oldError, double newError, double[] a);

  /**
   * Increment the current iteration and set the best coefficient values (using
   * {@link #copyCoefficients(double[]) }) if the fit was improved.
   *
   * <p>Sets the notSatisfied flag to false if the maximum number of iterations is reached.
   *
   * @param a The parameters
   * @param improved Flag to indicate if the parameters have improved the fit
   */
  protected void increment(double[] a, boolean improved) {
    iteration++;
    if (improved) {
      copyCoefficients(a);
    }

    // Stop if the maximum iterations have been reached
    if (iteration >= maximumIterations) {
      notSatisfied = false;
    }
  }

  /**
   * Called after each {@link #evaluate(double, double, double[]) } method to check if the fitting
   * should continue.
   *
   * @return True if the stopping criteria have not been met (i.e. fitting should continue)
   */
  public boolean areNotSatisfied() {
    return notSatisfied;
  }

  /**
   * Called once {@link #areNotSatisfied()} returns false to check if the stopping criteria were
   * successfully achieved or the fitting was terminated.
   *
   * @return True if the stopping criteria are achieved
   */
  public boolean areAchieved() {
    return areAchieved;
  }

  /**
   * Gets the iteration.
   *
   * @return the iteration.
   */
  public int getIteration() {
    return iteration;
  }

  /**
   * Sets the maximum iterations.
   *
   * @param maximumIterations the new maximum iterations
   */
  public void setMaximumIterations(int maximumIterations) {
    this.maximumIterations = maximumIterations;
  }

  /**
   * Gets the maximum iterations.
   *
   * @return the maximum iterations
   */
  public int getMaximumIterations() {
    return maximumIterations;
  }

  /**
   * Use this to pass in a recommended minimum number of iterations to the stopping criteria.
   *
   * <p>Note: Implementing classes may ignore this parameter
   *
   * @param minimumIterations the minimumIterations to set
   */
  public void setMinimumIterations(int minimumIterations) {
    this.minimumIterations = minimumIterations;
  }

  /**
   * Gets the minimum iterations.
   *
   * @return the minimum iterations
   */
  public int getMinimumIterations() {
    return minimumIterations;
  }

  /**
   * Sets the log. Used to output fit evaluations for each iteration.
   *
   * @param log the log to set
   */
  public void setLog(Logger log) {
    this.log = log;
  }

  /**
   * Gets the logger.
   *
   * @return the logger.
   */
  public Logger getLog() {
    return log;
  }
}
