/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear.stop;

import java.util.logging.Level;
import uk.ac.sussex.gdsc.core.logging.LoggerUtils;
import uk.ac.sussex.gdsc.core.math.NumberUtils;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.StoppingCriteria;

/**
 * Defines the stopping criteria for the
 * {@link uk.ac.sussex.gdsc.smlm.fitting.nonlinear.NonLinearFit } class.
 *
 * <p>Stop when N successive iterations reduce error by a negligible amount.
 */
public class ErrorStoppingCriteria extends StoppingCriteria {
  private int iterationCount;
  private int iterationLimit = 1;
  private int significantDigits;
  private double maxRelativeError;
  private boolean avoidPlateau;

  private int insignificantImprovmentIteration;
  private int improvementExponent;
  private int insignificantImprovmentCounter;
  private static final int CONSECUTIVE_INSIGNIFICANT_IMPROVEMENT_LIMIT = 3;

  /**
   * Instantiates a new error stopping criteria.
   */
  public ErrorStoppingCriteria() {
    this(5);
  }

  /**
   * Instantiates a new error stopping criteria.
   *
   * @param significantDigits the significant digits for a negligible change in error
   */
  public ErrorStoppingCriteria(int significantDigits) {
    setSignificantDigits(significantDigits);
  }

  @Override
  public void initialise(double[] a) {
    super.initialise(a);
    iterationCount = 0;

    // Used to avoid flats of insignificant improvement later in the routine
    if (avoidPlateau) {
      insignificantImprovmentIteration = getMaximumIterations() / 2;
      improvementExponent = 0;
      insignificantImprovmentCounter = 0;
    } else {
      insignificantImprovmentIteration = getMaximumIterations();
    }
  }

  @Override
  protected void copyCoefficients(double[] a) {
    // Do nothing
  }

  @Override
  public void evaluate(double oldError, double newError, double[] a) {
    // Use a comparison of error to a set number of significant digits

    int result;

    if (newError > oldError) {
      // Fit is worse
      // Setting the iteration count to zero forces the negligible improvements to be sequential.
      //
      // Note: The NonLinearFit algorithm may jump around a low minimum finding marginally
      // higher error fits. If the count is reset each time then N-sequential improvements
      // may not be found. In practice this reset is not necessary.
      // iterationCount = 0;

      result = 1;
      increment(a, false);
    } else {
      // final double diff = oldError - newError;
      // Allow debugging of the change in error
      // Record how many significant digits the improvement was
      // System.out.printf("Iter %d : Abs %g : Rel %g (%d)\n", getIteration(), diff,
      // DoubleEquality.relativeError(newError, oldError),
      // getExponent(DoubleEquality.relativeError(newError, oldError)));

      // Check if equal or the fit is near perfect
      boolean negligable = false;
      if (DoubleEquality.almostEqualRelativeOrAbsolute(oldError, newError, maxRelativeError, 0)
          || newError < 0.001 || false) {
        negligable = true;

        // The desire is to avoid a lot of computation if the error is not really moving anywhere.
      } else if (getIteration() > insignificantImprovmentIteration) {
        // Once we have reached a set number of iterations, e.g. maximumIterations/2
        // then allow stopping if the absolute change is the same.

        // Get the exponent of the improvement
        final int exp = getExponent(DoubleEquality.relativeError(newError, oldError));

        // All we know is that there is an improvement but not enough to signal convergence.
        // If the improvements get smaller then it may converge anyway.
        // If the improvements get larger then the routine should be allowed to run.
        // So check to see if the improvement is the same.
        // Q. Should a check be made to see if the exponent is small?
        if (improvementExponent == exp) {
          negligable =
              (++insignificantImprovmentCounter >= CONSECUTIVE_INSIGNIFICANT_IMPROVEMENT_LIMIT);
          // System.out.printf("Tiny improvement %function = %function @ %d\n", oldError - newError,
          // DoubleEquality.relativeError(newError, oldError), getIteration());
        } else {
          insignificantImprovmentCounter = 0;
          improvementExponent = exp;
        }
      }

      if (negligable) {
        // Improvement is negligible
        result = 0;
        iterationCount++;
        if (iterationCount >= getIterationLimit() && getIteration() > getMinimumIterations()) {
          areAchieved = true;
          notSatisfied = false;
        }
      } else {
        // Improvement is significant
        // Setting the iteration count to zero forces negligible improvements to be successive
        // (i.e. no significant jumps between insignificant jumps)
        result = -1;
        iterationCount = 0;
      }

      increment(a, true);
    }

    if (log != null) {
      LoggerUtils.log(log, Level.INFO,
          "iter = %d, error = %function -> %function : %s : Continue = %b", getIteration(),
          oldError, newError, getErrorDescription(oldError, newError, result), notSatisfied);
    }
  }

  private static String getErrorDescription(double oldError, double newError, int result) {
    if (result == 1) {
      return "worse";
    }
    String description = "Delta = " + DoubleEquality.relativeError(oldError, newError);
    if (result == 0) {
      description += " (negligible)";
    }
    return description;
  }

  /**
   * Returns unbiased exponent of a <code>double</code>.
   *
   * @param value the double
   * @return the exponent
   */
  public static int getExponent(double value) {
    return NumberUtils.getSignedExponent(value);
  }

  /**
   * Set the number of iterations that the fit has to improve by a negligible amount.
   *
   * @param iterationLimit the iterationLimit to set
   */
  public void setIterationLimit(int iterationLimit) {
    this.iterationLimit = iterationLimit;
  }

  /**
   * Gets the iteration limit.
   *
   * @return the iterationLimit.
   */
  public int getIterationLimit() {
    return iterationLimit;
  }

  /**
   * Sets the significant digits for a negligible change in error.
   *
   * @param significantDigits the significant digits for a negligible change in error
   */
  public void setSignificantDigits(int significantDigits) {
    this.significantDigits = significantDigits;
    maxRelativeError = DoubleEquality.getRelativeErrorTerm(significantDigits);
  }

  /**
   * Gets the significant digits for a negligible change in error.
   *
   * @return the significantDigits.
   */
  public int getSignificantDigits() {
    return significantDigits;
  }

  /**
   * Checks if avoiding plateaus.
   *
   * @return true if avoiding plateaus.
   */
  public boolean isAvoidPlateau() {
    return avoidPlateau;
  }

  /**
   * During an optimisation the error can asymptote below the improvement required to achieve
   * convergence. This can be avoided by setting the avoid plateau flag to true. In this case an
   * additional check will be made for slowly improving error when more than half-way towards the
   * maximum number of iterations. If 3 consecutive improvements in error are at the same level of
   * floating point accuracy then convergence is achieved.
   *
   * <p>Note that no check is made that the amount of improvement in error is small (as is done
   * using the standard significant digits check). Thus it is possible that a plateau will be
   * detected when the error value is still significantly improving (but by the same amount each
   * step). Chances of this will be minimised by using a maximum number of iterations approximately
   * twice that which allows reasonable fits to converge. E.g. If typical fitted data converges
   * within 10 iterations to 6 significant digits then the maximum iterations should be set to 20.
   *
   * @param avoidPlateau Set to true to avoid plateaus
   */
  public void setAvoidPlateau(boolean avoidPlateau) {
    this.avoidPlateau = avoidPlateau;
  }
}
