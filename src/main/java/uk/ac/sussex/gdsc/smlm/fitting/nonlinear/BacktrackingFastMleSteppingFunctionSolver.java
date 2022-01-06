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

import uk.ac.sussex.gdsc.core.utils.SortUtils;
import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.PoissonGradientProcedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient2Function;

/**
 * Uses the Fast MLE method to fit a gradient function with coefficients (a).
 *
 * <p>Calculates the Newton-Raphson update vector for a Poisson process using the first and second
 * partial derivatives.
 *
 * <p>Ref: Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically
 * minimum uncertainty. Nature Methods 7, 373-375 (supplementary note), Eq. 12.
 *
 * <p>Ref: Huang et al, (2015). Video-rate nanoscopy using sCMOS camera–specific single-molecule
 * localization algorithms. Nature Methods 10, 653–658.
 */
public class BacktrackingFastMleSteppingFunctionSolver extends FastMleSteppingFunctionSolver {
  /**
   * The minimum value between two doubles. ISO standard is 2^-52 = 2.220446049e-16. This computes
   * 2.220446049250313E-16.
   */
  private static double epsilon = Math.ulp(1.0);
  /** Line search constant. */
  private static final double ALF = 1.0e-4;

  private final LineStepSearch lineSearch = new LineStepSearch();
  /** Maximum step length used in line search. */
  private double[] maximumStepLength;
  private double maximumStepSize;

  /**
   * Create a new stepping function solver.
   *
   * @param function the function
   * @param maxRelativeError the max relative error
   * @param maxAbsoluteError the max absolute error
   * @throws NullPointerException if the function is null
   */
  public BacktrackingFastMleSteppingFunctionSolver(Gradient2Function function,
      double maxRelativeError, double maxAbsoluteError) {
    super(function, maxRelativeError, maxAbsoluteError);
  }

  /**
   * Create a new stepping function solver.
   *
   * @param function the function
   * @param tc the tolerance checker
   * @param bounds the bounds
   * @throws NullPointerException if the function or tolerance checker is null
   */
  public BacktrackingFastMleSteppingFunctionSolver(Gradient2Function function, ToleranceChecker tc,
      ParameterBounds bounds) {
    super(function, tc, bounds);
  }

  // Extend the Newton-Raphson method by implementing Line Search and Backtracking
  // (see Numerical Recipes in C++, 2nd Ed, page 388-389, function lnsrch).
  // This can still be done in the context of the stepping function solver.
  // Adjustments:
  // - Always ensure we compute the function value. This is needed to determine
  // if the step computed a better value.
  // - The first call to computeFitValue() computes derivatives. The value is obtained from
  // computePseudoLogLikelihood()
  // - All subsequent calls to computeFitValue() only call computeValue() and implement
  // backtracking if the value does not improve with the full Newton step.
  // - When a suitable value is achieved then this should be returned from computeFitValue()
  // - The next call to compute step must evaluate the derivatives.
  // - The function evaluations counter should be appropriately incremented

  @Override
  protected double[] prepareFitValue(double[] y, double[] a) {
    y = super.prepareFitValue(y, a);
    // We always compute the pseudolikelihood
    isPseudoLogLikelihood = true;

    // Configure maximum step length for each dimension using the bounds.
    // This is a simple check that can prevent wild Newton Raphson steps
    // that require a lot of backtracking.
    maximumStepLength = null;
    if (maximumStepSize > 0) {
      final double[] lower = bounds.getLower();
      final double[] upper = bounds.getUpper();
      if (lower != null && upper != null) {
        maximumStepLength = new double[lower.length];
        for (int i = 0; i < maximumStepLength.length; i++) {
          maximumStepLength[i] = (upper[i] - lower[i]) * maximumStepSize;
          if (maximumStepLength[i] <= 0) {
            maximumStepLength[i] = Double.POSITIVE_INFINITY;
          }
        }
      }
    }

    return y;
  }

  @Override
  protected double computeFitValue(double[] a) {
    if (firstEvaluation) {
      firstEvaluation = false;

      // The first call to computeFitValue() computes derivatives. The value is obtained from
      // computePseudoLogLikelihood()
      computeGradients(a);

      evaluations++;
      ll = gradientProcedure.computePseudoLogLikelihood();

      searchDirection = new double[a.length];
      oldA = a.clone();

      return ll;
    }

    // All subsequent calls to computeFitValue() only evaluate the value and implement
    // backtracking if the value does not improve with the full Newton step.
    for (int i = 0; i < searchDirection.length; i++) {
      // Configure the search direction with the full Newton step
      searchDirection[i] = a[i] - oldA[i];
    }

    oldA = lineSearch.lineSearch(oldA, ll, gradientProcedure.d1, searchDirection);
    ll = lineSearch.functionValue;

    // Update the parameters to reflect any backtracking
    System.arraycopy(oldA, 0, a, 0, a.length);

    return ll;
  }

  @Override
  protected void computeStep(double[] step) {
    if (tc.getIterations() > 0) {
      // After backtracking we must compute the derivatives.
      // Note we leave it to here (and not after the line search) so it
      // can be skipped if convergence is achieved.
      computeGradients(oldA);
    }
    super.computeStep(step);
  }

  /**
   * Gets the maximum step size.
   *
   * @return the maximum step size
   */
  public double getMaximumStepSize() {
    return maximumStepSize;
  }

  /**
   * Sets the maximum step size. This is expressed as a factor of the upper-lower bound. Steps
   * beyond this will cause the step to be truncated. Ignored if above 1. Set to below 1 to disable.
   *
   * <p>Note: Restriction of the steps is better controlled using the ParaneterBounds object.
   *
   * @param maximumStepSize the new maximum step size
   */
  public void setMaximumStepSize(double maximumStepSize) {
    this.maximumStepSize = (maximumStepSize >= 1) ? 0 : maximumStepSize;
  }

  @Override
  protected void initialiseAndRun(PoissonGradientProcedure procedure) {
    // Note: This overrides the FastMLESteppingFunctionSolver because the function
    // may not currently be initialised for gradients.
    // If backtracking was done then initialisation has only been done for the value.
    // Reinitialise using the current optimum.

    procedure.computeFisherInformation(oldA); // Use current optimum
  }

  /**
   * Internal class for a line search with backtracking.
   *
   * <p>Adapted from NR::lnsrch, as discussed in Numerical Recipes section 9.7. The algorithm has
   * been changed to find the maximum, support limits on the search direction in all dimensions and
   * check for bad function evaluations when backtracking.
   */
  private class LineStepSearch {
    /**
     * The function value at the new point.
     */
    double functionValue;

    /**
     * Given an n-dimension point, the function value and gradient at that point find a new point
     * along the given search direction so that the function value has decreased sufficiently.
     *
     * @param oldX The old point
     * @param oldF The old point function value
     * @param gradient The old point function gradient (only for the gradient indices)
     * @param searchDirection The search direction
     * @return The new point
     * @throws FunctionSolverException if the slope of the line search is positive
     */
    double[] lineSearch(double[] oldX, final double oldF, double[] gradient,
        double[] searchDirection) {
      double alam2 = 0.0;
      double f2 = 0.0;

      // New point
      final double[] x = new double[oldX.length];

      final int n = oldX.length;

      // Limit the search step size for each dimension
      if (maximumStepLength != null) {
        double scale = 1;
        for (int i = 0; i < n; i++) {
          if (Math.abs(searchDirection[i]) * scale > maximumStepLength[i]) {
            scale = maximumStepLength[i] / Math.abs(searchDirection[i]);
          }
        }
        if (scale < 1) {
          // Scale the entire search direction
          for (int i = 0; i < n; i++) {
            searchDirection[i] *= scale;
          }
        }
      }

      double slope = 0.0;
      final int[] gradientIndices =
          BacktrackingFastMleSteppingFunctionSolver.this.function.gradientIndices();
      for (int i = 0; i < gradient.length; i++) {
        slope += gradient[i] * searchDirection[gradientIndices[i]];
      }
      if (slope <= 0.0) {
        // The search direction for the NR step is the (first derivative / -second derivative).
        // If there are sign errors in the second derivative (it should be the same sign as the
        // first derivative) then the step will be in the 'wrong' direction.
        // Handle this with different options:
        switch (lineSearchMethod) {
          case NONE:
            throw new FunctionSolverException(FitStatus.LINE_SEARCH_ERROR,
                "Slope is negative: " + slope);

          case IGNORE:
            // Ignore any search direction that is in the opposite direction to the
            // first derivative gradient.
            slope = 0.0;
            for (int i = 0; i < gradient.length; i++) {
              final double slopeComponent = gradient[i] * searchDirection[gradientIndices[i]];
              if (slopeComponent < 0) {
                searchDirection[gradientIndices[i]] = 0;
              } else {
                slope += slopeComponent;
              }
            }
            if (slope == 0) {
              return setInsignificantStep(oldX, oldF);
            }
            break;

          case PARTIAL_IGNORE:
            // Progressively ignore any search direction that is in the opposite direction to
            // the first derivative gradient. Do this in order of the magnitude of the error
            final double[] slopeComponents = new double[gradient.length];
            final int[] indices = new int[slopeComponents.length];
            for (int i = 0; i < slopeComponents.length; i++) {
              slopeComponents[i] = gradient[i] * searchDirection[gradientIndices[i]];
              indices[i] = i;
            }
            SortUtils.sortIndices(indices, slopeComponents, false);
            int comp = 0;
            while (slope <= 0 && comp < slopeComponents.length
                && slopeComponents[indices[comp]] <= 0) {
              final int i = indices[comp];
              // Ignore this component
              searchDirection[gradientIndices[i]] = 0;
              comp++;
              // Recompute slope
              slope = 0;
              for (int k = comp; k < slopeComponents.length; k++) {
                slope += slopeComponents[indices[k]];
              }
            }
            if (slope == 0) {
              // All components have been removed so handle no slope
              return setInsignificantStep(oldX, oldF);
            }
            break;

          default:
            throw new IllegalStateException("Unknown line search method: " + lineSearchMethod);
        }
      }

      // Compute lambda min
      double test = 0.0;
      for (int i = 0; i < n; i++) {
        final double temp = Math.abs(searchDirection[i]) / max(Math.abs(oldX[i]), 1.0);
        if (temp > test) {
          test = temp;
        }
      }
      final double alamin = epsilon / test;

      // Always try the full step first
      double alam = 1.0;
      // Count the number of backtracking steps
      int backtracking = 0;
      for (;;) {
        if (alam < alamin) {
          // Convergence (insignificant step).
          return setInsignificantStep(oldX, oldF);
        }

        for (int i = 0; i < n; i++) {
          x[i] = oldX[i] + alam * searchDirection[i];
        }
        // Compute the pseudoLikelihood
        evaluations++;
        gradientProcedure.computeValue(x);
        functionValue = gradientProcedure.computePseudoLogLikelihood();
        if (functionValue >= oldF + ALF * alam * slope) {
          // Sufficient function decrease
          return x;
        }

        // Check for bad function evaluation
        if (!Double.isFinite(functionValue)) {
          // Reset backtracking
          backtracking = 0;

          alam *= 0.1;
          continue;
        }

        // Backtrack
        double tmplam;
        if (backtracking++ == 0) {
          // First backtrack iteration
          tmplam = -slope / (2.0 * (functionValue - oldF - slope));
          // Ensure the lambda is reduced, i.e. we take a step smaller than last time
          if (tmplam > 0.9 * alam) {
            tmplam = 0.9 * alam;
          }
        } else {
          // Subsequent backtracks
          final double rhs1 = functionValue - oldF - alam * slope;
          final double rhs2 = f2 - oldF - alam2 * slope;
          final double a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
          final double b =
              (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2);
          if (a == 0.0) {
            tmplam = -slope / (2.0 * b);
          } else {
            final double disc = b * b - 3.0 * a * slope;
            if (disc < 0.0) {
              tmplam = 0.5 * alam;
            } else if (b <= 0.0) {
              tmplam = (-b + Math.sqrt(disc)) / (3.0 * a);
            } else {
              tmplam = -slope / (b + Math.sqrt(disc));
            }
          }
          // Ensure the lambda is <= 0.5 lambda1, i.e. we take a step smaller than last time
          if (tmplam > 0.5 * alam) {
            tmplam = 0.5 * alam;
          }
        }

        alam2 = alam;
        f2 = functionValue;
        // Ensure the lambda is >= 0.1 lambda1, i.e. we take reasonable step
        alam = max(tmplam, 0.1 * alam);
      }
    }

    private double[] setInsignificantStep(double[] oldX, final double oldF) {
      // Since we use the old function and x then we do not need to compute the objective value
      tc.setConverged();
      functionValue = oldF;
      return oldX;
    }

    private double max(double v1, double v2) {
      return (v1 > v2) ? v1 : v2;
    }
  }
}
