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

import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.smlm.fitting.FunctionSolverType;
import uk.ac.sussex.gdsc.smlm.fitting.LseFunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.linear.EjmlLinearSolver;
import uk.ac.sussex.gdsc.smlm.function.GradientFunction;

/**
 * Abstract class with utility methods for the {@link LseFunctionSolver} interface.
 */
public abstract class LseBaseFunctionSolver extends BaseFunctionSolver
    implements LseFunctionSolver {
  /** The total sum of squares. */
  protected double totalSumOfSquares = Double.NaN;

  /**
   * Default constructor.
   *
   * @param function the function
   * @throws NullPointerException if the function is null
   */
  public LseBaseFunctionSolver(GradientFunction function) {
    super(FunctionSolverType.LSE, function);
  }

  @Override
  protected void preProcess() {
    totalSumOfSquares = Double.NaN;
  }

  /**
   * Gets the total sum of squares.
   *
   * @param y the y
   * @return the total sum of squares
   */
  public static double computeTotalSumOfSquares(double[] y) {
    double sx = 0;
    double ssx = 0;
    for (int i = y.length; i-- > 0;) {
      sx += y[i];
      ssx += y[i] * y[i];
    }
    return ssx - (sx * sx) / (y.length);
  }

  /**
   * Compute the error.
   *
   * @param value the value
   * @param noise the noise
   * @param numberOfFittedPoints the number of fitted points
   * @param numberOfFittedParameters the number of fitted parameters
   * @return the error
   */
  public static double getError(double value, double noise, int numberOfFittedPoints,
      int numberOfFittedParameters) {
    double error = value;

    // Divide by the uncertainty in the individual measurements yi to get the chi-squared
    if (noise > 0) {
      error /= numberOfFittedPoints * noise * noise;
    }

    // This updates the chi-squared value to the average error for a single fitted
    // point using the degrees of freedom (N-M)?
    // Note: This matches the mean squared error output from the MatLab fitting code.
    // If a noise estimate was provided for individual measurements then this will be the
    // reduced chi-square (see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2892436/)
    if (numberOfFittedPoints > numberOfFittedParameters) {
      error /= (numberOfFittedPoints - numberOfFittedParameters);
    } else {
      error = 0;
    }

    return error;
  }

  @Override
  public double getTotalSumOfSquares() {
    if (Double.isNaN(totalSumOfSquares) && lastY != null) {
      totalSumOfSquares = computeTotalSumOfSquares(lastY);
    }
    return totalSumOfSquares;
  }

  @Override
  public double getResidualSumOfSquares() {
    return value;
  }

  @Override
  public double getCoefficientOfDetermination() {
    return 1.0 - (value / getTotalSumOfSquares());
  }

  @Override
  public double getAdjustedCoefficientOfDetermination() {
    return MathUtils.getAdjustedCoefficientOfDetermination(getResidualSumOfSquares(),
        getTotalSumOfSquares(), getNumberOfFittedPoints(), getNumberOfFittedParameters());
  }

  @Override
  public double getMeanSquaredError() {
    return getResidualSumOfSquares() / (getNumberOfFittedPoints() - getNumberOfFittedParameters());
  }

  // Allow I and E as they have special meaning in the formula.
  // CHECKSTYLE.OFF: ParameterName

  /**
   * Compute the covariance matrix of the parameters of the function assuming a least squares fit of
   * a Poisson process.
   *
   * <p>Uses the Mortensen formula (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 25.
   *
   * <p>The method involves inversion of a matrix and may fail.
   *
   * <pre>
   * I = sum_i { Ei,a * Ei,b }
   * E = sum_i { Ei * Ei,a * Ei,b }
   * with
   * i the number of data points fit using least squares using a function of n variable parameters
   * Ei the expected value of the function at i
   * Ei,a the gradient the function at i with respect to parameter a
   * Ei,b the gradient the function at i with respect to parameter b
   * </pre>
   *
   * @param I the Iab matrix
   * @param E the Ei_Eia_Eib matrix
   * @return the covariance matrix (or null)
   */
  public static double[][] covariance(double[][] I, double[][] E) {
    final int n = I.length;

    // Invert the matrix
    final EjmlLinearSolver solver = EjmlLinearSolver.createForInversion(1e-2);
    if (!solver.invert(I)) {
      return null;
    }

    // Note that I now refers to I^-1 in the Mortensen notation

    final double[][] covar = new double[n][n];
    for (int a = 0; a < n; a++) {
      for (int b = 0; b < n; b++) {
        double var = 0;
        for (int ap = 0; ap < n; ap++) {
          for (int bp = 0; bp < n; bp++) {
            var += I[a][ap] * E[ap][bp] * I[bp][b];
          }
        }
        covar[a][b] = var;
      }
    }

    return covar;
  }

  /**
   * Compute the variance of the parameters of the function assuming a least squares fit of a
   * Poisson process.
   *
   * <p>Uses the Mortensen formula (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 25.
   *
   * <p>The method involves inversion of a matrix and may fail.
   *
   * <pre>
   * I = sum_i { Ei,a * Ei,b }
   * E = sum_i { Ei * Ei,a * Ei,b }
   * with
   * i the number of data points fit using least squares using a function of n variable parameters
   * Ei the expected value of the function at i
   * Ei,a the gradient the function at i with respect to parameter a
   * Ei,b the gradient the function at i with respect to parameter b
   * </pre>
   *
   * @param I the Iab matrix
   * @param E the Ei_Eia_Eib matrix
   * @return the variance (or null)
   */
  public static double[] variance(double[][] I, double[][] E) {
    final int n = I.length;

    // Invert the matrix
    final EjmlLinearSolver solver = EjmlLinearSolver.createForInversion(1e-2);
    if (!solver.invert(I)) {
      return null;
    }

    // Note that I now refers to I^-1 in the Mortensen notation

    final double[] covar = new double[n];
    for (int a = 0; a < n; a++) {
      // Note: b==a as we only do the diagonal
      double var = 0;
      for (int ap = 0; ap < n; ap++) {
        for (int bp = 0; bp < n; bp++) {
          var += I[a][ap] * E[ap][bp] * I[bp][a];
        }
      }
      covar[a] = var;
    }

    return covar;
  }
}
