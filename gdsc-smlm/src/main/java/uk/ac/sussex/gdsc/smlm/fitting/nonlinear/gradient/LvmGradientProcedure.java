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

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient;

import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Procedure;
import uk.ac.sussex.gdsc.smlm.function.ValueProcedure;

/**
 * Calculates the scaled Hessian matrix (the square matrix of second-order partial derivatives of a
 * function) and the scaled gradient vector of the function's partial first derivatives with respect
 * to the parameters. This is used within the Levenberg-Marquardt method to fit a nonlinear model
 * with coefficients (a) for a set of data points (x, y).
 *
 * <p>Note that the Hessian matrix is scaled by 1/2 and the gradient vector is scaled by -1/2 for
 * convenience in solving the non-linear model. See Numerical Recipes in C++, 2nd Ed. Equation
 * 15.5.8 for Nonlinear Models.
 */
public abstract class LvmGradientProcedure implements Gradient1Procedure, ValueProcedure {
  /** The y. */
  protected final double[] y;

  /** The function. */
  protected final Gradient1Function func;

  /**
   * The number of gradients.
   */
  public final int numberOfGradients;

  /**
   * The scaled gradient vector of the function's partial first derivatives with respect to the
   * parameters (size n).
   */
  public final double[] beta;
  /**
   * The value for the fit.
   */
  public double value;

  /** The y index counter. */
  protected int yi;

  /**
   * Instantiates a new procedure.
   *
   * @param y Data to fit
   * @param func Gradient function
   */
  public LvmGradientProcedure(final double[] y, final Gradient1Function func) {
    this.y = y;
    this.func = func;
    this.numberOfGradients = func.getNumberOfGradients();
    beta = new double[numberOfGradients];
  }

  /**
   * Instantiates a new procedure.
   *
   * @param y Data to fit
   * @param baseline Baseline pre-computed y-values
   * @param func Gradient function
   */
  public LvmGradientProcedure(final double[] y, final double[] baseline,
      final Gradient1Function func) {
    // Process the baseline. This could be part of the function that has been pre-evaluated
    if (baseline != null && baseline.length == y.length) {
      // For sum-of-squares we can just remove the baseline from the y-values
      this.y = new double[y.length];
      final int n = baseline.length;
      for (int i = 0; i < n; i++) {
        this.y[i] = y[i] - baseline[i];
      }
    } else {
      this.y = y;
    }
    this.func = func;
    this.numberOfGradients = func.getNumberOfGradients();
    beta = new double[numberOfGradients];
  }

  /**
   * Evaluate the function and compute the sum-of-squares and the curvature matrix.
   *
   * <p>A call to {@link #isNaNGradients()} will indicate if the gradients were invalid.
   *
   * @param a Set of coefficients for the function
   */
  public void gradient(final double[] a) {
    value = 0;
    yi = -1;
    initialiseGradient();
    func.initialise1(a);
    func.forEach((Gradient1Procedure) this);
    finishGradient();
  }

  /**
   * Initialise for the computation using first order gradients.
   */
  protected abstract void initialiseGradient();

  /**
   * Check gradients for NaN values.
   *
   * @return True if the current gradients contain NaN values
   */
  protected abstract boolean checkGradients();

  /**
   * Finish the computation using first order gradients.
   */
  protected abstract void finishGradient();

  /**
   * Evaluate the function and compute the sum-of-squares.
   *
   * @param a Set of coefficients for the function
   */
  public void value(final double[] a) {
    value = 0;
    yi = -1;
    initialiseValue();
    func.initialise0(a);
    func.forEach((ValueProcedure) this);
    finishValue();
  }

  /**
   * Initialise for the computation of the value.
   */
  protected abstract void initialiseValue();

  /**
   * Finish the computation of the value.
   */
  protected abstract void finishValue();

  /**
   * Get the scaled Hessian curvature matrix (size n*n).
   *
   * @return the alpha
   */
  public double[][] getAlphaMatrix() {
    final double[][] a = new double[numberOfGradients][numberOfGradients];
    getAlphaMatrix(a);
    return a;
  }

  /**
   * Get the scaled Hessian curvature matrix (size n*n) into the provided storage.
   *
   * @param alpha the alpha
   */
  public abstract void getAlphaMatrix(double[][] alpha);

  /**
   * Get the scaled Hessian curvature matrix (size n*n).
   *
   * @return the alpha
   */
  public double[] getAlphaLinear() {
    final double[] a = new double[numberOfGradients * numberOfGradients];
    getAlphaLinear(a);
    return a;
  }

  /**
   * Get the scaled Hessian curvature matrix (size n*n) into the provided storage.
   *
   * @param alpha the alpha
   */
  public abstract void getAlphaLinear(double[] alpha);

  /**
   * Get the scaled gradient vector (size n) into the provided storage.
   *
   * @param beta the beta
   */
  public void getBeta(double[] beta) {
    System.arraycopy(this.beta, 0, beta, 0, numberOfGradients);
  }

  /**
   * Checks if is na N gradients.
   *
   * @return True if the last calculation produced gradients with NaN values.
   */
  public boolean isNaNGradients() {
    return checkGradients();
  }

  /**
   * Convert linear data to a row/column format.
   *
   * @param data the data
   * @return the row/column format
   */
  protected double[][] toMatrix(double[] data) {
    return toMatrix(data, new double[numberOfGradients][numberOfGradients]);
  }

  /**
   * Convert linear data to a row/column format.
   *
   * @param data the data
   * @param out the row/column format
   * @return the row/column format
   */
  protected double[][] toMatrix(double[] data, double[][] out) {
    for (int i = 0, pos = 0; i < numberOfGradients; i++, pos += numberOfGradients) {
      System.arraycopy(data, pos, out[i], 0, numberOfGradients);
    }
    return out;
  }

  /**
   * Convert a the row/column format matrix to a linear data.
   *
   * @param data the data
   * @return the linear data
   */
  protected double[] toLinear(double[][] data) {
    return toLinear(data, new double[numberOfGradients * numberOfGradients]);
  }

  /**
   * Convert a the row/column format matrix to a linear data.
   *
   * @param data the data
   * @param out the linear data
   * @return the linear data
   */
  protected double[] toLinear(double[][] data, double[] out) {
    for (int i = 0, pos = 0; i < numberOfGradients; i++, pos += numberOfGradients) {
      System.arraycopy(data[i], 0, out, pos, numberOfGradients);
    }
    return out;
  }
}
