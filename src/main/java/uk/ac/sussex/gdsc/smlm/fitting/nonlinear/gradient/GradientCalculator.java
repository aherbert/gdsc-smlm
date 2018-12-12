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

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient;

import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.LSEBaseFunctionSolver;
import uk.ac.sussex.gdsc.smlm.function.NonLinearFunction;

/**
 * Calculates the Hessian matrix (the square matrix of second-order partial derivatives of a
 * function) and the scaled gradient vector of the function's partial first derivatives with respect
 * to the parameters. This is used within the Levenberg-Marquardt method to fit a nonlinear model
 * with coefficients (a) for a set of data points (x, y).
 *
 * <p>Note that the Hessian matrix is scaled by 1/2 and the gradient vector is scaled by -1/2 for
 * convenience in solving the non-linear model. See Numerical Recipes in C++, 2nd Ed. Equation
 * 15.5.8 for Nonlinear Models.
 */
public class GradientCalculator {
  /** The number of params. */
  public final int nparams;

  /** The bad gradients. */
  private boolean badGradients;

  /**
   * Instantiates a new gradient calculator.
   *
   * @param nparams The number of gradient parameters
   */
  public GradientCalculator(final int nparams) {
    this.nparams = nparams;
  }

  /**
   * Evaluate the function and compute the sum-of-squares and the curvature matrix.
   *
   * <p>A call to {@link #isNaNGradients()} will indicate if the gradients were invalid.
   *
   * @param x n observations
   * @param y Data to fit
   * @param a Set of m coefficients
   * @param alpha the scaled Hessian curvature matrix (size m*m)
   * @param beta the scaled gradient vector of the function's partial first derivatives with respect
   *        to the parameters (size m)
   * @param func Non-linear fitting function
   * @return The sum-of-squares value for the fit
   */
  public double findLinearised(final int[] x, final double[] y, final double[] a,
      final double[][] alpha, final double[] beta, final NonLinearFunction func) {
    double ssx = 0;
    final double[] dy_da = new double[a.length];

    for (int i = 0; i < nparams; i++) {
      beta[i] = 0;
      for (int j = 0; j <= i; j++) {
        alpha[i][j] = 0;
      }
    }

    func.initialise(a);

    if (func.canComputeWeights()) {
      final double[] w = new double[1];
      for (int i = 0; i < x.length; i++) {
        final double dy = y[i] - func.eval(x[i], dy_da, w);
        final double weight = getWeight(w[0]);

        // Compute:
        // - the scaled Hessian matrix (the square matrix of second-order partial derivatives of a
        // function;
        // that is, it describes the local curvature of a function of many variables.)
        // - the scaled gradient vector of the function's partial first derivatives with respect to
        // the parameters

        for (int j = 0; j < nparams; j++) {
          final double wgt = dy_da[j] * weight;

          for (int k = 0; k <= j; k++) {
            alpha[j][k] += wgt * dy_da[k];
          }

          beta[j] += wgt * dy;
        }

        ssx += dy * dy * weight;
      }
    } else {
      for (int i = 0; i < x.length; i++) {
        final double dy = y[i] - func.eval(x[i], dy_da);

        // Compute:
        // - the scaled Hessian matrix (the square matrix of second-order partial derivatives of a
        // function;
        // that is, it describes the local curvature of a function of many variables.)
        // - the scaled gradient vector of the function's partial first derivatives with respect to
        // the parameters

        for (int j = 0; j < nparams; j++) {
          final double wgt = dy_da[j];

          for (int k = 0; k <= j; k++) {
            alpha[j][k] += wgt * dy_da[k];
          }

          beta[j] += wgt * dy;
        }

        ssx += dy * dy;
      }
    }

    // Generate symmetric matrix
    for (int i = 0; i < nparams - 1; i++) {
      for (int j = i + 1; j < nparams; j++) {
        alpha[i][j] = alpha[j][i];
      }
    }

    return checkGradients(alpha, beta, nparams, ssx);
  }

  /**
   * Evaluate the function and compute the sum-of-squares and the curvature matrix.
   *
   * <p>A call to {@link #isNaNGradients()} will indicate if the gradients were invalid.
   *
   * <p>Allows disabling the use of gradients. The output alpha and beta will be reduced in size by
   * the number of indices
   *
   * @param x n observations
   * @param y Data to fit
   * @param a Set of m coefficients
   * @param alpha the scaled Hessian curvature matrix (size m*m)
   * @param beta the scaled gradient vector of the function's partial first derivatives with respect
   *        to the parameters (size m)
   * @param func Non-linear fitting function
   * @param ignore An array of size beta.length. Set the index to true to ignore the gradients.
   * @return The sum-of-squares value for the fit
   */
  public double findLinearised(final int[] x, final double[] y, final double[] a,
      final double[][] alpha, final double[] beta, final NonLinearFunction func, boolean[] ignore) {
    double ssx = 0;
    final double[] dy_da = new double[a.length];

    for (int i = 0; i < nparams; i++) {
      beta[i] = 0;
      for (int j = 0; j <= i; j++) {
        alpha[i][j] = 0;
      }
    }

    func.initialise(a);

    final int[] indices = new int[nparams];
    int nnparams = 0;
    for (int j = 0; j < nparams; j++) {
      if (ignore[j]) {
        continue;
      }
      indices[nnparams++] = j;
    }

    if (func.canComputeWeights()) {
      final double[] w = new double[1];
      for (int i = 0; i < x.length; i++) {
        final double dy = y[i] - func.eval(x[i], dy_da, w);
        final double weight = getWeight(w[0]);

        // Compute:
        // - the scaled Hessian matrix (the square matrix of second-order partial derivatives of a
        // function;
        // that is, it describes the local curvature of a function of many variables.)
        // - the scaled gradient vector of the function's partial first derivatives with respect to
        // the parameters

        for (int j = 0; j < nnparams; j++) {
          final double wgt = dy_da[indices[j]] * weight;

          for (int k = 0; k <= j; k++) {
            alpha[j][k] += wgt * dy_da[indices[k]];
          }

          beta[j] += wgt * dy;
        }

        ssx += dy * dy * weight;
      }
    } else {
      for (int i = 0; i < x.length; i++) {
        final double dy = y[i] - func.eval(x[i], dy_da);

        // Compute:
        // - the scaled Hessian matrix (the square matrix of second-order partial derivatives of a
        // function;
        // that is, it describes the local curvature of a function of many variables.)
        // - the scaled gradient vector of the function's partial first derivatives with respect to
        // the parameters

        for (int j = 0; j < nnparams; j++) {
          final double wgt = dy_da[indices[j]];

          for (int k = 0; k <= j; k++) {
            alpha[j][k] += wgt * dy_da[indices[k]];
          }

          beta[j] += wgt * dy;
        }

        ssx += dy * dy;
      }
    }

    // Generate symmetric matrix
    for (int i = 0; i < nparams - 1; i++) {
      for (int j = i + 1; j < nparams; j++) {
        alpha[i][j] = alpha[j][i];
      }
    }

    return checkGradients(alpha, beta, nparams, ssx);
  }

  /**
   * Evaluate the function and compute the sum-of-squares.
   *
   * @param x n observations
   * @param y The data
   * @param yFit The function data
   * @param a Set of m coefficients
   * @param func Non-linear fitting function
   * @return The sum-of-squares value for the fit
   */
  public double findLinearised(final int[] x, final double[] y, double[] yFit, final double[] a,
      final NonLinearFunction func) {
    double ssx = 0;

    func.initialise(a);

    if (yFit == null || yFit.length < x.length) {
      if (func.canComputeWeights()) {
        final double[] w = new double[1];
        for (int i = 0; i < x.length; i++) {
          final double dy = y[i] - func.evalw(x[i], w);
          final double weight = getWeight(w[0]);
          ssx += dy * dy * weight;
        }
      } else {
        for (int i = 0; i < x.length; i++) {
          final double dy = y[i] - func.eval(x[i]);
          ssx += dy * dy;
        }
      }
    } else if (func.canComputeWeights()) {
      final double[] w = new double[1];
      for (int i = 0; i < x.length; i++) {
        yFit[i] = func.evalw(x[i], w);
        final double dy = y[i] - yFit[i];
        final double weight = getWeight(w[0]);
        ssx += dy * dy * weight;
      }
    } else {
      for (int i = 0; i < x.length; i++) {
        yFit[i] = func.eval(x[i]);
        final double dy = y[i] - yFit[i];
        ssx += dy * dy;
      }
    }

    return ssx;
  }

  /**
   * Evaluate the function and compute the sum-of-squares and the curvature matrix. Assumes the n
   * observations (x) are sequential integers from 0.
   *
   * <p>If the function supports weights then these will be used to compute the SS and curvature
   * matrix.
   *
   * <p>A call to {@link #isNaNGradients()} will indicate if the gradients were invalid.
   *
   * @param n The number of data points
   * @param y Data to fit
   * @param a Set of m coefficients
   * @param alpha the scaled Hessian curvature matrix (size m*m)
   * @param beta the scaled gradient vector of the function's partial first derivatives with respect
   *        to the parameters (size m)
   * @param func Non-linear fitting function
   * @return The sum-of-squares value for the fit.
   * @see NonLinearFunction#eval(int, double[])
   * @see NonLinearFunction#eval(int, double[], double[])
   * @see NonLinearFunction#canComputeWeights()
   */
  public double findLinearised(final int n, final double[] y, final double[] a,
      final double[][] alpha, final double[] beta, final NonLinearFunction func) {
    double ssx = 0;
    final double[] dy_da = new double[a.length];

    for (int i = 0; i < nparams; i++) {
      beta[i] = 0;
      for (int j = 0; j <= i; j++) {
        alpha[i][j] = 0;
      }
    }

    func.initialise(a);

    if (func.canComputeWeights()) {
      final double[] w = new double[1];
      for (int i = 0; i < n; i++) {
        final double dy = y[i] - func.eval(i, dy_da, w);
        final double weight = getWeight(w[0]);

        // Compute:
        // - the scaled Hessian matrix (the square matrix of second-order partial derivatives of a
        // function;
        // that is, it describes the local curvature of a function of many variables.)
        // - the scaled gradient vector of the function's partial first derivatives with respect to
        // the parameters

        for (int j = 0; j < nparams; j++) {
          final double wgt = dy_da[j] * weight;

          for (int k = 0; k <= j; k++) {
            alpha[j][k] += wgt * dy_da[k];
          }

          beta[j] += wgt * dy;
        }

        ssx += dy * dy * weight;
      }
    } else {
      for (int i = 0; i < n; i++) {
        final double dy = y[i] - func.eval(i, dy_da);

        // Compute:
        // - the scaled Hessian matrix (the square matrix of second-order partial derivatives of a
        // function;
        // that is, it describes the local curvature of a function of many variables.)
        // - the scaled gradient vector of the function's partial first derivatives with respect to
        // the parameters

        for (int j = 0; j < nparams; j++) {
          final double wgt = dy_da[j];

          for (int k = 0; k <= j; k++) {
            alpha[j][k] += wgt * dy_da[k];
          }

          beta[j] += wgt * dy;
        }

        ssx += dy * dy;
      }
    }

    // Generate symmetric matrix
    for (int i = 0; i < nparams - 1; i++) {
      for (int j = i + 1; j < nparams; j++) {
        alpha[i][j] = alpha[j][i];
      }
    }

    return checkGradients(alpha, beta, nparams, ssx);
  }

  /**
   * Evaluate the function and compute the sum-of-squares and the curvature matrix. Assumes the n
   * observations (x) are sequential integers from 0.
   *
   * <p>If the function supports weights then these will be used to compute the SS and curvature
   * matrix.
   *
   * <p>A call to {@link #isNaNGradients()} will indicate if the gradients were invalid.
   *
   * <p>Allows disabling the use of gradients. The output alpha and beta will be reduced in size by
   * the number of indices that were disabled. The remaining positions will be zero filled.
   *
   * @param n The number of data points
   * @param y Data to fit
   * @param a Set of m coefficients
   * @param alpha the scaled Hessian curvature matrix (size m*m)
   * @param beta the scaled gradient vector of the function's partial first derivatives with respect
   *        to the parameters (size m)
   * @param func Non-linear fitting function
   * @param ignore An array of size beta.length. Set the index to true to ignore the gradients.
   * @return The sum-of-squares value for the fit.
   * @see NonLinearFunction#eval(int, double[])
   * @see NonLinearFunction#eval(int, double[], double[])
   * @see NonLinearFunction#canComputeWeights()
   */
  public double findLinearised(final int n, final double[] y, final double[] a,
      final double[][] alpha, final double[] beta, final NonLinearFunction func, boolean[] ignore) {
    double ssx = 0;
    final double[] dy_da = new double[a.length];

    for (int i = 0; i < nparams; i++) {
      beta[i] = 0;
      for (int j = 0; j <= i; j++) {
        alpha[i][j] = 0;
      }
    }

    func.initialise(a);

    final int[] indices = new int[nparams];
    int nnparams = 0;
    for (int j = 0; j < nparams; j++) {
      if (ignore[j]) {
        continue;
      }
      indices[nnparams++] = j;
    }

    if (func.canComputeWeights()) {
      final double[] w = new double[1];
      for (int i = 0; i < n; i++) {
        final double dy = y[i] - func.eval(i, dy_da, w);
        final double weight = getWeight(w[0]);

        // Compute:
        // - the scaled Hessian matrix (the square matrix of second-order partial derivatives of a
        // function;
        // that is, it describes the local curvature of a function of many variables.)
        // - the scaled gradient vector of the function's partial first derivatives with respect to
        // the parameters

        for (int j = 0; j < nnparams; j++) {
          final double wgt = dy_da[indices[j]] * weight;

          for (int k = 0; k <= j; k++) {
            alpha[j][k] += wgt * dy_da[indices[k]];
          }

          beta[j] += wgt * dy;
        }

        ssx += dy * dy * weight;
      }
    } else {
      for (int i = 0; i < n; i++) {
        final double dy = y[i] - func.eval(i, dy_da);

        // Compute:
        // - the scaled Hessian matrix (the square matrix of second-order partial derivatives of a
        // function. that is, it describes the local curvature of a function of many variables.)
        // - the scaled gradient vector of the function's partial first derivatives with respect to
        // the parameters

        for (int j = 0; j < nnparams; j++) {
          final double wgt = dy_da[indices[j]];

          for (int k = 0; k <= j; k++) {
            alpha[j][k] += wgt * dy_da[indices[k]];
          }

          beta[j] += wgt * dy;
        }

        ssx += dy * dy;
      }
    }

    // Generate symmetric matrix
    for (int i = 0; i < nparams - 1; i++) {
      for (int j = i + 1; j < nparams; j++) {
        alpha[i][j] = alpha[j][i];
      }
    }

    return checkGradients(alpha, beta, nparams, ssx);
  }

  /**
   * Evaluate the function and compute the sum-of-squares.
   *
   * @param n The number of data points
   * @param y The data
   * @param yFit The function data
   * @param a Set of m coefficients
   * @param func Non-linear fitting function
   * @return The sum-of-squares value for the fit
   */
  public double findLinearised(final int n, final double[] y, double[] yFit, final double[] a,
      final NonLinearFunction func) {
    double ssx = 0;

    func.initialise(a);

    if (yFit == null || yFit.length < n) {
      if (func.canComputeWeights()) {
        final double[] w = new double[1];
        for (int i = 0; i < n; i++) {
          final double dy = y[i] - func.evalw(i, w);
          final double weight = getWeight(w[0]);
          ssx += dy * dy * weight;
        }
      } else {
        for (int i = 0; i < n; i++) {
          final double dy = y[i] - func.eval(i);
          ssx += dy * dy;
        }
      }
    } else if (func.canComputeWeights()) {
      final double[] w = new double[1];
      for (int i = 0; i < n; i++) {
        yFit[i] = func.evalw(i, w);
        final double dy = y[i] - yFit[i];
        final double weight = getWeight(w[0]);
        ssx += dy * dy * weight;
      }
    } else {
      for (int i = 0; i < n; i++) {
        yFit[i] = func.eval(i);
        final double dy = y[i] - yFit[i];
        ssx += dy * dy;
      }
    }

    return ssx;
  }

  /**
   * Get the weight factor using the computed weight.
   *
   * <p>Check if the weight is below 1 and set to 1 to avoid excessive weights.
   *
   * @param weight The computed weight
   * @return The weight factor
   */
  protected double getWeight(final double weight) {
    // TODO - Check if there is a better way to smooth the weights rather than just truncating them
    // at 1
    return (weight < 1) ? 1 : 1.0 / weight;
  }

  /**
   * Check gradients for NaN values.
   *
   * @param alpha the alpha
   * @param beta the beta
   * @param nparams the number of params
   * @param ssx the sum of squares
   * @return the sum of squares
   */
  protected double checkGradients(final double[][] alpha, final double[] beta, int nparams,
      final double ssx) {
    badGradients = checkIsNaN(alpha, beta, nparams);
    return ssx;
  }

  /**
   * Check gradients for NaN values.
   *
   * @param alpha the alpha
   * @param nparams the number of params
   */
  protected void checkGradients(double[][] alpha, int nparams) {
    badGradients = checkIsNaN(alpha, nparams);
  }

  /**
   * Check gradients for NaN values.
   *
   * @param beta the beta
   * @param nparams the number of params
   */
  protected void checkGradients(final double[] beta, final int nparams) {
    badGradients = checkIsNaN(beta, nparams);
  }

  /**
   * Check is na N.
   *
   * @param alpha the alpha
   * @param beta the beta
   * @param nparams the nparams
   * @return true, if successful
   */
  private static boolean checkIsNaN(final double[][] alpha, final double[] beta,
      final int nparams) {
    for (int i = 0; i < nparams; i++) {
      if (Double.isNaN(beta[i])) {
        return true;
      }
      for (int j = 0; j <= i; j++) {
        if (Double.isNaN(alpha[i][j])) {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * Check gradients for NaN values.
   *
   * @param alpha the alpha
   * @param nparams the number of params
   * @return true, if successful
   */
  private static boolean checkIsNaN(final double[][] alpha, final int nparams) {
    for (int i = 0; i < nparams; i++) {
      for (int j = 0; j <= i; j++) {
        if (Double.isNaN(alpha[i][j])) {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * Check gradients for NaN values.
   *
   * @param beta the beta
   * @param nparams the number of params
   * @return true, if successful
   */
  private static boolean checkIsNaN(final double[] beta, final int nparams) {
    for (int i = 0; i < nparams; i++) {
      if (Double.isNaN(beta[i])) {
        return true;
      }
    }
    return false;
  }

  /**
   * Checks if the last calculation produced gradients with NaN values.
   *
   * @return True if the last calculation produced gradients with NaN values.
   */
  public boolean isNaNGradients() {
    return badGradients;
  }

  /**
   * Compute the Fisher's Information Matrix (I) assuming a Poisson process.
   *
   * <pre>
   * Iab = E [ ( d ln(L(x|p)) / da ) * ( d ln(L(x|p)) / db ) ]
   * p = parameters
   * x = observed values
   * L(x|p) = likelihood of X given p
   * E = expected value
   * </pre>
   *
   * <p>Note that this is only a true Fisher information diagonal if the function returns the
   * expected value for a Poisson process. In this case the equation reduces to:
   *
   * <pre>
   * Iaa = sum(i) (dYi da) * (dYi da) / Yi
   * </pre>
   *
   * <p>See Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically
   * minimum uncertainty. Nature Methods 7, 373-375 (supplementary note), Eq. 9.
   *
   * <p>A call to {@link #isNaNGradients()} will indicate if the gradients were invalid.
   *
   * @param x n observations
   * @param a Set of m coefficients (if null then the function must be pre-initialised)
   * @param func Non-linear fitting function. Must return positive values (other values are ignored)
   * @return I
   */
  public double[][] fisherInformationMatrix(int[] x, final double[] a,
      final NonLinearFunction func) {
    return fisherInformationMatrix(x.length, a, func);
  }

  /**
   * Compute the Fisher's Information Matrix (I) assuming a Poisson process.
   *
   * <pre>
   * Iab = E [ ( d ln(L(x|p)) / da ) * ( d ln(L(x|p)) / db ) ]
   * p = parameters
   * x = observed values
   * L(x|p) = likelihood of X given p
   * E = expected value
   * </pre>
   *
   * <p>Note that this is only a true Fisher information diagonal if the function returns the
   * expected value for a Poisson process. In this case the equation reduces to:
   *
   * <pre>
   * Iaa = sum(i) (dYi da) * (dYi da) / Yi
   * </pre>
   *
   * <p>See Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically
   * minimum uncertainty. Nature Methods 7, 373-375 (supplementary note), Eq. 9.
   *
   * <p>A call to {@link #isNaNGradients()} will indicate if the gradients were invalid.
   *
   * @param n The number of data points
   * @param a Set of m coefficients (if null then the function must be pre-initialised)
   * @param func Non-linear fitting function. Must return positive values (other values are ignored)
   * @return I
   */
  public double[][] fisherInformationMatrix(final int n, final double[] a,
      final NonLinearFunction func) {
    final double[] dy_da = new double[nparams];

    final double[][] alpha = new double[nparams][nparams];

    if (a != null) {
      func.initialise(a);
    }

    for (int i = 0; i < n; i++) {
      final double yi = func.eval(i, dy_da);
      if (yi > 0) {
        for (int j = 0; j < nparams; j++) {
          final double dy_db = dy_da[j] / yi;

          for (int k = 0; k <= j; k++) {
            alpha[j][k] += dy_db * dy_da[k];
          }
        }
      }
    }

    // Generate symmetric matrix
    for (int i = 0; i < nparams - 1; i++) {
      for (int j = i + 1; j < nparams; j++) {
        alpha[i][j] = alpha[j][i];
      }
    }

    checkGradients(alpha, nparams);
    return alpha;
  }

  /**
   * Evaluate the function and compute the sum-of-squares and the gradient with respect to the model
   * parameters.
   *
   * <p>A call to {@link #isNaNGradients()} will indicate if the gradients were invalid.
   *
   * @param x n observations
   * @param y Data to fit
   * @param a Set of m coefficients
   * @param df_da the gradient vector of the function's partial first derivatives with respect to
   *        the parameters (size m)
   * @param func Non-linear fitting function
   * @return The sum-of-squares value for the fit
   */
  public double evaluate(final int[] x, final double[] y, final double[] a, final double[] df_da,
      final NonLinearFunction func) {
    double ssx = 0;
    final double[] dy_da = new double[nparams];

    zero(df_da);

    func.initialise(a);

    for (int i = 0; i < x.length; i++) {
      final double dy = y[i] - func.eval(x[i], dy_da);

      // Compute:
      // - the gradient vector of the function's partial first derivatives with respect to the
      // parameters

      for (int j = 0; j < nparams; j++) {
        df_da[j] += dy_da[j] * dy;
      }

      ssx += dy * dy;
    }

    checkGradients(df_da, nparams);

    // Apply a factor of -2 to compute the actual gradients:
    // See Numerical Recipes in C++, 2nd Ed. Equation 15.5.6 for Nonlinear Models
    for (int j = 0; j < nparams; j++) {
      df_da[j] *= -2;
    }

    return ssx;
  }

  /**
   * Zero the working region of the input matrix alpha and vector beta.
   *
   * @param beta the beta
   */
  protected void zero(final double[] beta) {
    for (int i = 0; i < nparams; i++) {
      beta[i] = 0;
    }
  }

  /**
   * Compute the covariance matrix for the parameters of the function assuming a least squares fit
   * of a Poisson process.
   *
   * <p>Uses the Mortensen formula (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 25.
   *
   * <p>The method involves inversion of a matrix and may fail.
   *
   * @param x n observations
   * @param a Set of m coefficients (if null then the function must be pre-initialised)
   * @param func Non-linear fitting function
   * @return the covariance matrix (or null)
   */
  public double[][] covarianceMatrix(int[] x, final double[] a, final NonLinearFunction func) {
    return covarianceMatrix(x.length, a, func);
  }

  /**
   * Compute the covariance matrix for the parameters of the function assuming a least squares fit
   * of a Poisson process.
   *
   * <p>Uses the Mortensen formula (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 25.
   *
   * <p>The method involves inversion of a matrix and may fail.
   *
   * @param n The number of data points
   * @param theta Set of m coefficients (if null then the function must be pre-initialised)
   * @param func Non-linear fitting function
   * @return the covariance matrix (or null)
   */
  public double[][] covarianceMatrix(final int n, final double[] theta,
      final NonLinearFunction func) {
    // Same notation as Mortensen
    final double[] Eix = new double[nparams];

    final double[][] I = new double[nparams][nparams];
    final double[][] Ei_Eia_Eib = new double[nparams][nparams];

    if (theta != null) {
      func.initialise(theta);
    }

    for (int i = 0; i < n; i++) {
      final double Ei = func.eval(i, Eix);
      for (int a = 0; a < nparams; a++) {
        for (int b = 0; b <= a; b++) {
          final double v = Eix[a] * Eix[b];
          I[a][b] += v;
          Ei_Eia_Eib[a][b] += Ei * v;
        }
      }
    }

    // Generate symmetric matrix
    for (int i = 0; i < nparams - 1; i++) {
      for (int j = i + 1; j < nparams; j++) {
        I[i][j] = I[j][i];
        Ei_Eia_Eib[i][j] = I[j][i];
      }
    }

    checkGradients(I, nparams);
    if (isNaNGradients()) {
      return null;
    }

    return LSEBaseFunctionSolver.covariance(I, Ei_Eia_Eib);
  }

  /**
   * Compute the variance of the parameters of the function assuming a least squares fit of a
   * Poisson process.
   *
   * <p>Uses the Mortensen formula (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 25.
   *
   * <p>The method involves inversion of a matrix and may fail.
   *
   * @param x n observations
   * @param a Set of m coefficients (if null then the function must be pre-initialised)
   * @param func Non-linear fitting function
   * @return the variance (or null)
   */
  public double[] variance(int[] x, final double[] a, final NonLinearFunction func) {
    return variance(x.length, a, func);
  }

  /**
   * Compute the variance of the parameters of the function assuming a least squares fit of a
   * Poisson process.
   *
   * <p>Uses the Mortensen formula (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 25.
   *
   * <p>The method involves inversion of a matrix and may fail.
   *
   * @param n The number of data points
   * @param theta Set of m coefficients (if null then the function must be pre-initialised)
   * @param func Non-linear fitting function
   * @return the variance (or null)
   */
  public double[] variance(final int n, final double[] theta, final NonLinearFunction func) {
    // Same notation as Mortensen
    final double[] Eix = new double[nparams];

    final double[][] I = new double[nparams][nparams];
    final double[][] Ei_Eia_Eib = new double[nparams][nparams];

    if (theta != null) {
      func.initialise(theta);
    }

    for (int i = 0; i < n; i++) {
      final double Ei = func.eval(i, Eix);
      for (int a = 0; a < nparams; a++) {
        for (int b = 0; b <= a; b++) {
          final double v = Eix[a] * Eix[b];
          I[a][b] += v;
          Ei_Eia_Eib[a][b] += Ei * v;
        }
      }
    }

    // Generate symmetric matrix
    for (int i = 0; i < nparams - 1; i++) {
      for (int j = i + 1; j < nparams; j++) {
        I[i][j] = I[j][i];
        Ei_Eia_Eib[i][j] = Ei_Eia_Eib[j][i];
      }
    }

    checkGradients(I, nparams);
    if (isNaNGradients()) {
      return null;
    }

    // System.out.println("Gradient calc");
    // System.out.println(new DenseMatrix64F(I).toString());
    // System.out.println(new DenseMatrix64F(Ei_Eia_Eib).toString());

    return LSEBaseFunctionSolver.variance(I, Ei_Eia_Eib);
  }
}
