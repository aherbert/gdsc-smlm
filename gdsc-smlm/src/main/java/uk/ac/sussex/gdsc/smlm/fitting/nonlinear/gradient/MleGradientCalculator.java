/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

import uk.ac.sussex.gdsc.smlm.function.NonLinearFunction;
import uk.ac.sussex.gdsc.smlm.function.PoissonCalculator;

/**
 * Calculates the Hessian matrix (the square matrix of second-order partial derivatives of a
 * function) and the gradient vector of the function's partial first derivatives with respect to the
 * parameters. This is used within the Levenberg-Marquardt method to fit a nonlinear model with
 * coefficients (a) for a set of data points (x, y).
 *
 * <p>This calculator computes a modified Chi-squared expression to perform Maximum Likelihood
 * Estimation assuming Poisson model. See Laurence &amp; Chromy (2010) Efficient maximum likelihood
 * estimator. Nature Methods 7, 338-339. The input data must be Poisson distributed for this to be
 * relevant.
 */
public class MleGradientCalculator extends GradientCalculator {

  /** The log for min. */
  private static final double LOG_FOR_MIN = Math.log(Double.MIN_VALUE);

  /**
   * Instantiates a new MLE gradient calculator.
   *
   * @param nparams The number of gradient parameters
   */
  public MleGradientCalculator(final int nparams) {
    super(nparams);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Note: if the function returns a negative value then it is set to zero.
   *
   * @param y Data to fit (must be strictly positive Poisson data)
   * @return The MLE chi-squared value
   */
  @Override
  public double findLinearised(final int[] x, final double[] y, final double[] a,
      final double[][] alpha, final double[] beta, final NonLinearFunction func) {
    double chisq = 0;
    final double[] dfiDa = new double[nparams];

    zero(alpha, beta);

    func.initialise(a);

    for (int i = 0; i < x.length; i++) {
      // Function must produce a positive output.
      final double xi = y[i];

      // The code provided in Laurence & Chromy (2010) Nature Methods 7, 338-339, SI
      // effectively ignores any function value below zero. This could lead to a
      // situation where the best chisq value can be achieved by setting the output
      // function to produce 0 for all evaluations. To cope with this we heavily
      // penalise the chisq value.
      // Optimally the function should be bounded to always produce a positive number.
      final double fi = func.eval(i, dfiDa);

      if (fi <= 0) {
        // We assume xi is positive
        if (xi != 0) {
          // Penalise the chi-squared value by assuming fi is a very small positive value
          chisq += (-xi - xi * LOG_FOR_MIN);
        }

        // We ignore this contribution to the gradient for stability
        // compute(alpha, beta, dfiDa, Double.MIN_VALUE, xi);

        // We assume y[i] is positive but must handle zero
      } else if (xi <= 0.0) {
        chisq += fi;
        compute0(beta, dfiDa, fi);
      } else {
        chisq += (fi - xi - xi * Math.log(fi / xi));
        compute(alpha, beta, dfiDa, fi, xi);
      }
    }

    symmetric(alpha);

    // Move the factor of 2 to the end
    return checkGradients(alpha, beta, nparams, chisq * 2);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Note: if the function returns a negative value then it is set to zero.
   *
   * @param y Data to fit (must be strictly positive Poisson data)
   * @return The MLE chi-squared value
   */
  @Override
  public double findLinearised(final int[] x, final double[] y, final double[] a,
      final double[][] alpha, final double[] beta, final NonLinearFunction func, boolean[] ignore) {
    double chisq = 0;
    final double[] dfiDa = new double[nparams];

    zero(alpha, beta);

    func.initialise(a);

    final int[] indices = new int[nparams];
    int nnparams = 0;
    for (int j = 0; j < nparams; j++) {
      if (ignore[j]) {
        continue;
      }
      indices[nnparams++] = j;
    }

    for (int i = 0; i < x.length; i++) {
      // Function must produce a positive output.
      final double xi = y[i];

      // The code provided in Laurence & Chromy (2010) Nature Methods 7, 338-339, SI
      // effectively ignores the function value below zero. This could lead to a
      // situation where the best chisq value can be achieved by setting the output
      // function to produce 0 for all evaluations. To cope with this we heavily
      // penalise the chisq value.
      // Optimally the function should be bounded to always produce a positive number.
      final double fi = func.eval(i, dfiDa);

      if (fi <= 0) {
        // We assume xi is positive
        if (xi != 0) {
          // Penalise the chi-squared value by assuming fi is a very small positive value
          chisq += (-xi - xi * LOG_FOR_MIN);
        }

        // We ignore this contribution to the gradient for stability
        // compute(alpha, beta, dfiDa, Double.MIN_VALUE, xi);
      } else {
        // We assume y[i] is positive but must handle zero
        if (xi <= 0.0) {
          chisq += fi;
        } else {
          chisq += (fi - xi - xi * Math.log(fi / xi));
        }

        compute(alpha, beta, dfiDa, fi, xi, indices, nnparams);
      }
    }

    symmetric(alpha);

    // Move the factor of 2 to the end
    return checkGradients(alpha, beta, nparams, chisq * 2);
  }

  /**
   * {@inheritDoc}
   *
   * @param y Data to fit (must be strictly positive Poisson data)
   * @return The MLE chi-squared value
   */
  @Override
  public double findLinearised(final int[] x, final double[] y, double[] fx, final double[] a,
      final NonLinearFunction func) {
    double chisq = 0;

    func.initialise(a);

    if (fx == null || fx.length < x.length) {
      for (int i = 0; i < x.length; i++) {
        // Function must produce a positive output.
        final double xi = y[i];

        // The code provided in Laurence & Chromy (2010) Nature Methods 7, 338-339, SI
        // effectively ignores the function value below zero. This could lead to a
        // situation where the best chisq value can be achieved by setting the output
        // function to produce 0 for all evaluations. To cope with this we heavily
        // penalise the chisq value.
        // Optimally the function should be bounded to always produce a positive number.
        final double fi = func.eval(i);

        if (fi <= 0) {
          // We assume xi is positive
          if (xi != 0) {
            // Penalise the chi-squared value by assuming fi is a very small positive value
            chisq += (-xi - xi * LOG_FOR_MIN);
          }

          // We assume y[i] is positive but must handle zero
        } else if (xi <= 0.0) {
          chisq += fi;
        } else {
          chisq += (fi - xi - xi * Math.log(fi / xi));
        }
      }
    } else {
      for (int i = 0; i < x.length; i++) {
        // Function must produce a positive output.
        final double xi = y[i];

        // The code provided in Laurence & Chromy (2010) Nature Methods 7, 338-339, SI
        // effectively ignores the function value below zero. This could lead to a
        // situation where the best chisq value can be achieved by setting the output
        // function to produce 0 for all evaluations. To cope with this we heavily
        // penalise the chisq value.
        // Optimally the function should be bounded to always produce a positive number.
        final double fi = func.eval(i);
        fx[i] = fi;

        if (fi <= 0) {
          // We assume xi is positive
          if (xi != 0) {
            // Penalise the chi-squared value by assuming fi is a very small positive value
            chisq += (-xi - xi * LOG_FOR_MIN);
          }

          // We assume y[i] is positive but must handle zero
        } else if (xi <= 0.0) {
          chisq += fi;
        } else {
          chisq += (fi - xi - xi * Math.log(fi / xi));
        }
      }
    }

    // Move the factor of 2 to the end
    return chisq * 2;
  }

  /**
   * {@inheritDoc}
   *
   * @param y Data to fit (must be strictly positive Poisson data)
   * @return The MLE chi-squared value
   */
  @Override
  public double findLinearised(final int n, final double[] y, final double[] a,
      final double[][] alpha, final double[] beta, final NonLinearFunction func) {
    double chisq = 0;
    final double[] dfiDa = new double[nparams];

    zero(alpha, beta);

    func.initialise(a);

    for (int i = 0; i < n; i++) {
      // Function must produce a positive output.
      final double xi = y[i];

      // The code provided in Laurence & Chromy (2010) Nature Methods 7, 338-339, SI
      // effectively ignores the function value below zero. This could lead to a
      // situation where the best chisq value can be achieved by setting the output
      // function to produce 0 for all evaluations. To cope with this we heavily
      // penalise the chisq value.
      // Optimally the function should be bounded to always produce a positive number.
      final double fi = func.eval(i, dfiDa);

      if (fi <= 0) {
        // We assume xi is positive
        if (xi != 0) {
          // Penalise the chi-squared value by assuming fi is a very small positive value
          chisq += (-xi - xi * LOG_FOR_MIN);
        }

        // We ignore this contribution to the gradient for stability
        // compute(alpha, beta, dfiDa, Double.MIN_VALUE, xi);

        // We assume y[i] is positive but must handle zero
      } else if (xi <= 0.0) {
        chisq += fi;
        compute0(beta, dfiDa, fi);
      } else {
        chisq += (fi - xi - xi * Math.log(fi / xi));
        compute(alpha, beta, dfiDa, fi, xi);
      }
    }

    symmetric(alpha);

    // Move the factor of 2 to the end
    return checkGradients(alpha, beta, nparams, chisq * 2);
  }

  /**
   * {@inheritDoc}
   *
   * @param y Data to fit (must be strictly positive Poisson data)
   * @return The MLE chi-squared value
   */
  @Override
  public double findLinearised(final int n, final double[] y, final double[] a,
      final double[][] alpha, final double[] beta, final NonLinearFunction func, boolean[] ignore) {
    double chisq = 0;
    final double[] dfiDa = new double[nparams];

    zero(alpha, beta);

    func.initialise(a);

    final int[] indices = new int[nparams];
    int nnparams = 0;
    for (int j = 0; j < nparams; j++) {
      if (ignore[j]) {
        continue;
      }
      indices[nnparams++] = j;
    }

    for (int i = 0; i < n; i++) {
      // Function must produce a positive output.
      final double xi = y[i];

      // The code provided in Laurence & Chromy (2010) Nature Methods 7, 338-339, SI
      // effectively ignores the function value below zero. This could lead to a
      // situation where the best chisq value can be achieved by setting the output
      // function to produce 0 for all evaluations. To cope with this we heavily
      // penalise the chisq value.
      // Optimally the function should be bounded to always produce a positive number.
      final double fi = func.eval(i, dfiDa);

      if (fi <= 0) {
        // We assume xi is positive
        if (xi != 0) {
          // Penalise the chi-squared value by assuming fi is a very small positive value
          chisq += (-xi - xi * LOG_FOR_MIN);
        }

        // We ignore this contribution to the gradient for stability
        // compute(alpha, beta, dfiDa, Double.MIN_VALUE, xi);
      } else {
        // We assume y[i] is positive but must handle zero
        if (xi <= 0.0) {
          chisq += fi;
        } else {
          chisq += (fi - xi - xi * Math.log(fi / xi));
        }

        compute(alpha, beta, dfiDa, fi, xi, indices, nnparams);
      }
    }

    symmetric(alpha);

    // Move the factor of 2 to the end
    return checkGradients(alpha, beta, nparams, chisq * 2);
  }

  /**
   * Find linearised.
   *
   * @param n the n
   * @param y Data to fit (must be strictly positive Poisson data)
   * @param fx the y fit
   * @param a the a
   * @param func the func
   * @return The MLE chi-squared value
   */
  @Override
  public double findLinearised(final int n, final double[] y, double[] fx, final double[] a,
      final NonLinearFunction func) {
    double chisq = 0;

    func.initialise(a);

    if (fx == null || fx.length < n) {
      for (int i = 0; i < n; i++) {
        // Function must produce a positive output.
        final double xi = y[i];

        // The code provided in Laurence & Chromy (2010) Nature Methods 7, 338-339, SI
        // effectively ignores the function value below zero. This could lead to a
        // situation where the best chisq value can be achieved by setting the output
        // function to produce 0 for all evaluations. To cope with this we heavily
        // penalise the chisq value.
        // Optimally the function should be bounded to always produce a positive number.
        final double fi = func.eval(i);

        if (fi <= 0) {
          // We assume xi is positive
          if (xi != 0) {
            // Penalise the chi-squared value by assuming fi is a very small positive value
            chisq += (-xi - xi * LOG_FOR_MIN);
          }

          // We assume y[i] is positive but must handle zero
        } else if (xi <= 0.0) {
          chisq += fi;
        } else {
          chisq += (fi - xi - xi * Math.log(fi / xi));
        }
      }
    } else {
      for (int i = 0; i < n; i++) {
        // Function must produce a positive output.
        final double xi = y[i];

        // The code provided in Laurence & Chromy (2010) Nature Methods 7, 338-339, SI
        // effectively ignores the function value below zero. This could lead to a
        // situation where the best chisq value can be achieved by setting the output
        // function to produce 0 for all evaluations. To cope with this we heavily
        // penalise the chisq value.
        // Optimally the function should be bounded to always produce a positive number.
        final double fi = func.eval(i);
        fx[i] = fi;

        if (fi <= 0) {
          // We assume xi is positive
          if (xi != 0) {
            // Penalise the chi-squared value by assuming fi is a very small positive value
            chisq += (-xi - xi * LOG_FOR_MIN);
          }

          // We assume y[i] is positive but must handle zero
        } else if (xi <= 0.0) {
          chisq += fi;
        } else {
          chisq += (fi - xi - xi * Math.log(fi / xi));
        }
      }
    }

    // Move the factor of 2 to the end
    return chisq * 2;
  }

  /**
   * Zero the working region of the input matrix alpha and vector beta.
   *
   * @param alpha the alpha
   * @param beta the beta
   */
  protected void zero(final double[][] alpha, final double[] beta) {
    for (int i = 0; i < nparams; i++) {
      beta[i] = 0;
      for (int j = 0; j <= i; j++) {
        alpha[i][j] = 0;
      }
    }
  }

  /**
   * Compute the matrix alpha and vector beta.
   *
   * @param alpha the alpha
   * @param beta the beta
   * @param dfiDa the gradient of the function with respect to each parameter a
   * @param fi the function value at index i
   * @param xi the data value at index i
   */
  protected void compute(final double[][] alpha, final double[] beta, final double[] dfiDa,
      final double fi, final double xi) {
    final double xi_fi = xi / fi;
    final double xi_fi2 = xi_fi / fi;
    final double e = 1 - (xi_fi);
    // For + summation:
    // final double e = (xi_fi) - 1;

    // Compute:
    // Laurence &amp; Chromy (2010) Nature Methods 7, 338-339, SI
    // alpha - the Hessian matrix (the square matrix of second-order partial derivatives of a
    // function;
    // that is, it describes the local curvature of a function of many variables.)
    // beta - the gradient vector of the function's partial first derivatives with respect to the
    // parameters

    for (int k = 0; k < nparams; k++) {
      final double wgt = dfiDa[k] * xi_fi2;

      for (int l = 0; l <= k; l++) {
        // This is the non-optimised version:
        // alpha[k][l] += dfiDa[k] * dfiDa[l] * xi / (fi * fi);
        alpha[k][l] += wgt * dfiDa[l];
      }

      // This is the non-optimised version:
      // beta[k] -= (1 - xi / fi) * dfiDa[k];
      beta[k] -= e * dfiDa[k];
      // + summation:
      // beta[k] += e * dfiDa[k];
    }
  }

  /**
   * Compute the matrix alpha and vector beta.
   *
   * @param alpha the alpha
   * @param beta the beta
   * @param dfiDa the gradient of the function with respect to each parameter a
   * @param fi the function value at index i
   * @param xi the data value at index i
   * @param indices the indices
   * @param nnparams the nnparams
   */
  protected void compute(final double[][] alpha, final double[] beta, final double[] dfiDa,
      final double fi, final double xi, int[] indices, int nnparams) {
    final double xi_fi = xi / fi;
    final double xi_fi2 = xi_fi / fi;
    final double e = 1 - (xi_fi);

    // Compute:
    // Laurence &amp; Chromy (2010) Nature Methods 7, 338-339, SI
    // alpha - the Hessian matrix (the square matrix of second-order partial derivatives of a
    // function;
    // that is, it describes the local curvature of a function of many variables.)
    // beta - the gradient vector of the function's partial first derivatives with respect to the
    // parameters

    for (int k = 0; k < nnparams; k++) {
      final double wgt = dfiDa[indices[k]] * xi_fi2;

      for (int l = 0; l <= k; l++) {
        alpha[k][l] += wgt * dfiDa[indices[l]];
      }

      beta[k] -= e * dfiDa[k];
    }
  }

  /**
   * Compute the matrix alpha and vector beta when the data value is zero.
   *
   * @param beta the beta
   * @param dfiDa the gradient of the function with respect to each parameter a
   * @param fi the function value at index i
   */
  protected void compute0(final double[] beta, final double[] dfiDa, final double fi) {
    // Assume xi is zero. This removes most of the computation

    for (int k = 0; k < nparams; k++) {
      beta[k] -= dfiDa[k];
    }
  }

  /**
   * Generate a symmetric matrix alpha.
   *
   * @param alpha the alpha
   */
  protected void symmetric(final double[][] alpha) {
    for (int i = 0; i < nparams - 1; i++) {
      for (int j = i + 1; j < nparams; j++) {
        alpha[i][j] = alpha[j][i];
      }
    }
  }

  /**
   * Evaluate the function and compute the MLE chi-squared value and the gradient with respect to
   * the model parameters.
   *
   * <p>A call to {@link #isNaNGradients()} will indicate if the gradients were invalid.
   *
   * @param x n observations
   * @param y Data to fit
   * @param a Set of m coefficients
   * @param dfDa the gradient vector of the function's partial first derivatives with respect to the
   *        parameters (size m)
   * @param func Non-linear fitting function
   * @return The MLE chi-squared value
   */
  @Override
  public double evaluate(final int[] x, final double[] y, final double[] a, final double[] dfDa,
      final NonLinearFunction func) {
    double chisq = 0;
    final double[] dfiDa = new double[nparams];

    zero(dfDa);

    func.initialise(a);

    for (int i = 0; i < x.length; i++) {
      // Function must produce a positive output.
      final double xi = y[i];

      // The code provided in Laurence & Chromy (2010) Nature Methods 7, 338-339, SI
      // effectively ignores the function value below zero. This could lead to a
      // situation where the best chisq value can be achieved by setting the output
      // function to produce 0 for all evaluations. To cope with this we heavily
      // penalise the chisq value.
      // Optimally the function should be bounded to always produce a positive number.
      final double fi = func.eval(i, dfiDa);

      if (fi <= 0) {
        // We assume xi is positive
        if (xi != 0) {
          // Penalise the chi-squared value by assuming fi is a very small positive value
          chisq += (-xi - xi * LOG_FOR_MIN);
        }

        //// Note: We ignore this contribution to the gradient for stability
        // final double e = 1 - (xi / Double.MIN_VALUE);
        // for (int k = 0; k < nparams; k++)
        // {
        // dfDa[k] += e * dfiDa[k];
        // }
      } else {
        // We assume y[i] is positive but must handle zero
        if (xi <= 0.0) {
          chisq += fi;
        } else {
          chisq += (fi - xi - xi * Math.log(fi / xi));
        }

        final double e = 1 - (xi / fi);

        // Compute:
        // Laurence &amp; Chromy (2010) Nature Methods 7, 338-339, SI, Equation 6
        // dfDa - the gradient vector of the function's partial first derivatives with respect to
        // the parameters

        for (int k = 0; k < nparams; k++) {
          dfDa[k] += e * dfiDa[k];
        }
      }
    }

    checkGradients(dfDa, nparams);

    // Move the factor of 2 to the end
    for (int j = 0; j < nparams; j++) {
      dfDa[j] *= 2;
    }
    return chisq * 2;
  }

  /**
   * Compute the Poisson log likelihood for the data.
   *
   * @param x the data
   * @param a the parameters for the function
   * @param func the function (must evaluate to strictly positive for all the data points)
   * @return the Poisson log likelihood
   */
  public double logLikelihood(final double[] x, final double[] a, final NonLinearFunction func) {
    double ll = 0;
    func.initialise(a);
    for (int i = 0; i < x.length; i++) {
      ll += PoissonCalculator.logLikelihood(func.eval(i), x[i]);
    }
    return ll;
  }
}
