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

package uk.ac.sussex.gdsc.smlm.function;

import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.data.NotImplementedException;
import uk.ac.sussex.gdsc.smlm.fitting.FisherInformationMatrix;

import java.util.Arrays;

/**
 * This is a wrapper for any function to compute the negative log-likelihood.
 */
public abstract class LikelihoodWrapper {
  /** The function. */
  protected final NonLinearFunction f;

  /** The parameters. */
  protected final double[] a;

  /** The observed values. */
  protected final double[] data;

  /** The number of observed values. */
  protected final int n;

  /** The number of variables. */
  protected final int nVariables;

  private double lastScore;
  private double[] lastVariables;

  /**
   * Initialise the function.
   *
   * <p>The input parameters must be the full parameters for the non-linear function. Only those
   * parameters with gradient indices should be passed in to the functions to obtain the value (and
   * gradient).
   *
   * @param f The function to be used to calculated the expected values
   * @param a The initial parameters for the function
   * @param k The observed values
   * @param n The number of observed values
   */
  public LikelihoodWrapper(NonLinearFunction f, double[] a, double[] k, int n) {
    this.f = f;
    this.a = Arrays.copyOf(a, a.length);
    this.data = k;
    this.n = n;
    nVariables = f.gradientIndices().length;
  }

  /**
   * Copy the variables into the appropriate parameter positions for the NonLinearFunction.
   *
   * @param variables the variables
   */
  protected void initialiseFunction(double[] variables) {
    final int[] gradientIndices = f.gradientIndices();
    for (int i = 0; i < gradientIndices.length; i++) {
      a[gradientIndices[i]] = variables[i];
    }
    f.initialise(a);
  }

  /**
   * Check if the variable match those last used for computation of the value.
   *
   * @param variables the variables
   * @return True if the variables are the same
   */
  private boolean sameVariables(double[] variables) {
    if (lastVariables != null) {
      for (int i = 0; i < variables.length; i++) {
        if (variables[i] != lastVariables[i]) {
          return false;
        }
      }
      return true;
    }
    return false;
  }

  /**
   * Compute the negative log likelihood. Returns positive infinity if the likelihood is zero at any
   * point in the observed values.
   *
   * @param variables The variables of the function
   * @return The negative log likelihood
   */
  public double likelihood(double[] variables) {
    // Check if we have a cached score
    if (sameVariables(variables)) {
      return lastScore;
    }

    initialiseFunction(variables);
    lastVariables = variables.clone();
    return lastScore = computeLikelihood();
  }

  /**
   * Compute the negative log likelihood. Returns positive infinity if the likelihood is zero at any
   * point in the observed values.
   *
   * <p>The wrapped NonLinearFunction will be correctly initialised before this function is called.
   *
   * @return The negative log likelihood
   */
  protected abstract double computeLikelihood();

  /**
   * Compute the negative log likelihood and the gradient. Returns positive infinity if the
   * likelihood is zero at any point in the observed values. In this case the gradient computed is
   * invalid.
   *
   * @param variables The variables of the function
   * @param gradient The gradient (must be equal length to the variables array)
   * @return The negative log likelihood
   * @throws NotImplementedException If the sub-class cannot provide gradients
   */
  public double likelihood(double[] variables, double[] gradient) throws NotImplementedException {
    initialiseFunction(variables);
    lastVariables = variables.clone();
    return lastScore = computeLikelihood(gradient);
  }

  /**
   * Compute the negative log likelihood and the gradient. Returns positive infinity if the
   * likelihood is zero at any point in the observed values. In this case the gradient computed is
   * invalid.
   *
   * <p>The wrapped NonLinearFunction will be correctly initialised before this function is called
   *
   * @param gradient The gradient (must be equal length to the variables array)
   * @return The negative log likelihood
   * @throws NotImplementedException If the sub-class cannot provide gradients
   */
  protected double computeLikelihood(double[] gradient) throws NotImplementedException {
    throw new NotImplementedException();
  }

  /**
   * Compute the negative log likelihood at observed value i. Returns positive infinity if the
   * likelihood is zero at the observed value.
   *
   * @param variables The variables of the function
   * @param i Observed value i
   * @return The negative log likelihood
   */
  public double likelihood(double[] variables, int i) {
    initialiseFunction(variables);
    return computeLikelihood(i);
  }

  /**
   * Compute the negative log likelihood at observed value i. Returns positive infinity if the
   * likelihood is zero at the observed value.
   *
   * <p>The wrapped NonLinearFunction will be correctly initialised before this function is called
   *
   * @param i Observed value i
   * @return The negative log likelihood
   */
  protected abstract double computeLikelihood(int i);

  /**
   * Compute the negative log likelihood and gradient of the function at observed value i. Returns
   * positive infinity if the likelihood is zero at the observed value. In this case the gradient
   * computed will be invalid.
   *
   * @param variables The variables of the function
   * @param gradient The gradient (must be equal length to the variables array)
   * @param i Observed value i
   * @return The negative log likelihood
   * @throws NotImplementedException If the sub-class cannot provide gradients
   */
  public double likelihood(double[] variables, double[] gradient, int i)
      throws NotImplementedException {
    initialiseFunction(variables);
    return computeLikelihood(gradient, i);
  }

  /**
   * Compute the negative log likelihood and gradient of the function at observed value i. Returns
   * positive infinity if the likelihood is zero at the observed value. In this case the gradient
   * computed will be invalid.
   *
   * <p>The wrapped NonLinearFunction will be correctly initialised before this function is called
   *
   * @param gradient The gradient (must be equal length to the variables array)
   * @param i Observed value i
   * @return The negative log likelihood
   * @throws NotImplementedException If the sub-class cannot provide gradients
   */
  protected double computeLikelihood(double[] gradient, int i) throws NotImplementedException {
    throw new NotImplementedException();
  }

  /**
   * Specify if the likelihood function can compute gradients. If false then the calls to the
   * likelihood functions to compute the gradient will throw a
   * {@link uk.ac.sussex.gdsc.core.data.NotImplementedException }
   *
   * @return True if the likelihood function can compute gradients
   */
  public abstract boolean canComputeGradient();

  /**
   * Compute the Fisher's Information Matrix (I) for fitted variables.
   *
   * <p>Note that this is only a true Fisher information matrix if the function returns the expected
   * value for a Poisson process. In this case the equation reduces to:
   *
   * <pre>
   * Iaa = sum(i) (dYi da) * (dYi da) / Yi
   * </pre>
   *
   * <p>See Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically
   * minimum uncertainty. Nature Methods 7, 373-375 (supplementary note), Eq. 9.
   *
   * @param variables The variables of the function
   * @return Fisher's Information Matrix (I)
   */
  public double[][] fisherInformation(final double[] variables) {
    initialiseFunction(variables);

    final double[] du_da = new double[nVariables];

    final double[][] I = new double[nVariables][nVariables];

    for (int k = 0; k < n; k++) {
      final double uk = f.eval(k, du_da);
      final double yk = 1 / uk;
      for (int i = 0; i < nVariables; i++) {
        final double du_dai = yk * du_da[i];
        for (int j = 0; j <= i; j++) {
          I[i][j] += du_dai * du_da[j];
        }
      }
    }

    // Generate symmetric matrix
    for (int i = 0; i < nVariables - 1; i++) {
      for (int j = i + 1; j < nVariables; j++) {
        if (Double.isNaN(I[j][i])) {
          throw new DataException("Invalid gradients");
        }
        I[i][j] = I[j][i];
      }
    }

    return I;
  }

  /**
   * Compute the Cramer-Rao Lower Bound (CRLB) variance for fitted variables using the central
   * diagonal of the inverted Fisher's Information Matrix (I).
   *
   * <p>The information matrix is inverted and the central diagonal returned.
   *
   * @param variables The variables of the function
   * @return CRLB (or null if inversion failed)
   */
  public double[] crlb(final double[] variables) {
    return crlb(variables, false);
  }

  /**
   * Compute the Cramer-Rao Lower Bound (CRLB) variance for fitted variables using the central
   * diagonal of the inverted Fisher's Information Matrix (I).
   *
   * <p>The information matrix is inverted and the central diagonal returned. If the inversion fails
   * then the routine optionally returns the reciprocal of the diagonal element to find a (possibly
   * loose) lower bound.
   *
   * @param variables The variables of the function
   * @param allowReciprocal the allow reciprocal flag
   * @return CRLB (or null if inversion failed)
   */
  public double[] crlb(final double[] variables, boolean allowReciprocal) {
    final double[][] I = fisherInformation(variables);
    return new FisherInformationMatrix(I).crlb(allowReciprocal);
  }
}
