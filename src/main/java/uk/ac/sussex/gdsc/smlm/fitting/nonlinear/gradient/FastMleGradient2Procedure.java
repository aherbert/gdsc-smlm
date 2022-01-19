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

import java.util.Arrays;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Procedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient2Function;
import uk.ac.sussex.gdsc.smlm.function.Gradient2Procedure;
import uk.ac.sussex.gdsc.smlm.function.PoissonCalculator;
import uk.ac.sussex.gdsc.smlm.function.ValueProcedure;

/**
 * Calculates the Newton-Raphson update vector for a Poisson process using the first and second
 * partial derivatives.
 *
 * <p>Ref: Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically
 * minimum uncertainty. Nature Methods 7, 373-375 (supplementary note), Eq. 12.
 */
public class FastMleGradient2Procedure
    implements ValueProcedure, Gradient1Procedure, Gradient2Procedure {
  /** The data to fit (must be positive, i.e. the value of a Poisson process). */
  protected final double[] x;
  /** The function. */
  protected final Gradient2Function func;
  /** The poisson calculator. */
  protected PoissonCalculator poissonCalculator;

  /**
   * The number of gradients.
   */
  public final int numberOfGradients;
  /**
   * The first derivative of the Poisson log likelihood with respect to each parameter.
   */
  public final double[] d1;
  /**
   * The second derivative of the Poisson log likelihood with respect to each parameter.
   */
  public final double[] d2;

  // Names u and k match the uk nomenclature of Smith, et al.
  // CHECKSTYLE.OFF: MemberName

  /** Counter. */
  protected int k;

  /**
   * The value of the function. This is updated by calls to {@link #computeValue(double[])},
   * {@link #computeFirstDerivative(double[])}, {@link #computeSecondDerivative(double[])}
   */
  public final double[] u;

  // CHECKSTYLE.ON: MemberName

  /**
   * Instantiates a new procedure.
   *
   * @param x Data to fit (must be positive, i.e. the value of a Poisson process)
   * @param func Gradient function (must produce a strictly positive value, i.e. the mean of a
   *        Poisson process)
   */
  public FastMleGradient2Procedure(final double[] x, final Gradient2Function func) {
    this.x = x;
    this.u = new double[x.length];
    this.func = func;
    this.numberOfGradients = func.getNumberOfGradients();
    d1 = new double[numberOfGradients];
    d2 = new double[numberOfGradients];
  }

  /**
   * Calculates the first and second derivative of the Poisson log likelihood with respect to each
   * parameter.
   *
   * @param a Set of coefficients for the function
   */
  public void computeSecondDerivative(final double[] a) {
    k = 0;
    reset2();
    func.initialise2(a);
    func.forEach((Gradient2Procedure) this);
  }

  /**
   * Reset the first and second derivative vectors.
   */
  protected void reset2() {
    Arrays.fill(d1, 0);
    Arrays.fill(d2, 0);
  }

  /**
   * Calculates the first derivative of the Poisson log likelihood with respect to each parameter.
   *
   * @param a Set of coefficients for the function
   * @return the first derivative of the Poisson log likelihood with respect to each parameter
   */
  public double[] computeFirstDerivative(final double[] a) {
    k = 0;
    reset1();
    func.initialise1(a);
    func.forEach((Gradient1Procedure) this);
    return d1;
  }

  /**
   * Reset the first derivative vector.
   */
  protected void reset1() {
    Arrays.fill(d1, 0);
  }

  /**
   * Compute the value of the function.
   *
   * @param a the a
   * @return the double[]
   */
  public double[] computeValue(final double[] a) {
    k = 0;
    func.initialise0(a);
    func.forEach((ValueProcedure) this);
    return u;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Calculates the Newton-Raphson update vector for a Poisson process. Variables are named as
   * per the Smith, et al (2010) paper.
   */
  @Override
  public void execute(double uk, double[] dukDt, double[] d2ukDt2) {
    u[k] = uk;
    final double xk = x[k++];
    if (xk == 0) {
      for (int i = 0; i < numberOfGradients; i++) {
        d1[i] -= dukDt[i];
        d2[i] -= d2ukDt2[i];
      }
    } else {
      final double xk_uk_minus1 = xk / uk - 1.0;
      final double xk_uk2 = xk / (uk * uk);
      for (int i = 0; i < numberOfGradients; i++) {
        d1[i] += dukDt[i] * xk_uk_minus1;
        d2[i] += d2ukDt2[i] * xk_uk_minus1 - dukDt[i] * dukDt[i] * xk_uk2;
      }
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>Variables are named as per the Smith, et al (2010) paper.
   */
  @Override
  public void execute(double uk, double[] dukDt) {
    u[k] = uk;
    final double xk = x[k++];
    if (xk == 0) {
      for (int i = 0; i < numberOfGradients; i++) {
        d1[i] -= dukDt[i];
      }
    } else {
      final double xk_uk_minus1 = xk / uk - 1.0;
      for (int i = 0; i < numberOfGradients; i++) {
        d1[i] += dukDt[i] * xk_uk_minus1;
      }
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>Variables are named as per the Smith, et al (2010) paper.
   */
  @Override
  public void execute(double uk) {
    u[k++] = uk;
  }

  /**
   * Calculates the Poisson log likelihood.
   *
   * @param a Set of coefficients for the function
   * @return the Poisson log likelihood
   */
  public double computeLogLikelihood(final double[] a) {
    computeValue(a);
    return computeLogLikelihood();
  }

  /**
   * Calculates the Poisson log likelihood using the last value of the function.
   *
   * @return the Poisson log likelihood
   */
  public double computeLogLikelihood() {
    return getPoissonCalculator().logLikelihood(u);
  }

  /**
   * Calculates the pseudo Poisson log likelihood using the last value of the function.
   *
   * <p>The pseudo log-likelihood is equivalent to the log-likelihood without subtracting the
   * log(x!) term. It can be converted to the log-likelihood by subtracting
   * {@link #computeLogXFactorialTerm()}.
   *
   * <p>This term is suitable for use in maximum likelihood routines.
   *
   * <pre>
   * pseudo ll = x * log(u) - u
   * </pre>
   *
   * @return the pseudo Poisson log likelihood
   */
  public double computePseudoLogLikelihood() {
    return getPoissonCalculator().pseudoLogLikelihood(u);
  }

  /**
   * Computes the log X factorial term to convert the pseudo log-likelihood to the log-likelihood.
   *
   * <pre>
   * ll = pseudo ll - log(x!)
   * </pre>
   *
   * @return the log X factorial term
   */
  public double computeLogXFactorialTerm() {
    return getPoissonCalculator().getLogXFactorialTerm();
  }

  /**
   * Gets the poisson calculator, creating using the values x if necessary.
   *
   * @return the poisson calculator
   */
  private PoissonCalculator getPoissonCalculator() {
    PoissonCalculator c = poissonCalculator;
    if (c == null) {
      c = new PoissonCalculator(x);
      poissonCalculator = c;
    }
    return c;
  }

  /**
   * Calculates the Poisson log likelihood ratio.
   *
   * @param a Set of coefficients for the function
   * @return the Poisson log likelihood ratio
   */
  public double computeLogLikelihoodRatio(final double[] a) {
    computeValue(a);
    return computeLogLikelihoodRatio();
  }

  /**
   * Calculates the Poisson log likelihood ratio using the last value of the function.
   *
   * @return the Poisson log likelihood ratio
   */
  public double computeLogLikelihoodRatio() {
    return computeLogLikelihoodRatio(computeLogLikelihood());
  }

  /**
   * Calculates the Poisson log likelihood ratio using the given log-likelihood.
   *
   * @param logLikelihood the log likelihood
   * @return the Poisson log likelihood ratio
   */
  public double computeLogLikelihoodRatio(double logLikelihood) {
    return getPoissonCalculator().getLogLikelihoodRatio(logLikelihood);
  }

  /**
   * Checks if is na N gradients.
   *
   * @return True if the last update calculation produced gradients with NaN values.
   */
  public boolean isNaNGradients() {
    for (int i = numberOfGradients; i-- > 0;) {
      if (Double.isNaN(d1[i])) {
        return true;
      }
      if (Double.isNaN(d2[i])) {
        return true;
      }
    }
    return false;
  }
}
