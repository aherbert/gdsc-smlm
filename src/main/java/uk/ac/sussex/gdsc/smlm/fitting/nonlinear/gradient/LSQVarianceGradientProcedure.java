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

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient;

import uk.ac.sussex.gdsc.smlm.fitting.linear.EJMLLinearSolver;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Procedure;

import java.util.Arrays;

/**
 * Compute the variance of the parameters of the function assuming a least squares fit of a Poisson
 * process.
 *
 * <p>Uses the Mortensen formula (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 25.
 *
 * <p>Note: If the fit is for a Poisson process convolved with a Gamma distribution, e.g. in the
 * case of an EM-CCD camera, and the function describes a 2D point-spread function (PSF) with long
 * tails then the variance for position parameters will be scaled by a factor of 2 (See Mortensen,
 * SI 4.3 for assumptions and proof using MLE). The estimated variance of other function parameters
 * (e.g. background, total intensity) will be incorrect.
 */
public class LSQVarianceGradientProcedure implements Gradient1Procedure {
  /** Returned when computation was OK. */
  public static final int STATUS_OK = 0;
  /** Returned when computation failed due to NaN values in the gradients. */
  public static final int STATUS_BAD_GRADIENTS = 1;
  /** Returned when computation failed to invert the information matrix. */
  public static final int STATUS_FAILED_INVERSION = 2;

  /** The function. */
  protected final Gradient1Function func;

  /**
   * The number of gradients.
   */
  public final int n;

  /** Working space for I = sum_i { Ei,a * Ei,b }. */
  protected final double[] I;

  /** Working space for E = sum_i { Ei * Ei,a * Ei,b }. */
  protected final double[] E;

  /** The solver. */
  protected final EJMLLinearSolver solver;

  /** The variance. */
  public final double[] variance;

  /**
   * Instantiates a new procedure.
   *
   * @param func Gradient function
   */
  public LSQVarianceGradientProcedure(final Gradient1Function func) {
    this(func, EJMLLinearSolver.createForInversion(1e-2));
  }

  /**
   * Instantiates a new procedure.
   *
   * @param func Gradient function
   * @param solver The solver used to invert the Fisher information matrix to find the Cramér–Rao
   *        lower bound (CRLB).
   * @throws IllegalArgumentException if the solver is null
   */
  public LSQVarianceGradientProcedure(final Gradient1Function func, EJMLLinearSolver solver) {
    if (solver == null) {
      throw new IllegalArgumentException("The solver cannot be null");
    }

    this.func = func;
    this.n = func.getNumberOfGradients();

    I = new double[n * n];
    E = new double[I.length];
    variance = new double[n];

    this.solver = solver;
  }

  /**
   * Evaluate the function and compute the variance of the parameters. The result will be zero if
   * OK. Otherwise the failure reason can be checked using the status constants.
   *
   * @param a Set of coefficients for the function (if null the function must be pre-initialised)
   * @return the result
   */
  public int variance(final double[] a) {
    initialise();
    if (a != null) {
      func.initialise1(a);
    }
    func.forEach(this);
    if (finish()) {
      return STATUS_BAD_GRADIENTS;
    }
    if (!solver.invert(I, n)) {
      return STATUS_FAILED_INVERSION;
    }
    computeVariance();
    return STATUS_OK;
  }

  @Override
  public void execute(final double Ei, double[] Eix) {
    for (int a = 0; a < n; a++) {
      for (int b = 0, j = a * n; b <= a; b++, j++) {
        final double v = Eix[a] * Eix[b];
        I[j] += v;
        E[j] += Ei * v;
      }
    }
  }

  /**
   * Initialise for the computation using the first order gradients.
   */
  protected void initialise() {
    for (int a = 0; a < n; a++) {
      for (int b = 0, j = a * n; b <= a; b++, j++) {
        // System.out.printf("I[%d] = 0; E[%d] = 0;\n", j, j);
        I[j] = 0;
        E[j] = 0;
      }
    }
    Arrays.fill(variance, 0);
  }

  /**
   * Finish the computation using the first order gradients. Check the gradients are OK then
   * generate symmetric matrix for I and E.
   *
   * @return true, if the gradient computation failed (e.g. NaN gradients
   */
  protected boolean finish() {
    // System.out.printf("if (");
    for (int a = 0; a < n; a++) {
      for (int b = 0, j = a * n; b <= a; b++, j++) {
        // System.out.printf("%s I[%d]!=I[%d] ", (a+b==0)?"":"||", j, j);
        if (Double.isNaN(I[j])) {
          return true;
        }
      }
    }
    // System.out.printf(") return true;\n");
    // Generate symmetric matrix
    for (int a = 0; a < n; a++) {
      for (int b = 0; b < a; b++) {
        final int j = a * n + b;
        final int k = b * n + a;
        // System.out.printf("I[%d] = I[%d];\n", k, j);
        // System.out.printf("E[%d] = E[%d];\n", k, j);
        I[k] = I[j];
        E[k] = E[j];
      }
    }

    // System.out.println("Procedure");
    // System.out.println(DenseMatrix64F.wrap(n, n, I).toString());
    // System.out.println(DenseMatrix64F.wrap(n, n, E).toString());

    return false;
  }

  /**
   * Compute the variance using the inverted I and E matrices.
   */
  protected void computeVariance() {
    for (int a = 0; a < n; a++) {
      // Note: b==a as we only do the diagonal
      // System.out.printf("variance[%d] = \n", a);
      double v = 0;
      for (int ap = 0; ap < n; ap++) {
        for (int bp = 0; bp < n; bp++) {
          // System.out.printf("%s I[%d]*E[%d]*I[%d]\n", (ap+bp==0)?"":"+", a*n+ap, ap*n+bp,
          // bp*n+a);
          v += I[a * n + ap] * E[ap * n + bp] * I[bp * n + a];
        }
      }
      // System.out.printf(";\n");
      variance[a] = v;
    }
  }
}
