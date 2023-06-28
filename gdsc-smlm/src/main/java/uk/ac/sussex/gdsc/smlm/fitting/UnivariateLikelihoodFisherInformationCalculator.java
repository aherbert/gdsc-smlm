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

package uk.ac.sussex.gdsc.smlm.fitting;

import java.util.Arrays;
import java.util.Objects;
import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.smlm.function.FisherInformation;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Procedure;

/**
 * Calculator for the Fisher information, a symmetric positive definite matrix containing the amount
 * of information that an observable random variable X carries about an unknown parameter θ of a
 * distribution that models X.
 *
 * <pre>
 * Iij = E [ (d log f(X;θ) / dθi) (d log f(X;θ) / dθj) | θ ]
 * E = Expected value
 * f(X;θ) = Likelihood function for data X given parameters θ
 * </pre>
 *
 * <p>The calculator assumes a gradient function outputs a single value (v) that can be input into a
 * univariate probability distribution. The Fisher information is derived using using the equation
 * of Chao, et al (2013) Nature Methods, 10, 335-338, SI Eq S3.
 *
 * <pre>
 * Iij = sum(k) (d v(θ,k) / dθi) . (d v(θ,k) / dθj) . E [ (d ln(p(z|v(θ,k)) / d v(θ,k) )^2 ]
 * k = the number of points over which the function is evaluated
 * v(θ,k) = the function value at point k given parameters θ
 * p(z|v(θ,k) = the likelihood function for data z given value v
 * </pre>
 *
 * <p>This reduces to:
 *
 * <pre>
 * Iij = sum(k) (d v(θ,k) / dθi) . (d v(θ,k) / dθj) . I(v(θ,k))
 * I(v(θ,k)) = the Fisher information of the likelihood function for value v at point k
 * </pre>
 *
 */
public class UnivariateLikelihoodFisherInformationCalculator
    implements FisherInformationCalculator {
  /** The gradient function. */
  protected final Gradient1Function gf;

  /** The Fisher information for each evaluated point in the function. */
  protected final FisherInformation[] fi;

  /**
   * Instantiates a new univariate likelihood fisher information calculator.
   *
   * @param gf the gradient function
   * @param fi the fisher information of the output value of the function
   */
  public UnivariateLikelihoodFisherInformationCalculator(Gradient1Function gf,
      FisherInformation fi) {
    Objects.requireNonNull(fi, "fi");
    this.gf = Objects.requireNonNull(gf, "gf");
    this.fi = new FisherInformation[gf.size()];
    Arrays.fill(this.fi, fi);
  }

  /**
   * Instantiates a new univariate likelihood fisher information calculator.
   *
   * @param gf the gradient function
   * @param fi the fisher information of the output of each value of the function
   */
  public UnivariateLikelihoodFisherInformationCalculator(Gradient1Function gf,
      FisherInformation[] fi) {
    this.gf = Objects.requireNonNull(gf, "gf");
    this.fi = Objects.requireNonNull(fi, "fi");
    if (fi.length != gf.size()) {
      throw new IllegalArgumentException(
          "Fisher information must be provided for each function value");
    }
  }

  /**
   * {@inheritDoc}
   *
   * @throws DataException If the Fisher information cannot be computed for a function value
   * @throws DataException If the Fisher information is infinite for a function value
   */
  @Override
  public FisherInformationMatrix compute(double[] parameters) {
    final int n = gf.getNumberOfGradients();
    final double[] data = new double[n * (n + 1) / 2];
    gf.initialise1(parameters);
    gf.forEach(new Gradient1Procedure() {
      // CHECKSTYLE.OFF: MemberName
      int k = -1;
      // CHECKSTYLE.ON: MemberName

      @Override
      public void execute(double value, double[] dvDt) {
        k++;
        if (!fi[k].isValid(value)) {
          return;
        }
        // Get the Fisher information of the value
        final double f = fi[k].getFisherInformation(value);
        if (f == 0) {
          // No summation
          return;
        }
        if (f == Double.POSITIVE_INFINITY) {
          throw new DataException("Fisher information is infinite at f(" + k + ")");
        }

        // Compute the actual matrix data
        for (int i = 0, c = 0; i < n; i++) {
          final double wgt = f * dvDt[i];
          for (int j = 0; j <= i; j++) {
            data[c++] += wgt * dvDt[j];
          }
        }
      }
    });
    // Generate symmetric matrix
    final double[] matrix = new double[n * n];
    for (int i = 0, c = 0; i < n; i++) {
      for (int j = 0; j <= i; j++) {
        matrix[i * n + j] = matrix[j * n + i] = data[c++];
      }
    }
    return new FisherInformationMatrix(matrix, n);
  }
}
