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
package uk.ac.sussex.gdsc.smlm.fitting;

/**
 * Calculator for the Fisher information, a symmetric positive definite matrix containing the amount
 * of information that an observable random variable X carries about an unknown parameter θ of a
 * distribution that models X. <p> The calculator will compute the Fisher Information Matrix (I)
 * using numerical gradients:
 *
 * <pre>
 * Iij = E [ (d log f(X;θ) / dθi) (d log f(X;θ) / dθj) | θ ]
 * E = Expected value
 * f(X;θ) = Likelihood function for data X given parameters θ
 * </pre>
 */
public interface FisherInformationCalculator {
  /**
   * Compute the Fisher information, a symmetric positive definite matrix containing the amount of
   * information that an observable random variable X carries about an unknown parameter θ of a
   * distribution that models X.
   *
   * @param parameters the parameters (θ)
   * @return the fisher information matrix
   */
  public FisherInformationMatrix compute(double[] parameters);
}
