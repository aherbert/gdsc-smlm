/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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

/**
 * Provides functions to compute the gradient of the likelihood for a univariate distribution.
 */
public interface GradientLikelihoodFunction extends LikelihoodFunction {
  /**
   * Compute the likelihood of an observation x given a parameter value theta.
   *
   * <p>This is the probability mass function P(X=x|θ) or the probability density function f(x|θ)
   * depending on parameter θ.
   *
   * @param x The observed value (x)
   * @param theta The parameter value (θ)
   * @param dpDt the gradient d P(X=x|θ) dθ
   * @return The likelihood
   */
  double likelihood(double x, double theta, double[] dpDt);
}
