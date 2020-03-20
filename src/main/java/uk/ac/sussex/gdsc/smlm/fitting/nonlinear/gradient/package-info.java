/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

/**
 * Provides computation of the problem for a generic (non-linear) function that models the given the
 * data.
 *
 * <p>Types of problem supported are:
 *
 * <ul>
 *
 * <li>The Hessian-type matrix of partial first-order derivatives used in the <a
 * href="https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm"> Levenberg–Marquardt
 * algorithm</a>.
 *
 * <li>The Hessian-type matrix of partial first-order derivatives used in Laurence &amp; Chromy
 * (2010) Maximum Likelihood Estimation variant of the Levenberg–Marquardt algorithm.
 *
 * <li>The Hessian-type matrix of partial first-order derivatives used in Lin, et al (2017) weighted
 * variant of the Levenberg–Marquardt algorithm.
 *
 * <li>The Fisher information matrix for a Poisson process.
 *
 * <li>The Smith et al, (2010) Newton-Raphson update vector for a Poisson process using the first
 * and second partial derivatives.
 *
 * </ul>
 *
 * @since 1.0.0
 */

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient;
