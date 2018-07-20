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
/**
 * Provides interfaces and classes for functions that compute values and gradients.
 * <p>
 * Examples of supported function are:
 * <ul>
 * <li>Fast computation of the <a href="https://en.wikipedia.org/wiki/Error_function">Error function</a>.
 * <li>Poisson probability functions with Fisher information.
 * <li>Poisson-Gaussian probability functions with Fisher information (for CCD/sCMOS cameras).
 * <li>Poisson-Gamma-Gaussian probability functions with Fisher information (for EM-CCD cameras).
 * <li>Fast log computation using tabulated values.
 * <li>Camera noise models (e.g. CCD, EM-CCD).
 * </ul>
 * 
 * @since 1.0.0
 */
package uk.ac.sussex.gdsc.smlm.function;
