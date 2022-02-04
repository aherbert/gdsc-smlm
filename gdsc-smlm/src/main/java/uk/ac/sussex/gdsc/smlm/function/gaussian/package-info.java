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

/**
 * Provides classes for computation of 2D Gaussian function.
 *
 * <p>The 2D Gaussian function is used as a model of the Point Spread Function (PSF) of a microscope
 * and can be fit to Single Molecule Localisation Microscopy data.
 *
 * <p>The functions in this package can compute the Gaussian function over a 2D range that
 * represents for example the pixels in an image. This package evaluates the function as a single
 * point for each (X,Y) position. That is the function integral over the 1x1 area surrounding point
 * (X,Y) is approximated using the function value at (X,Y).
 *
 * <p>The functions can support rotated elliptical 2D Gaussians.
 *
 * @since 1.0
 */

package uk.ac.sussex.gdsc.smlm.function.gaussian;
