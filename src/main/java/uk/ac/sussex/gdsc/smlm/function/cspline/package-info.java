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

/**
 * Provides classes for computation of cubic spline (CSpline) functions.
 *
 * <p>A cubic spline function can model the Point Spread Function (PSF) of a microscope and can be
 * fit to Single Molecule Localisation Microscopy data.
 *
 * <p>A cubic spline is defined using (X,Y,Z) dimensions. To model a PSF the Z coordinate is used to
 * define the XY slice of the spline. The functions in this package can compute the spline function
 * over a 2D (X,Y) range from the Z-slice that represents for example the pixels in an image. This
 * package evaluates the function as a single point for each (X,Y) position. That is the function
 * integral over the 1x1 area surrounding point (X,Y) is approximated using the function value at
 * (X,Y).
 *
 * @since 1.0.0
 */

package uk.ac.sussex.gdsc.smlm.function.cspline;
