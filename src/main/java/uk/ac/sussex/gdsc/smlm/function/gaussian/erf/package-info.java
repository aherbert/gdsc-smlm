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
 * Provides classes for computation of 2D Gaussian function using the Error function.
 *
 * <p>The 2D Gaussian function is used as a model of the Point Spread Function (PSF) of a microscope
 * and can be fit to Single Molecule Localisation Microscopy data.
 *
 * <p>The functions in this package can compute the Gaussian function over a 2D range that
 * represents for example the pixels in an image. This package evaluates the function as an integral
 * over the 1x1 area surrounding point (X,Y).
 *
 * <p>The function is evaluated using the XY separability of the Gaussian. Each dimension must align
 * with the XY axes and so the functions <strong>cannot</strong> support rotated elliptical 2D
 * Gaussians.
 *
 * <p>This work is based on the paper:<br> Smith et al, (2010). Fast, single-molecule localisation
 * that achieves theoretically minimum uncertainty. Nature Methods 7, 373-375 (supplementary note).
 *
 * @since 1.0.0
 */

package uk.ac.sussex.gdsc.smlm.function.gaussian.erf;
