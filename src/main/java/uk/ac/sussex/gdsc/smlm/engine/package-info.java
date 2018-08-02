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
 * Provides the engine to perform filtering and fitting of Single Molecule
 * Localisation Microscopy (SMLM) data.
 * <p>
 * Data is processed as 2D frames from an image source. The engine distributes the
 * input image frames to workers for analysis. Each worker filters the frames to
 * identify candidate localisations and then fits the candidates using a Gaussian 2D function.
 * The fits are assessed using criteria such as localisation precision, Signal-to-Noise Ratio (SNR),
 * drift, and spot width. Fitting may not process all candidates due to the use of smart
 * stopping criteria based on the acceptance rate.
 *
 * @since 1.0.0
 */
package uk.ac.sussex.gdsc.smlm.engine;
