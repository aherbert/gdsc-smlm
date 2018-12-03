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
 * Provides sort functions.
 *
 * @since 1.0.0
 */
package uk.ac.sussex.gdsc.smlm.results.sort;

/*
 * Java 1.7 Arrays.sort can throw an exception when Object1.compareTo(Object2) does not equal
 * -Object2.compareTo(Object1). This was silently ignored in previous JVMs. <p> The exception is
 * throw in java.utils.ComparableTimSort: throw new
 * IllegalArgumentException("Comparison method violates its general contract!"); <p> This can occur
 * when sorting mixed lists of PeakResult objects using their compareTo method as sub-classes may
 * have different a compareTo method. To work around this the compareTo method has been removed from
 * the PeakResult object and comparators added to this package for common comparisons.
 */
