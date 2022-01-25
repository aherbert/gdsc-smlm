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

package uk.ac.sussex.gdsc.smlm.utils;

/**
 * Class for computing standard math functions.
 *
 * <p>Provides an entry point to routines, or their equivalent, found in {@link java.lang.Math}.
 */
public final class StdMath {
  /** No public constructor. */
  private StdMath() {}

  /**
   * Returns Euler's number <i>e</i> raised to the power of {@code x}.
   *
   * @param x the exponent to raise <i>e</i> to.
   * @return the value <i>e</i><sup>{@code x}</sup>, where <i>e</i> is the base of the natural
   *         logarithms.
   * @see Math#exp(double)
   */
  public static double exp(double x) {
    return org.apache.commons.math3.util.FastMath.exp(x);
    //return Math.exp(x);
  }
}
