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

package uk.ac.sussex.gdsc.smlm.ij.plugins.pcpalm;

/**
 * Used to store all the information required for the PC-PALM analysis.
 */
public class Molecule {
  /** The x. */
  public double x;
  /** The y. */
  public double y;
  /** The precision. */
  public double precision;
  /** The photons. */
  public double photons;

  /** Used to construct a single linked list of molecules. */
  public Molecule next;

  /**
   * Instantiates a new molecule.
   *
   * @param x the x
   * @param y the y
   * @param precision the precision
   * @param photons the photons
   */
  public Molecule(double x, double y, double precision, double photons) {
    this.x = x;
    this.y = y;
    this.precision = precision;
    this.photons = photons;
  }

  /**
   * Get the distance.
   *
   * @param other the other
   * @return the distance
   */
  public double distance(Molecule other) {
    final double dx = x - other.x;
    final double dy = y - other.y;
    return Math.sqrt(dx * dx + dy * dy);
  }

  /**
   * Get the squared distance.
   *
   * @param other the other
   * @return the squared distance
   */
  public double distance2(Molecule other) {
    final double dx = x - other.x;
    final double dy = y - other.y;
    return dx * dx + dy * dy;
  }
}
