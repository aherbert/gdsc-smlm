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

package uk.ac.sussex.gdsc.smlm.model;

/**
 * Contains methods for sampling the spatial position of a molecule.
 *
 * <p>The centre of the distribution is [0,0,0]. Therefore the coordinates can be negative or
 * positive.
 */
public interface SpatialDistribution {
  /**
   * Get the next position. Note the centre of the distribution is [0,0,0].
   *
   * <p>Note: May return null if no more positions are available.
   *
   * @return The next position [x,y,z]
   */
  public double[] next();

  /**
   * Check if the coordinates are within the distribution bounds.
   *
   * @param xyz the xyz
   * @return True if the coordinates are within the distribution bounds
   */
  public boolean isWithin(double[] xyz);

  /**
   * Check if the coordinates are within the distribution bounds in the XY dimensions. If the
   * distribution is dependent on the Z-dimension (e.g. 3D objects) then this can return the same as
   * the {@link #isWithin(double[])} method.
   *
   * @param xyz the xyz
   * @return True if the coordinates are within the distribution bounds in the XY dimensions
   */
  public boolean isWithinXy(double[] xyz);

  /**
   * Initialise the distribution with a set of coordinates. This can be used before calls to
   * {@link #isWithin(double[])} or {@link #isWithinXy(double[])} if the implementation depends on
   * knowing the original coordinate location.
   *
   * @param xyz the xyz
   */
  public void initialise(double[] xyz);
}
