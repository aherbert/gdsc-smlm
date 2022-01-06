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

package uk.ac.sussex.gdsc.smlm.function.gaussian;

/**
 * Implements a no astigmatism model of a 2D Gaussian function. The width is constant.
 */
public class NullAstigmatismZModel implements AstigmatismZModel {
  /** The width in x. */
  public final double sx;
  /** The width in y. */
  public final double sy;

  /**
   * Instantiates a new null astigmatism Z model.
   *
   * @param sx the sx
   * @param sy the sy
   * @throws IllegalArgumentException if the widths are not positive
   */
  public NullAstigmatismZModel(double sx, double sy) {
    if (!(sx > 0 && sy > 0)) {
      throw new IllegalArgumentException("Width must be positive");
    }
    this.sx = sx;
    this.sy = sy;
  }

  @Override
  public double getSx(double z) {
    return sx;
  }

  @Override
  public double getSx(double z, double[] dsdz) {
    dsdz[0] = 0;
    return sx;
  }

  @Override
  public double getSx2(double z, double[] dsdz) {
    dsdz[0] = 0;
    dsdz[1] = 0;
    return sx;
  }

  @Override
  public double getSy(double z) {
    return sy;
  }

  @Override
  public double getSy(double z, double[] dsdz) {
    dsdz[0] = 0;
    return sy;
  }

  @Override
  public double getSy2(double z, double[] dsdz) {
    dsdz[0] = 0;
    dsdz[1] = 0;
    return sy;
  }
}
