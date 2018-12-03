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
package uk.ac.sussex.gdsc.smlm.ij.ij3d;

import org.scijava.vecmath.Color4f;

/**
 * Interface for shape objects that represent a set of items to have transparent colours.
 */
public interface TransparentItemShape extends ItemShape {
  /**
   * Sets the color for each item.
   *
   * @param color the new color
   * @throws IllegalArgumentException if the number of colours is incorrect
   */
  public void setItemColor4(final Color4f[] color) throws IllegalArgumentException;

  /**
   * Sets the alpha for each item. 1 is opaque and 0 is fully transparent.
   *
   * @param alpha the new alpha
   * @throws IllegalArgumentException if the number of colours is incorrect
   */
  public void setItemAlpha(final float[] alpha) throws IllegalArgumentException;

  /**
   * Sets the alpha for each item. 1 is opaque and 0 is fully transparent.
   *
   * @param alpha the new alpha
   * @throws IllegalArgumentException if the number of colours is incorrect
   */
  public void setItemAlpha(final float alpha) throws IllegalArgumentException;

  /**
   * Gets the alpha for each item.
   *
   * @param alpha the new alpha
   * @throws IllegalArgumentException if the number of colours is incorrect
   */
  public void getItemAlpha(final float[] alpha) throws IllegalArgumentException;
}
