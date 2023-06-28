/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.results;

/**
 * Represent a named results source. Does not support data provision.
 */
public class NullSource extends ImageSource {

  /**
   * Create a new image source.
   *
   * @param name the name
   */
  public NullSource(String name) {
    super(name);
  }

  @Override
  protected boolean openSource() {
    return false;
  }

  @Override
  protected void closeSource() {
    // Nothing to do
  }

  @Override
  protected boolean initialiseSequentialRead() {
    return false;
  }

  @Override
  protected float[] nextRawFrame() {
    return null;
  }

  @Override
  protected float[] getRawFrame(int frame) {
    return null;
  }

  @Override
  public boolean isValid(int frame) {
    return false;
  }
}
