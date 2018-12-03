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
package uk.ac.sussex.gdsc.smlm.results.filter;

/**
 * Null implementation of the CoordinateStore interface.
 */
public class NullCoordinateStore implements CoordinateStore {
  /** An instance to ignore calls to the CoordinateStore interface. */
  public static final NullCoordinateStore INSTANCE = new NullCoordinateStore();

  /**
   * Instantiates a new null coordinate store.
   */
  private NullCoordinateStore() {}

  /**
   * Creates an instance if the argument is null, else return the argument.
   *
   * @param coordinateStore the coordinate store (may be null)
   * @return the coordinate store (not null)
   */
  public static CoordinateStore replaceIfNull(CoordinateStore coordinateStore) {
    return (coordinateStore == null) ? INSTANCE : coordinateStore;
  }

  /** {@inheritDoc} */
  @Override
  public double getXYResolution() {
    return 0;
  }

  /** {@inheritDoc} */
  @Override
  public double getZResolution() {
    return 0;
  }

  /** {@inheritDoc} */
  @Override
  public void addToQueue(double x, double y, double z) {
    // Do nothing
  }

  /** {@inheritDoc} */
  @Override
  public void flush() {
    // Do nothing
  }

  /** {@inheritDoc} */
  @Override
  public void add(double x, double y, double z) {
    // Do nothing
  }

  /** {@inheritDoc} */
  @Override
  public void clear() {
    // Do nothing
  }

  /** {@inheritDoc} */
  @Override
  public boolean contains(double x, double y, double z) {
    return false;
  }

  /** {@inheritDoc} */
  @Override
  public double[] find(double x, double y, double z) {
    return null;
  }

  /** {@inheritDoc} */
  @Override
  public CoordinateStore newInstance() {
    return this;
  }

  /** {@inheritDoc} */
  @Override
  public CoordinateStore resize(int minx, int miny, int maxx, int maxy) {
    return this;
  }

  /** {@inheritDoc} */
  @Override
  public int getMinX() {
    return 0;
  }

  /** {@inheritDoc} */
  @Override
  public int getMinY() {
    return 0;
  }

  /** {@inheritDoc} */
  @Override
  public int getWidth() {
    return 0;
  }

  /** {@inheritDoc} */
  @Override
  public int getHeight() {
    return 0;
  }
}
