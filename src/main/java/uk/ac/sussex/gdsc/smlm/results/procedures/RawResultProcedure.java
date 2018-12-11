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

package uk.ac.sussex.gdsc.smlm.results.procedures;

import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;

/**
 * Contains functionality to obtain the standard data for results.
 */
//@formatter:off
public class RawResultProcedure extends AbstractResultProcedure implements
        BIXYZResultProcedure,
    IResultProcedure,
    BResultProcedure,
    XYZResultProcedure
//@formatter:on
{
  /** The background. */
  public float[] background;

  /** The intensity. */
  public float[] intensity;

  /** The x. */
  public float[] x;

  /** The y. */
  public float[] y;

  /** The z. */
  public float[] z;

  /**
   * Instantiates a new standard result procedure.
   *
   * @param results the results
   */
  public RawResultProcedure(MemoryPeakResults results) {
    super(results);
  }

  /**
   * Gets the BIXYZ data in the configured units.
   *
   * @throws DataException if conversion to the required units is not possible
   */
  public void getBIXYZ() throws DataException {
    i = 0;
    this.background = allocate(this.background);
    this.intensity = allocate(this.intensity);
    this.x = allocate(this.x);
    this.y = allocate(this.y);
    this.z = allocate(this.z);
    results.forEachNative((BIXYZResultProcedure) this);
  }

  @Override
  public void executeBIXYZ(float background, float intensity, float x, float y, float z) {
    this.background[i] = background;
    this.intensity[i] = intensity;
    this.x[i] = x;
    this.y[i] = y;
    this.z[i] = z;
    i++;
  }

  /**
   * Gets the I data in the configured units.
   *
   * @throws DataException if conversion to the required units is not possible
   */
  public void getI() throws DataException {
    i = 0;
    this.intensity = allocate(this.intensity);
    results.forEachNative((IResultProcedure) this);
  }

  @Override
  public void executeI(float intensity) {
    this.intensity[i] = intensity;
    i++;
  }

  /**
   * Gets the B data in the configured units.
   *
   * @throws DataException if conversion to the required units is not possible
   */
  public void getB() throws DataException {
    i = 0;
    this.background = allocate(this.background);
    results.forEachNative((BResultProcedure) this);
  }

  @Override
  public void executeB(float background) {
    this.background[i] = background;
    i++;
  }

  /**
   * Gets the XYZ data in the configured units.
   *
   * @throws DataException if conversion to the required units is not possible
   */
  public void getXyz() throws DataException {
    i = 0;
    this.x = allocate(this.x);
    this.y = allocate(this.y);
    this.z = allocate(this.z);
    results.forEachNative((XYZResultProcedure) this);
  }

  @Override
  public void executeXYZ(float x, float y, float z) {
    this.x[i] = x;
    this.y[i] = y;
    this.z[i] = z;
    i++;
  }
}
