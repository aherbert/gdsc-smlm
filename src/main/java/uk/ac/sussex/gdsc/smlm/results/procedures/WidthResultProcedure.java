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
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;

/**
 * Contains functionality to obtain the standard calibrated data for results.
 */
//@formatter:off
public class WidthResultProcedure extends UnitResultProcedure implements
    WResultProcedure,
    WxWyResultProcedure {
  //@formatter:on

  /** The x width. */
  public float[] wx;

  /** The y width. */
  public float[] wy;

  /**
   * Instantiates a new width result procedure.
   *
   * @param results the results
   * @param distanceUnit the distance unit
   */
  public WidthResultProcedure(MemoryPeakResults results, DistanceUnit distanceUnit) {
    super(results, distanceUnit);
  }

  /**
   * Instantiates a new width result procedure.
   *
   * @param results the results
   */
  public WidthResultProcedure(MemoryPeakResults results) {
    super(results);
  }

  /**
   * Gets the W data in the configured units.
   *
   * @throws DataException if conversion to the required units is not possible
   */
  public void getW() {
    counter = 0;
    this.wx = allocate(this.wx);
    results.forEach(getDistanceUnit(), (WResultProcedure) this);
  }

  @Override
  public void executeW(float width) {
    this.wx[counter++] = width;
  }

  /**
   * Gets the WxWy data in the configured units.
   *
   * @throws DataException if conversion to the required units is not possible
   */
  public void getWxWy() {
    counter = 0;
    this.wx = allocate(this.wx);
    this.wy = allocate(this.wy);
    results.forEach(getDistanceUnit(), (WxWyResultProcedure) this);
  }

  @Override
  public void executeWxWy(float wx, float wy) {
    this.wx[counter] = wx;
    this.wy[counter] = wy;
    counter++;
  }
}
