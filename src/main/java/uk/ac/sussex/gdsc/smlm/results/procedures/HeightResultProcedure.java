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
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;

/**
 * Contains functionality to obtain the standard calibrated data for results.
 */
public class HeightResultProcedure extends UnitResultProcedure implements HResultProcedure {
  /** The height. */
  public float[] h;

  /**
   * Instantiates a new width result procedure.
   *
   * @param results the results
   * @param intensityUnit the intensity unit
   */
  public HeightResultProcedure(MemoryPeakResults results, IntensityUnit intensityUnit) {
    super(results, intensityUnit);
  }

  /**
   * Instantiates a new width result procedure.
   *
   * @param results the results
   */
  public HeightResultProcedure(MemoryPeakResults results) {
    super(results);
  }

  /**
   * Gets the height data in the configured units.
   *
   * @throws DataException if conversion to the required units is not possible
   */
  public void getH() throws DataException {
    i = 0;
    this.h = allocate(this.h);
    results.forEach(getIntensityUnit(), this);
  }

  @Override
  public void executeH(float a) {
    this.h[i++] = a;
  }
}
