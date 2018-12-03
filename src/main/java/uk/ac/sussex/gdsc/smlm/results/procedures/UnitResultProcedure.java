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

import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;

/**
 * Contains core functionality to for result procedures.
 */
public abstract class UnitResultProcedure extends AbstractResultProcedure {
  /** The distance unit. */
  private DistanceUnit distanceUnit;

  /** The intensity unit. */
  private IntensityUnit intensityUnit;

  /**
   * Instantiates a new abstract result procedure.
   *
   * @param results the results
   * @param distanceUnit the distance unit
   * @param intensityUnit the intensity unit
   */
  public UnitResultProcedure(MemoryPeakResults results, DistanceUnit distanceUnit,
      IntensityUnit intensityUnit) {
    super(results);
    this.setDistanceUnit(distanceUnit);
    this.setIntensityUnit(intensityUnit);
  }

  /**
   * Instantiates a new abstract result procedure.
   *
   * @param results the results
   * @param distanceUnit the distance unit
   */
  public UnitResultProcedure(MemoryPeakResults results, DistanceUnit distanceUnit) {
    this(results, distanceUnit, IntensityUnit.PHOTON);
  }

  /**
   * Instantiates a new abstract result procedure.
   *
   * @param results the results
   * @param intensityUnit the intensity unit
   */
  public UnitResultProcedure(MemoryPeakResults results, IntensityUnit intensityUnit) {
    this(results, DistanceUnit.PIXEL, intensityUnit);
  }

  /**
   * Instantiates a new abstract result procedure.
   *
   * @param results the results
   */
  public UnitResultProcedure(MemoryPeakResults results) {
    this(results, DistanceUnit.PIXEL, IntensityUnit.PHOTON);
  }

  /**
   * Gets the distance unit.
   *
   * @return the distance unit
   */
  public DistanceUnit getDistanceUnit() {
    return distanceUnit;
  }

  /**
   * Sets the distance unit.
   *
   * @param distanceUnit the new distance unit
   */
  public void setDistanceUnit(DistanceUnit distanceUnit) {
    if (distanceUnit == null) {
      throw new IllegalArgumentException("unit must not be null");
    }
    this.distanceUnit = distanceUnit;
  }

  /**
   * Gets the intensity unit.
   *
   * @return the intensity unit
   */
  public IntensityUnit getIntensityUnit() {
    return intensityUnit;
  }

  /**
   * Sets the intensity unit.
   *
   * @param intensityUnit the new intensity unit
   */
  public void setIntensityUnit(IntensityUnit intensityUnit) {
    if (intensityUnit == null) {
      throw new IllegalArgumentException("unit must not be null");
    }
    this.intensityUnit = intensityUnit;
  }
}
