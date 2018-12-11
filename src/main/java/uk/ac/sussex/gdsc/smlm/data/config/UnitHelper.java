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

package uk.ac.sussex.gdsc.smlm.data.config;

import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;

/**
 * Contains helper functions for the units. Adds functionality to the protocol buffer unit enums.
 */
public class UnitHelper {
  /**
   * Gets the name.
   *
   * @param unit the unit
   * @return the name
   */
  public static String getName(DistanceUnit unit) {
    switch (unit) {
      case NM:
        return "nanometre";
      case PIXEL:
        return "pixel";
      case UM:
        return "micrometre";
      default:
        return "unknown";
    }
  }

  /**
   * Gets the short name.
   *
   * @param unit the unit
   * @return the short name
   */
  public static String getShortName(DistanceUnit unit) {
    switch (unit) {
      case NM:
        return "nm";
      case PIXEL:
        return "px";
      case UM:
        return "μm";
      default:
        return "na";
    }
  }

  /**
   * Guess distance unit from short name.
   *
   * @param name the name
   * @return the distance unit
   */
  public static DistanceUnit guessDistanceUnitFromShortName(String name) {
    if (name != null) {
      if (name.equalsIgnoreCase("px")) {
        return DistanceUnit.PIXEL;
      }
      // There is no uppercase μm
      if (name.equals("μm") || name.equalsIgnoreCase("um")) {
        return DistanceUnit.UM;
      }
      if (name.equalsIgnoreCase("nm")) {
        return DistanceUnit.NM;
      }
    }
    return null;
  }

  /**
   * Gets the name.
   *
   * @param unit the unit
   * @return the name
   */
  public static String getName(IntensityUnit unit) {
    switch (unit) {
      case COUNT:
        return "count";
      case PHOTON:
        return "photon";
      default:
        return "unknown";
    }
  }

  /**
   * Gets the short name.
   *
   * @param unit the unit
   * @return the short name
   */
  public static String getShortName(IntensityUnit unit) {
    switch (unit) {
      case COUNT:
        return "count";
      case PHOTON:
        return "photon";
      default:
        return "na";
    }
  }

  /**
   * Guess intensity unit from short name.
   *
   * @param name the name
   * @return the intensity unit
   */
  public static IntensityUnit guessIntensityUnitFromShortName(String name) {
    if (name != null) {
      if (name.equalsIgnoreCase("photon")) {
        return IntensityUnit.PHOTON;
      }
      if (name.equalsIgnoreCase("count")) {
        return IntensityUnit.COUNT;
      }
    }
    return null;
  }

  /**
   * Gets the name.
   *
   * @param unit the unit
   * @return the name
   */
  public static String getName(AngleUnit unit) {
    switch (unit) {
      case DEGREE:
        return "degree";
      case RADIAN:
        return "radian";
      default:
        return "unknown";
    }
  }

  /**
   * Gets the short name.
   *
   * @param unit the unit
   * @return the short name
   */
  public static String getShortName(AngleUnit unit) {
    switch (unit) {
      case DEGREE:
        return "°";
      case RADIAN:
        return "rad";
      default:
        return "na";
    }
  }

  /**
   * Guess angle unit from short name.
   *
   * @param name the name
   * @return the angle unit
   */
  public static AngleUnit guessAngleUnitFromShortName(String name) {
    if (name != null) {
      if (name.equalsIgnoreCase("°")) {
        return AngleUnit.DEGREE;
      }
      if (name.equalsIgnoreCase("rad")) {
        return AngleUnit.RADIAN;
      }
    }
    return null;
  }

  /**
   * Gets the name.
   *
   * @param unit the unit
   * @return the name
   */
  public static String getName(TimeUnit unit) {
    switch (unit) {
      case FRAME:
        return "frame";
      case MILLISECOND:
        return "millisecond";
      case SECOND:
        return "second";
      default:
        return "unknown";
    }
  }

  /**
   * Gets the short name.
   *
   * @param unit the unit
   * @return the short name
   */
  public static String getShortName(TimeUnit unit) {
    switch (unit) {
      case FRAME:
        return "frame";
      case MILLISECOND:
        return "ms";
      case SECOND:
        return "s";
      default:
        return "na";
    }
  }

  /**
   * Guess time unit from short name.
   *
   * @param name the name
   * @return the time unit
   */
  public static TimeUnit guessTimeUnitFromShortName(String name) {
    if (name != null) {
      if (name.equalsIgnoreCase("frame")) {
        return TimeUnit.FRAME;
      }
      if (name.equalsIgnoreCase("s")) {
        return TimeUnit.SECOND;
      }
      if (name.equalsIgnoreCase("ms")) {
        return TimeUnit.MILLISECOND;
      }
    }
    return null;
  }
}
