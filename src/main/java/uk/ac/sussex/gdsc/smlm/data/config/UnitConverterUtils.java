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

import uk.ac.sussex.gdsc.core.data.utils.AddMultiplyTypeConverter;
import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.data.utils.IdentityTypeConverter;
import uk.ac.sussex.gdsc.core.data.utils.MultiplyAddTypeConverter;
import uk.ac.sussex.gdsc.core.data.utils.MultiplyTypeConverter;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;

/**
 * Factory for creating unit converters.
 */
public final class UnitConverterUtils {

  /** No public constructor. */
  private UnitConverterUtils() {}

  /**
   * Creates the {@link DistanceUnit} converter.
   *
   * @param from the unit to convert from
   * @param to the unit to convert to
   * @param nmPerPixel the nm per pixel
   * @return the unit converter
   * @throws ConversionException if a converter cannot be created
   */
  public static TypeConverter<DistanceUnit> createConverter(DistanceUnit from, DistanceUnit to,
      double nmPerPixel) {
    if (from == to) {
      return new IdentityTypeConverter<>(from);
    }
    checkNotNull(from, to);

    switch (from) {
      case NM:
        switch (to) {
          case PIXEL:
            return new MultiplyTypeConverter<>(from, to, 1.0 / checkNmPerPixel(nmPerPixel));
          case UM:
            return new MultiplyTypeConverter<>(from, to, 1e-3);
          default:
            throw new ConversionException(from + " to " + to);
        }

      case PIXEL:
        switch (to) {
          case NM:
            return new MultiplyTypeConverter<>(from, to, checkNmPerPixel(nmPerPixel));
          case UM:
            return new MultiplyTypeConverter<>(from, to, checkNmPerPixel(nmPerPixel) / 1e3);
          default:
            throw new ConversionException(from + " to " + to);
        }

      case UM:
        switch (to) {
          case NM:
            return new MultiplyTypeConverter<>(from, to, 1e3);
          case PIXEL:
            return new MultiplyTypeConverter<>(from, to, 1e3 / checkNmPerPixel(nmPerPixel));
          default:
            throw new ConversionException(from + " to " + to);
        }

      default:
        throw new ConversionException(from + " to " + to);
    }
  }

  /**
   * Creates the {@link TimeUnit} converter.
   *
   * @param from the unit to convert from
   * @param to the unit to convert to
   * @param msPerFrame the ms per frame
   * @return the unit converter
   * @throws ConversionException if a converter cannot be created
   */
  public static TypeConverter<TimeUnit> createConverter(TimeUnit from, TimeUnit to,
      double msPerFrame) {
    if (from == to) {
      return new IdentityTypeConverter<>(from);
    }
    checkNotNull(from, to);

    switch (from) {
      case MILLISECOND:
        switch (to) {
          case FRAME:
            return new MultiplyTypeConverter<>(from, to, 1.0 / checkMsPerFrame(msPerFrame));
          case SECOND:
            return new MultiplyTypeConverter<>(from, to, 1e-3);
          default:
            throw new ConversionException(from + " to " + to);
        }

      case FRAME:
        switch (to) {
          case MILLISECOND:
            return new MultiplyTypeConverter<>(from, to, checkMsPerFrame(msPerFrame));
          case SECOND:
            return new MultiplyTypeConverter<>(from, to, checkMsPerFrame(msPerFrame) / 1e3);
          default:
            throw new ConversionException(from + " to " + to);
        }

      case SECOND:
        switch (to) {
          case MILLISECOND:
            return new MultiplyTypeConverter<>(from, to, 1e3);
          case FRAME:
            return new MultiplyTypeConverter<>(from, to, 1e3 / checkMsPerFrame(msPerFrame));
          default:
            throw new ConversionException(from + " to " + to);
        }

      default:
        throw new ConversionException(from + " to " + to);
    }
  }

  /**
   * Creates the {@link AngleUnit} converter.
   *
   * @param from the unit to convert from
   * @param to the unit to convert to
   * @return the unit converter
   * @throws ConversionException if a converter cannot be created
   */
  public static TypeConverter<AngleUnit> createConverter(AngleUnit from, AngleUnit to) {
    if (from == to) {
      return new IdentityTypeConverter<>(from);
    }
    checkNotNull(from, to);

    switch (from) {
      case RADIAN:
        if (to == AngleUnit.DEGREE) {
          return new MultiplyTypeConverter<>(from, to, 180.0 / Math.PI);
        }
        throw new ConversionException(from + " to " + to);

      case DEGREE:
        if (to == AngleUnit.RADIAN) {
          return new MultiplyTypeConverter<>(from, to, Math.PI / 180.0);
        }
        throw new ConversionException(from + " to " + to);

      default:
        throw new ConversionException(from + " to " + to);
    }
  }

  /**
   * Creates the {@link IntensityUnit} converter.
   *
   * @param from the unit to convert from
   * @param to the unit to convert to
   * @param offset the offset
   * @param countPerPhoton the count per photon
   * @return the unit converter
   * @throws ConversionException if a converter cannot be created
   */
  public static TypeConverter<IntensityUnit> createConverter(IntensityUnit from, IntensityUnit to,
      double offset, double countPerPhoton) {
    if (from == to) {
      return new IdentityTypeConverter<>(from);
    }
    checkNotNull(from, to);

    switch (from) {
      case COUNT:
        if (to == IntensityUnit.PHOTON) {
          return new AddMultiplyTypeConverter<>(from, to, checkOffset(offset, -1),
              1.0 / checkCountPerPhoton(countPerPhoton));
        }
        throw new ConversionException(from + " to " + to);

      case PHOTON:
        if (to == IntensityUnit.COUNT) {
          return new MultiplyAddTypeConverter<>(from, to, checkCountPerPhoton(countPerPhoton),
              checkOffset(offset, 1));
        }
        throw new ConversionException(from + " to " + to);

      default:
        throw new ConversionException(from + " to " + to);
    }
  }

  /**
   * Creates the {@link IntensityUnit} converter.
   *
   * @param from the unit to convert from
   * @param to the unit to convert to
   * @param countPerPhoton the count per photon
   * @return the unit converter
   * @throws ConversionException if a converter cannot be created
   */
  public static TypeConverter<IntensityUnit> createConverter(IntensityUnit from, IntensityUnit to,
      double countPerPhoton) {
    if (from == to) {
      return new IdentityTypeConverter<>(from);
    }
    checkNotNull(from, to);

    switch (from) {
      case COUNT:
        if (to == IntensityUnit.PHOTON) {
          return new MultiplyTypeConverter<>(from, to, 1.0 / checkCountPerPhoton(countPerPhoton));
        }
        throw new ConversionException(from + " to " + to);

      case PHOTON:
        if (to == IntensityUnit.COUNT) {
          return new MultiplyTypeConverter<>(from, to, checkCountPerPhoton(countPerPhoton));
        }
        throw new ConversionException(from + " to " + to);

      default:
        throw new ConversionException(from + " to " + to);
    }
  }

  /**
   * Check the from and to units are not null.
   *
   * @param from the unit to convert from
   * @param to the unit to convert to
   * @throws ConversionException if any arguments are null
   */
  private static void checkNotNull(Object from, Object to) {
    if (from == null) {
      throw new ConversionException("from must not be null");
    }
    if (to == null) {
      throw new ConversionException("to must not be null");
    }
  }

  private static double checkNmPerPixel(double nmPerPixel) {
    if (!(nmPerPixel > 0 && nmPerPixel <= java.lang.Double.MAX_VALUE)) {
      throw new ConversionException("nm/pixel must be positive");
    }
    return nmPerPixel;
  }

  private static double checkMsPerFrame(double msPerFrame) {
    if (!(msPerFrame > 0 && msPerFrame <= java.lang.Double.MAX_VALUE)) {
      throw new ConversionException("ms/frame must be positive");
    }
    return msPerFrame;
  }

  private static double checkOffset(double offset, int sign) {
    if (!(offset >= 0 && offset <= java.lang.Double.MAX_VALUE)) {
      throw new ConversionException("count/photon must be positive");
    }
    return offset * sign;
  }

  private static double checkCountPerPhoton(double countPerPhoton) {
    if (!(countPerPhoton > 0 && countPerPhoton <= java.lang.Double.MAX_VALUE)) {
      throw new ConversionException("count/photon must be positive");
    }
    return countPerPhoton;
  }
}
