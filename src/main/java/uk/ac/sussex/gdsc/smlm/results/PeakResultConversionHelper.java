/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

import uk.ac.sussex.gdsc.core.data.utils.Converter;
import uk.ac.sussex.gdsc.core.data.utils.IdentityTypeConverter;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationHelper;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.ConfigurationException;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFParameter;
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;

/**
 * Helper class for converting peak results.
 */
public class PeakResultConversionHelper {
  private final Calibration calibration;
  private final PSF psf;

  /**
   * Instantiates a new peak result conversion helper.
   *
   * @param calibration the calibration
   * @param psf the psf
   */
  public PeakResultConversionHelper(Calibration calibration, PSF psf) {
    this.calibration = calibration;
    this.psf = psf;
  }

  /** The intensity unit. */
  private IntensityUnit intensityUnit;

  /** The intensity converter. */
  private TypeConverter<IntensityUnit> intensityConverter;

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
    intensityConverter = null;
    this.intensityUnit = intensityUnit;
  }

  /**
   * Checks for intensity unit.
   *
   * @return true, if successful
   */
  public boolean hasIntensityUnit() {
    return intensityUnit != null;
  }

  /**
   * Checks for intensity converter.
   *
   * @return true, if successful
   */
  public boolean hasIntensityConverter() {
    return intensityConverter != null && intensityConverter.to() != null;
  }

  /**
   * Gets the intensity converter for the configured units. If the calibration is null then an
   * identity converter is returned.
   *
   * @return the intensity converter
   */
  public TypeConverter<IntensityUnit> getIntensityConverter() {
    if (intensityConverter == null) {
      intensityConverter = (calibration == null) ? new IdentityTypeConverter<>(null)
          : CalibrationHelper.getIntensityConverterSafe(calibration, intensityUnit);
    }
    return intensityConverter;
  }

  /** The distance unit. */
  private DistanceUnit distanceUnit;

  /** The distance converter. */
  private TypeConverter<DistanceUnit> distanceConverter;

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
    distanceConverter = null;
    this.distanceUnit = distanceUnit;
  }

  /**
   * Checks for distance unit.
   *
   * @return true, if successful
   */
  public boolean hasDistanceUnit() {
    return distanceUnit != null;
  }

  /**
   * Checks for distance converter.
   *
   * @return true, if successful
   */
  public boolean hasDistanceConverter() {
    return distanceConverter != null && distanceConverter.to() != null;
  }

  /**
   * Gets the distance converter for the configured units. If the calibration is null then an
   * identity converter is returned.
   *
   * @return the distance converter
   */
  public TypeConverter<DistanceUnit> getDistanceConverter() {
    if (distanceConverter == null) {
      distanceConverter = (calibration == null) ? new IdentityTypeConverter<>(null)
          : CalibrationHelper.getDistanceConverterSafe(calibration, distanceUnit);
    }
    return distanceConverter;
  }

  /** The angle unit. */
  private AngleUnit angleUnit;

  /** The angle converter. */
  private TypeConverter<AngleUnit> angleConverter;

  /**
   * Gets the angle unit.
   *
   * @return the angle unit
   */
  public AngleUnit getAngleUnit() {
    return angleUnit;
  }

  /**
   * Sets the angle unit.
   *
   * @param angleUnit the new angle unit
   */
  public void setAngleUnit(AngleUnit angleUnit) {
    angleConverter = null;
    this.angleUnit = angleUnit;
  }

  /**
   * Checks for angle unit.
   *
   * @return true, if successful
   */
  public boolean hasAngleUnit() {
    return angleUnit != null;
  }

  /**
   * Checks for angle converter.
   *
   * @return true, if successful
   */
  public boolean hasAngleConverter() {
    return angleConverter != null && angleConverter.to() != null;
  }

  /**
   * Gets the angle converter for the configured units. If the calibration is null then an identity
   * converter is returned.
   *
   * @return the angle converter
   */
  public TypeConverter<AngleUnit> getAngleConverter() {
    if (angleConverter == null) {
      angleConverter = (calibration == null) ? new IdentityTypeConverter<>(null)
          : CalibrationHelper.getAngleConverterSafe(calibration, angleUnit);
    }
    return angleConverter;
  }

  /**
   * Gets the converters for the peak results parameters. This includes the standard parameters and
   * any additional parameters defined in the PSF. If a parameter unit type is undefined then an
   * identity converter is created.
   *
   * @return the converters
   */
  public Converter[] getConverters() {
    final LocalList<Converter> list = new LocalList<>(5);

    getIntensityConverter();
    getDistanceConverter();

    list.add(intensityConverter);
    list.add(intensityConverter);
    list.add(distanceConverter);
    list.add(distanceConverter);
    list.add(distanceConverter);
    if (psf != null) {
      try {
        for (final PSFParameter p : PsfHelper.getParameters(psf)) {
          switch (p.getUnit()) {
            case DISTANCE:
              list.add(distanceConverter);
              break;
            case INTENSITY:
              list.add(intensityConverter);
              break;
            case ANGLE:
              list.add(getAngleConverter());
              break;
            default:
              list.add(new IdentityTypeConverter<>(p.getUnit()));
          }
        }
      } catch (final ConfigurationException ex) {
        // Ignore
      }
    }
    return list.toArray(new Converter[0]);
  }

  /**
   * Gets the names for the peak results parameters. This includes the standard parameters and any
   * additional parameters defined in the PSF. If a parameter name is undefined then unknown is
   * returned.
   *
   * @return the converters
   */
  public String[] getNames() {
    final LocalList<String> list = new LocalList<>(5);

    list.add("Background");
    list.add("Intensity");
    list.add("X");
    list.add("Y");
    list.add("Z");
    if (psf != null) {
      try {
        for (final PSFParameter p : PsfHelper.getParameters(psf)) {
          final String name = p.getName();
          list.add(TextUtils.isNullOrEmpty(name) ? "unknown" : name);
        }
      } catch (final ConfigurationException ex) {
        // Ignore
      }
    }
    return list.toArray(new String[0]);
  }

  /**
   * Gets the unit names for the peak results parameters. This includes the standard parameters and
   * any additional parameters defined in the PSF. If a parameter unit is undefined then an empty
   * string is returned.
   *
   * @return the converters
   */
  public String[] getUnitNames() {
    final LocalList<String> list = new LocalList<>(5);

    getIntensityConverter();
    getDistanceConverter();

    final String safeIntensityUnit =
        (intensityConverter.to() != null) ? UnitHelper.getShortName(intensityConverter.to()) : "";
    final String safeDistanceUnit =
        (distanceConverter.to() != null) ? UnitHelper.getShortName(distanceConverter.to()) : "";
    String safeAngleUnit = null;

    list.add(safeIntensityUnit);
    list.add(safeIntensityUnit);
    list.add(safeDistanceUnit);
    list.add(safeDistanceUnit);
    list.add(safeDistanceUnit);
    if (psf != null) {
      try {
        for (final PSFParameter p : PsfHelper.getParameters(psf)) {
          switch (p.getUnit()) {
            case DISTANCE:
              list.add(safeDistanceUnit);
              break;
            case INTENSITY:
              list.add(safeIntensityUnit);
              break;
            case ANGLE:
              safeAngleUnit = getOrCreateAngleUnit(safeAngleUnit);
              list.add(safeAngleUnit);
              break;
            default:
              list.add("");
          }
        }
      } catch (final ConfigurationException ex) {
        // Ignore
      }
    }
    return list.toArray(new String[0]);
  }

  private String getOrCreateAngleUnit(String angleUnit) {
    if (angleUnit == null) {
      getAngleConverter();
      angleUnit = (angleConverter.to() != null) ? UnitHelper.getShortName(angleConverter.to()) : "";
    }
    return angleUnit;
  }

  /**
   * Test if the current converters will cause the calibration to be changed, i.e. new units.
   *
   * @return true, if successful
   */
  public boolean isCalibrationChanged() {
    if (calibration == null) {
      return false;
    }

    if (hasIntensityConverter() && getIntensityConverter().from() != getIntensityConverter().to()) {
      return true;
    }
    if (hasDistanceConverter() && getDistanceConverter().from() != getDistanceConverter().to()) {
      return true;
    }
    return (hasAngleConverter() && getAngleConverter().from() != getAngleConverter().to());
  }

  /**
   * Test if the current converters will changed to desired units.
   *
   * @return true, if successful
   */
  public boolean isValidConversion() {
    boolean bad = false;

    // The conversion is bad if the output unit is specified and either:
    // there is no converter; or the converter will output the wrong units.

    //@formatter:off
    bad |= (hasIntensityUnit()
        && (!hasIntensityConverter() || intensityUnit != getIntensityConverter().to()));
    bad |= (hasDistanceUnit()
        && (!hasDistanceConverter() || distanceUnit != getDistanceConverter().to()));
    bad |= (hasAngleUnit()
        && (!hasAngleConverter() || angleUnit != getAngleConverter().to()));
    //@formatter:on

    return !bad;
  }

  /**
   * Gets the calibration, updated with the current output units of the converters.
   *
   * @return the calibration
   */
  public Calibration getCalibration() {
    if (calibration == null) {
      return null;
    }

    final CalibrationWriter calibrationWriter = new CalibrationWriter(this.calibration);
    if (hasIntensityConverter()) {
      calibrationWriter.setIntensityUnit(getIntensityConverter().to());
    }
    if (hasDistanceConverter()) {
      calibrationWriter.setDistanceUnit(getDistanceConverter().to());
    }
    if (hasAngleConverter()) {
      calibrationWriter.setAngleUnit(getAngleConverter().to());
    }

    return calibrationWriter.getCalibration();
  }
}
