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

import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.AngleCalibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.DistanceCalibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.IntensityCalibration;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

@SuppressWarnings({"javadoc"})
public class CalibrationReaderTest {
  double nmPerPixel = 104.5;
  double gain = 45;
  double bias = 100;

  @Test
  public void canGetDistanceConverter() {
    final Calibration.Builder builder = Calibration.newBuilder();
    final DistanceCalibration.Builder distanceBuilder = builder.getDistanceCalibrationBuilder();

    distanceBuilder.setNmPerPixel(nmPerPixel);
    distanceBuilder.setDistanceUnit(DistanceUnit.PIXEL);

    final Calibration c = builder.build();

    final CalibrationReader reader = new CalibrationReader(c);

    final TypeConverter<DistanceUnit> distanceConverter =
        reader.getDistanceConverter(DistanceUnit.NM);
    Assertions.assertEquals(distanceConverter.from(), DistanceUnit.PIXEL);
    Assertions.assertEquals(distanceConverter.to(), DistanceUnit.NM);

    final TypeConverter<DistanceUnit> distanceConverter2 =
        CalibrationHelper.getDistanceConverter(c, DistanceUnit.NM);
    Assertions.assertEquals(distanceConverter2.from(), DistanceUnit.PIXEL);
    Assertions.assertEquals(distanceConverter2.to(), DistanceUnit.NM);

    Assertions.assertEquals(distanceConverter.getFunction(), distanceConverter2.getFunction());
  }

  @Test
  public void canGetIntensityConverter() {
    final Calibration.Builder builder = Calibration.newBuilder();
    final IntensityCalibration.Builder intensityBuilder = builder.getIntensityCalibrationBuilder();

    intensityBuilder.setIntensityUnit(IntensityUnit.PHOTON);
    intensityBuilder.setCountPerPhoton(gain);

    final Calibration c = builder.build();

    final CalibrationReader reader = new CalibrationReader(c);

    final TypeConverter<IntensityUnit> intensityConverter =
        reader.getIntensityConverter(IntensityUnit.COUNT);
    Assertions.assertEquals(intensityConverter.from(), IntensityUnit.PHOTON);
    Assertions.assertEquals(intensityConverter.to(), IntensityUnit.COUNT);

    final TypeConverter<IntensityUnit> intensityConverter2 =
        CalibrationHelper.getIntensityConverter(c, IntensityUnit.COUNT);
    Assertions.assertEquals(intensityConverter2.from(), IntensityUnit.PHOTON);
    Assertions.assertEquals(intensityConverter2.to(), IntensityUnit.COUNT);

    Assertions.assertEquals(intensityConverter.getFunction(), intensityConverter2.getFunction());
  }

  @Test
  public void canGetAngleConverter() {
    final Calibration.Builder builder = Calibration.newBuilder();
    final AngleCalibration.Builder psfBuilder = builder.getAngleCalibrationBuilder();

    psfBuilder.setAngleUnit(AngleUnit.RADIAN);

    final Calibration c = builder.build();

    final CalibrationReader reader = new CalibrationReader(c);

    final TypeConverter<AngleUnit> angleConverter = reader.getAngleConverter(AngleUnit.DEGREE);
    Assertions.assertEquals(angleConverter.from(), AngleUnit.RADIAN);
    Assertions.assertEquals(angleConverter.to(), AngleUnit.DEGREE);

    final TypeConverter<AngleUnit> angleConverter2 =
        CalibrationHelper.getAngleConverter(c, AngleUnit.DEGREE);
    Assertions.assertEquals(angleConverter2.from(), AngleUnit.RADIAN);
    Assertions.assertEquals(angleConverter2.to(), AngleUnit.DEGREE);

    Assertions.assertEquals(angleConverter.getFunction(), angleConverter2.getFunction());
  }
}
