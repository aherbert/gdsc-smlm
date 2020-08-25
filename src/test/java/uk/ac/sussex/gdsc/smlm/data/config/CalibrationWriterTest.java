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

package uk.ac.sussex.gdsc.smlm.data.config;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;

@SuppressWarnings({"javadoc"})
class CalibrationWriterTest {
  @SeededTest
  void canWrite(RandomSeed seed) {
    final UniformRandomProvider rng = RngUtils.create(seed.getSeed());
    for (int i = 0; i < 100; i++) {
      canWrite(rng);
    }
  }

  private static void canWrite(UniformRandomProvider rng) {
    final double qe = rng.nextDouble();
    final double bias = 1 + rng.nextDouble();
    final double exposureTime = 1 + rng.nextDouble();
    final double gain = 1 + rng.nextDouble();
    final double nmPerPixel = 1 + rng.nextDouble();
    final double readNoise = 1 + rng.nextDouble();
    final AngleUnit angleUnit = AngleUnit.values()[rng.nextInt(AngleUnit.values().length - 1)];
    final CameraType cameraType = CameraType.values()[rng.nextInt(CameraType.values().length - 1)];
    final DistanceUnit distanceUnit =
        DistanceUnit.values()[rng.nextInt(DistanceUnit.values().length - 1)];
    final IntensityUnit intensityUnit =
        IntensityUnit.values()[rng.nextInt(IntensityUnit.values().length - 1)];

    final CalibrationWriter writer = new CalibrationWriter();

    Assertions.assertEquals(writer.getQuantumEfficiency(), 0);
    Assertions.assertEquals(writer.getBias(), 0);
    Assertions.assertEquals(writer.getExposureTime(), 0);
    Assertions.assertEquals(writer.getCountPerPhoton(), 0);
    Assertions.assertEquals(writer.getNmPerPixel(), 0);
    Assertions.assertEquals(writer.getReadNoise(), 0);
    Assertions.assertFalse(writer.hasQuantumEfficiency());
    Assertions.assertFalse(writer.hasBias());
    Assertions.assertFalse(writer.hasExposureTime());
    Assertions.assertFalse(writer.hasCountPerPhoton());
    Assertions.assertFalse(writer.hasNmPerPixel());
    Assertions.assertFalse(writer.hasReadNoise());
    Assertions.assertEquals(writer.getAngleUnit(), AngleUnit.ANGLE_UNIT_NA);
    Assertions.assertEquals(writer.getCameraType(), CameraType.CAMERA_TYPE_NA);
    Assertions.assertEquals(writer.getDistanceUnit(), DistanceUnit.DISTANCE_UNIT_NA);
    Assertions.assertEquals(writer.getIntensityUnit(), IntensityUnit.INTENSITY_UNIT_NA);

    writer.setQuantumEfficiency(qe);
    writer.setBias(bias);
    writer.setExposureTime(exposureTime);
    writer.setCountPerPhoton(gain);
    writer.setNmPerPixel(nmPerPixel);
    writer.setReadNoise(readNoise);
    writer.setAngleUnit(angleUnit);
    writer.setCameraType(cameraType);
    writer.setDistanceUnit(distanceUnit);
    writer.setIntensityUnit(intensityUnit);

    Assertions.assertEquals(writer.getQuantumEfficiency(), qe);
    Assertions.assertEquals(writer.getBias(), bias);
    Assertions.assertEquals(writer.getExposureTime(), exposureTime);
    Assertions.assertEquals(writer.getCountPerPhoton(), gain);
    Assertions.assertEquals(writer.getNmPerPixel(), nmPerPixel);
    Assertions.assertEquals(writer.getReadNoise(), readNoise);
    Assertions.assertTrue(writer.hasQuantumEfficiency());
    Assertions.assertTrue(writer.hasBias());
    Assertions.assertTrue(writer.hasExposureTime());
    Assertions.assertTrue(writer.hasCountPerPhoton());
    Assertions.assertTrue(writer.hasNmPerPixel());
    Assertions.assertTrue(writer.hasReadNoise());
    Assertions.assertEquals(writer.getAngleUnit(), angleUnit);
    Assertions.assertEquals(writer.getCameraType(), cameraType);
    Assertions.assertEquals(writer.getDistanceUnit(), distanceUnit);
    Assertions.assertEquals(writer.getIntensityUnit(), intensityUnit);

    final CalibrationReader reader = new CalibrationReader(writer.getCalibration());

    Assertions.assertEquals(reader.getQuantumEfficiency(), qe);
    Assertions.assertEquals(reader.getBias(), bias);
    Assertions.assertEquals(reader.getExposureTime(), exposureTime);
    Assertions.assertEquals(reader.getCountPerPhoton(), gain);
    Assertions.assertEquals(reader.getNmPerPixel(), nmPerPixel);
    Assertions.assertEquals(reader.getReadNoise(), readNoise);
    Assertions.assertTrue(reader.hasQuantumEfficiency());
    Assertions.assertTrue(reader.hasBias());
    Assertions.assertTrue(reader.hasExposureTime());
    Assertions.assertTrue(reader.hasCountPerPhoton());
    Assertions.assertTrue(reader.hasNmPerPixel());
    Assertions.assertTrue(reader.hasReadNoise());
    Assertions.assertEquals(reader.getAngleUnit(), angleUnit);
    Assertions.assertEquals(reader.getCameraType(), cameraType);
    Assertions.assertEquals(reader.getDistanceUnit(), distanceUnit);
    Assertions.assertEquals(reader.getIntensityUnit(), intensityUnit);
  }
}
