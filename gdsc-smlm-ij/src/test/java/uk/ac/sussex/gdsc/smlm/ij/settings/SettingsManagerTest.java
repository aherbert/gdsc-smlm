/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.ij.settings;

import java.io.File;
import java.io.IOException;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;

@SuppressWarnings({"javadoc"})
class SettingsManagerTest {
  @SeededTest
  void canReadWriteConfiguration(RandomSeed seed) throws IOException {
    final UniformRandomProvider rand = RngUtils.create(seed.get());

    final Calibration.Builder builder = Calibration.newBuilder();
    builder.getCameraCalibrationBuilder().setBias(rand.nextDouble());
    final Calibration exp = builder.build();

    final String dir = SettingsManager.getSettingsDirectory();
    // System.out.println(dir);
    final File tmp = createTempDirectory(false);
    try {
      SettingsManager.setSettingsDirectory(tmp.getPath());
      Calibration obs = SettingsManager
          .readCalibration(SettingsManager.FLAG_SILENT | SettingsManager.FLAG_NO_DEFAULT);
      Assertions.assertTrue(obs == null, "Failed to read null");
      obs = SettingsManager.readCalibration(SettingsManager.FLAG_SILENT);
      Assertions.assertTrue(exp.getDefaultInstanceForType().equals(obs), "Failed to read default");
      Assertions.assertTrue(SettingsManager.writeSettings(exp), "Failed to write");
      obs = SettingsManager.readCalibration(0);
      Assertions.assertTrue(exp.equals(obs), "Not equal");
      SettingsManager.clearSettings(exp.getClass());
      obs = SettingsManager
          .readCalibration(SettingsManager.FLAG_SILENT | SettingsManager.FLAG_NO_DEFAULT);
      Assertions.assertTrue(obs == null, "Failed to clear");
    } finally {
      // Reset
      tmp.delete(); // Will work if the directory is empty
      SettingsManager.setSettingsDirectory(dir);
    }
  }

  private static File createTempDirectory(boolean create) throws IOException {
    final File temp;

    temp = File.createTempFile(SettingsManagerTest.class.getSimpleName(), ".tmp");
    temp.deleteOnExit();

    if (!(temp.delete())) {
      throw new IOException("Could not delete temp file: " + temp.getAbsolutePath());
    }

    if (create && !(temp.mkdir())) {
      throw new IOException("Could not create temp directory: " + temp.getAbsolutePath());
    }

    // System.out.println(temp.getPath());

    return (temp);
  }
}
