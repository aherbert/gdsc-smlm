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
package uk.ac.sussex.gdsc.smlm.ij.settings;

import java.io.File;
import java.io.IOException;

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.test.TestSettings;

@SuppressWarnings({ "javadoc" })
public class SettingsManagerTest
{
	@Test
	public void canReadWriteConfiguration() throws IOException
	{
		final RandomGenerator rand = TestSettings.getRandomGenerator();

		final Calibration.Builder builder = Calibration.newBuilder();
		builder.getCameraCalibrationBuilder().setBias(rand.nextDouble());
		final Calibration e = builder.build();
		Calibration o;

		final String dir = SettingsManager.getSettingsDirectory();
		//System.out.println(dir);
		final File tmp = createTempDirectory(false);
		try
		{
			SettingsManager.setSettingsDirectory(tmp.getPath());
			o = SettingsManager.readCalibration(SettingsManager.FLAG_SILENT | SettingsManager.FLAG_NO_DEFAULT);
			Assertions.assertTrue(o == null, "Failed to read null");
			o = SettingsManager.readCalibration(SettingsManager.FLAG_SILENT);
			Assertions.assertTrue(e.getDefaultInstanceForType().equals(o), "Failed to read default");
			Assertions.assertTrue(SettingsManager.writeSettings(e), "Failed to write");
			o = SettingsManager.readCalibration(0);
			Assertions.assertTrue(e.equals(o), "Not equal");
			SettingsManager.clearSettings(e.getClass());
			o = SettingsManager.readCalibration(SettingsManager.FLAG_SILENT | SettingsManager.FLAG_NO_DEFAULT);
			Assertions.assertTrue(o == null, "Failed to clear");
		}
		finally
		{
			// Reset
			tmp.delete(); // Will work if the directory is empty
			SettingsManager.setSettingsDirectory(dir);
		}
	}

	public static File createTempDirectory(boolean create) throws IOException
	{
		final File temp;

		temp = File.createTempFile(SettingsManagerTest.class.getSimpleName(), ".tmp");
		temp.deleteOnExit();

		if (!(temp.delete()))
			throw new IOException("Could not delete temp file: " + temp.getAbsolutePath());

		if (create && !(temp.mkdir()))
			throw new IOException("Could not create temp directory: " + temp.getAbsolutePath());

		//System.out.println(temp.getPath());

		return (temp);
	}
}
