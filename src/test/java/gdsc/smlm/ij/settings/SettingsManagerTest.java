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
package gdsc.smlm.ij.settings;

import java.io.File;
import java.io.IOException;

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import gdsc.smlm.data.config.CalibrationProtos.Calibration;
import gdsc.test.TestSettings;

@SuppressWarnings({ "javadoc" })
public class SettingsManagerTest
{
	@Test
	public void canReadWriteConfiguration() throws IOException
	{
		RandomGenerator rand = TestSettings.getRandomGenerator();

		Calibration.Builder builder = Calibration.newBuilder();
		builder.getCameraCalibrationBuilder().setBias(rand.nextDouble());
		Calibration e = builder.build();
		Calibration o;

		String dir = SettingsManager.getSettingsDirectory();
		//System.out.println(dir);
		File tmp = createTempDirectory(false);
		try
		{
			SettingsManager.setSettingsDirectory(tmp.getPath());
			o = SettingsManager.readCalibration(SettingsManager.FLAG_SILENT | SettingsManager.FLAG_NO_DEFAULT);
			Assert.assertTrue("Failed to read null", o == null);
			o = SettingsManager.readCalibration(SettingsManager.FLAG_SILENT);
			Assert.assertTrue("Failed to read default", e.getDefaultInstanceForType().equals(o));
			Assert.assertTrue("Failed to write", SettingsManager.writeSettings(e));
			o = SettingsManager.readCalibration(0);
			Assert.assertTrue("Not equal", e.equals(o));
			SettingsManager.clearSettings(e.getClass());
			o = SettingsManager.readCalibration(SettingsManager.FLAG_SILENT | SettingsManager.FLAG_NO_DEFAULT);
			Assert.assertTrue("Failed to clear", o == null);
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
		{
			throw new IOException("Could not delete temp file: " + temp.getAbsolutePath());
		}

		if (create && !(temp.mkdir()))
		{
			throw new IOException("Could not create temp directory: " + temp.getAbsolutePath());
		}

		//System.out.println(temp.getPath());

		return (temp);
	}
}
