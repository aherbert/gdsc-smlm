package gdsc.smlm.ij.settings;

import java.io.File;
import java.io.IOException;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.smlm.data.config.CalibrationConfig.Calibration;

public class SettingsManagerTest
{
	@Test
	public void canReadWriteConfiguration() throws IOException
	{
		RandomGenerator rand = new Well19937c(30051977);

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
