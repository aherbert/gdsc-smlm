package gdsc.smlm.ij.plugins;

import java.awt.Color;
import java.awt.Component;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Label;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.regex.Pattern;

import org.apache.commons.math3.util.FastMath;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.core.ij.Utils;
import gdsc.core.utils.UnicodeReader;
import gdsc.smlm.data.config.SMLMSettings.DistanceUnit;
import gdsc.smlm.data.config.SMLMSettings.IntensityUnit;
import gdsc.smlm.data.config.SMLMSettings.TimeUnit;
import gdsc.smlm.data.config.UnitConverterFactory;
import gdsc.smlm.data.utils.TypeConverter;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.ij.settings.CreateDataSettings;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.AttributePeakResult;
import gdsc.smlm.results.Calibration;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

/**
 * Loads generic localisation files into memory
 */
public class LoadLocalisations implements PlugIn
{
	// Time units for the exposure time cannot be in frames as this makes no sense
	private static EnumSet<TimeUnit> set = EnumSet.allOf(TimeUnit.class);
	static
	{
		set.remove(TimeUnit.FRAME);
	}

	public static class Localisation
	{
		int t, id;
		float x, y, z, intensity, sx = -1, sy = -1, precision = -1;
	}

	public static class LocalisationList extends ArrayList<Localisation>
	{
		private static final long serialVersionUID = 6616011992365324247L;

		public final TimeUnit timeUnit;
		public final DistanceUnit distanceUnit;
		public final IntensityUnit intensityUnit;
		public final double gain;
		public final double pixelPitch;
		public final double exposureTime;

		public LocalisationList(TimeUnit timeUnit, DistanceUnit distanceUnit, IntensityUnit intensityUnit, double gain,
				double pixelPitch, double exposureTime)
		{
			this.timeUnit = timeUnit;
			this.distanceUnit = distanceUnit;
			this.intensityUnit = intensityUnit;
			this.gain = gain;
			this.pixelPitch = pixelPitch;
			this.exposureTime = exposureTime;
		}

		private LocalisationList(int timeUnit, int distanceUnit, int intensityUnit, double gain, double pixelPitch,
				double exposureTime)
		{
			this((TimeUnit) set.toArray()[timeUnit], DistanceUnit.values()[distanceUnit],
					IntensityUnit.values()[intensityUnit], gain, pixelPitch, exposureTime);
		}

		public MemoryPeakResults toPeakResults()
		{
			// Convert exposure time to milliseconds
			TypeConverter<TimeUnit> timeConverter = UnitConverterFactory.createConverter(timeUnit, TimeUnit.MILLISECOND,
					1);
			// Convert precision to nm
			TypeConverter<DistanceUnit> distanceConverter = UnitConverterFactory.createConverter(distanceUnit,
					DistanceUnit.NM, pixelPitch);

			MemoryPeakResults results = new MemoryPeakResults();
			results.setName(name);
			Calibration calibration = new Calibration(pixelPitch, gain, timeConverter.convert(exposureTime));
			calibration.setDistanceUnit(distanceUnit);
			calibration.setIntensityUnit(intensityUnit);
			results.setCalibration(calibration);

			for (int i = 0; i < size(); i++)
			{
				final Localisation l = get(i);
				final float[] params = new float[7];
				if (l.intensity <= 0)
					params[Gaussian2DFunction.SIGNAL] = 1;
				else
					params[Gaussian2DFunction.SIGNAL] = (float) (l.intensity);
				params[Gaussian2DFunction.X_POSITION] = (float) (l.x);
				params[Gaussian2DFunction.Y_POSITION] = (float) (l.y);
				// We may not have read in the widths
				if (l.sx == -1)
					params[Gaussian2DFunction.X_SD] = 1;
				else
					params[Gaussian2DFunction.X_SD] = (float) (l.sx);
				if (l.sy == -1)
					params[Gaussian2DFunction.Y_SD] = 1;
				else
					params[Gaussian2DFunction.Y_SD] = (float) (l.sy);
				// Store the z-position in the error field.
				// Q. Should this be converted?
				AttributePeakResult peakResult = new AttributePeakResult(l.t,
						(int) params[Gaussian2DFunction.X_POSITION], (int) params[Gaussian2DFunction.Y_POSITION], 0,
						l.z, 0, params, null);
				peakResult.setId(l.id);
				// Convert to nm
				peakResult.setPrecision(distanceConverter.convert(l.precision));
				results.add(peakResult);
			}

			// Convert to preferred units
			results.convertToPreferenceUnits();

			return results;
		}
	}

	private static final String TITLE = "Load Localisations";
	private static boolean limitZ = false;
	private boolean myLimitZ = false;
	private static double minz = -5;
	private static double maxz = 5;

	private static int it = 0;
	private static int iid = -1;
	private static int ix = 1;
	private static int iy = 2;
	private static int iz = -1;
	private static int ii = 3;
	private static int isx = -1;
	private static int isy = -1;
	private static int ip = -1;
	private static int header = 1;
	private static String comment = "#";
	private static String delimiter = "\\s+";
	private static String name = "Localisations";
	private static int timeUnit = 0;
	private static int distanceUnit = 0;
	private static int intensityUnit = 0;
	private static double gain;
	private static double pixelPitch;
	private static double exposureTime;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		GlobalSettings globalSettings = SettingsManager.loadSettings();
		CreateDataSettings settings = globalSettings.getCreateDataSettings();

		String[] path = Utils.decodePath(settings.localisationsFilename);
		OpenDialog chooser = new OpenDialog("Localisations_File", path[0], path[1]);
		if (chooser.getFileName() == null)
			return;

		settings.localisationsFilename = chooser.getDirectory() + chooser.getFileName();
		SettingsManager.saveSettings(globalSettings);

		LocalisationList localisations = loadLocalisations(settings.localisationsFilename);

		if (localisations == null)
			// Cancelled
			return;

		if (localisations.isEmpty())
		{
			IJ.error(TITLE, "No localisations could be loaded");
			return;
		}

		MemoryPeakResults results = localisations.toPeakResults();

		// Ask the user what depth to use to create the in-memory results
		if (!getZDepth(results))
			return;
		if (myLimitZ)
		{
			MemoryPeakResults results2 = new MemoryPeakResults(results.size());
			results2.setName(name);
			results2.copySettings(results);

			for (PeakResult peak : results.getResults())
			{
				if (peak.error < minz || peak.error > maxz)
					continue;
				results2.add(peak);
			}
			results = results2;
		}

		// Create the in-memory results
		if (results.size() > 0)
		{
			MemoryPeakResults.addResults(results);
		}

		IJ.showStatus(String.format("Loaded %d localisations", results.size()));
		if (myLimitZ)
			Utils.log("Loaded %d localisations, z between %.2f - %.2f", results.size(), minz, maxz);
		else
			Utils.log("Loaded %d localisations", results.size());
	}

	private boolean getZDepth(MemoryPeakResults results)
	{
		// The z-depth is stored in pixels in the error field
		double min = results.getHead().error;
		double max = min;
		for (PeakResult peak : results.getResults())
		{
			if (min > peak.error)
				min = peak.error;
			else if (max < peak.error)
				max = peak.error;
		}

		// No z-depth
		if (min == max && min == 0)
			return true;

		maxz = FastMath.min(maxz, max);
		minz = FastMath.max(minz, min);

		// Display in nm
		final double pp = results.getNmPerPixel();
		min *= pp;
		max *= pp;

		String msg = String.format("%d localisations with %.2f <= z <= %.2f", results.size(), min, max);

		min = Math.floor(min);
		max = Math.ceil(max);

		GenericDialog gd = new GenericDialog(TITLE);
		gd.addMessage(msg);
		gd.addCheckbox("Limit Z-depth", limitZ);
		gd.addSlider("minZ", min, max, minz * pp);
		gd.addSlider("maxZ", min, max, maxz * pp);
		gd.showDialog();
		if (gd.wasCanceled() || gd.invalidNumber())
		{
			return false;
		}
		myLimitZ = limitZ = gd.getNextBoolean();
		minz = gd.getNextNumber() / pp;
		maxz = gd.getNextNumber() / pp;
		return true;
	}

	static LocalisationList loadLocalisations(String filename)
	{
		if (!getFields())
			return null;

		LocalisationList localisations = new LocalisationList(timeUnit, distanceUnit, intensityUnit, gain, pixelPitch,
				exposureTime);

		final boolean hasComment = !Utils.isNullOrEmpty(comment);
		int errors = 0;
		int count = 0;
		int h = Math.max(0, header);

		BufferedReader input = null;
		try
		{
			FileInputStream fis = new FileInputStream(filename);
			input = new BufferedReader(new UnicodeReader(fis, null));
			Pattern p = Pattern.compile(delimiter);

			String line;
			while ((line = input.readLine()) != null)
			{
				// Skip header
				if (h-- > 0)
					continue;
				// Skip empty lines
				if (line.length() == 0)
					continue;
				// Skip comments
				if (hasComment && line.startsWith(comment))
					continue;

				count++;
				final String[] fields = p.split(line);

				Localisation l = new Localisation();
				try
				{
					if (it >= 0)
						l.t = Integer.parseInt(fields[it]);
					if (iid >= 0)
						l.id = Integer.parseInt(fields[iid]);
					l.x = Float.parseFloat(fields[ix]);
					l.y = Float.parseFloat(fields[iy]);
					if (iz >= 0)
						l.z = Float.parseFloat(fields[iz]);
					if (ii >= 0)
						l.intensity = Float.parseFloat(fields[ii]);
					if (isx >= 0)
						l.sy = l.sx = Integer.parseInt(fields[isx]);
					if (isy >= 0)
						l.sy = Integer.parseInt(fields[isy]);
					if (ip >= 0)
						l.precision = Float.parseFloat(fields[ip]);

					localisations.add(l);
				}
				catch (NumberFormatException e)
				{
					if (errors++ == 0)
						Utils.log("%s error on record %d: %s", TITLE, count, e.getMessage());
				}
				catch (IndexOutOfBoundsException e)
				{
					if (errors++ == 0)
						Utils.log("%s error on record %d: %s", TITLE, count, e.getMessage());
				}
			}
		}
		catch (IOException e)
		{
			Utils.log("%s IO error: %s", TITLE, e.getMessage());
		}
		finally
		{
			try
			{
				if (input != null)
					input.close();
			}
			catch (IOException e)
			{
				// Ignore
			}
			if (errors != 0)
				Utils.log("%s has %d / %d error lines", TITLE, errors, count);
		}

		return localisations;
	}

	private static boolean getFields()
	{
		GenericDialog gd = new GenericDialog(TITLE);

		gd.addMessage("Load delimited localisations");
		gd.addStringField("Dataset_name", name, 30);

		gd.addMessage("Calibration:");
		gd.addNumericField("Pixel_size", pixelPitch, 3, 8, "nm");
		gd.addNumericField("Gain", gain, 3, 8, "Count/photon");
		gd.addNumericField("Exposure_time", exposureTime, 3, 8, "");

		String[] tUnits = SettingsManager.getNames(set.toArray());
		gd.addChoice("Time_unit", tUnits, tUnits[timeUnit]);

		gd.addMessage("Records:");
		gd.addNumericField("Header_lines", header, 0);
		gd.addStringField("Comment", comment);
		gd.addStringField("Delimiter", delimiter);
		String[] dUnits = SettingsManager.distanceUnitNames;
		gd.addChoice("Distance_unit", dUnits, dUnits[distanceUnit]);
		String[] iUnits = SettingsManager.intensityUnitNames;
		gd.addChoice("Intensity_unit", iUnits, iUnits[intensityUnit]);

		gd.addMessage("Define the fields:");
		Label l = (Label) gd.getMessage();
		gd.addNumericField("T", it, 0);
		gd.addNumericField("ID", iid, 0);
		gd.addNumericField("X", ix, 0);
		gd.addNumericField("Y", iy, 0);
		gd.addNumericField("Z", iz, 0);
		gd.addNumericField("Intensity", ii, 0);
		gd.addNumericField("Sx", isx, 0);
		gd.addNumericField("Sy", isy, 0);
		gd.addNumericField("Precision", ip, 0);

		// Rearrange
		if (gd.getLayout() != null)
		{
			GridBagLayout grid = (GridBagLayout) gd.getLayout();

			int xOffset = 0, yOffset = 0;
			int lastY = -1, rowCount = 0;
			for (Component comp : gd.getComponents())
			{
				// Check if this should be the second major column
				if (comp == l)
				{
					xOffset += 2;
					yOffset = yOffset - rowCount + 1; // Skip title row
				}
				// Reposition the field
				GridBagConstraints c = grid.getConstraints(comp);
				if (lastY != c.gridy)
					rowCount++;
				lastY = c.gridy;
				c.gridx = c.gridx + xOffset;
				c.gridy = c.gridy + yOffset;
				c.insets.left = c.insets.left + 10 * xOffset;
				c.insets.top = 0;
				c.insets.bottom = 0;
				grid.setConstraints(comp, c);
			}

			if (IJ.isLinux())
				gd.setBackground(new Color(238, 238, 238));
		}

		gd.showDialog();
		if (gd.wasCanceled())
		{
			return false;
		}

		name = getNextString(gd, name);
		pixelPitch = gd.getNextNumber();
		gain = gd.getNextNumber();
		exposureTime = gd.getNextNumber();
		timeUnit = gd.getNextChoiceIndex();

		header = (int) gd.getNextNumber();
		comment = gd.getNextString();
		delimiter = getNextString(gd, delimiter);

		distanceUnit = gd.getNextChoiceIndex();
		intensityUnit = gd.getNextChoiceIndex();

		int[] columns = new int[9];
		for (int i = 0; i < columns.length; i++)
			columns[i] = (int) gd.getNextNumber();

		{
			int i = 0;
			it = columns[i++];
			iid = columns[i++];
			ix = columns[i++];
			iy = columns[i++];
			iz = columns[i++];
			ii = columns[i++];
			isx = columns[i++];
			isy = columns[i++];
			ip = columns[i++];
		}

		// Validate after reading the dialog (so the static fields store the last entered values)

		if (gd.invalidNumber())
		{
			IJ.error(TITLE, "Invalid number in input fields");
			return false;
		}
		for (int i = 0; i < columns.length; i++)
		{
			if (columns[i] < 0)
				continue;
			for (int j = i + 1; j < columns.length; j++)
			{
				if (columns[j] < 0)
					continue;
				if (columns[i] == columns[j])
				{
					IJ.error(TITLE, "Duplicate indicies: " + columns[i]);
					return false;
				}
			}
		}
		if (gain <= 0 || pixelPitch <= 0)
		{
			IJ.error(TITLE, "Require positive gain and pixel pitch");
			return false;
		}
		if (ix < 0 || iy < 0)
		{
			IJ.error(TITLE, "Require valid X and Y indices");
			return false;
		}
		return true;
	}

	private static String getNextString(GenericDialog gd, String defaultValue)
	{
		String value = gd.getNextString();
		if (Utils.isNullOrEmpty(value))
			return defaultValue;
		return value;
	}
}
