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

import gdsc.core.data.utils.ConversionException;
import gdsc.core.data.utils.IdentityTypeConverter;
import gdsc.core.data.utils.TypeConverter;

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
import gdsc.core.utils.NotImplementedException;
import gdsc.core.utils.TextUtils;
import gdsc.core.utils.UnicodeReader;
import gdsc.smlm.data.config.CalibrationHelper;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import gdsc.smlm.data.config.PSFProtos.PSFType;
import gdsc.smlm.data.config.UnitProtos.TimeUnit;
import gdsc.smlm.data.config.UnitConverterFactory;
import gdsc.smlm.data.config.UnitHelper;
import gdsc.smlm.data.config.GUIProtos.LoadLocalisationsSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.AttributePeakResult;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.procedures.BIXYZResultProcedure;
import gdsc.smlm.results.procedures.PeakResultProcedure;
import ij.IJ;
import ij.gui.ExtendedGenericDialog;
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
			CalibrationWriter calibration = new CalibrationWriter();
			calibration.setNmPerPixel(pixelPitch);
			calibration.setCountPerPhoton(gain);
			calibration.setExposureTime(timeConverter.convert(exposureTime));
			calibration.setDistanceUnit(distanceUnit);
			calibration.setIntensityUnit(intensityUnit);
			results.setCalibration(calibration.getCalibration());

			if (size() > 0)
			{
				// Guess the PSF type from the first localisation
				Localisation l = get(0);
				PSFType psfType = PSFType.CUSTOM;
				if (l.sx != -1)
				{
					psfType = PSFType.ONE_AXIS_GAUSSIAN_2D;
					if (l.sy != -1)
					{
						psfType = PSFType.TWO_AXIS_GAUSSIAN_2D;
					}
				}
				results.setPSF(PSFHelper.create(psfType));

				float intensity = (l.intensity <= 0) ? 1 : (float) (l.intensity);
				float x = (float) (l.x);
				float y = (float) (l.y);
				float z = (float) (l.z);

				for (int i = 0; i < size(); i++)
				{
					l = get(i);
					float[] params;
					switch (psfType)
					{
						case CUSTOM:
							params = PeakResult.createParams(0, intensity, x, y, z);
							break;
						case ONE_AXIS_GAUSSIAN_2D:
							params = Gaussian2DPeakResultHelper.createOneAxisParams(0, intensity, x, y, z,
									(float) l.sx);
							break;
						case TWO_AXIS_GAUSSIAN_2D:
							params = Gaussian2DPeakResultHelper.createTwoAxisParams(0, intensity, x, y, z, (float) l.sx,
									(float) l.sy);
							break;
						default:
							throw new NotImplementedException("Unsupported PSF type: " + psfType);
					}
					AttributePeakResult peakResult = new AttributePeakResult(l.t, (int) x, (int) y, 0, 0, 0, params,
							null);
					peakResult.setId(l.id);
					// Convert to nm
					peakResult.setPrecision(distanceConverter.convert(l.precision));
					results.add(peakResult);
				}
			}

			// Convert to preferred units. This can be done even if the results are empty.
			results.convertToPreferredUnits();

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

		LoadLocalisationsSettings.Builder settings = SettingsManager.readLoadLocalisationsSettings(0).toBuilder();

		String[] path = Utils.decodePath(settings.getLocalisationsFilename());
		OpenDialog chooser = new OpenDialog("Localisations_File", path[0], path[1]);
		if (chooser.getFileName() == null)
			return;

		settings.setLocalisationsFilename(chooser.getDirectory() + chooser.getFileName());
		SettingsManager.writeSettings(settings.build());

		LocalisationList localisations = loadLocalisations(settings.getLocalisationsFilename());

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
			final MemoryPeakResults results2 = new MemoryPeakResults(results.size());
			results2.setName(name);
			results2.copySettings(results);

			results.forEach(new PeakResultProcedure()
			{
				public void execute(PeakResult peak)
				{
					if (peak.error >= minz && peak.error <= maxz)
						results2.add(peak);
				}
			});
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

	private class ZResultProcedure implements BIXYZResultProcedure
	{
		float min, max;

		public void executeBIXYZ(float background, float intensity, float x, float y, float z)
		{
			if (min > z)
				min = z;
			if (max < z)
				max = z;
		}
	}

	private boolean getZDepth(MemoryPeakResults results)
	{
		final ZResultProcedure p = new ZResultProcedure();
		results.forEachNative(p);

		double min = p.min;
		double max = p.max;

		// No z-depth
		if (min == max && min == 0)
			return true;

		maxz = FastMath.min(maxz, max);
		minz = FastMath.max(minz, min);

		// Display in nm
		TypeConverter<DistanceUnit> c = new IdentityTypeConverter<DistanceUnit>(null);
		String unit = "";

		DistanceUnit nativeUnit = results.getDistanceUnit();
		if (nativeUnit != null)
		{
			unit = UnitHelper.getShortName(nativeUnit);
			try
			{
				c = CalibrationHelper.getDistanceConverter(results.getCalibration(), DistanceUnit.NM);
				unit = UnitHelper.getShortName(DistanceUnit.NM);
			}
			catch (ConversionException e)
			{

			}
		}
		min = c.convert(min);
		max = c.convert(max);

		String msg = String.format("%d localisations with %.2f <= z <= %.2f (%s)", results.size(), min, max, unit);

		min = Math.floor(min);
		max = Math.ceil(max);

		GenericDialog gd = new GenericDialog(TITLE);
		gd.addMessage(msg);
		gd.addCheckbox("Limit Z-depth", limitZ);
		gd.addSlider("minZ", min, max, c.convert(minz));
		gd.addSlider("maxZ", min, max, c.convert(maxz));
		gd.showDialog();
		if (gd.wasCanceled() || gd.invalidNumber())
		{
			return false;
		}
		myLimitZ = limitZ = gd.getNextBoolean();
		minz = c.convertBack(gd.getNextNumber());
		maxz = c.convertBack(gd.getNextNumber());
		return true;
	}

	static LocalisationList loadLocalisations(String filename)
	{
		if (!getFields())
			return null;

		LocalisationList localisations = new LocalisationList(timeUnit, distanceUnit, intensityUnit, gain, pixelPitch,
				exposureTime);

		final boolean hasComment = !TextUtils.isNullOrEmpty(comment);
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
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

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
		String[] dUnits = SettingsManager.getDistanceUnitNames();
		gd.addChoice("Distance_unit", dUnits, distanceUnit);
		String[] iUnits = SettingsManager.getIntensityUnitNames();
		gd.addChoice("Intensity_unit", iUnits, intensityUnit);

		gd.addMessage("Define the fields:");
		gd.addNumericField("T", it, 0);
		gd.addNumericField("ID", iid, 0);
		gd.addNumericField("X", ix, 0);
		gd.addNumericField("Y", iy, 0);
		gd.addNumericField("Z", iz, 0);
		gd.addNumericField("Intensity", ii, 0);
		gd.addNumericField("Sx", isx, 0);
		gd.addNumericField("Sy", isy, 0);
		gd.addNumericField("Precision", ip, 0);

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
		if (TextUtils.isNullOrEmpty(value))
			return defaultValue;
		return value;
	}
}
