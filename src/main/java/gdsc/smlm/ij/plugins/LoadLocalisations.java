package gdsc.smlm.ij.plugins;

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

import uk.ac.sussex.gdsc.core.ij.ImageJUtils; import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils; import uk.ac.sussex.gdsc.core.utils.TextUtils; import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.UnicodeReader;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.ij.settings.CreateDataSettings;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.ExtendedPeakResult;
import gdsc.smlm.results.MemoryPeakResults;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import org.apache.commons.math3.util.FastMath;

/**
 * Loads generic localisation files into memory
 */
public class LoadLocalisations implements PlugIn
{
	static public class Localisation
	{
		int t, id;
		float x, y, z, intensity, sx = 1, sy = 1;
	}

	private static final String TITLE = "Load Localisations";
	private static boolean limitZ = false;
	private static double minz = -5;
	private static double maxz = 5;

	private static int it = 0;
	private static int iid = 1;
	private static int ix = 2;
	private static int iy = 3;
	private static int iz = 4;
	private static int ii = 5;
	private static int isx = -1;
	private static int isy = -1;
	private static int header = 1;
	private static String comment = "#";
	private static String delimiter = "\\t";
	private static String name = "Localisations";

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

		String[] path = ImageJUtils.decodePath(settings.localisationsFilename);
		OpenDialog chooser = new OpenDialog("Localisations_File", path[0], path[1]);
		if (chooser.getFileName() == null)
			return;

		settings.localisationsFilename = chooser.getDirectory() + chooser.getFileName();
		SettingsManager.saveSettings(globalSettings);

		List<Localisation> localisations = loadLocalisations(settings.localisationsFilename);

		if (localisations.isEmpty())
		{
			IJ.error(TITLE, "No localisations could be loaded");
			return;
		}

		// Ask the user what depth to use to create the in-memory results
		if (!getZDepth(localisations))
			return;

		// Create the in-memory results
		MemoryPeakResults results = new MemoryPeakResults();
		results.setName(name);

		for (Localisation l : localisations)
		{
			if (limitZ)
			{
				if (l.z < minz || l.z > maxz)
					continue;
			}

			float[] params = new float[7];
			params[Gaussian2DFunction.SIGNAL] = l.intensity;
			params[Gaussian2DFunction.X_POSITION] = l.x;
			params[Gaussian2DFunction.Y_POSITION] = l.y;
			params[Gaussian2DFunction.X_SD] = l.sx;
			params[Gaussian2DFunction.Y_SD] = l.sy;
			results.add(new ExtendedPeakResult(l.t, (int) l.x, (int) l.y, 0, 0, 0, params, null, l.t, l.id));
		}

		if (results.size() > 0)
		{
			ResultsManager.checkCalibration(results);
			MemoryPeakResults.addResults(results);
		}

		IJ.showStatus(String.format("Loaded %d localisations", results.size()));
		if (limitZ)
			ImageJUtils.log("Loaded %d localisations, z between %.2f - %.2f", results.size(), minz, maxz);
		else
			ImageJUtils.log("Loaded %d localisations", results.size());
	}

	private boolean getZDepth(List<Localisation> localisations)
	{
		double min = localisations.get(0).z;
		double max = min;
		for (Localisation l : localisations)
		{
			if (min > l.z)
				min = l.z;
			if (max < l.z)
				max = l.z;
		}

		maxz = FastMath.min(maxz, max);
		minz = FastMath.max(minz, min);

		String msg = String.format("%d localisations with %.2f <= z <= %.2f", localisations.size(), min, max);

		min = Math.floor(min);
		max = Math.ceil(max);

		GenericDialog gd = new GenericDialog(TITLE);
		gd.addMessage(msg);
		gd.addCheckbox("Limit Z-depth", limitZ);
		gd.addSlider("minZ", min, max, minz);
		gd.addSlider("maxZ", min, max, maxz);
		gd.showDialog();
		if (gd.wasCanceled() || gd.invalidNumber())
		{
			return false;
		}
		limitZ = gd.getNextBoolean();
		minz = gd.getNextNumber();
		maxz = gd.getNextNumber();
		return true;
	}

	private static List<Localisation> loadLocalisations(String filename)
	{
		List<Localisation> localisations = new ArrayList<Localisation>();

		getFields();

		final boolean hasComment = !TextUtils.isNullOrEmpty(comment);
		int errors = 0;
		int count = 0;
		int h = Math.abs(header);

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

					localisations.add(l);
				}
				catch (NumberFormatException e)
				{
					if (errors++ == 0)
						ImageJUtils.log("%s error on record %d: %s", TITLE, count, e.getMessage());
				}
				catch (IndexOutOfBoundsException e)
				{
					if (errors++ == 0)
						ImageJUtils.log("%s error on record %d: %s", TITLE, count, e.getMessage());
				}
			}
		}
		catch (IOException e)
		{
			// ignore
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
				ImageJUtils.log("%s has %d / %d error lines", TITLE, errors, count);
		}

		return localisations;
	}

	private static boolean getFields()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addMessage("Define the fields");
		gd.addNumericField("T", it, 0);
		gd.addNumericField("ID", iid, 0);
		gd.addNumericField("X", ix, 0);
		gd.addNumericField("Y", iy, 0);
		gd.addNumericField("Z", iz, 0);
		gd.addNumericField("Intensity", ii, 0);
		gd.addNumericField("Sx", isx, 0);
		gd.addNumericField("Sy", isy, 0);
		gd.addNumericField("Header_lines", header, 0);
		gd.addStringField("Comment", comment);
		gd.addStringField("Delimiter", delimiter);
		gd.addStringField("Name", name);
		gd.showDialog();
		if (gd.wasCanceled())
		{
			return false;
		}
		int[] columns = new int[8];
		for (int i = 0; i < columns.length; i++)
			columns[i] = (int) gd.getNextNumber();
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
		int i = 0;
		it = columns[i++];
		iid = columns[i++];
		ix = columns[i++];
		iy = columns[i++];
		iz = columns[i++];
		ii = columns[i++];
		isx = columns[i++];
		isy = columns[i++];
		header = (int) gd.getNextNumber();
		comment = gd.getNextString();
		delimiter = getNextString(gd, delimiter);
		name = getNextString(gd, name);
		if (ix < 0 || iy < 0 || ix == iy)
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
