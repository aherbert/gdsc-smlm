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

import gdsc.smlm.fitting.function.Gaussian2DFunction;
import gdsc.smlm.ij.settings.CreateDataSettings;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.model.LocalisationModel;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.utils.UnicodeReader;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Loads the localisation files created by Create Data into memory
 */
public class LoadLocalisations implements PlugIn
{
	private static final String TITLE = "Load Localisations";
	private static boolean limitZ = false;
	private static double minz = -5;
	private static double maxz = 5;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		GlobalSettings globalSettings = SettingsManager.loadSettings();
		CreateDataSettings settings = globalSettings.getCreateDataSettings();

		String[] path = Utils.decodePath(settings.localisationsFilename);
		OpenDialog chooser = new OpenDialog("Localisations_File", path[0], path[1]);
		if (chooser.getFileName() == null)
			return;

		settings.localisationsFilename = chooser.getDirectory() + chooser.getFileName();
		SettingsManager.saveSettings(globalSettings);

		List<LocalisationModel> localisations = loadLocalisations(settings.localisationsFilename);

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
		results.setName("Localisations");

		for (LocalisationModel l : localisations)
		{
			if (limitZ)
			{
				if (l.getZ() < minz || l.getZ() > maxz)
					continue;
			}

			float[] params = new float[7];
			params[Gaussian2DFunction.AMPLITUDE] = (float) (l.getIntensity() / (2 * Math.PI));
			params[Gaussian2DFunction.X_POSITION] = (float) l.getX();
			params[Gaussian2DFunction.Y_POSITION] = (float) l.getY();
			params[Gaussian2DFunction.X_WIDTH] = 1;
			params[Gaussian2DFunction.Y_WIDTH] = 1;
			results.add(new PeakResult(l.getTime(), (int) l.getX(), (int) l.getY(), 0, 0, 0, params, null));
		}

		if (results.size() > 0)
			MemoryPeakResults.addResults(results);

		IJ.showStatus(String.format("Loaded %d localisations", results.size()));
		if (limitZ)
			Utils.log("Loaded %d localisations, z between %.2f - %.2f", results.size(), minz, maxz);
		else
			Utils.log("Loaded %d localisations", results.size());
	}

	private boolean getZDepth(List<LocalisationModel> localisations)
	{
		double min = localisations.get(0).getZ();
		double max = min;
		for (LocalisationModel l : localisations)
		{
			if (min > l.getZ())
				min = l.getZ();
			if (max < l.getZ())
				max = l.getZ();
		}

		maxz = Math.min(maxz, max);
		minz = Math.max(minz, min);

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

	public static List<LocalisationModel> loadLocalisations(String filename)
	{
		List<LocalisationModel> localisations = new ArrayList<LocalisationModel>();

		BufferedReader input = null;
		try
		{
			FileInputStream fis = new FileInputStream(filename);
			input = new BufferedReader(new UnicodeReader(fis, null));

			String line;
			while ((line = input.readLine()) != null)
			{
				if (line.length() == 0)
					continue;
				if (line.charAt(0) == '#')
					continue;

				String[] fields = line.split("\t");
				if (fields.length >= 6)
				{
					int t = Integer.parseInt(fields[0]);
					int id = Integer.parseInt(fields[1]);
					float x = Float.parseFloat(fields[2]);
					float y = Float.parseFloat(fields[3]);
					float z = Float.parseFloat(fields[4]);
					float intensity = Float.parseFloat(fields[5]);

					localisations.add(new LocalisationModel(id, t, x, y, z, intensity, LocalisationModel.SINGLE));
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
		}

		return localisations;
	}
}
