package gdsc.smlm.ij.settings;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;

import com.thoughtworks.xstream.XStream;
import com.thoughtworks.xstream.XStreamException;
import com.thoughtworks.xstream.io.xml.DomDriver;

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

import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.utils.CodeReporter;
import ij.Prefs;

/**
 * Manage the settings for the gdsc.fitting package
 */
public class SettingsManager
{
	private static final String DEFAULT_FILENAME = System.getProperty("user.home") +
			System.getProperty("file.separator") + "gdsc.smlm.settings.xml";

	private static XStream xs = null;

	/**
	 * @return The settings filename (from the ImageJ preferences or the default in the home directory)
	 */
	public static String getSettingsFilename()
	{
		String filename = Prefs.get(Constants.settingsFilename, DEFAULT_FILENAME);
		return filename;
	}

	/**
	 * Save settings filename.
	 *
	 * @param filename the filename
	 */
	public static void saveSettingsFilename(String filename)
	{
		Prefs.set(Constants.settingsFilename, filename);
	}

	/**
	 * Convert a list of objects into names (e.g. pass in (Object[])enum.getValues()). The first letter is capitalised.
	 * The rest of the name is converted to lowercase if it is all uppercase. Remaining mixed case names are left
	 * unchanged.
	 * <p>
	 * Used to convert the settings enumerations into names used with dialogs.
	 * 
	 * @param objects
	 * @return the names
	 */
	public static String[] getNames(Object... objects)
	{
		String[] names = new String[objects.length];
		for (int i = 0; i < names.length; i++)
		{
			String name = objects[i].toString();

			if (name.length() > 0)
			{
				// Capitalise first letter
				if (Character.isLowerCase(name.charAt(0)))
					name = Character.toUpperCase(name.charAt(0)) + name.substring(1);

				//// Check if all upper-case
				//boolean isUpper = true;
				//for (int j = 0; j < name.length(); j++)
				//	if (Character.isLetter(name.charAt(j)) && !Character.isUpperCase(name.charAt(j)))
				//		isUpper = false;
				//
				//if (isUpper) // Use sentence case
				//	name = name.charAt(0) + name.substring(1).toLowerCase();
			}
			names[i] = name;
		}
		return names;
	}

	/**
	 * Save the configuration to file
	 * 
	 * @param config
	 * @param filename
	 * @return True if saved
	 */
	public static boolean saveFitEngineConfiguration(FitEngineConfiguration config, String filename)
	{
		XStream xs = createXStream();
		FileOutputStream fs = null;
		try
		{
			fs = new FileOutputStream(filename);
			xs.toXML(config, fs);
			return true;
		}
		catch (FileNotFoundException ex)
		{
			//ex.printStackTrace();
		}
		catch (XStreamException ex)
		{
			ex.printStackTrace();
		}
		finally
		{
			if (fs != null)
			{
				try
				{
					fs.close();
				}
				catch (IOException e)
				{
					e.printStackTrace();
				}
			}
		}
		return false;
	}

	private static XStream createXStream()
	{
		if (xs == null)
		{
			xs = new XStream(new DomDriver());

			// Map the object names/fields for a nicer configuration file
			xs.alias("gdsc.fitting.settings", GlobalSettings.class);
			xs.alias("gdsc.fitting.configuration", FitEngineConfiguration.class);

			xs.aliasField("gdsc.fitting.configuration", GlobalSettings.class, "fitEngineConfiguration");

			xs.aliasField("peakParameters", FitEngineConfiguration.class, "fitConfiguration");
			xs.aliasField("smoothing", FitEngineConfiguration.class, "smooth");
			//xs.aliasField("width0", FitConfiguration.class, "initialPeakWidth0");
			//xs.aliasField("width1", FitConfiguration.class, "initialPeakWidth1");
			//xs.aliasField("angle", FitConfiguration.class, "initialAngle");
			xs.aliasField("enableValidation", FitConfiguration.class, "fitValidation");

			xs.omitField(FitConfiguration.class, "flags");
			xs.omitField(FitConfiguration.class, "signalThreshold");
			//xs.omitField(FitConfiguration.class, "noise");
			xs.omitField(FitConfiguration.class, "enableValidation");
			xs.omitField(FitConfiguration.class, "computeResiduals");
			
			// Smart filter settings
			xs.omitField(FitConfiguration.class, "directFilter");
			xs.omitField(FitConfiguration.class, "dynamicPeakResult");
			xs.omitField(FitConfiguration.class, "filterResult");
			xs.omitField(FitConfiguration.class, "widthEnabled");
			xs.omitField(FitConfiguration.class, "offset");
		}

		return xs;
	}

	/**
	 * Load the configuration within the specified file
	 * 
	 * @param filename
	 * @return The configuration (or null)
	 */
	public static FitEngineConfiguration unsafeLoadFitEngineConfiguration(String filename)
	{
		XStream xs = createXStream();
		FitEngineConfiguration config = null;

		FileInputStream fs = null;
		try
		{
			fs = new FileInputStream(filename);
			config = (FitEngineConfiguration) xs.fromXML(fs);
		}
		catch (ClassCastException ex)
		{
			//ex.printStackTrace();
		}
		catch (FileNotFoundException ex)
		{
			//ex.printStackTrace();
		}
		catch (XStreamException ex)
		{
			ex.printStackTrace();
		}
		finally
		{
			if (fs != null)
			{
				try
				{
					fs.close();
				}
				catch (IOException e)
				{
					e.printStackTrace();
				}
			}
		}

		return config;
	}

	/**
	 * Load the configuration within the specified file
	 * 
	 * @param filename
	 * @return The configuration (or a default instance)
	 */
	public static FitEngineConfiguration loadFitEngineConfiguration(String filename)
	{
		FitEngineConfiguration config = unsafeLoadFitEngineConfiguration(filename);
		if (config == null)
			config = new FitEngineConfiguration(new FitConfiguration());
		return config;
	}

	/**
	 * Save the settings to file
	 * 
	 * @param settings
	 * @param filename
	 * @return True if saved
	 */
	public static boolean saveSettings(GlobalSettings settings, String filename)
	{
		CodeReporter.debug(new Throwable());
		XStream xs = createXStream();
		CodeReporter.debug(new Throwable());
		FileOutputStream fs = null;
		try
		{
			CodeReporter.debug(new Throwable());
			fs = new FileOutputStream(filename);
			CodeReporter.debug(new Throwable());
			xs.toXML(settings, fs);
			CodeReporter.debug(new Throwable());
			return true;
		}
		catch (FileNotFoundException ex)
		{
			ex.printStackTrace();
		}
		catch (XStreamException ex)
		{
			ex.printStackTrace();
		}
		finally
		{
			if (fs != null)
			{
				try
				{
					fs.close();
				}
				catch (IOException e)
				{
					e.printStackTrace();
				}
			}
		}
		CodeReporter.debug(new Throwable());
		return false;
	}

	/**
	 * Save the settings to the default file
	 * 
	 * @param settings
	 * @return True if saved
	 */
	public static boolean saveSettings(GlobalSettings settings)
	{
		return saveSettings(settings, getSettingsFilename());
	}

	/**
	 * Load the settings within the specified file
	 * 
	 * @param filename
	 * @return The settings (or null)
	 */
	public static GlobalSettings unsafeLoadSettings(String filename)
	{
		XStream xs = createXStream();
		GlobalSettings config = null;

		FileInputStream fs = null;
		try
		{
			fs = new FileInputStream(filename);
			config = (GlobalSettings) xs.fromXML(fs);
		}
		catch (ClassCastException ex)
		{
			//ex.printStackTrace();
		}
		catch (FileNotFoundException ex)
		{
			//ex.printStackTrace();
		}
		catch (XStreamException ex)
		{
			ex.printStackTrace();
		}
		finally
		{
			if (fs != null)
			{
				try
				{
					fs.close();
				}
				catch (IOException e)
				{
					e.printStackTrace();
				}
			}
		}

		return config;
	}

	/**
	 * Load the settings within the specified file
	 * 
	 * @param filename
	 * @return The settings (or a default instance)
	 */
	public static GlobalSettings loadSettings(String filename)
	{
		GlobalSettings config = unsafeLoadSettings(filename);
		if (config == null)
			config = new GlobalSettings();
		// This should not be null
		config.getFitEngineConfiguration().initialiseState();
		config.getCreateDataSettings().initialiseState();
		return config;
	}

	/**
	 * Load the settings from the default file
	 * 
	 * @return The settings (or a default instance)
	 */
	public static GlobalSettings loadSettings()
	{
		return loadSettings(getSettingsFilename());
	}
}
