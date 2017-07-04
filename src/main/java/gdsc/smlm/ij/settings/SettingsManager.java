package gdsc.smlm.ij.settings;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.EnumSet;

import com.google.protobuf.InvalidProtocolBufferException;
import com.google.protobuf.Message;
import com.google.protobuf.Parser;
import com.thoughtworks.xstream.XStream;
import com.thoughtworks.xstream.XStreamException;
import com.thoughtworks.xstream.io.xml.DomDriver;

import gdsc.core.clustering.ClusteringAlgorithm;
import gdsc.core.utils.BitFlags;
import gdsc.smlm.data.NamedObject;
import gdsc.smlm.data.config.CalibrationConfig.Calibration;
import gdsc.smlm.data.config.CalibrationConfigHelper;
import gdsc.smlm.data.config.FitConfig.DataFilterMethod;
import gdsc.smlm.data.config.FitConfig.DataFilterType;
import gdsc.smlm.data.config.FitConfig.FitEngineSettings;
import gdsc.smlm.data.config.FitConfig.FitSolver;
import gdsc.smlm.data.config.FitConfig.NoiseEstimatorMethod;
import gdsc.smlm.data.config.FitConfigHelper;
import gdsc.smlm.data.config.ResultsConfig.ResultsFileFormat;
import gdsc.smlm.data.config.ResultsConfig.ResultsImageType;
import gdsc.smlm.data.config.ResultsConfig.ResultsSettings;
import gdsc.smlm.data.config.ResultsConfigHelper;
import gdsc.smlm.data.config.UnitConfig.AngleUnit;
import gdsc.smlm.data.config.UnitConfig.DistanceUnit;
import gdsc.smlm.data.config.UnitConfig.IntensityUnit;
import gdsc.smlm.data.config.UnitConfig.TimeUnit;
import gdsc.smlm.data.config.UnitHelper;

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
import ij.IJ;
import ij.Prefs;

/**
 * Manage the settings for ImageJ plugins
 */
public class SettingsManager
{
	private static final String DEFAULT_FILENAME = System.getProperty("user.home") +
			System.getProperty("file.separator") + "gdsc.smlm.settings.xml";
	private static final String DEFAULT_DIRECTORY = System.getProperty("user.home") +
			System.getProperty("file.separator") + ".gdsc.smlm";

	/** Use this to suppress warnings. */
	public static final int FLAG_SILENT = 0x00000001;
	/** Use this flag to suppress returning a default instance. */
	public static final int FLAG_NO_DEFAULT = 0x00000002;

	/** The settings directory. */
	private static File settingsDirectory;
	static
	{
		setSettingsDirectory(Prefs.get(Constants.settingsDirectory, DEFAULT_DIRECTORY));
	}

	/**
	 * Gets the settings directory.
	 *
	 * @return The settings directory (from the ImageJ preferences or the default under the home directory)
	 */
	public static String getSettingsDirectory()
	{
		return settingsDirectory.getPath();
	}

	/**
	 * Save the settings directory.
	 *
	 * @param directory
	 *            the directory
	 */
	public static void setSettingsDirectory(String directory)
	{
		Prefs.set(Constants.settingsDirectory, directory);
		settingsDirectory = new File(directory);
		try
		{
			if (!settingsDirectory.exists())
				settingsDirectory.mkdirs();
		}
		catch (Exception e)
		{
			IJ.log("Unable create settings directory: " + e.getMessage());
		}
	}

	private static XStream xs = null;

	// Encapsulate the values and names of enums for lazy loading

	private static DistanceUnit[] _DistanceUnitValues;

	public static DistanceUnit[] getDistanceUnitValues()
	{
		if (_DistanceUnitValues == null)
			initDistanceUnit();
		return _DistanceUnitValues;
	}

	private static String[] _DistanceUnitNames;

	public static String[] getDistanceUnitNames()
	{
		if (_DistanceUnitNames == null)
			initDistanceUnit();
		return _DistanceUnitNames;
	}

	private static void initDistanceUnit()
	{
		EnumSet<DistanceUnit> d = EnumSet.allOf(DistanceUnit.class);
		d.remove(DistanceUnit.UNRECOGNIZED);
		//d.remove(DistanceUnit.DISTANCE_UNIT_NA);
		_DistanceUnitValues = d.toArray(new DistanceUnit[d.size()]);
		_DistanceUnitNames = new String[_DistanceUnitValues.length];
		for (int i = 0; i < _DistanceUnitValues.length; i++)
		{
			_DistanceUnitNames[i] = getName(UnitHelper.getName(_DistanceUnitValues[i]),
					UnitHelper.getShortName(_DistanceUnitValues[i]));
		}
	}

	private static IntensityUnit[] _IntensityUnitValues;

	public static IntensityUnit[] getIntensityUnitValues()
	{
		if (_IntensityUnitValues == null)
			initIntensityUnit();
		return _IntensityUnitValues;
	}

	private static String[] _IntensityUnitNames;

	public static String[] getIntensityUnitNames()
	{
		if (_IntensityUnitNames == null)
			initIntensityUnit();
		return _IntensityUnitNames;
	}

	private static void initIntensityUnit()
	{
		EnumSet<IntensityUnit> d = EnumSet.allOf(IntensityUnit.class);
		d.remove(IntensityUnit.UNRECOGNIZED);
		//d.remove(IntensityUnit.INTENSITY_UNIT_NA);
		_IntensityUnitValues = d.toArray(new IntensityUnit[d.size()]);
		_IntensityUnitNames = new String[_IntensityUnitValues.length];
		for (int i = 0; i < _IntensityUnitValues.length; i++)
		{
			_IntensityUnitNames[i] = getName(UnitHelper.getName(_IntensityUnitValues[i]),
					UnitHelper.getShortName(_IntensityUnitValues[i]));
		}
	}

	private static AngleUnit[] _AngleUnitValues;

	public static AngleUnit[] getAngleUnitValues()
	{
		if (_AngleUnitValues == null)
			initAngleUnit();
		return _AngleUnitValues;
	}

	private static String[] _AngleUnitNames;

	public static String[] getAngleUnitNames()
	{
		if (_AngleUnitNames == null)
			initAngleUnit();
		return _AngleUnitNames;
	}

	private static void initAngleUnit()
	{
		EnumSet<AngleUnit> d = EnumSet.allOf(AngleUnit.class);
		d.remove(AngleUnit.UNRECOGNIZED);
		//d.remove(AngleUnit.ANGLE_UNIT_NA);
		_AngleUnitValues = d.toArray(new AngleUnit[d.size()]);
		_AngleUnitNames = new String[_AngleUnitValues.length];
		for (int i = 0; i < _AngleUnitValues.length; i++)
		{
			_AngleUnitNames[i] = getName(UnitHelper.getName(_AngleUnitValues[i]),
					UnitHelper.getShortName(_AngleUnitValues[i]));
		}
	}

	private static TimeUnit[] _TimeUnitValues;

	public static TimeUnit[] getTimeUnitValues()
	{
		if (_TimeUnitValues == null)
			initTimeUnit();
		return _TimeUnitValues;
	}

	private static String[] _TimeUnitNames;

	public static String[] getTimeUnitNames()
	{
		if (_TimeUnitNames == null)
			initTimeUnit();
		return _TimeUnitNames;
	}

	private static void initTimeUnit()
	{
		EnumSet<TimeUnit> d = EnumSet.allOf(TimeUnit.class);
		d.remove(TimeUnit.UNRECOGNIZED);
		//d.remove(TimeUnit.TIME_UNIT_NA);
		_TimeUnitValues = d.toArray(new TimeUnit[d.size()]);
		_TimeUnitNames = new String[_TimeUnitValues.length];
		for (int i = 0; i < _TimeUnitValues.length; i++)
		{
			_TimeUnitNames[i] = getName(UnitHelper.getName(_TimeUnitValues[i]),
					UnitHelper.getShortName(_TimeUnitValues[i]));
		}
	}

	private static ResultsImageType[] _ResultsImageTypeValues;

	public static ResultsImageType[] getResultsImageTypeValues()
	{
		if (_ResultsImageTypeValues == null)
			initResultsImageType();
		return _ResultsImageTypeValues;
	}

	private static String[] _ResultsImageTypeNames;

	public static String[] getResultsImageTypeNames()
	{
		if (_ResultsImageTypeNames == null)
			initResultsImageType();
		return _ResultsImageTypeNames;
	}

	private static void initResultsImageType()
	{
		EnumSet<ResultsImageType> d = EnumSet.allOf(ResultsImageType.class);
		d.remove(ResultsImageType.UNRECOGNIZED);
		_ResultsImageTypeValues = d.toArray(new ResultsImageType[d.size()]);
		_ResultsImageTypeNames = new String[_ResultsImageTypeValues.length];
		for (int i = 0; i < _ResultsImageTypeValues.length; i++)
		{
			_ResultsImageTypeNames[i] = ResultsConfigHelper.getName(_ResultsImageTypeValues[i]);
		}
	}

	private static ResultsFileFormat[] _ResultsFileFormatValues;

	public static ResultsFileFormat[] getResultsFileFormatValues()
	{
		if (_ResultsFileFormatValues == null)
			initResultsFileFormat();
		return _ResultsFileFormatValues;
	}

	private static String[] _ResultsFileFormatNames;

	public static String[] getResultsFileFormatNames()
	{
		if (_ResultsFileFormatNames == null)
			initResultsFileFormat();
		return _ResultsFileFormatNames;
	}

	private static void initResultsFileFormat()
	{
		EnumSet<ResultsFileFormat> d = EnumSet.allOf(ResultsFileFormat.class);
		d.remove(ResultsFileFormat.UNRECOGNIZED);
		_ResultsFileFormatValues = d.toArray(new ResultsFileFormat[d.size()]);
		_ResultsFileFormatNames = new String[_ResultsFileFormatValues.length];
		for (int i = 0; i < _ResultsFileFormatValues.length; i++)
		{
			_ResultsFileFormatNames[i] = ResultsConfigHelper.getName(_ResultsFileFormatValues[i]);
		}
	}

	private static DataFilterType[] _DataFilterTypeValues;

	public static DataFilterType[] getDataFilterTypeValues()
	{
		if (_DataFilterTypeValues == null)
			initDataFilterType();
		return _DataFilterTypeValues;
	}

	private static String[] _DataFilterTypeNames;

	public static String[] getDataFilterTypeNames()
	{
		if (_DataFilterTypeNames == null)
			initDataFilterType();
		return _DataFilterTypeNames;
	}

	private static void initDataFilterType()
	{
		EnumSet<DataFilterType> d = EnumSet.allOf(DataFilterType.class);
		d.remove(DataFilterType.UNRECOGNIZED);
		_DataFilterTypeValues = d.toArray(new DataFilterType[d.size()]);
		_DataFilterTypeNames = new String[_DataFilterTypeValues.length];
		for (int i = 0; i < _DataFilterTypeValues.length; i++)
		{
			_DataFilterTypeNames[i] = FitConfigHelper.getName(_DataFilterTypeValues[i]);
		}
	}

	private static DataFilterMethod[] _DataFilterMethodValues;

	public static DataFilterMethod[] getDataFilterMethodValues()
	{
		if (_DataFilterMethodValues == null)
			initDataFilterMethod();
		return _DataFilterMethodValues;
	}

	private static String[] _DataFilterMethodNames;

	public static String[] getDataFilterMethodNames()
	{
		if (_DataFilterMethodNames == null)
			initDataFilterMethod();
		return _DataFilterMethodNames;
	}

	private static void initDataFilterMethod()
	{
		EnumSet<DataFilterMethod> d = EnumSet.allOf(DataFilterMethod.class);
		d.remove(DataFilterMethod.UNRECOGNIZED);
		_DataFilterMethodValues = d.toArray(new DataFilterMethod[d.size()]);
		_DataFilterMethodNames = new String[_DataFilterMethodValues.length];
		for (int i = 0; i < _DataFilterMethodValues.length; i++)
		{
			_DataFilterMethodNames[i] = FitConfigHelper.getName(_DataFilterMethodValues[i]);
		}
	}

	private static FitSolver[] _FitSolverValues;

	public static FitSolver[] getFitSolverValues()
	{
		if (_FitSolverValues == null)
			initFitSolver();
		return _FitSolverValues;
	}

	private static String[] _FitSolverNames;

	public static String[] getFitSolverNames()
	{
		if (_FitSolverNames == null)
			initFitSolver();
		return _FitSolverNames;
	}

	private static void initFitSolver()
	{
		EnumSet<FitSolver> d = EnumSet.allOf(FitSolver.class);
		d.remove(FitSolver.UNRECOGNIZED);
		_FitSolverValues = d.toArray(new FitSolver[d.size()]);
		_FitSolverNames = new String[_FitSolverValues.length];
		for (int i = 0; i < _FitSolverValues.length; i++)
		{
			_FitSolverNames[i] = FitConfigHelper.getName(_FitSolverValues[i]);
		}
	}

	private static NoiseEstimatorMethod[] _NoiseEstimatorMethodValues;

	public static NoiseEstimatorMethod[] getNoiseEstimatorMethodValues()
	{
		if (_NoiseEstimatorMethodValues == null)
			initNoiseEstimatorMethod();
		return _NoiseEstimatorMethodValues;
	}

	private static String[] _NoiseEstimatorMethodNames;

	public static String[] getNoiseEstimatorMethodNames()
	{
		if (_NoiseEstimatorMethodNames == null)
			initNoiseEstimatorMethod();
		return _NoiseEstimatorMethodNames;
	}

	private static void initNoiseEstimatorMethod()
	{
		EnumSet<NoiseEstimatorMethod> d = EnumSet.allOf(NoiseEstimatorMethod.class);
		d.remove(NoiseEstimatorMethod.UNRECOGNIZED);
		_NoiseEstimatorMethodValues = d.toArray(new NoiseEstimatorMethod[d.size()]);
		_NoiseEstimatorMethodNames = new String[_NoiseEstimatorMethodValues.length];
		for (int i = 0; i < _NoiseEstimatorMethodValues.length; i++)
		{
			_NoiseEstimatorMethodNames[i] = FitConfigHelper.getName(_NoiseEstimatorMethodValues[i]);
		}
	}

	public final static String[] clusteringAlgorithmNames;

	static
	{
		// TODO Update these as the configuration objects are auto-generated

		clusteringAlgorithmNames = SettingsManager.getNames((Object[]) ClusteringAlgorithm.values());
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
			String name;
			if (objects[i] instanceof NamedObject)
			{
				NamedObject o = (NamedObject) objects[i];
				name = getName(o.getName(), o.getShortName());
			}
			else
			{
				name = objects[i].toString();
			}

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

	private static String getName(String name, String shortName)
	{
		if (!name.equals(shortName))
			name += " (" + shortName + ")";
		return name;
	}

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
	 * @param filename
	 *            the filename
	 */
	public static void saveSettingsFilename(String filename)
	{
		Prefs.set(Constants.settingsFilename, filename);
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
			config = new FitEngineConfiguration();
		return config;
	}

	/**
	 * Save the settings to file.
	 * <p>
	 * If this fails then an error message is written to the ImageJ log
	 *
	 * @param settings
	 *            the settings
	 * @param filename
	 *            the filename
	 * @return True if saved
	 */
	public static boolean saveSettings(GlobalSettings settings, String filename)
	{
		return saveSettings(settings, filename, false);
	}

	/**
	 * Save the settings to file.
	 * <p>
	 * If this fails then an error message is written to the ImageJ log
	 *
	 * @param settings
	 *            the settings
	 * @param filename
	 *            the filename
	 * @param silent
	 *            Set to true to suppress writing an error message to the ImageJ log
	 * @return True if saved
	 */
	public static boolean saveSettings(GlobalSettings settings, String filename, boolean silent)
	{
		XStream xs = createXStream();
		FileOutputStream fs = null;
		try
		{
			fs = new FileOutputStream(filename);
			xs.toXML(settings, fs);
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
		if (!silent)
			IJ.log("Unable to save settings to: " + filename);
		return false;
	}

	/**
	 * Save the settings to the default file.
	 * <p>
	 * If this fails then an error message is written to the ImageJ log
	 *
	 * @param settings
	 *            the settings
	 * @return True if saved
	 */
	public static boolean saveSettings(GlobalSettings settings)
	{
		return saveSettings(settings, getSettingsFilename(), false);
	}

	/**
	 * Load the settings within the specified file
	 * <p>
	 * If this fails then an error message is written to the ImageJ log
	 * 
	 * @param filename
	 * @param silent
	 *            Set to true to suppress writing an error message to the ImageJ log
	 * @return The settings (or null)
	 */
	public static GlobalSettings unsafeLoadSettings(String filename, boolean silent)
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

		if (config == null)
			if (!silent)
				IJ.log("Unable to load settings from: " + filename);

		return config;
	}

	/**
	 * Load the settings from the input stream. The stream will be closed.
	 * <p>
	 * If this fails then an error message is written to the ImageJ log
	 *
	 * @param inputStream
	 *            the input stream
	 * @param silent
	 *            Set to true to suppress writing an error message to the ImageJ log
	 * @return The settings (or null)
	 */
	public static GlobalSettings unsafeLoadSettings(InputStream inputStream, boolean silent)
	{
		XStream xs = createXStream();
		GlobalSettings config = null;

		try
		{
			config = (GlobalSettings) xs.fromXML(inputStream);
		}
		catch (ClassCastException ex)
		{
			//ex.printStackTrace();
		}
		catch (XStreamException ex)
		{
			ex.printStackTrace();
		}
		finally
		{
			try
			{
				inputStream.close();
			}
			catch (IOException e)
			{
				e.printStackTrace();
			}
		}

		if (config == null)
			if (!silent)
				IJ.log("Unable to load settings from input stream");

		return config;
	}

	/**
	 * Load the settings within the specified file
	 * <p>
	 * If this fails then an error message is written to the ImageJ log
	 * 
	 * @param filename
	 * @return The settings (or a default instance)
	 */
	public static GlobalSettings loadSettings(String filename)
	{
		return loadSettings(filename, false);
	}

	/**
	 * Load the settings within the specified file
	 * <p>
	 * If this fails then an error message is written to the ImageJ log
	 * 
	 * @param filename
	 * @param silent
	 *            Set to true to suppress writing an error message to the ImageJ log
	 * @return The settings (or a default instance)
	 */
	public static GlobalSettings loadSettings(String filename, boolean silent)
	{
		GlobalSettings config = unsafeLoadSettings(filename, silent);
		if (config == null)
		{
			config = new GlobalSettings();
		}
		else
		{
			config.getFitEngineConfiguration().initialiseState();
			config.getCreateDataSettings().initialiseState();
		}
		return config;
	}

	/**
	 * Load the settings from the default file
	 * <p>
	 * If this fails then an error message is written to the ImageJ log
	 * 
	 * @return The settings (or a default instance)
	 */
	public static GlobalSettings loadSettings()
	{
		return loadSettings(getSettingsFilename(), false);
	}

	/**
	 * Creates the settings file using the class name to create the filename in the settings directory.
	 *
	 * @param clazz
	 *            the clazz
	 * @return the file
	 */
	private static File createSettingsFile(Class<?> clazz)
	{
		return new File(settingsDirectory, clazz.getSimpleName().toLowerCase() + ".settings");
	}

	/**
	 * Write a message to a settings file in the settings directory.
	 *
	 * @param message
	 *            the message
	 * @return true, if successful
	 */
	public static boolean writeSettings(Message message)
	{
		return writeSettings(message, 0);
	}

	/**
	 * Write a message to a settings file in the settings directory.
	 *
	 * @param message
	 *            the message
	 * @param flags
	 *            the flags
	 * @return true, if successful
	 */
	public static boolean writeSettings(Message message, int flags)
	{
		return writeMessage(message, createSettingsFile(message.getClass()), BitFlags.anySet(flags, FLAG_SILENT));
	}

	/**
	 * Clear the settings file for the given class.
	 *
	 * @param clazz
	 *            the class
	 * @return true, if the settings are cleared
	 */
	public static boolean clearSettings(Class<?> clazz)
	{
		File file = createSettingsFile(clazz);
		try
		{
			if (file.exists())
			{
				return file.delete();
			}
			return true; // Already clear
		}
		catch (SecurityException e)
		{
			IJ.log("Unable to clear the settings: " + e.getMessage());
		}
		return false;
	}

	/**
	 * Simple class to allow generic message reading.
	 *
	 * @param <T>
	 *            the generic message type
	 */
	public static class ConfigurationReader<T extends Message>
	{
		/** the default instance of the message type */
		private T t;

		/**
		 * Instantiates a new configuration reader.
		 *
		 * @param t
		 *            the default instance of the message type
		 */
		public ConfigurationReader(T t)
		{
			this.t = t;
		}

		/**
		 * Read the message.
		 *
		 * @return the message
		 */
		public T read()
		{
			return read(0);
		}

		/**
		 * Read the message.
		 *
		 * @param flags
		 *            the flags
		 * @return the message
		 */
		@SuppressWarnings("unchecked")
		public T read(int flags)
		{
			T c = (T) readMessage(t.getParserForType(), createSettingsFile(t.getClass()),
					BitFlags.anySet(flags, FLAG_SILENT));
			if (c == null && BitFlags.anyNotSet(flags, FLAG_NO_DEFAULT))
				c = t;
			return c;
		}
	}

	/**
	 * Read the Calibration from the settings file in the settings directory.
	 *
	 * @return the Calibration
	 */
	public static Calibration readCalibration()
	{
		return readCalibration(0);
	}

	/**
	 * Read the Calibration from the settings file in the settings directory.
	 *
	 * @param flags
	 *            the flags
	 * @return the Calibration
	 */
	public static Calibration readCalibration(int flags)
	{
		return new ConfigurationReader<Calibration>(CalibrationConfigHelper.defaultCalibration).read(flags);
	}

	/**
	 * Read the ResultsSettings from the settings file in the settings directory.
	 *
	 * @return the ResultsSettings
	 */
	public static ResultsSettings readResultsSettings()
	{
		return readResultsSettings(0);
	}

	/**
	 * Read the ResultsSettings from the settings file in the settings directory.
	 *
	 * @param flags
	 *            the flags
	 * @return the ResultsSettings
	 */
	public static ResultsSettings readResultsSettings(int flags)
	{
		return new ConfigurationReader<ResultsSettings>(ResultsConfigHelper.defaultResultsSettings).read(flags);
	}

	/**
	 * Read the FitEngineSettings from the settings file in the settings directory.
	 *
	 * @return the FitEngineSettings
	 */
	public static FitEngineSettings readFitEngineSettings()
	{
		return readFitEngineSettings(0);
	}

	/**
	 * Read the FitEngineSettings from the settings file in the settings directory.
	 *
	 * @param flags
	 *            the flags
	 * @return the FitEngineSettings
	 */
	public static FitEngineSettings readFitEngineSettings(int flags)
	{
		return new ConfigurationReader<FitEngineSettings>(FitConfigHelper.defaultFitEngineSettings).read(flags);
	}

	/**
	 * Write the message to file.
	 * <p>
	 * If this fails then an error message is written to the ImageJ log
	 *
	 * @param message
	 *            the message
	 * @param filename
	 *            the filename
	 * @param silent
	 *            Set to true to suppress writing an error message to the ImageJ log
	 * @return True if written
	 */
	public static boolean writeMessage(Message message, String filename, boolean silent)
	{
		return writeMessage(message, new File(filename), silent);
	}

	/**
	 * Write the message to file.
	 * <p>
	 * If this fails then an error message is written to the ImageJ log
	 *
	 * @param message
	 *            the message
	 * @param file
	 *            the file
	 * @param silent
	 *            Set to true to suppress writing an error message to the ImageJ log
	 * @return True if written
	 */
	public static boolean writeMessage(Message message, File file, boolean silent)
	{
		FileOutputStream fs = null;
		try
		{
			fs = new FileOutputStream(file);
			return writeMessage(message, fs, silent);
		}
		catch (FileNotFoundException e)
		{
			//e.printStackTrace();
			if (!silent)
				IJ.log("Unable to write message: " + e.getMessage());
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

	/**
	 * Write the message to file.
	 * <p>
	 * If this fails then an error message is written to the ImageJ log
	 *
	 * @param message
	 *            the message
	 * @param output
	 *            the output
	 * @param silent
	 *            Set to true to suppress writing an error message to the ImageJ log
	 * @return True if saved
	 */
	public static boolean writeMessage(Message message, OutputStream output, boolean silent)
	{
		try
		{
			message.writeDelimitedTo(output);
			return true;
		}
		catch (IOException e)
		{
			if (!silent)
				IJ.log("Unable to write message: " + e.getMessage());
		}
		return false;
	}

	/**
	 * Read the message from file.
	 * <p>
	 * If this fails then an error message is written to the ImageJ log
	 *
	 * @param parser
	 *            the parser
	 * @param filename
	 *            the filename
	 * @param silent
	 *            Set to true to suppress writing an error message to the ImageJ log
	 * @return the message
	 */
	public static Message readMessage(Parser<? extends Message> parser, String filename, boolean silent)
	{
		return readMessage(parser, new File(filename), silent);
	}

	/**
	 * Read the message from file.
	 * <p>
	 * If this fails then an error message is written to the ImageJ log
	 *
	 * @param parser
	 *            the parser
	 * @param file
	 *            the file
	 * @param silent
	 *            Set to true to suppress writing an error message to the ImageJ log
	 * @return the message
	 */
	public static Message readMessage(Parser<? extends Message> parser, File file, boolean silent)
	{
		FileInputStream fs = null;
		try
		{
			fs = new FileInputStream(file);
			return readMessage(parser, fs, silent);
		}
		catch (FileNotFoundException e)
		{
			//e.printStackTrace();
			if (!silent)
				IJ.log("Unable to read message to: " + e.getMessage());
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
		return null;
	}

	/**
	 * Read the message from file.
	 * <p>
	 * If this fails then an error message is written to the ImageJ log
	 *
	 * @param parser
	 *            the parser
	 * @param input
	 *            the input
	 * @param silent
	 *            Set to true to suppress writing an error message to the ImageJ log
	 * @return the message
	 */
	public static Message readMessage(Parser<? extends Message> parser, InputStream input, boolean silent)
	{
		try
		{
			return parser.parseDelimitedFrom(input);
		}
		catch (InvalidProtocolBufferException e)
		{
			//e.printStackTrace();
			if (!silent)
				IJ.log("Unable to read message to: " + e.getMessage());
		}
		return null;
	}
}
