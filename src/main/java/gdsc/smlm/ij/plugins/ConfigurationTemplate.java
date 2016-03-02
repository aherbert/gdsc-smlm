package gdsc.smlm.ij.plugins;

import gdsc.smlm.engine.DataFilter;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.fitting.FitSolver;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.Utils;
import ij.IJ;
import ij.plugin.PlugIn;

import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * This plugin loads configuration templates for the localisation fitting settings
 */
public class ConfigurationTemplate implements PlugIn
{
	private class Template
	{
		final GlobalSettings settings;
		final boolean custom;

		public Template(GlobalSettings settings, boolean custom)
		{
			this.settings = settings;
			this.custom = custom;

		}
	}

	private static HashMap<String, Template> map;
	private static ArrayList<String> names;
	private static String configurationDirectory;

	static
	{
		map = new HashMap<String, ConfigurationTemplate.Template>();
		names = new ArrayList<String>();

		String currentUsersHomeDir = System.getProperty("user.home");
		configurationDirectory = currentUsersHomeDir + File.separator + "gdsc.smlm";

		// Q. What settings should be in the template?

		FitConfiguration fitConfig = new FitConfiguration();
		FitEngineConfiguration config = new FitEngineConfiguration(fitConfig);

		fitConfig.setPrecisionUsingBackground(true);
		config.setFailuresLimit(1);

		// LSE
		fitConfig.setFitSolver(FitSolver.LVM);
		config.setDataFilter(DataFilter.MEAN, 1.2, 0);
		fitConfig.setCoordinateShiftFactor(1.2);
		fitConfig.setSignalStrength(35);
		fitConfig.setMinPhotons(30);
		fitConfig.setWidthFactor(1.8);
		fitConfig.setPrecisionThreshold(45);
		addTemplate("PALM LSE", config, false);

		// Change settings for different fit engines
		fitConfig.setFitSolver(FitSolver.MLE);
		config.setDataFilter(DataFilter.GAUSSIAN, 1.2, 0);
		fitConfig.setCoordinateShiftFactor(1.2);
		fitConfig.setSignalStrength(32);
		fitConfig.setMinPhotons(30);
		fitConfig.setWidthFactor(1.8);
		fitConfig.setPrecisionThreshold(47);
		addTemplate("PALM MLE", config, false);

		fitConfig.setModelCamera(true);
		fitConfig.setCoordinateShiftFactor(1.5);
		fitConfig.setSignalStrength(30);
		fitConfig.setMinPhotons(30);
		fitConfig.setWidthFactor(1.8);
		fitConfig.setPrecisionThreshold(50);
		addTemplate("PALM MLE Camera", config, false);

		// TODO: Add settings for STORM ...
	}

	private static void addTemplate(String name, FitEngineConfiguration config, boolean custom)
	{
		GlobalSettings settings = new GlobalSettings();
		settings.setFitEngineConfiguration(config.clone());
		addTemplate(name, settings, custom);
	}

	private static void addTemplate(String name, GlobalSettings settings, boolean custom)
	{
		// Maintain the names in the order they are added
		if (!map.containsKey(name))
			names.add(name);
		map.put(name, new ConfigurationTemplate().new Template(settings, custom));
	}

	/**
	 * Get the template configuration
	 * 
	 * @param name
	 *            The name of the template
	 * @return The template
	 */
	public static GlobalSettings getTemplate(String name)
	{
		Template template = map.get(name);
		return (template == null) ? null : template.settings;
	}

	/**
	 * Check if this is a custom template, i.e. not a standard GDSC SMLM template
	 * 
	 * @param name
	 *            The name of the template
	 * @return True if a custom template
	 */
	public static boolean isCustomTemplate(String name)
	{
		Template template = map.get(name);
		return (template == null) ? null : template.custom;
	}

	/**
	 * Get the names of the available templates
	 * 
	 * @param includeNone
	 *            Set to true to include [None] in the list of names
	 * @return The template names
	 */
	public static String[] getTemplateNames(boolean includeNone)
	{
		int length = (includeNone) ? names.size() + 1 : names.size();
		String[] templateNames = new String[length];
		int i = 0;
		if (includeNone)
			templateNames[i++] = "[None]";
		for (String name : names)
			templateNames[i++] = name;
		return templateNames;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);
		
		// Allow the user to specify a configuration directory
		String newDirectory = Utils.getDirectory("Template_directory", configurationDirectory);

		if (newDirectory == null)
			return;

		configurationDirectory = newDirectory;

		// Search the configuration directory and add any custom templates that can be deserialised from XML files
		File[] fileList = (new File(configurationDirectory)).listFiles(new FilenameFilter()
		{
			public boolean accept(File arg0, String arg1)
			{
				return arg1.toLowerCase().endsWith("xml");
			}
		});
		if (fileList == null)
			return;

		int count = 0;
		for (File file : fileList)
		{
			if (file.isFile())
			{
				GlobalSettings settings = SettingsManager.unsafeLoadSettings(file.getPath());
				if (settings != null)
				{
					count++;
					String name = Utils.removeExtension(file.getName());
					addTemplate(name, settings, true);
				}
			}
		}
		IJ.showMessage("Loaded " + Utils.pleural(count, "result"));
	}
}
