package gdsc.smlm.ij.plugins;

import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.HashMap;

import gdsc.core.ij.Utils;
import gdsc.smlm.engine.DataFilter;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.fitting.FitSolver;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import ij.IJ;
import ij.plugin.PlugIn;
import ij.util.StringSorter;

/**
 * This plugin loads configuration templates for the localisation fitting settings
 */
public class ConfigurationTemplate implements PlugIn
{
	private class Template
	{
		GlobalSettings settings;
		final boolean custom;
		final File file;
		long timestamp;

		public Template(GlobalSettings settings, boolean custom, File file)
		{
			this.settings = settings;
			this.custom = custom;
			this.file = file;
			timestamp = (file != null) ? file.lastModified() : 0;
		}

		public void update()
		{
			// Check if we can update from the file
			if (file != null)
			{
				if (file.lastModified() != timestamp)
				{
					GlobalSettings settings = SettingsManager.unsafeLoadSettings(file.getPath());
					if (settings != null)
					{
						this.settings = settings;
						timestamp = file.lastModified();
					}
				}
			}
		}

		/**
		 * Save the settings to file
		 * 
		 * @param file
		 *
		 * @return true, if successful, False if failed (or no file to save to)
		 */
		public boolean save(File file)
		{
			boolean result = false;
			if (file != null)
			{
				result = SettingsManager.saveSettings(settings, file.getPath());
			}
			return result;
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
		fitConfig.setMinWidthFactor(fitConfig.getMinWidthFactor()); // TODO
		fitConfig.setWidthFactor(1.8);
		fitConfig.setPrecisionThreshold(45);
		addTemplate("PALM LSE", config, false);

		// Add settings for STORM ...
		config.setResidualsThreshold(0.4);
		addTemplate("STORM LSE", config, false);
		config.setResidualsThreshold(1);

		// Change settings for different fit engines
		fitConfig.setFitSolver(FitSolver.MLE);
		config.setDataFilter(DataFilter.GAUSSIAN, 1.2, 0);
		fitConfig.setCoordinateShiftFactor(1.2);
		fitConfig.setSignalStrength(32);
		fitConfig.setMinPhotons(30);
		fitConfig.setMinWidthFactor(fitConfig.getMinWidthFactor()); // TODO
		fitConfig.setWidthFactor(1.8);
		fitConfig.setPrecisionThreshold(47);
		addTemplate("PALM MLE", config, false);

		// Add settings for STORM ...
		config.setResidualsThreshold(0.4);
		addTemplate("STORM MLE", config, false);
		config.setResidualsThreshold(1);

		fitConfig.setModelCamera(true);
		fitConfig.setCoordinateShiftFactor(1.5);
		fitConfig.setSignalStrength(30);
		fitConfig.setMinPhotons(30);
		fitConfig.setMinWidthFactor(fitConfig.getMinWidthFactor()); // TODO
		fitConfig.setWidthFactor(1.8);
		fitConfig.setPrecisionThreshold(50);
		addTemplate("PALM MLE Camera", config, false);

		// Add settings for STORM ...
		config.setResidualsThreshold(0.4);
		addTemplate("STORM MLE Camera", config, false);
		config.setResidualsThreshold(1);
	}

	private static void addTemplate(String name, FitEngineConfiguration config, boolean custom)
	{
		GlobalSettings settings = new GlobalSettings();
		settings.setFitEngineConfiguration(config.clone());
		addTemplate(name, settings, custom, null);
	}

	private static void addTemplate(String name, GlobalSettings settings, boolean custom, File file)
	{
		// Maintain the names in the order they are added
		if (!map.containsKey(name))
			names.add(name);
		map.put(name, new ConfigurationTemplate().new Template(settings, custom, file));
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
		if (template == null)
			return null;

		template.update();

		return template.settings;
	}

	/**
	 * Save template configuration. If an existing template exists with the same name it will be over-written. If an
	 * existing template was loaded from file it will be saved back to the same file, or optionally a different file.
	 *
	 * @param name
	 *            The name of the template
	 * @param settings
	 *            The template settings
	 * @param file
	 *            The file to save the template (over-riding the file the template was loaded from)
	 * @return true, if successful
	 */
	public static boolean saveTemplate(String name, GlobalSettings settings, File file)
	{
		Template template = map.get(name);
		if (template == null)
		{
			addTemplate(name, settings, true, file);
			return true;
		}
		template.settings = settings;
		if (file != null)
			return template.save(file);
		if (template.file != null)
			return template.save(template.file);
		return true;
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

		// Sort partially numerically
		String[] list = new String[fileList.length];
		int n = 0;
		for (File file : fileList)
		{
			if (file.isFile())
			{
				list[n++] = file.getPath();
			}
		}
		list = StringSorter.sortNumerically(list);

		int count = 0;
		for (String path : list)
		{
			GlobalSettings settings = SettingsManager.unsafeLoadSettings(path);
			if (settings != null)
			{
				count++;
				File file = new File(path);
				String name = Utils.removeExtension(file.getName());
				addTemplate(name, settings, true, file);
			}
		}
		IJ.showMessage("Loaded " + Utils.pleural(count, "result"));
	}
}
