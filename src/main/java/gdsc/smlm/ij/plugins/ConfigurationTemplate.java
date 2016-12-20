package gdsc.smlm.ij.plugins;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Pattern;

import gdsc.core.ij.Utils;
import gdsc.core.utils.TurboList;
import gdsc.core.utils.TurboList.SimplePredicate;
import gdsc.smlm.engine.DataFilter;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.fitting.FitSolver;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.util.StringSorter;

/**
 * This plugin loads configuration templates for the localisation fitting settings
 */
public class ConfigurationTemplate implements PlugIn
{
	static class TemplateResource
	{
		final String path;
		final String name;
		final boolean optional;

		TemplateResource(String path, String name, boolean optional)
		{
			this.path = path;
			this.name = name;
			this.optional = optional;
		}
	}

	private static class Template
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
					GlobalSettings settings = SettingsManager.unsafeLoadSettings(file.getPath(), false);
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
	private static boolean selectStandardTemplates = false;
	private static boolean chooseDirectory = true;
	private static String configurationDirectory;
	// Used for the multiMode option 
	private static ArrayList<String> selected;

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
		fitConfig.setMinWidthFactor(1 / 1.8); // Original code used the reciprocal
		fitConfig.setWidthFactor(1.8);
		fitConfig.setPrecisionThreshold(45);
		addTemplate("PALM LSE", config, false);

		// Add settings for STORM ...
		config.setResidualsThreshold(0.4);
		config.setFailuresLimit(3);
		addTemplate("STORM LSE", config, false);
		config.setResidualsThreshold(1);
		config.setFailuresLimit(1);

		// Change settings for different fit engines
		fitConfig.setFitSolver(FitSolver.MLE);
		config.setDataFilter(DataFilter.GAUSSIAN, 1.2, 0);
		fitConfig.setCoordinateShiftFactor(1.2);
		fitConfig.setSignalStrength(32);
		fitConfig.setMinPhotons(30);
		fitConfig.setMinWidthFactor(1 / 1.8); // Original code used the reciprocal
		fitConfig.setWidthFactor(1.8);
		fitConfig.setPrecisionThreshold(47);
		addTemplate("PALM MLE", config, false);

		// Add settings for STORM ...
		config.setResidualsThreshold(0.4);
		config.setFailuresLimit(3);
		addTemplate("STORM MLE", config, false);
		config.setResidualsThreshold(1);
		config.setFailuresLimit(1);

		fitConfig.setModelCamera(true);
		fitConfig.setCoordinateShiftFactor(1.5);
		fitConfig.setSignalStrength(30);
		fitConfig.setMinPhotons(30);
		fitConfig.setMinWidthFactor(1 / 1.8); // Original code used the reciprocal
		fitConfig.setWidthFactor(1.8);
		fitConfig.setPrecisionThreshold(50);
		addTemplate("PALM MLE Camera", config, false);

		// Add settings for STORM ...
		config.setResidualsThreshold(0.4);
		config.setFailuresLimit(3);
		addTemplate("STORM MLE Camera", config, false);

		loadStandardTemplates();
	}

	private static void addTemplate(String name, FitEngineConfiguration config, boolean custom)
	{
		GlobalSettings settings = new GlobalSettings();
		settings.setFitEngineConfiguration(config.clone());
		addTemplate(name, settings, custom, null);
	}

	/**
	 * Load standard templates (those not marked as optional).
	 */
	private static void loadStandardTemplates()
	{
		TemplateResource[] templates = listTemplates(true, false);
		loadTemplates(templates);
	}

	/**
	 * List the templates from package resources.
	 *
	 * @param loadMandatory
	 *            Set to true to list the mandatory templates
	 * @param loadOptional
	 *            Set to true to list the optional templates
	 * @return the templates
	 */
	static TemplateResource[] listTemplates(boolean loadMandatory, boolean loadOptional)
	{
		// Load templates from package resources
		String templateDir = "/gdsc/smlm/templates/";
		Class<ConfigurationTemplate> resourceClass = ConfigurationTemplate.class;
		InputStream templateListStream = resourceClass.getResourceAsStream(templateDir + "list.txt");
		if (templateListStream == null)
			return new TemplateResource[0];
		BufferedReader input = new BufferedReader(new InputStreamReader(templateListStream));
		String line;
		ArrayList<TemplateResource> list = new ArrayList<TemplateResource>();
		Pattern p = Pattern.compile("\\*");
		try
		{
			while ((line = input.readLine()) != null)
			{
				String template = line;
				boolean optional = true;
				int index = template.indexOf('*');
				if (index >= 0)
				{
					template = p.matcher(template).replaceAll("");
					optional = false;
				}
				if (optional)
				{
					if (!loadOptional)
						continue;
				}
				else
				{
					if (!loadMandatory)
						continue;
				}
				String name = Utils.removeExtension(template);
				list.add(new TemplateResource(templateDir + template, name, optional));
			}
		}
		catch (IOException e)
		{

		}
		return list.toArray(new TemplateResource[list.size()]);
	}

	/**
	 * Load templates from package resources.
	 *
	 * @param templates
	 *            the templates
	 */
	static void loadTemplates(TemplateResource[] templates)
	{
		if (templates == null || templates.length == 0)
			return;
		Class<ConfigurationTemplate> resourceClass = ConfigurationTemplate.class;
		for (TemplateResource template : templates)
		{
			// Skip those already done
			if (map.containsKey(template.name))
				continue;
			
			InputStream templateStream = resourceClass.getResourceAsStream(template.path);
			if (templateStream == null)
				continue;
			GlobalSettings settings = SettingsManager.unsafeLoadSettings(templateStream, true);
			if (settings != null)
			{
				addTemplate(template.name, settings, false, null);
			}
		}
	}

	private static void addTemplate(String name, GlobalSettings settings, boolean custom, File file)
	{
		// Maintain the names in the order they are added
		if (!map.containsKey(name))
			names.add(name);
		map.put(name, new Template(settings, custom, file));
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

		GenericDialog gd = new GenericDialog("Template Configuration");
		gd.addCheckbox("Select_templates", selectStandardTemplates);
		gd.addCheckbox("Choose_directory", chooseDirectory);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		selectStandardTemplates = gd.getNextBoolean();
		chooseDirectory = gd.getNextBoolean();

		if (selectStandardTemplates)
			loadSelectedStandardTemplates();
		if (chooseDirectory)
			loadTemplatesFromDirectory();
	}

	private void loadSelectedStandardTemplates()
	{
		final TemplateResource[] templates = listTemplates(false, true);
		if (templates.length == 0)
			return;

		MultiDialog md = new MultiDialog("Select Templates", new MultiDialog.BaseItems()
		{
			public int size()
			{
				return templates.length;
			}

			public String getFormattedName(int i)
			{
				return templates[i].name;
			}
		});
		md.addSelected(selected);

		md.showDialog();

		if (md.wasCanceled())
			return;

		selected = md.getSelectedResults();
		if (selected.isEmpty())
			return;
		
		// Use list filtering to get the selected templates
		TurboList<TemplateResource> list = new TurboList<TemplateResource>(Arrays.asList(templates));
		list.removeIf(new SimplePredicate<TemplateResource>()
		{
			public boolean test(TemplateResource t)
			{
				return !(selected.contains(t.name));
			}
		});
		loadTemplates(list.toArray(new TemplateResource[list.size()]));
	}

	private void loadTemplatesFromDirectory()
	{
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
			GlobalSettings settings = SettingsManager.unsafeLoadSettings(path, false);
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
