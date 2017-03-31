package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Choice;

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
import java.util.LinkedHashMap;
import java.util.Map;
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
import ij.ImagePlus;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.ImageWindow;
import ij.gui.NonBlockingGenericDialog;
import ij.io.Opener;
import ij.plugin.PlugIn;
import ij.util.StringSorter;

/**
 * This plugin loads configuration templates for the localisation fitting settings
 */
public class ConfigurationTemplate implements PlugIn, DialogListener
{
	private String TITLE;
	private static String template = "";
	private static boolean close = true;
	private ImagePlus imp;

	/**
	 * Describes the details of a template that can be loaded from the JAR resources folder
	 */
	static class TemplateResource
	{
		final String path;
		final String tifPath;
		final String name;
		final boolean optional;

		TemplateResource(String path, String name, boolean optional, String tifPath)
		{
			this.path = path;
			this.name = name;
			this.optional = optional;
			this.tifPath = tifPath;
		}

		@Override
		public String toString()
		{
			String text = String.format("path=%s, name=%s, optional=%b", path, name, optional);
			if (tifPath != null)
				text += ", tifPath=" + tifPath;
			return text;
		}
	}

	private static class Template
	{
		GlobalSettings settings;
		final boolean custom;
		final File file;
		long timestamp;
		// An example image from the data used to build the template
		String tifPath;

		public Template(GlobalSettings settings, boolean custom, File file, String tifPath)
		{
			this.settings = settings;
			this.custom = custom;
			this.file = file;
			timestamp = (file != null) ? file.lastModified() : 0;

			// Resource templates may have a tif image as a resource
			if (tifPath != null)
			{
				this.tifPath = tifPath;
			}
			// Templates with a file may have a corresponding tif image
			else if (file != null)
			{
				tifPath = Utils.replaceExtension(file.getPath(), ".tif");
				if (new File(tifPath).exists())
					this.tifPath = tifPath;
			}
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

		public boolean hasImage()
		{
			return tifPath != null;
		}

		public ImagePlus loadImage()
		{
			if (!hasImage())
				return null;

			Opener opener = new Opener();
			opener.setSilentMode(true);

			// The tifPath may be a system resource or it may be a file 
			File file = new File(tifPath);
			if (file.exists())
			{
				// Load directly from a file path
				return opener.openImage(tifPath);
			}

			// IJ has support for loading TIFs from an InputStream
			Class<ConfigurationTemplate> resourceClass = ConfigurationTemplate.class;
			InputStream inputStream = resourceClass.getResourceAsStream(tifPath);
			if (inputStream != null)
			{
				return opener.openTiff(inputStream, Utils.removeExtension(file.getName()));
			}

			return null;
		}
	}

	private static LinkedHashMap<String, Template> map;
	private static boolean selectStandardTemplates = false;
	private static boolean selectCustomDirectory = true;
	private static String configurationDirectory;
	// Used for the multiMode option 
	private static ArrayList<String> selected;

	static
	{
		// Maintain the names in the order they are added
		map = new LinkedHashMap<String, ConfigurationTemplate.Template>();

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
		addTemplate("PALM LSE", config);

		// Add settings for STORM ...
		config.setResidualsThreshold(0.4);
		config.setFailuresLimit(3);
		addTemplate("STORM LSE", config);
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
		addTemplate("PALM MLE", config);

		// Add settings for STORM ...
		config.setResidualsThreshold(0.4);
		config.setFailuresLimit(3);
		addTemplate("STORM MLE", config);
		config.setResidualsThreshold(1);
		config.setFailuresLimit(1);

		fitConfig.setModelCamera(true);
		fitConfig.setCoordinateShiftFactor(1.5);
		fitConfig.setSignalStrength(30);
		fitConfig.setMinPhotons(30);
		fitConfig.setMinWidthFactor(1 / 1.8); // Original code used the reciprocal
		fitConfig.setWidthFactor(1.8);
		fitConfig.setPrecisionThreshold(50);
		addTemplate("PALM MLE Camera", config);

		// Add settings for STORM ...
		config.setResidualsThreshold(0.4);
		config.setFailuresLimit(3);
		addTemplate("STORM MLE Camera", config);

		loadStandardTemplates();
	}

	/**
	 * Adds the template using the configuration. This should be used to add templates that have not been produced using
	 * benchmarking on a specific image. Those can be added to the /gdsc/smlm/templates/ resources directory.
	 *
	 * @param name
	 *            the name
	 * @param config
	 *            the config
	 */
	private static void addTemplate(String name, FitEngineConfiguration config)
	{
		GlobalSettings settings = new GlobalSettings();
		settings.setFitEngineConfiguration(config.clone());
		addTemplate(name, settings, false, null, null);
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
				// Skip comment character
				if (line.length() == 0 || line.charAt(0) == '#')
					continue;

				String template = line;
				boolean optional = true;
				// Mandatory templates have a '*' suffix 
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
				// Check the resource exists
				String path = templateDir + template;
				InputStream templateStream = resourceClass.getResourceAsStream(path);
				if (templateStream == null)
					continue;

				// Create a simple name
				String name = Utils.removeExtension(template);

				// Check if an example TIF file exists for the template
				String tifPath = templateDir + name + ".tif";
				InputStream tifStream = resourceClass.getResourceAsStream(tifPath);
				if (tifStream == null)
					tifPath = null;

				list.add(new TemplateResource(path, name, optional, tifPath));
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
	static int loadTemplates(TemplateResource[] templates)
	{
		if (templates == null || templates.length == 0)
			return 0;
		int count = 0;
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
				count++;
				addTemplate(template.name, settings, false, null, template.tifPath);
			}
		}
		return count;
	}

	private static void addTemplate(String name, GlobalSettings settings, boolean custom, File file, String tifPath)
	{
		map.put(name, new Template(settings, custom, file, tifPath));
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

	public static ImagePlus getTemplateImage(String name)
	{
		Template template = map.get(name);
		if (template == null)
			return null;
		return template.loadImage();
	}

	static void clearTemplates()
	{
		map.clear();
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
			addTemplate(name, settings, true, file, null);
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
	 * Get the names of the available templates.
	 *
	 * @return The template names
	 */
	public static String[] getTemplateNames()
	{
		return getTemplateNames(false);
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
		int length = (includeNone) ? map.size() + 1 : map.size();
		String[] templateNames = new String[length];
		int i = 0;
		if (includeNone)
			templateNames[i++] = "[None]";
		for (String name : map.keySet())
			templateNames[i++] = name;
		return templateNames;
	}

	/**
	 * Get the names of the available templates that have an example image.
	 *
	 * @return The template names
	 */
	public static String[] getTemplateNamesWithImage()
	{
		TurboList<String> templateNames = new TurboList<String>(map.size());
		for (Map.Entry<String, Template> entry : map.entrySet())
			if (entry.getValue().hasImage())
				templateNames.add(entry.getKey());
		return templateNames.toArray(new String[templateNames.size()]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if ("images".equals(arg))
		{
			showTemplateImages();
			return;
		}

		TITLE = "Template Configuration";
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addCheckbox("Select_standard_templates", selectStandardTemplates);
		gd.addCheckbox("Select_custom_directory", selectCustomDirectory);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		selectStandardTemplates = gd.getNextBoolean();
		selectCustomDirectory = gd.getNextBoolean();

		if (selectStandardTemplates)
			loadSelectedStandardTemplates();
		if (selectCustomDirectory)
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
		int count = loadTemplates(list.toArray(new TemplateResource[list.size()]));
		IJ.showMessage("Loaded " + Utils.pleural(count, "standard template"));
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
				addTemplate(name, settings, true, file, null);
			}
		}
		IJ.showMessage("Loaded " + Utils.pleural(count, "custom template"));
	}

	private void showTemplateImages()
	{
		TITLE = "Template Example Images";

		String[] names = getTemplateNamesWithImage();
		if (names.length == 0)
		{
			IJ.error(TITLE, "No templates with example images");
			return;
		}

		NonBlockingGenericDialog gd = new NonBlockingGenericDialog(TITLE);
		gd.addMessage("View the example source image");
		gd.addChoice("Template", names, template);
		gd.addCheckbox("Close_image_on_exit", close);
		gd.hideCancelButton();
		gd.addDialogListener(this);

		// Show the first template
		showTemplateImage(((Choice) (gd.getChoices().get(0))).getSelectedItem());

		gd.showDialog();
		template = gd.getNextChoice();
		close = gd.getNextBoolean();

		if (close)
			closeImage();
	}

	private void showTemplateImage(String name)
	{
		ImagePlus imp = getTemplateImage(name);
		if (imp == null)
		{
			IJ.error(TITLE, "Failed to load example image for template: " + name);
		} else {
			this.imp = Utils.display(TITLE, imp.getStack());
			if (Utils.isNewWindow())
			{
				// Zoom a bit
				ImageWindow iw = this.imp.getWindow();
				for (int i = 7; i-- > 0 && Math.max(iw.getWidth(), iw.getHeight()) < 512;)
				{
					iw.getCanvas().zoomIn(0, 0);
				}
			}
		}
	}

	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
	{
		if (e != null && e.getSource() instanceof Choice)
		{
			showTemplateImage(((Choice) (e.getSource())).getSelectedItem());
		}
		return true;
	}

	private void closeImage()
	{
		if (imp != null)
			imp.close();
	}
}
