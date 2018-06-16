/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 * 
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Choice;
import java.awt.Point;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;

import gdsc.core.ij.Utils;
import gdsc.core.utils.TextUtils;
import gdsc.core.utils.TurboList;
import gdsc.smlm.data.NamedObject;
import gdsc.smlm.data.config.FitProtos.DataFilterMethod;
import gdsc.smlm.data.config.FitProtos.FitSolver;
import gdsc.smlm.data.config.FitProtos.PrecisionMethod;
import gdsc.smlm.data.config.GUIProtos;
import gdsc.smlm.data.config.GUIProtos.ConfigurationTemplateSettings;
import gdsc.smlm.data.config.GUIProtos.ConfigurationTemplateSettings.Builder;
import gdsc.smlm.data.config.GUIProtos.DefaultTemplate;
import gdsc.smlm.data.config.GUIProtos.DefaultTemplateSettings;
import gdsc.smlm.data.config.TemplateProtos.TemplateSettings;
import gdsc.smlm.engine.FitConfiguration;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.ij.settings.SettingsManager;
import gnu.trove.map.hash.TIntObjectHashMap;
import ij.IJ;
import ij.ImageListener;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.ImageWindow;
import ij.gui.NonBlockingGenericDialog;
import ij.io.Opener;
import ij.plugin.PlugIn;
import ij.text.TextWindow;
import ij.util.StringSorter;

/**
 * This plugin loads configuration templates for the localisation fitting settings.
 */
public class ConfigurationTemplate implements PlugIn, DialogListener, ImageListener
{
	/**
	 * Describes the details of a template that can be loaded from the JAR resources folder.
	 */
	static class TemplateResource
	{
		/** The path. */
		final String path;

		/** The tif path. */
		final String tifPath;

		/** The name. */
		final String name;

		/**
		 * Instantiates a new template resource.
		 *
		 * @param path
		 *            the path
		 * @param name
		 *            the name
		 * @param tifPath
		 *            the tif path
		 */
		TemplateResource(String path, String name, String tifPath)
		{
			this.path = path;
			this.name = name;
			this.tifPath = tifPath;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Object#toString()
		 */
		@Override
		public String toString()
		{
			String text = String.format("path=%s, name=%s", path, name);
			if (tifPath != null)
				text += ", tifPath=" + tifPath;
			return text;
		}
	}

	/**
	 * The template type.
	 */
	private enum TemplateType
	{
		/** A template that was create using inline code. */
		INLINE,
		/** A template loaded from a jar resource. */
		RESOURCE,
		/** A custom template, e.g. loaded from file or saved from another plugin. */
		CUSTOM
	}

	/**
	 * The Class Template.
	 */
	private static class Template
	{
		/** The settings. */
		TemplateSettings settings;

		/** The template type. */
		TemplateType templateType;

		/** The file. */
		final File file;

		/** The timestamp. */
		long timestamp;

		/** The tif path. */
		// An example image from the data used to build the template
		String tifPath;

		/**
		 * Instantiates a new template.
		 *
		 * @param settings
		 *            the settings
		 * @param templateType
		 *            the template type
		 * @param file
		 *            the file
		 * @param tifPath
		 *            the tif path
		 */
		public Template(TemplateSettings settings, TemplateType templateType, File file, String tifPath)
		{
			this.settings = settings;
			this.templateType = templateType;
			this.file = file;
			timestamp = (file != null) ? file.lastModified() : 0;

			// Resource templates may have a tif image as a resource
			if (!TextUtils.isNullOrEmpty(tifPath))
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

		/**
		 * Update.
		 */
		public void update()
		{
			// Check if we can update from the file
			if (file != null)
			{
				if (file.lastModified() != timestamp)
				{
					TemplateSettings.Builder builder = TemplateSettings.newBuilder();
					if (SettingsManager.fromJSON(file, builder, 0))
					{
						this.settings = builder.build();
						timestamp = file.lastModified();
					}
				}
			}
		}

		/**
		 * Save the settings to file.
		 *
		 * @param file
		 *            the file
		 * @return true, if successful, False if failed (or no file to save to)
		 */
		public boolean save(File file)
		{
			boolean result = false;
			if (file != null)
			{
				result = SettingsManager.toJSON(settings, file, SettingsManager.FLAG_JSON_WHITESPACE);
				timestamp = file.lastModified();
			}
			return result;
		}

		/**
		 * Checks for image.
		 *
		 * @return true, if successful
		 */
		public boolean hasImage()
		{
			return tifPath != null;
		}

		/**
		 * Load image.
		 *
		 * @return the image plus
		 */
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

	/** A set of inline templates. These can be loaded. */
	private static LinkedHashMap<String, Template> inlineTemplates;
	/** The current set of templates that will be listed as loaded. */
	private static LinkedHashMap<String, Template> map;
	private String TITLE;
	private ImagePlus imp;
	private int currentSlice = 0;
	private TextWindow resultsWindow, infoWindow;
	private int templateId;
	private String headings;
	private TIntObjectHashMap<String> text;
	private boolean templateImage = false;

	static
	{
		// Maintain the names in the order they are added
		map = new LinkedHashMap<String, ConfigurationTemplate.Template>();

		createInlineTemplates();

		restoreLoadedTemplates();
	}

	private static void createInlineTemplates()
	{
		// Q. What settings should be in the template?
		inlineTemplates = new LinkedHashMap<String, ConfigurationTemplate.Template>();

		FitEngineConfiguration config = new FitEngineConfiguration();
		FitConfiguration fitConfig = config.getFitConfiguration();

		fitConfig.setPrecisionMethod(PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND);
		config.setFailuresLimit(1);

		// LSE
		fitConfig.setFitSolver(FitSolver.LVM_LSE);
		config.setDataFilter(DataFilterMethod.MEAN, 1.2, false, 0);
		fitConfig.setCoordinateShiftFactor(1.2);
		fitConfig.setSignalStrength(5);
		fitConfig.setMinPhotons(30);
		fitConfig.setMinWidthFactor(1 / 1.8); // Original code used the reciprocal
		fitConfig.setWidthFactor(1.8);
		fitConfig.setPrecisionThreshold(45);
		addInlineTemplate("PALM LSE", config);

		// Add settings for STORM ...
		config.setResidualsThreshold(0.4);
		config.setFailuresLimit(3);
		addInlineTemplate("STORM LSE", config);
		config.setResidualsThreshold(1);
		config.setFailuresLimit(1);

		// Change settings for different fit engines
		fitConfig.setFitSolver(FitSolver.MLE);
		config.setDataFilter(DataFilterMethod.GAUSSIAN, 1.2, false, 0);
		fitConfig.setCoordinateShiftFactor(1.2);
		fitConfig.setSignalStrength(4.5);
		fitConfig.setMinPhotons(30);
		fitConfig.setMinWidthFactor(1 / 1.8); // Original code used the reciprocal
		fitConfig.setWidthFactor(1.8);
		fitConfig.setPrecisionThreshold(47);
		addInlineTemplate("PALM MLE", config);

		// Add settings for STORM ...
		config.setResidualsThreshold(0.4);
		config.setFailuresLimit(3);
		addInlineTemplate("STORM MLE", config);
		config.setResidualsThreshold(1);
		config.setFailuresLimit(1);

		fitConfig.setModelCamera(true);
		fitConfig.setCoordinateShiftFactor(1.5);
		fitConfig.setSignalStrength(4.5);
		fitConfig.setMinPhotons(30);
		fitConfig.setMinWidthFactor(1 / 1.8); // Original code used the reciprocal
		fitConfig.setWidthFactor(1.8);
		fitConfig.setPrecisionThreshold(50);
		addInlineTemplate("PALM MLE Camera", config);

		// Add settings for STORM ...
		config.setResidualsThreshold(0.4);
		config.setFailuresLimit(3);
		addInlineTemplate("STORM MLE Camera", config);
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
	private static void addInlineTemplate(String name, FitEngineConfiguration config)
	{
		TemplateSettings.Builder builder = TemplateSettings.newBuilder();
		builder.setFitEngineSettings(config.getFitEngineSettings());
		Template template = new Template(builder.build(), TemplateType.INLINE, null, null);
		inlineTemplates.put(name, template);
	}

	private static String[] listInlineTemplates()
	{
		// Turn the keys into an array 
		ArrayList<String> names = new ArrayList<String>(inlineTemplates.keySet());
		return names.toArray(new String[inlineTemplates.size()]);
	}

	/**
	 * Restore the templates that were loaded.
	 * <p>
	 * Given the list of standard templates is manipulated only by this plugin this should
	 * be the same set of templates as that used last time by the user.
	 */
	private static void restoreLoadedTemplates()
	{
		// Allow this to fail silently
		DefaultTemplateSettings settings = SettingsManager.readDefaultTemplateSettings(SettingsManager.FLAG_SILENT);
		if (settings.getDefaultTemplatesCount() == 0)
			return;

		HashMap<String, TemplateResource> templateResources = null;

		// Process in order so that the order is preserved, i.e. do not bulk load each type
		for (DefaultTemplate d : settings.getDefaultTemplatesList())
		{
			switch (d.getTemplateType())
			{
				case CUSTOM_TEMPLATE:
					loadCustomTemplate(d.getName(), d.getFilename(), d.getTifFilename());
					break;

				case INLINE_TEMPLATE:
					Template t = inlineTemplates.get(d.getName());
					if (t != null)
						map.put(d.getName(), t);
					break;

				case RESOURCE_TEMPLATE:
					if (templateResources == null)
					{
						TemplateResource[] list = listTemplateResources();
						templateResources = new HashMap<String, TemplateResource>(list.length);
						for (TemplateResource r : list)
							templateResources.put(r.name, r);
					}
					loadTemplateResource(templateResources.get(d.getName()));
					break;

				default:
					break;
			}
		}

		if (map.size() != settings.getDefaultTemplatesCount())
		{
			// This occurs if we cannot reload some of the templates. 
			// Prevent this from happening again.
			saveLoadedTemplates();
		}
	}

	private static void loadCustomTemplate(String name, String path, String tifPath)
	{
		TemplateSettings.Builder builder = TemplateSettings.newBuilder();
		File file = new File(path);
		if (SettingsManager.fromJSON(file, builder, 0))
		{
			addTemplate(name, builder.build(), TemplateType.CUSTOM, file, tifPath);
		}
	}

	private static void loadTemplateResource(TemplateResource template)
	{
		if (template == null)
			return;
		Class<ConfigurationTemplate> resourceClass = ConfigurationTemplate.class;
		InputStream templateStream = resourceClass.getResourceAsStream(template.path);
		if (templateStream == null)
			return;

		InputStreamReader reader = new InputStreamReader(templateStream);
		TemplateSettings.Builder builder = TemplateSettings.newBuilder();
		if (SettingsManager.fromJSON(reader, builder, 0
		//SettingsManager.FLAG_SILENT
		))
		{
			addTemplate(template.name, builder.build(), TemplateType.RESOURCE, null, template.tifPath);
		}
	}

	/**
	 * Save the templates currently available in memory as the default templates to load on start-up.
	 */
	private static void saveLoadedTemplates()
	{
		DefaultTemplateSettings.Builder settings = DefaultTemplateSettings.newBuilder();
		DefaultTemplate.Builder defaultTemplate = DefaultTemplate.newBuilder();
		for (Entry<String, Template> entry : map.entrySet())
		{
			defaultTemplate.clear();
			defaultTemplate.setName(entry.getKey());
			Template t = entry.getValue();
			switch (t.templateType)
			{
				case CUSTOM:
					defaultTemplate.setTemplateType(GUIProtos.TemplateType.CUSTOM_TEMPLATE);
					break;
				case INLINE:
					defaultTemplate.setTemplateType(GUIProtos.TemplateType.INLINE_TEMPLATE);
					break;
				case RESOURCE:
					defaultTemplate.setTemplateType(GUIProtos.TemplateType.RESOURCE_TEMPLATE);
					break;
			}
			if (t.file != null)
				defaultTemplate.setFilename(t.file.getPath());
			if (t.tifPath != null)
				defaultTemplate.setTifFilename(t.tifPath);
			settings.addDefaultTemplates(defaultTemplate.build());
		}
		SettingsManager.writeSettings(settings);
	}

	/**
	 * List the templates from package resources.
	 *
	 * @return the templates
	 */
	static TemplateResource[] listTemplateResources()
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
		try
		{
			while ((line = input.readLine()) != null)
			{
				// Skip comment character
				if (line.length() == 0 || line.charAt(0) == '#')
					continue;

				//System.out.println(line);

				// Check the resource exists
				String path = templateDir + line;
				InputStream templateStream = resourceClass.getResourceAsStream(path);
				if (templateStream == null)
					continue;

				// Create a simple name
				String name = Utils.removeExtension(line);

				// Check if an example TIF file exists for the template
				String tifPath = templateDir + name + ".tif";
				InputStream tifStream = resourceClass.getResourceAsStream(tifPath);
				if (tifStream == null)
					tifPath = null;

				list.add(new TemplateResource(path, name, tifPath));
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
	 * @return the int
	 */
	static int loadTemplateResources(TemplateResource[] templates)
	{
		if (templates == null || templates.length == 0)
			return 0;
		int count = 0;
		Class<ConfigurationTemplate> resourceClass = ConfigurationTemplate.class;
		TemplateSettings.Builder builder = TemplateSettings.newBuilder();
		for (TemplateResource template : templates)
		{
			// Skip those already done
			if (map.containsKey(template.name))
				continue;

			InputStream templateStream = resourceClass.getResourceAsStream(template.path);
			if (templateStream == null)
				continue;

			//// Debug by printing the entire resource
			//{
			//	InputStreamReader reader = new InputStreamReader(resourceClass.getResourceAsStream(template.path));
			//	BufferedReader input = new BufferedReader(reader);
			//	String line;
			//	try
			//	{
			//		while ((line = input.readLine()) != null)
			//			System.out.println(line);
			//	}
			//	catch (Exception e)
			//	{
			//	}
			//}

			InputStreamReader reader = new InputStreamReader(templateStream);
			builder.clear();
			if (SettingsManager.fromJSON(reader, builder, 0
			//SettingsManager.FLAG_SILENT
			))
			{
				count++;
				addTemplate(template.name, builder.build(), TemplateType.RESOURCE, null, template.tifPath);
			}
		}
		return count;
	}

	/**
	 * Adds the template.
	 *
	 * @param name
	 *            the name
	 * @param settings
	 *            the settings
	 * @param templateType
	 *            the template type
	 * @param file
	 *            the file
	 * @param tifPath
	 *            the tif path
	 * @return the template
	 */
	private static Template addTemplate(String name, TemplateSettings settings, TemplateType templateType, File file,
			String tifPath)
	{
		Template template = new Template(settings, templateType, file, tifPath);
		map.put(name, template);
		return template;
	}

	/**
	 * Get the template configuration.
	 *
	 * @param name
	 *            The name of the template
	 * @return The template
	 */
	public static TemplateSettings getTemplate(String name)
	{
		Template template = map.get(name);
		if (template == null)
			return null;

		template.update();

		return template.settings;
	}

	/**
	 * Gets the template file.
	 *
	 * @param name
	 *            the name
	 * @return the template file
	 */
	public static File getTemplateFile(String name)
	{
		Template template = map.get(name);
		if (template == null)
			return null;
		return template.file;
	}

	/**
	 * Gets the template image.
	 *
	 * @param name
	 *            the name
	 * @return the template image
	 */
	public static ImagePlus getTemplateImage(String name)
	{
		Template template = map.get(name);
		if (template == null)
			return null;
		return template.loadImage();
	}

	/**
	 * Clear templates. Used for testing so made package level.
	 */
	static void clearTemplates()
	{
		map.clear();
	}

	/**
	 * Save template configuration. If an existing template exists with the same name it will be over-written. If an
	 * existing template was loaded from file it will be saved back to the same file, or optionally a new file.
	 *
	 * @param name
	 *            The name of the template
	 * @param settings
	 *            The template settings
	 * @param file
	 *            The file to save the template (over-riding the file the template was loaded from)
	 * @return true, if successful
	 */
	public static boolean saveTemplate(String name, TemplateSettings settings, File file)
	{
		Template template = map.get(name);
		if (template != null)
		{
			// Keep the file to allow it to be loaded on start-up
			if (file == null)
				file = template.file;
		}

		// Replace any existing template with a new one
		template = new Template(settings, TemplateType.CUSTOM, file, null);

		boolean result = true;
		if (file != null)
			result = template.save(file);
		if (result)
		{
			// Update the loaded templates
			map.put(name, template);
			saveLoadedTemplates();
		}
		return result;
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
		return (template == null) ? null : template.templateType == TemplateType.CUSTOM;
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
	 * Get the names of the available templates.
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

	//@formatter:off
	public enum TemplateOption implements NamedObject
	{
		LOAD_STANDARD_TEMPLATES("Load standard templates"),
		LOAD_CUSTOM_TEMPLATES("Load custom templates"),
		REMOVE_LOADED_TEMPLATES("Remove loaded templates"),
		VIEW_TEMPLATE("View template"),
		VIEW_TEMPLATE_IMAGE("View image example for template");

		private String name;
		TemplateOption(String name)
		{
			this.name = name;
		}
		@Override
		public String getName()
		{
			return name;
		}

		@Override
		public String getShortName()
		{
			return name;
		}
		
		static TemplateOption forNumber(int i)
		{
			TemplateOption[] values = TemplateOption.values();
			if (i < 0 || i >= values.length)
				i = 0;
			return values[i];
		}
	}
	//@formatter:on

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	@Override
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		ConfigurationTemplateSettings.Builder settings = SettingsManager.readConfigurationTemplateSettings(0)
				.toBuilder();

		TITLE = "Template Manager";
		GenericDialog gd = new GenericDialog(TITLE);
		String[] options = SettingsManager.getNames((Object[]) TemplateOption.values());
		gd.addChoice("Option", options, options[settings.getOption()]);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		int option = gd.getNextChoiceIndex();
		settings.setOption(option);
		switch (TemplateOption.forNumber(option))
		{
			case LOAD_CUSTOM_TEMPLATES:
				loadSelectedCustomTemplatesFromDirectory(settings);
				break;
			case LOAD_STANDARD_TEMPLATES:
				loadSelectedStandardTemplates(settings);
				break;
			case REMOVE_LOADED_TEMPLATES:
				removeLoadedTemplates(settings);
				break;
			case VIEW_TEMPLATE:
				showTemplate(settings);
				break;
			case VIEW_TEMPLATE_IMAGE:
				showTemplateImages(settings);
				break;
			default:
				break;

		}
		SettingsManager.writeSettings(settings);
	}

	/**
	 * Load selected standard templates.
	 *
	 * @param settings
	 *            the settings
	 */
	private void loadSelectedStandardTemplates(ConfigurationTemplateSettings.Builder settings)
	{
		final String[] inlineNames = listInlineTemplates();
		final TemplateResource[] templates = listTemplateResources();
		if (templates.length + inlineNames.length == 0)
			return;

		MultiDialog md = new MultiDialog("Select templates", new MultiDialog.BaseItems()
		{
			@Override
			public int size()
			{
				return templates.length + inlineNames.length;
			}

			@Override
			public String getFormattedName(int i)
			{
				if (i < inlineNames.length)
					return inlineNames[i];
				return templates[i - inlineNames.length].name;
			}
		});
		md.addSelected(settings.getSelectedStandardTemplatesList());

		md.showDialog();

		if (md.wasCanceled())
			return;

		ArrayList<String> selected = md.getSelectedResults();
		if (selected.isEmpty())
			return;

		// Save
		settings.clearSelectedStandardTemplates();
		settings.addAllSelectedStandardTemplates(selected);

		int count = 0;

		// Keep a hash of those not loaded from inline resources
		final HashSet<String> remaining = new HashSet<String>(selected.size());
		for (int i = 0; i < selected.size(); i++)
		{
			String name = selected.get(i);
			// Try and get the template from inline resources
			Template t = inlineTemplates.get(name);
			if (t != null)
			{
				count++;
				map.put(name, t);
			}
			else
			{
				remaining.add(name);
			}
		}

		if (!remaining.isEmpty())
		{
			// Build a list of resources to load
			TurboList<TemplateResource> list = new TurboList<TemplateResource>(remaining.size());
			for (TemplateResource t : templates)
			{
				if (remaining.contains(t.name))
					list.add(t);
			}
			count += loadTemplateResources(list.toArray(new TemplateResource[list.size()]));
		}

		if (count > 0)
			saveLoadedTemplates();
		IJ.showMessage("Loaded " + TextUtils.pleural(count, "standard template"));
	}

	/**
	 * Load templates from directory.
	 *
	 * @param settings
	 *            the settings
	 */
	private void loadSelectedCustomTemplatesFromDirectory(ConfigurationTemplateSettings.Builder settings)
	{
		// Allow the user to specify a configuration directory
		String newDirectory = Utils.getDirectory("Template_directory", settings.getConfigurationDirectory());

		if (newDirectory == null)
			return;

		settings.setConfigurationDirectory(newDirectory);

		// Search the configuration directory and add any custom templates that can be deserialised from XML files
		File[] fileList = (new File(newDirectory)).listFiles(new FileFilter()
		{
			@Override
			public boolean accept(File file)
			{
				// We can try and deserialise everything that is not a tif image
				// (which may be the template source image example)
				return file.isFile() && !file.getName().toLowerCase().endsWith("tif");
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
		final String[] sortedList = StringSorter.sortNumerically(list);

		// Select
		MultiDialog md = new MultiDialog("Select templates", new MultiDialog.BaseItems()
		{
			@Override
			public int size()
			{
				return sortedList.length;
			}

			@Override
			public String getFormattedName(int i)
			{
				String[] path = Utils.decodePath(sortedList[i]);
				return path[1];
			}
		});
		md.addSelected(settings.getSelectedCustomTemplatesList());

		md.showDialog();

		if (md.wasCanceled())
			return;

		ArrayList<String> selected = md.getSelectedResults();
		if (selected.isEmpty())
			return;

		// Save
		settings.clearSelectedCustomTemplates();
		settings.addAllSelectedCustomTemplates(selected);

		int count = 0;
		TemplateSettings.Builder builder = TemplateSettings.newBuilder();
		for (String path : selected)
		{
			builder.clear();
			File file = new File(newDirectory, path);
			if (SettingsManager.fromJSON(file, builder, 0))
			{
				count++;
				String name = Utils.removeExtension(file.getName());
				// Assume the Tif image will be detected automatically
				addTemplate(name, builder.build(), TemplateType.CUSTOM, file, null);
			}
		}

		if (count > 0)
			saveLoadedTemplates();
		IJ.showMessage("Loaded " + TextUtils.pleural(count, "custom template"));
	}

	private void removeLoadedTemplates(Builder settings)
	{
		if (map.isEmpty())
		{
			IJ.error(TITLE, "No templates are currently loaded");
			return;
		}
		final String[] names = getTemplateNames();
		MultiDialog md = new MultiDialog("Select templates to remove", new MultiDialog.BaseItems()
		{
			@Override
			public int size()
			{
				return names.length;
			}

			@Override
			public String getFormattedName(int i)
			{
				return names[i];
			}
		});

		md.showDialog();

		if (md.wasCanceled())
			return;

		ArrayList<String> selected = md.getSelectedResults();
		if (selected.isEmpty())
			// Nothing to do
			return;

		if (selected.size() == map.size())
		{
			clearTemplates();
		}
		else
		{
			for (String name : selected)
			{
				map.remove(name);
			}
		}
		saveLoadedTemplates();
	}

	/**
	 * Show template.
	 *
	 * @param settings
	 *            the settings
	 */
	private void showTemplate(ConfigurationTemplateSettings.Builder settings)
	{
		if (map.isEmpty())
		{
			IJ.error(TITLE, "No templates are currently loaded");
			return;
		}
		final String[] names = getTemplateNames();

		NonBlockingGenericDialog gd = new NonBlockingGenericDialog(TITLE);
		gd.addMessage("View the template");
		gd.addChoice("Template", names, settings.getTemplate());
		gd.addCheckbox("Close_on_exit", settings.getClose());
		gd.hideCancelButton();
		gd.addDialogListener(this);

		// Show the first template
		String template = ((Choice) (gd.getChoices().get(0))).getSelectedItem();
		showTemplate(template);

		gd.showDialog();

		// There is no cancel so read the settings.
		settings.setTemplate(gd.getNextChoice());
		settings.setClose(gd.getNextBoolean());

		if (settings.getClose())
		{
			closeInfo();
		}
	}

	/**
	 * Show template image.
	 *
	 * @param name
	 *            the name
	 */
	private void showTemplate(String name)
	{
		Template template = map.get(name);
		if (template == null)
		{
			IJ.error(TITLE, "Failed to load template: " + name);
			return;
		}

		template.update();

		if (infoWindow == null || !infoWindow.isVisible())
		{
			infoWindow = new TextWindow(TITLE + " Info", "", "", 450, 600);
		}

		infoWindow.getTextPanel().clear();
		add("Template", name);
		add("Type", (template.templateType == TemplateType.CUSTOM) ? "Custom" : "Standard");
		add("File", (template.file == null) ? null : template.file.getPath());
		add("Tif Image", (template.tifPath == null) ? null : template.tifPath);
		infoWindow.append("");
		infoWindow.append(template.settings.toString());
		infoWindow.getTextPanel().scrollToTop();
	}

	private void add(String key, String value)
	{
		if (value == null)
			return;
		infoWindow.append(key + " : " + value);
	}

	/**
	 * Show template images.
	 *
	 * @param settings
	 *            the settings
	 */
	private void showTemplateImages(ConfigurationTemplateSettings.Builder settings)
	{
		TITLE = "Template Example Images";

		String[] names = getTemplateNamesWithImage();
		if (names.length == 0)
		{
			IJ.error(TITLE, "No templates with example images");
			return;
		}

		// Follow when the image slice is changed
		ImagePlus.addImageListener(this);

		NonBlockingGenericDialog gd = new NonBlockingGenericDialog(TITLE);
		gd.addMessage("View the example source image");
		gd.addChoice("Template", names, settings.getTemplate());
		gd.addCheckbox("Close_on_exit", settings.getClose());
		gd.hideCancelButton();
		templateImage = true;
		gd.addDialogListener(this);

		// Show the first template
		String template = ((Choice) (gd.getChoices().get(0))).getSelectedItem();
		showTemplateImage(template);

		gd.showDialog();

		// There is no cancel so read the settings.
		settings.setTemplate(gd.getNextChoice());
		settings.setClose(gd.getNextBoolean());

		ImagePlus.removeImageListener(this);

		if (settings.getClose())
		{
			if (imp != null)
				imp.close();
			closeResults();
			closeInfo();
		}
	}

	/**
	 * Show template image.
	 *
	 * @param name
	 *            the name
	 */
	private void showTemplateImage(String name)
	{
		ImagePlus imp = getTemplateImage(name);
		if (imp == null)
		{
			IJ.error(TITLE, "Failed to load example image for template: " + name);
		}
		else
		{
			this.imp = displayTemplate(TITLE, imp);
			if (Utils.isNewWindow())
			{
				// Zoom a bit
				ImageWindow iw = this.imp.getWindow();
				for (int i = 7; i-- > 0 && Math.max(iw.getWidth(), iw.getHeight()) < 512;)
				{
					iw.getCanvas().zoomIn(0, 0);
				}
			}
			createResults(this.imp);

			showTemplateInfo(name);
		}
	}

	/**
	 * Display the template image in an image window with the specified title. If the window exists it will be reused
	 * and the appropriate properties updated.
	 *
	 * @param title
	 *            the title
	 * @param templateImp
	 *            the template image
	 * @return the image plus
	 */
	public static ImagePlus displayTemplate(String title, ImagePlus templateImp)
	{
		ImagePlus imp = Utils.display(title, templateImp.getStack());
		imp.setOverlay(templateImp.getOverlay());
		imp.setProperty("Info", templateImp.getProperty("Info"));
		imp.setCalibration(templateImp.getCalibration());
		return imp;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.gui.DialogListener#dialogItemChanged(ij.gui.GenericDialog, java.awt.AWTEvent)
	 */
	@Override
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
	{
		if (e != null && e.getSource() instanceof Choice)
		{
			String template = ((Choice) (e.getSource())).getSelectedItem();
			if (templateImage)
				showTemplateImage(template);
			else
				showTemplate(template);

		}
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.ImageListener#imageOpened(ij.ImagePlus)
	 */
	@Override
	public void imageOpened(ImagePlus imp)
	{

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.ImageListener#imageClosed(ij.ImagePlus)
	 */
	@Override
	public void imageClosed(ImagePlus imp)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.ImageListener#imageUpdated(ij.ImagePlus)
	 */
	@Override
	public void imageUpdated(ImagePlus imp)
	{
		if (imp != null && imp == this.imp)
		{
			updateResults(imp.getCurrentSlice());
		}
	}

	/**
	 * Creates a results window showing the localisation results from a template image. This will be positioned next to
	 * the input template image plus if it is currently displayed.
	 *
	 * @param templateImp
	 *            the template image
	 * @return the text window
	 */
	public TextWindow createResults(ImagePlus templateImp)
	{
		if (TITLE == null)
			TITLE = templateImp.getTitle();
		templateId = templateImp.getID();
		currentSlice = 0;
		headings = "";
		text = new TIntObjectHashMap<String>();
		Object info = templateImp.getProperty("Info");
		if (info != null)
		{
			// First line is the headings
			String[] lines = info.toString().split("\n");
			headings = lines[0].replace(' ', '\t');

			// The remaining lines are the data for each stack position
			StringBuilder sb = new StringBuilder();
			int last = 0;
			for (int i = 1; i < lines.length; i++)
			{
				// Get the position
				String[] data = lines[i].split(" ");
				int slice = Integer.parseInt(data[0]);
				if (last != slice)
				{
					text.put(last, sb.toString());
					last = slice;
					sb.setLength(0);
				}
				sb.append(slice);
				for (int j = 1; j < data.length; j++)
				{
					sb.append('\t').append(data[j]);
				}
				sb.append('\n');
			}
			text.put(last, sb.toString());
		}
		return updateResults(templateImp.getCurrentSlice());
	}

	/**
	 * Update the results window using the current selected slice from the template image.
	 *
	 * @param slice
	 *            the slice
	 * @return the text window
	 */
	public TextWindow updateResults(int slice)
	{
		if (slice == currentSlice || text == null)
			return resultsWindow;
		currentSlice = slice;

		if (resultsWindow == null || !resultsWindow.isVisible())
		{
			resultsWindow = new TextWindow(TITLE + " Results", headings, "", 450, 250);
			// Put next to the image
			ImagePlus imp = WindowManager.getImage(templateId);
			if (imp != null && imp.getWindow() != null)
			{
				ImageWindow iw = imp.getWindow();
				Point p = iw.getLocation();
				p.x += iw.getWidth();
				resultsWindow.setLocation(p);
			}
		}

		resultsWindow.getTextPanel().clear();
		String data = text.get(slice);
		if (!TextUtils.isNullOrEmpty(data))
			resultsWindow.append(data);
		return resultsWindow;
	}

	/**
	 * Close results.
	 */
	public void closeResults()
	{
		if (resultsWindow != null)
			resultsWindow.close();
	}

	/**
	 * Show the info from the template.
	 *
	 * @param name
	 *            the name
	 */
	private void showTemplateInfo(String name)
	{
		TemplateSettings settings = getTemplate(name);
		if (settings == null || settings.getNotesCount() == 0)
			return;
		if (infoWindow == null || !infoWindow.isVisible())
		{
			infoWindow = new TextWindow(TITLE + " Info", "", "", 450, 250);

			// Put underneath the results window
			if (resultsWindow != null)
			{
				Point p = resultsWindow.getLocation();
				p.y += resultsWindow.getHeight();
				infoWindow.setLocation(p);
			}
		}

		infoWindow.getTextPanel().clear();
		for (String note : settings.getNotesList())
			// Text window cannot show tabs
			infoWindow.append(note.replace('\t', ','));
	}

	/**
	 * Close info.
	 */
	private void closeInfo()
	{
		if (infoWindow != null)
			infoWindow.close();
	}
}
