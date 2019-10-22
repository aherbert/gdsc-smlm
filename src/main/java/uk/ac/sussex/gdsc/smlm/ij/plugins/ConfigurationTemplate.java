/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

import uk.ac.sussex.gdsc.core.ij.ImageAdapter;
import uk.ac.sussex.gdsc.core.ij.ImageJPluginLoggerHelper;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.MultiDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.utils.AlphaNumericComparator;
import uk.ac.sussex.gdsc.core.utils.FileUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.TurboList;
import uk.ac.sussex.gdsc.core.utils.ValidationUtils;
import uk.ac.sussex.gdsc.smlm.data.NamedObject;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.DataFilterMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitSolver;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.PrecisionMethod;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.ConfigurationTemplateSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.DefaultTemplate;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.DefaultTemplateSettings;
import uk.ac.sussex.gdsc.smlm.data.config.TemplateProtos.TemplateSettings;
import uk.ac.sussex.gdsc.smlm.engine.FitConfiguration;
import uk.ac.sussex.gdsc.smlm.engine.FitEngineConfiguration;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;

import gnu.trove.map.hash.TIntObjectHashMap;

import ij.IJ;
import ij.ImageListener;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.ImageWindow;
import ij.gui.NonBlockingGenericDialog;
import ij.io.Opener;
import ij.plugin.PlugIn;
import ij.text.TextWindow;

import org.apache.commons.lang3.ArrayUtils;

import java.awt.AWTEvent;
import java.awt.Choice;
import java.awt.Point;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Map.Entry;
import java.util.stream.Collectors;

/**
 * This plugin loads configuration templates for the localisation fitting settings.
 */
public class ConfigurationTemplate implements PlugIn {
  /** The title for the template manager plugin. */
  private static final String TITLE = "Template Manager";

  /** A set of inline templates. These can be loaded. */
  private static Map<String, Template> inlineTemplates;
  /** The current set of templates that will be listed as loaded. */
  private static Map<String, Template> templates;

  private String title;
  private ImagePlus imp;
  private int currentSlice;
  private TextWindow resultsWindow;
  private TextWindow infoWindow;
  private int templateId;
  private String headings;
  private TIntObjectHashMap<String> text;
  private boolean templateImage;

  static {
    final Map<String, Template> localInlineTemplates = createInlineTemplates();

    // Make maps synchronized
    inlineTemplates = Collections.synchronizedMap(localInlineTemplates);
    templates = Collections.synchronizedMap(restoreLoadedTemplates(localInlineTemplates));
  }

  /**
   * Describes the details of a template that can be loaded from the JAR resources folder.
   */
  static class TemplateResource {
    /** The path. */
    final String path;

    /** The tif path. */
    final String tifPath;

    /** The name. */
    final String name;

    /**
     * Instantiates a new template resource.
     *
     * @param path the path
     * @param name the name
     * @param tifPath the tif path
     */
    TemplateResource(String path, String name, String tifPath) {
      this.path = path;
      this.name = name;
      this.tifPath = tifPath;
    }

    @Override
    public String toString() {
      String text = String.format("path=%s, name=%s", path, name);
      if (tifPath != null) {
        text += ", tifPath=" + tifPath;
      }
      return text;
    }
  }

  /**
   * The template type.
   */
  private enum TemplateType {
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
  private static class Template {
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
     * @param settings the settings
     * @param templateType the template type
     * @param file the file
     * @param tifPath the tif path
     */
    public Template(TemplateSettings settings, TemplateType templateType, File file,
        String tifPath) {
      this.settings = settings;
      this.templateType = templateType;
      this.file = file;
      timestamp = (file != null) ? file.lastModified() : 0;

      // Resource templates may have a tif image as a resource
      if (!TextUtils.isNullOrEmpty(tifPath)) {
        this.tifPath = tifPath;
      } else if (file != null) {
        tifPath = FileUtils.replaceExtension(file.getPath(), ".tif");
        if (new File(tifPath).exists()) {
          this.tifPath = tifPath;
        }
      }
    }

    /**
     * Update.
     */
    public void update() {
      // Check if we can update from the file
      if (file != null && file.lastModified() != timestamp) {
        final TemplateSettings.Builder builder = TemplateSettings.newBuilder();
        if (SettingsManager.fromJson(file, builder, 0)) {
          this.settings = builder.build();
          timestamp = file.lastModified();
        }
      }
    }

    /**
     * Save the settings to file.
     *
     * @param file the file
     * @return true, if successful, False if failed (or no file to save to)
     */
    public boolean save(File file) {
      boolean result = false;
      if (file != null) {
        result = SettingsManager.toJson(settings, file, SettingsManager.FLAG_JSON_WHITESPACE);
        timestamp = file.lastModified();
      }
      return result;
    }

    /**
     * Checks for image.
     *
     * @return true, if successful
     */
    public boolean hasImage() {
      return tifPath != null;
    }

    /**
     * Load image.
     *
     * @return the image plus
     */
    public ImagePlus loadImage() {
      if (!hasImage()) {
        return null;
      }

      final Opener opener = new Opener();
      opener.setSilentMode(true);

      // The tifPath may be a system resource or it may be a file
      final File tifFile = new File(tifPath);
      if (tifFile.exists()) {
        // Load directly from a file path
        return opener.openImage(tifPath);
      }

      // IJ has support for loading TIFs from an InputStream
      final Class<ConfigurationTemplate> resourceClass = ConfigurationTemplate.class;
      ImagePlus imp = null;
      try (InputStream inputStream = resourceClass.getResourceAsStream(tifPath)) {
        if (inputStream != null) {
          imp = opener.openTiff(inputStream, FileUtils.removeExtension(tifFile.getName()));
        }
      } catch (final IOException ex) {
        // Ignore
      }

      return imp;
    }
  }

  private static Map<String, Template> createInlineTemplates() {
    final Map<String, Template> inlineTemplates = new LinkedHashMap<>();

    // Q. What settings should be in the template?
    final FitEngineConfiguration config = new FitEngineConfiguration();
    final FitConfiguration fitConfig = config.getFitConfiguration();

    fitConfig.setPrecisionMethod(PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND);
    config.setFailuresLimit(1);

    // LSE
    fitConfig.setFitSolver(FitSolver.LVM_LSE);
    config.setDataFilter(DataFilterMethod.MEAN, 1.2, false, 0);
    fitConfig.setCoordinateShiftFactor(1.2);
    fitConfig.setSignalStrength(5);
    fitConfig.setMinPhotons(30);
    fitConfig.setMinWidthFactor(1 / 1.8); // Original code used the reciprocal
    fitConfig.setMaxWidthFactor(1.8);
    fitConfig.setPrecisionThreshold(45);
    addInlineTemplate(inlineTemplates, "PALM LSE", config);

    // Add settings for STORM ...
    config.setResidualsThreshold(0.4);
    config.setFailuresLimit(3);
    addInlineTemplate(inlineTemplates, "STORM LSE", config);
    config.setResidualsThreshold(1);
    config.setFailuresLimit(1);

    // Change settings for different fit engines
    fitConfig.setFitSolver(FitSolver.MLE);
    config.setDataFilter(DataFilterMethod.GAUSSIAN, 1.2, false, 0);
    fitConfig.setCoordinateShiftFactor(1.2);
    fitConfig.setSignalStrength(4.5);
    fitConfig.setMinPhotons(30);
    fitConfig.setMinWidthFactor(1 / 1.8); // Original code used the reciprocal
    fitConfig.setMaxWidthFactor(1.8);
    fitConfig.setPrecisionThreshold(47);
    addInlineTemplate(inlineTemplates, "PALM MLE", config);

    // Add settings for STORM ...
    config.setResidualsThreshold(0.4);
    config.setFailuresLimit(3);
    addInlineTemplate(inlineTemplates, "STORM MLE", config);
    config.setResidualsThreshold(1);
    config.setFailuresLimit(1);

    fitConfig.setModelCamera(true);
    fitConfig.setCoordinateShiftFactor(1.5);
    fitConfig.setSignalStrength(4.5);
    fitConfig.setMinPhotons(30);
    fitConfig.setMinWidthFactor(1 / 1.8); // Original code used the reciprocal
    fitConfig.setMaxWidthFactor(1.8);
    fitConfig.setPrecisionThreshold(50);
    addInlineTemplate(inlineTemplates, "PALM MLE Camera", config);

    // Add settings for STORM ...
    config.setResidualsThreshold(0.4);
    config.setFailuresLimit(3);
    addInlineTemplate(inlineTemplates, "STORM MLE Camera", config);

    return inlineTemplates;
  }

  /**
   * Adds the template using the configuration. This should be used to add templates that have not
   * been produced using benchmarking on a specific image. Those can be added to the
   * {@code /uk/ac/sussex/gdsc/smlm/templates/} resources directory.
   *
   * @param inlineTemplates the inline templates
   * @param name the name
   * @param config the config
   */
  private static void addInlineTemplate(Map<String, Template> inlineTemplates, String name,
      FitEngineConfiguration config) {
    final TemplateSettings.Builder builder = TemplateSettings.newBuilder();
    builder.setFitEngineSettings(config.getFitEngineSettings());
    final Template template = new Template(builder.build(), TemplateType.INLINE, null, null);
    inlineTemplates.put(name, template);
  }

  private static String[] listInlineTemplates() {
    // Turn the keys into an array
    return inlineTemplates.keySet().toArray(new String[0]);
  }

  /**
   * Restore the templates that were loaded.
   *
   * <p>Given the list of standard templates is manipulated only by this plugin this should be the
   * same set of templates used last time by the user.
   *
   * @param inlineTemplates the inline templates
   * @return the templates
   */
  private static Map<String, Template>
      restoreLoadedTemplates(Map<String, Template> inlineTemplates) {
    Map<String, Template> map = new LinkedHashMap<>();

    // Allow this to fail silently
    final DefaultTemplateSettings settings =
        SettingsManager.readDefaultTemplateSettings(SettingsManager.FLAG_SILENT);
    if (settings.getDefaultTemplatesCount() == 0) {
      return map;
    }

    HashMap<String, TemplateResource> templateResources = null;

    // Process in order so that the order is preserved, i.e. do not bulk load each type
    for (final DefaultTemplate d : settings.getDefaultTemplatesList()) {
      switch (d.getTemplateType()) {
        case CUSTOM_TEMPLATE:
          loadCustomTemplate(map, d.getName(), d.getFilename(), d.getTifFilename());
          break;

        case INLINE_TEMPLATE:
          final Template t = inlineTemplates.get(d.getName());
          if (t != null) {
            map.put(d.getName(), t);
          }
          break;

        case RESOURCE_TEMPLATE:
          if (templateResources == null) {
            final TemplateResource[] list = listTemplateResources();
            templateResources = new HashMap<>(list.length);
            for (final TemplateResource r : list) {
              templateResources.put(r.name, r);
            }
          }
          loadTemplateResource(map, templateResources.get(d.getName()));
          break;

        default:
          break;
      }
    }

    if (map.size() != settings.getDefaultTemplatesCount()) {
      // This occurs if we cannot reload some of the templates.
      // Prevent this from happening again.
      saveLoadedTemplates(map);
    }

    return map;
  }

  private static void loadCustomTemplate(Map<String, Template> map, String name, String path,
      String tifPath) {
    if (TextUtils.isNullOrEmpty(path)) {
      return;
    }
    final TemplateSettings.Builder builder = TemplateSettings.newBuilder();
    final File file = new File(path);
    if (SettingsManager.fromJson(file, builder, 0)) {
      addTemplate(map, name, builder.build(), TemplateType.CUSTOM, file, tifPath);
    }
  }

  private static void loadTemplateResource(Map<String, Template> map, TemplateResource template) {
    if (template == null) {
      return;
    }
    final Class<ConfigurationTemplate> resourceClass = ConfigurationTemplate.class;
    final InputStream templateStream = resourceClass.getResourceAsStream(template.path);
    if (templateStream == null) {
      return;
    }

    try (InputStreamReader reader = new InputStreamReader(templateStream, StandardCharsets.UTF_8)) {
      final TemplateSettings.Builder builder = TemplateSettings.newBuilder();
      if (SettingsManager.fromJson(reader, builder, 0
      // SettingsManager.FLAG_SILENT
      )) {
        addTemplate(map, template.name, builder.build(), TemplateType.RESOURCE, null,
            template.tifPath);
      }
    } catch (final IOException ex) {
      // Ignore
    }
  }

  /**
   * Save the provided templates as the default templates to load on start-up.
   *
   * @param map the template map
   */
  private static void saveLoadedTemplates(Map<String, Template> map) {
    final DefaultTemplateSettings.Builder settings = DefaultTemplateSettings.newBuilder();
    final DefaultTemplate.Builder defaultTemplate = DefaultTemplate.newBuilder();
    for (final Entry<String, Template> entry : map.entrySet()) {
      defaultTemplate.clear();
      defaultTemplate.setName(entry.getKey());
      final Template t = entry.getValue();
      switch (t.templateType) {
        case CUSTOM:
          defaultTemplate.setTemplateType(GUIProtos.TemplateType.CUSTOM_TEMPLATE);
          break;
        case INLINE:
          defaultTemplate.setTemplateType(GUIProtos.TemplateType.INLINE_TEMPLATE);
          break;
        default:
          ValidationUtils.checkArgument(t.templateType == TemplateType.RESOURCE, t.templateType);
          defaultTemplate.setTemplateType(GUIProtos.TemplateType.RESOURCE_TEMPLATE);
          break;
      }
      if (t.file != null) {
        defaultTemplate.setFilename(t.file.getPath());
      }
      if (t.tifPath != null) {
        defaultTemplate.setTifFilename(t.tifPath);
      }
      settings.addDefaultTemplates(defaultTemplate.build());
    }
    SettingsManager.writeSettings(settings);
  }

  /**
   * List the templates from package resources.
   *
   * @return the templates
   */
  static TemplateResource[] listTemplateResources() {
    // Load templates from package resources
    final String templateDir = "/uk/ac/sussex/gdsc/smlm/templates/";
    final Class<ConfigurationTemplate> resourceClass = ConfigurationTemplate.class;
    try (InputStream templateListStream =
        resourceClass.getResourceAsStream(templateDir + "list.txt")) {
      if (templateListStream == null) {
        return new TemplateResource[0];
      }
      try (final BufferedReader input =
          new BufferedReader(new InputStreamReader(templateListStream, StandardCharsets.UTF_8))) {
        String line;
        final ArrayList<TemplateResource> list = new ArrayList<>();
        while ((line = input.readLine()) != null) {
          // Skip comment character
          if (line.length() == 0 || line.charAt(0) == '#') {
            continue;
          }

          // Check the resource exists
          final String path = templateDir + line;
          try (InputStream templateStream = resourceClass.getResourceAsStream(path)) {
            if (templateStream == null) {
              continue;
            }
          }

          // Create a simple name
          final String name = FileUtils.removeExtension(line);

          // Check if an example TIF file exists for the template
          String tifPath = templateDir + name + ".tif";
          try (InputStream tifStream = resourceClass.getResourceAsStream(tifPath)) {
            if (tifStream == null) {
              tifPath = null;
            }
          }

          list.add(new TemplateResource(path, name, tifPath));
        }
        return list.toArray(new TemplateResource[0]);
      }
    } catch (final IOException ex) {
      // Ignore
    }
    return new TemplateResource[0];
  }

  /**
   * Load templates from package resources.
   *
   * @param templateResources the template resources
   * @return the number loaded
   */
  static int loadTemplateResources(TemplateResource[] templateResources) {
    if (ArrayUtils.getLength(templateResources) == 0) {
      return 0;
    }
    int count = 0;
    final Class<ConfigurationTemplate> resourceClass = ConfigurationTemplate.class;
    final TemplateSettings.Builder builder = TemplateSettings.newBuilder();
    for (final TemplateResource template : templateResources) {
      // Skip those already done
      if (templates.containsKey(template.name)) {
        continue;
      }

      final InputStream templateStream = resourceClass.getResourceAsStream(template.path);
      if (templateStream == null) {
        continue;
      }

      try (InputStreamReader reader =
          new InputStreamReader(templateStream, StandardCharsets.UTF_8)) {
        builder.clear();
        if (SettingsManager.fromJson(reader, builder, 0
        // SettingsManager.FLAG_SILENT
        )) {
          count++;
          addTemplate(templates, template.name, builder.build(), TemplateType.RESOURCE, null,
              template.tifPath);
        }
      } catch (final IOException ex) {
        // Ignore
      }
    }
    return count;
  }

  /**
   * Adds the template.
   *
   * @param map the map
   * @param name the name
   * @param settings the settings
   * @param templateType the template type
   * @param file the file
   * @param tifPath the tif path
   * @return the template
   */
  private static Template addTemplate(Map<String, Template> map, String name,
      TemplateSettings settings, TemplateType templateType, File file, String tifPath) {
    final Template template = new Template(settings, templateType, file, tifPath);
    map.put(name, template);
    return template;
  }

  /**
   * Get the template configuration.
   *
   * @param name The name of the template
   * @return The template
   */
  public static TemplateSettings getTemplate(String name) {
    final Template template = templates.get(name);
    if (template == null) {
      return null;
    }

    template.update();

    return template.settings;
  }

  /**
   * Gets the template file.
   *
   * @param name the name
   * @return the template file
   */
  public static File getTemplateFile(String name) {
    final Template template = templates.get(name);
    if (template == null) {
      return null;
    }
    return template.file;
  }

  /**
   * Gets the template image.
   *
   * @param name the name
   * @return the template image
   */
  public static ImagePlus getTemplateImage(String name) {
    final Template template = templates.get(name);
    if (template == null) {
      return null;
    }
    return template.loadImage();
  }

  /**
   * Clear templates. Used for testing so made package level.
   */
  static void clearTemplates() {
    templates.clear();
  }

  /**
   * Save template configuration. If an existing template exists with the same name it will be
   * over-written. If an existing template was loaded from file it will be saved back to the same
   * file, or optionally a new file.
   *
   * @param name The name of the template
   * @param settings The template settings
   * @param file The file to save the template (over-riding the file the template was loaded from)
   * @return true, if successful
   */
  public static boolean saveTemplate(String name, TemplateSettings settings, File file) {
    Template template = templates.get(name);
    if (template != null
        // Keep the file to allow it to be loaded on start-up
        && file == null) {
      file = template.file;
    }

    // Replace any existing template with a new one
    template = new Template(settings, TemplateType.CUSTOM, file, null);

    boolean result = true;
    if (file != null) {
      result = template.save(file);
    }
    if (result) {
      // Update the loaded templates
      templates.put(name, template);
      // If the template came from a file ensure it can be restored
      if (file != null) {
        saveLoadedTemplates(templates);
      }
    }
    return result;
  }

  /**
   * Check if this is a custom template, i.e. not a standard GDSC SMLM template.
   *
   * @param name The name of the template
   * @return True if a custom template
   */
  public static boolean isCustomTemplate(String name) {
    final Template template = templates.get(name);
    return template != null && template.templateType == TemplateType.CUSTOM;
  }

  /**
   * Get the names of the available templates.
   *
   * @return The template names
   */
  public static String[] getTemplateNames() {
    return getTemplateNames(false);
  }

  /**
   * Get the names of the available templates.
   *
   * @param includeNone Set to true to include [None] in the list of names
   * @return The template names
   */
  public static String[] getTemplateNames(boolean includeNone) {
    final TurboList<String> names = new TurboList<>(templates.size() + 1);
    if (includeNone) {
      names.add("[None]");
    }
    names.addAll(templates.keySet());
    return names.toArray(new String[0]);
  }

  /**
   * Get the names of the available templates.
   *
   * @return The template names
   */
  private static List<String> getTemplateNamesAsList() {
    return new TurboList<>(templates.keySet());
  }

  /**
   * Get the names of the available templates that have an example image.
   *
   * @return The template names
   */
  public static String[] getTemplateNamesWithImage() {
    return templates.entrySet().stream().filter(e -> e.getValue().hasImage()).map(Map.Entry::getKey)
        .toArray(String[]::new);
  }

  /**
   * The template option.
   */
  //@formatter:off
  public enum TemplateOption implements NamedObject
  {
    /** Load standard templates. */
    LOAD_STANDARD_TEMPLATES{ @Override public String getName() { return "Load standard templates"; }},
    /** Load custom templates. */
    LOAD_CUSTOM_TEMPLATES{ @Override public String getName() { return "Load custom templates"; }},
    /** Remove loaded templates. */
    REMOVE_LOADED_TEMPLATES{ @Override public String getName() { return "Remove loaded templates"; }},
    /** View template. */
    VIEW_TEMPLATE{ @Override public String getName() { return "View template"; }},
    /** View image example for template. */
    VIEW_TEMPLATE_IMAGE{ @Override public String getName() { return "View image example for template"; }};

    @Override
    public String getShortName()
    {
      return getName();
    }

    /**
     * Get the template option for the number.
     *
     * @param number
     *            the number
     * @return the template option
     */
    public static TemplateOption forNumber(int number)
    {
      final TemplateOption[] values = TemplateOption.values();
      if (number < 0 || number >= values.length) {
        number = 0;
      }
      return values[number];
    }
  }
  //@formatter:on

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    final ConfigurationTemplateSettings.Builder settings =
        SettingsManager.readConfigurationTemplateSettings(0).toBuilder();

    title = "Template Manager";
    final GenericDialog gd = new GenericDialog(title);
    final String[] options = SettingsManager.getNames((Object[]) TemplateOption.values());
    gd.addChoice("Option", options, options[settings.getOption()]);
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    final int option = gd.getNextChoiceIndex();
    settings.setOption(option);
    switch (TemplateOption.forNumber(option)) {
      case LOAD_CUSTOM_TEMPLATES:
        loadSelectedCustomTemplatesFromDirectory(settings);
        break;
      case LOAD_STANDARD_TEMPLATES:
        loadSelectedStandardTemplates(settings);
        break;
      case REMOVE_LOADED_TEMPLATES:
        removeLoadedTemplates();
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
   * @param settings the settings
   */
  private static void
      loadSelectedStandardTemplates(ConfigurationTemplateSettings.Builder settings) {
    final String[] inlineNames = listInlineTemplates();
    final TemplateResource[] templateResources = listTemplateResources();
    if (templateResources.length + inlineNames.length == 0) {
      return;
    }

    final List<String> items = new TurboList<>();
    Arrays.stream(inlineNames).forEach(items::add);
    Arrays.stream(templateResources).forEach(t -> items.add(t.name));
    final MultiDialog md = new MultiDialog("Select templates", items);
    md.setSelected(settings.getSelectedStandardTemplatesList());

    md.showDialog();

    if (md.wasCancelled()) {
      return;
    }

    final List<String> selected = md.getSelectedResults();
    if (selected.isEmpty()) {
      return;
    }

    // Save
    settings.clearSelectedStandardTemplates();
    settings.addAllSelectedStandardTemplates(selected);

    int count = 0;

    // Keep a hash of those not loaded from inline resources
    final HashSet<String> remaining = new HashSet<>(selected.size());
    for (int i = 0; i < selected.size(); i++) {
      final String name = selected.get(i);
      // Try and get the template from inline resources
      final Template t = inlineTemplates.get(name);
      if (t != null) {
        count++;
        templates.put(name, t);
      } else {
        remaining.add(name);
      }
    }

    if (!remaining.isEmpty()) {
      // Build a list of resources to load
      final TurboList<TemplateResource> list = new TurboList<>(remaining.size());
      for (final TemplateResource t : templateResources) {
        if (remaining.contains(t.name)) {
          list.add(t);
        }
      }
      count += loadTemplateResources(list.toArray(new TemplateResource[list.size()]));
    }

    if (count > 0) {
      saveLoadedTemplates(templates);
    }
    IJ.showMessage("Loaded " + TextUtils.pleural(count, "standard template"));
  }

  /**
   * Load templates from directory.
   *
   * @param settings the settings
   */
  private static void
      loadSelectedCustomTemplatesFromDirectory(ConfigurationTemplateSettings.Builder settings) {
    // Allow the user to specify a configuration directory
    final String newDirectory =
        ImageJUtils.getDirectory("Template_directory", settings.getConfigurationDirectory());

    if (newDirectory == null) {
      // Cancelled dialog
      return;
    }

    settings.setConfigurationDirectory(newDirectory);

    // Search the configuration directory and add everything that is not a tif image
    // (which may be the template source image example).
    final File[] fileList = (new File(newDirectory))
        .listFiles(file -> file.isFile() && !file.getName().toLowerCase(Locale.US).endsWith("tif"));
    if (fileList == null) {
      IJ.error(TITLE, "No files in template directory: " + newDirectory);
      return;
    }

    // Sort partially numerically
    final List<String> list = Arrays.stream(fileList).map(File::getPath)
        .sorted(AlphaNumericComparator.NULL_IS_MORE_INSTANCE).collect(Collectors.toList());

    // Select
    final MultiDialog md = new MultiDialog("Select templates", list);
    md.setDisplayConverter(path -> ImageJUtils.decodePath(path)[1]);
    md.setSelected(settings.getSelectedCustomTemplatesList());

    md.showDialog();

    if (md.wasCancelled()) {
      return;
    }

    final List<String> selected = md.getSelectedResults();
    if (selected.isEmpty()) {
      return;
    }

    // Save
    settings.clearSelectedCustomTemplates();
    settings.addAllSelectedCustomTemplates(selected);

    int count = 0;
    final TemplateSettings.Builder builder = TemplateSettings.newBuilder();
    for (final String path : selected) {
      builder.clear();
      final File file = new File(path);
      if (SettingsManager.fromJson(file, builder, 0)) {
        count++;
        final String name = FileUtils.removeExtension(file.getName());
        // Assume the tif image will be detected automatically
        addTemplate(templates, name, builder.build(), TemplateType.CUSTOM, file, null);
      } else {
        ImageJPluginLoggerHelper.getDefaultLogger()
            .info(() -> "Failed to load template file: " + file);
      }
    }

    if (count > 0) {
      saveLoadedTemplates(templates);
    }
    IJ.showMessage(TITLE, "Loaded " + TextUtils.pleural(count, "custom template"));
  }

  private void removeLoadedTemplates() {
    if (templates.isEmpty()) {
      IJ.error(title, "No templates are currently loaded");
      return;
    }

    final MultiDialog md = new MultiDialog("Select templates to remove", getTemplateNamesAsList());
    md.showDialog();

    if (md.wasCancelled()) {
      return;
    }

    final List<String> selected = md.getSelectedResults();
    if (selected.isEmpty()) {
      // Nothing to do
      return;
    }

    if (selected.size() == templates.size()) {
      clearTemplates();
    } else {
      for (final String name : selected) {
        templates.remove(name);
      }
    }
    saveLoadedTemplates(templates);
  }

  /**
   * Show template.
   *
   * @param settings the settings
   */
  private void showTemplate(ConfigurationTemplateSettings.Builder settings) {
    if (templates.isEmpty()) {
      IJ.error(title, "No templates are currently loaded");
      return;
    }
    final String[] names = getTemplateNames();

    final NonBlockingGenericDialog gd = new NonBlockingGenericDialog(title);
    gd.addMessage("View the template");
    gd.addChoice("Template", names, settings.getTemplate());
    gd.addCheckbox("Close_on_exit", settings.getClose());
    gd.hideCancelButton();
    gd.addDialogListener(this::dialogItemChanged);

    // Show the first template
    final String template = ((Choice) (gd.getChoices().get(0))).getSelectedItem();
    showTemplate(template);

    gd.showDialog();

    // There is no cancel so read the settings.
    settings.setTemplate(gd.getNextChoice());
    settings.setClose(gd.getNextBoolean());

    if (settings.getClose()) {
      closeInfo();
    }
  }

  /**
   * Show template image.
   *
   * @param name the name
   */
  private void showTemplate(String name) {
    final Template template = templates.get(name);
    if (template == null) {
      IJ.error(title, "Failed to load template: " + name);
      return;
    }

    template.update();

    if (infoWindow == null || !infoWindow.isVisible()) {
      infoWindow = new TextWindow(title + " Info", "", "", 450, 600);
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

  private void add(String key, String value) {
    if (value == null) {
      return;
    }
    infoWindow.append(key + " : " + value);
  }

  /**
   * Show template images.
   *
   * @param settings the settings
   */
  private void showTemplateImages(ConfigurationTemplateSettings.Builder settings) {
    title = "Template Example Images";

    final String[] names = getTemplateNamesWithImage();
    if (names.length == 0) {
      IJ.error(title, "No templates with example images");
      return;
    }

    // Follow when the image slice is changed
    ImageListener listener;
    if (imp != null) {
      listener = new ImageAdapter() {
        @Override
        public void imageUpdated(ImagePlus imp) {
          if (imp == ConfigurationTemplate.this.imp) {
            updateResults(imp.getCurrentSlice());
          }
        }
      };
      ImagePlus.addImageListener(listener);
    } else {
      listener = null;
    }

    final NonBlockingGenericDialog gd = new NonBlockingGenericDialog(title);
    gd.addMessage("View the example source image");
    gd.addChoice("Template", names, settings.getTemplate());
    gd.addCheckbox("Close_on_exit", settings.getClose());
    gd.hideCancelButton();
    templateImage = true;
    gd.addDialogListener(this::dialogItemChanged);

    // Show the first template
    final String template = ((Choice) (gd.getChoices().get(0))).getSelectedItem();
    showTemplateImage(template);

    gd.showDialog();

    // There is no cancel so read the settings.
    settings.setTemplate(gd.getNextChoice());
    settings.setClose(gd.getNextBoolean());

    // This is null safe so always do this
    ImagePlus.removeImageListener(listener);

    if (settings.getClose()) {
      if (imp != null) {
        imp.close();
      }
      closeResults();
      closeInfo();
    }
  }

  /**
   * Show template image.
   *
   * @param name the name
   */
  private void showTemplateImage(String name) {
    final ImagePlus tmpImp = getTemplateImage(name);
    if (tmpImp == null) {
      IJ.error(title, "Failed to load example image for template: " + name);
    } else {
      final WindowOrganiser windowOrganiser = new WindowOrganiser();
      this.imp = displayTemplate(title, tmpImp, null);
      if (windowOrganiser.isNotEmpty()) {
        // Zoom a bit
        final ImageWindow iw = this.imp.getWindow();
        for (int i = 7; i-- > 0 && Math.max(iw.getWidth(), iw.getHeight()) < 512;) {
          iw.getCanvas().zoomIn(0, 0);
        }
      }
      createResults(this.imp);

      showTemplateInfo(name);
    }
  }

  /**
   * Display the template image in an image window with the specified title. If the window exists it
   * will be reused and the appropriate properties updated.
   *
   * @param title the title
   * @param templateImp the template image
   * @param windowOrganiser the window organiser (new windows will be added to this if not null)
   * @return the image plus
   */
  public static ImagePlus displayTemplate(String title, ImagePlus templateImp,
      WindowOrganiser windowOrganiser) {
    final ImagePlus imp = ImageJUtils.display(title, templateImp.getStack(), windowOrganiser);
    imp.setOverlay(templateImp.getOverlay());
    imp.setProperty("Info", templateImp.getProperty("Info"));
    imp.setCalibration(templateImp.getCalibration());
    return imp;
  }

  private boolean dialogItemChanged(@SuppressWarnings("unused") GenericDialog gd, AWTEvent event) {
    if (event != null && event.getSource() instanceof Choice) {
      final String template = ((Choice) (event.getSource())).getSelectedItem();
      if (templateImage) {
        showTemplateImage(template);
      } else {
        showTemplate(template);
      }
    }
    return true;
  }

  /**
   * Creates a results window showing the localisation results from a template image. This will be
   * positioned next to the input template image plus if it is currently displayed.
   *
   * @param templateImp the template image
   * @return the text window
   */
  public TextWindow createResults(ImagePlus templateImp) {
    if (title == null) {
      title = templateImp.getTitle();
    }
    templateId = templateImp.getID();
    currentSlice = 0;
    headings = "";
    text = new TIntObjectHashMap<>();
    final Object info = templateImp.getProperty("Info");
    if (info != null) {
      // First line is the headings
      final String[] lines = info.toString().split("\n");
      headings = lines[0].replace(' ', '\t');

      // The remaining lines are the data for each stack position
      final StringBuilder sb = new StringBuilder();
      int last = 0;
      for (int i = 1; i < lines.length; i++) {
        // Get the position
        final String[] data = lines[i].split(" ");
        final int slice = Integer.parseInt(data[0]);
        if (last != slice) {
          text.put(last, sb.toString());
          last = slice;
          sb.setLength(0);
        }
        sb.append(slice);
        for (int j = 1; j < data.length; j++) {
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
   * @param slice the slice
   * @return the text window
   */
  public TextWindow updateResults(int slice) {
    if (slice == currentSlice || text == null) {
      return resultsWindow;
    }
    currentSlice = slice;

    if (resultsWindow == null || !resultsWindow.isVisible()) {
      resultsWindow = new TextWindow(title + " Results", headings, "", 450, 250);
      // Put next to the image
      final ImagePlus tmpImp = WindowManager.getImage(templateId);
      if (tmpImp != null && tmpImp.getWindow() != null) {
        final ImageWindow iw = tmpImp.getWindow();
        final Point p = iw.getLocation();
        p.x += iw.getWidth();
        resultsWindow.setLocation(p);
      }
    }

    resultsWindow.getTextPanel().clear();
    final String data = text.get(slice);
    if (!TextUtils.isNullOrEmpty(data)) {
      resultsWindow.append(data);
    }
    return resultsWindow;
  }

  /**
   * Close results.
   */
  public void closeResults() {
    if (resultsWindow != null) {
      resultsWindow.close();
    }
  }

  /**
   * Show the info from the template.
   *
   * @param name the name
   */
  private void showTemplateInfo(String name) {
    final TemplateSettings settings = getTemplate(name);
    if (settings == null || settings.getNotesCount() == 0) {
      return;
    }
    if (infoWindow == null || !infoWindow.isVisible()) {
      infoWindow = new TextWindow(title + " Info", "", "", 450, 250);

      // Put underneath the results window
      if (resultsWindow != null) {
        final Point p = resultsWindow.getLocation();
        p.y += resultsWindow.getHeight();
        infoWindow.setLocation(p);
      }
    }

    infoWindow.getTextPanel().clear();
    for (final String note : settings.getNotesList()) {
      // Text window cannot show tabs
      infoWindow.append(note.replace('\t', ','));
    }
  }

  /**
   * Close info.
   */
  private void closeInfo() {
    if (infoWindow != null) {
      infoWindow.close();
    }
  }
}
