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

import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.UnicodeReader;
import uk.ac.sussex.gdsc.smlm.Version;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.macro.ExtensionDescriptor;
import ij.macro.Functions;
import ij.macro.MacroExtension;
import ij.plugin.PlugIn;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.LinkedList;

/**
 * Contains help dialogs for the GDSC ImageJ plugins.
 */
public class About implements PlugIn, MacroExtension {
  private static final String TITLE = "GDSC SMLM ImageJ Plugins";
  /** The help url for the SMLM plugins. */
  public static final String HELP_URL =
      "http://www.sussex.ac.uk/gdsc/intranet/microscopy/imagej/smlm_plugins";
  private static final String YEAR = "2019";

  /**
   * The configure option.
   */
  enum ConfigureOption {
    //@formatter:off
    /** Install. */
    INSTALL{ @Override
    public String getName() { return "Install"; }},
    /** Remove. */
    REMOVE{ @Override
    public String getName() { return "Remove"; }},
    /** Edit. */
    EDIT{ @Override
    public String getName() { return "Edit & Install"; }};
    //@formatter:on

    @Override
    public String toString() {
      return getName();
    }

    /**
     * Gets the name.
     *
     * @return the name
     */
    public abstract String getName();
  }

  @SuppressWarnings("unused")
  @Override
  public void run(String arg) {
    SMLMUsageTracker.recordPlugin(this.getClass(), arg);

    if (arg.equals("about")) {
      showAbout();
      return;
    }

    if (arg.equals("uninstall")) {
      showUninstallDialog();
      return;
    }

    if (arg.equals("toolset")) {
      installResource("/macros/toolsets/SMLM Tools.txt", "macros",
          "toolsets" + File.separator + "SMLM Tools.txt", "SMLM toolset",
          "Select the toolset from the ImageJ 'More Tools' menu to load buttons on to the "
              + "ImageJ menu bar.",
          ConfigureOption.INSTALL, ConfigureOption.REMOVE);
      return;
    }

    if (arg.equals("config")) {
      final int result = installResource("/uk/ac/sussex/gdsc/smlm/plugins.config", "plugins",
          "smlm.config", "SMLM Tools Configuration",
          "The configuration file is used to specify which plugins to display on the SMLM Tools "
              + "window. Creating a custom file will need to be repeated when the available "
              + "plugins change.",
          ConfigureOption.INSTALL, ConfigureOption.EDIT, ConfigureOption.REMOVE);
      // If install/remove was successful then reload the GDSC SMLM Panel if it is showing.
      if (result != -1 && SMLMTools.isFrameVisible()) {
        SMLMTools.closeFrame();
        new SMLMTools();
      }
      return;
    }

    if (arg.equals("ext")) {
      setupExtensions();
      return;
    }

    showAbout();
  }

  /**
   * Show uninstall dialog.
   */
  public static void showUninstallDialog() {
    IJ.showMessage(TITLE, "To uninstall this plugin, move the SMLM jar out\n"
        + "of the plugins folder and restart ImageJ.");
  }

  /**
   * Show about dialog.
   */
  public static void showAbout() {
    // Locate the README.txt file and load that into the dialog. Include revision
    final Class<About> resourceClass = About.class;

    StringBuilder msg = new StringBuilder();
    String helpURL = HELP_URL;

    try (BufferedReader input = new BufferedReader(new UnicodeReader(
        resourceClass.getResourceAsStream("/uk/ac/sussex/gdsc/smlm/README.txt"), null))) {
      // Read the contents of the README file
      String line;
      while ((line = input.readLine()) != null) {
        if (line.contains("http:")) {
          helpURL = line;
        } else {
          if (line.equals("")) {
            line = " "; // Required to insert a line in the GenericDialog
          }
          msg.append(line).append('\n');
        }
      }
    } catch (final IOException ex) {
      // Default message
      msg.append("GDSC SMLM Plugins for ImageJ\n");
      msg.append(" \n");
      msg.append("Copyright (C) ").append(YEAR).append(" Alex Herbert\n");
      msg.append("MRC Genome Damage and Stability Centre\n");
      msg.append("University of Sussex, UK\n");
    }

    // Build final message
    msg = new StringBuilder(msg.toString().trim());
    addVersion(msg, "GDSC-SMLM", Version.getVersion(), Version.getBuildDate(),
        Version.getBuildNumber());
    addVersion(msg, "GDSC-Core", uk.ac.sussex.gdsc.core.VersionUtils.getVersion(),
        uk.ac.sussex.gdsc.core.VersionUtils.getBuildDate(),
        uk.ac.sussex.gdsc.core.VersionUtils.getBuildNumber());

    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addMessage(msg.toString());
    gd.addHelp(helpURL);
    gd.hideCancelButton();
    gd.showDialog();
  }

  private static void addVersion(StringBuilder msg, String name, String version, String buildDate,
      String buildNumber) {
    final boolean hasVersion = !TextUtils.isNullOrEmpty(version);
    final boolean hasBuildDate = !TextUtils.isNullOrEmpty(buildDate);
    final boolean hasBuildNumber = !TextUtils.isNullOrEmpty(buildNumber);
    if (hasVersion || hasBuildDate || hasBuildNumber) {
      msg.append("\n \n").append(name).append("\n");
    }
    if (hasVersion) {
      msg.append("Version : ").append(version).append("\n");
    }
    if (hasBuildDate) {
      msg.append("Build Date : ").append(buildDate).append("\n");
    }
    if (hasBuildNumber) {
      msg.append("Build Number : ").append(buildNumber).append("\n");
    }
  }

  /**
   * Install resource.
   *
   * @param resource the resource
   * @param ijDirectory the ij directory
   * @param destinationName the destination name
   * @param resourceTitle the resource title
   * @param notes the notes
   * @param options the options
   * @return -1 on error, 0 if installed, 1 if removed
   */
  private static int installResource(String resource, String ijDirectory, String destinationName,
      String resourceTitle, String notes, ConfigureOption... options) {
    final Class<About> resourceClass = About.class;

    final String dir = IJ.getDirectory(ijDirectory);
    if (dir == null) {
      IJ.error("Unable to locate " + ijDirectory + " directory");
      return -1;
    }

    final EnumSet<ConfigureOption> opt = EnumSet.of(options[0], options);

    GenericDialog gd = new GenericDialog(TITLE);
    final String filename = dir + destinationName;
    final boolean fileExists = new File(filename).exists();
    final StringBuilder sb = new StringBuilder();
    sb.append("Configure resource '").append(resourceTitle).append("' at:\n \n").append(filename);
    if (notes != null) {
      sb.append("\n \n").append(uk.ac.sussex.gdsc.core.utils.XmlUtils.lineWrap(notes, 80, 0, null));
    }

    gd.addMessage(sb.toString());

    // Configure the options
    String[] choices = new String[3];
    final ConfigureOption[] optChoices = new ConfigureOption[choices.length];
    int count = 0;
    if (opt.contains(ConfigureOption.INSTALL)) {
      choices[count] = ConfigureOption.INSTALL.toString();
      if (fileExists) {
        choices[count] += " (overwrite)";
      }
      optChoices[count] = ConfigureOption.INSTALL;
      count++;
    }
    if (opt.contains(ConfigureOption.EDIT)) {
      choices[count] = ConfigureOption.EDIT.toString();
      if (fileExists) {
        choices[count] += " (overwrite)";
      }
      optChoices[count] = ConfigureOption.EDIT;
      count++;
    }
    if (opt.contains(ConfigureOption.REMOVE) && fileExists) {
      choices[count] = ConfigureOption.REMOVE.toString();
      optChoices[count] = ConfigureOption.REMOVE;
      count++;
    }

    if (count == 0) {
      return -1;
    }
    choices = Arrays.copyOf(choices, count);
    gd.addChoice("Option", choices, choices[0]);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return -1;
    }

    final ConfigureOption choice = optChoices[gd.getNextChoiceIndex()];

    if (choice == ConfigureOption.REMOVE) {
      try {
        Files.delete(Paths.get(filename));
        return 1;
      } catch (final SecurityException | IOException ex) {
        IJ.error("Unable to remove existing file");
      }
      return -1;
    }

    // Read the file
    final LinkedList<String> contents = new LinkedList<>();
    try (BufferedReader input =
        new BufferedReader(new UnicodeReader(resourceClass.getResourceAsStream(resource), null))) {
      String line;
      while ((line = input.readLine()) != null) {
        contents.add(line);
      }
    } catch (final IOException ex) {
      IJ.error("Unable to install " + resourceTitle + ".\n \n" + ex.getMessage());
      return -1;
    }

    if (choice == ConfigureOption.EDIT) {
      // Allow the user to edit the file contents
      gd = new GenericDialog(TITLE);
      gd.addMessage("Edit the file contents before install:");
      sb.setLength(0);
      for (final String line : contents) {
        sb.append(line).append('\n');
      }
      gd.addTextAreas(sb.toString(), null, 20, 80);
      gd.showDialog();
      if (gd.wasOKed()) {
        contents.clear();
        final String text = gd.getNextText();
        for (final String line : text.split("\n")) {
          contents.add(line);
        }
      }
    }

    // Install the file
    try (BufferedWriter output = Files.newBufferedWriter(Paths.get(filename))) {
      for (final String content : contents) {
        output.write(content);
        output.newLine();
      }
    } catch (final IOException ex) {
      IJ.error("Unable to install " + resourceTitle + ".\n \n" + ex.getMessage());
    }
    return 0;
  }

  private void setupExtensions() {
    Functions.registerExtensions(this);
  }

  @Override
  public String handleExtension(String name, Object[] args) {
    if (name == null) {
      return "";
    }
    if (name.equals("getNumberOfSpecies")) {
      return TraceDiffusion.getNumberOfSpecies(args);
    }
    if (name.equals("getD")) {
      return TraceDiffusion.getD(args);
    }
    if (name.equals("getF")) {
      return TraceDiffusion.getF(args);
    }
    if (name.equals("getSpecies")) {
      return TraceDiffusion.getSpecies(args);
    }
    return "";
  }

  @Override
  public ExtensionDescriptor[] getExtensionFunctions() {
    final ArrayList<ExtensionDescriptor> list = new ArrayList<>(3);
    list.add(ExtensionDescriptor.newDescriptor("getNumberOfSpecies", this,
        MacroExtension.ARG_NUMBER + MacroExtension.ARG_OUTPUT));
    list.add(ExtensionDescriptor.newDescriptor("getD", this, MacroExtension.ARG_NUMBER,
        MacroExtension.ARG_NUMBER + MacroExtension.ARG_OUTPUT));
    list.add(ExtensionDescriptor.newDescriptor("getF", this, MacroExtension.ARG_NUMBER,
        MacroExtension.ARG_NUMBER + MacroExtension.ARG_OUTPUT));
    list.add(ExtensionDescriptor.newDescriptor("getSpecies", this, MacroExtension.ARG_NUMBER,
        MacroExtension.ARG_NUMBER + MacroExtension.ARG_OUTPUT,
        MacroExtension.ARG_NUMBER + MacroExtension.ARG_OUTPUT));
    return list.toArray(new ExtensionDescriptor[list.size()]);
  }
}
