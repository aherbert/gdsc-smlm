/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

import ij.Executer;
import ij.IJ;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GUI;
import ij.plugin.frame.PlugInFrame;
import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Panel;
import java.awt.Point;
import java.awt.ScrollPane;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicReference;
import uk.ac.sussex.gdsc.core.utils.UnicodeReader;

/**
 * Build a frame window to run all the GDSC SMLM ImageJ plugins defined in
 * /uk/ac/sussex/gdsc/smlm/plugins.config. Also add these commands to the plugins menu.
 */
public class SmlmTools extends PlugInFrame {
  private static final long serialVersionUID = -5457127382849923056L;
  private static final String TITLE = "GDSC SMLM ImageJ Plugins";
  private static final String OPT_LOCATION = "SMLM_Plugins.location";

  private static final AtomicReference<PlugInFrame> instance = new AtomicReference<>();

  // Store the screen dimension
  private static Dimension screenDimension;

  static {
    screenDimension = IJ.getScreenSize();
  }

  private final HashMap<String, String[]> plugins = new HashMap<>();
  private boolean addSpacer;

  /**
   * Constructor.
   *
   * <p>Create a frame showing all the available plugins within the user
   * [ImageJ]/plugins/smlm.config file or the default /uk/ac/sussex/gdsc/smlm/plugins.config file.
   */
  public SmlmTools() {
    super(TITLE);

    // Only allow one instance to run
    final Frame frame = instance.get();

    if (frame != null) {
      frame.toFront();
      return;
    }

    if (!createFrame()) {
      return;
    }

    instance.set(this);
    WindowManager.addWindow(this);

    pack();
    final Point loc = Prefs.getLocation(OPT_LOCATION);
    if (loc != null) {
      setLocation(loc);
    } else {
      GUI.center(this);
    }
    if (IJ.isMacOSX()) {
      setResizable(false);
    }
    setVisible(true);
  }

  @Override
  public void close() {
    Prefs.saveLocation(OPT_LOCATION, getLocation());
    instance.compareAndSet(this, null);
    super.close();
  }

  /**
   * Checks if the instance of the SMLM Tools Frame is visible.
   *
   * @return True if the instance of the SMLM Tools Frame is visible.
   */
  public static boolean isFrameVisible() {
    final Frame frame = instance.get();
    return (frame != null && frame.isVisible());
  }

  /**
   * Close the instance of the SMLM Tools Frame.
   */
  public static void closeFrame() {
    final PlugInFrame frame = instance.getAndUpdate(fr -> null);
    if (frame != null) {
      frame.close();
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    // Do nothing. The frame has been created and the buttons run the plugins.
  }

  private boolean createFrame() {
    final ArrayList<String[]> configPlugins = new ArrayList<>();

    // Locate all the GDSC SMLM plugins using the plugins.config:
    try (InputStream readmeStream = getToolsPluginsConfig()) {
      // Read into memory
      int gaps = 0;
      try (BufferedReader input = new BufferedReader(new UnicodeReader(readmeStream, null))) {
        String line;
        while ((line = input.readLine()) != null) {
          if (line.startsWith("#")) {
            continue;
          }
          final String[] tokens = line.split(",");
          if (tokens.length == 3) {
            // Only copy the entries from the Plugins menu
            if (!ignore(tokens)) {
              if (!configPlugins.isEmpty()
                  // Multiple gaps indicates a new column
                  && gaps > 1) {
                configPlugins.add(new String[] {"next", ""});
              }
              gaps = 0;
              configPlugins.add(new String[] {tokens[1].trim(), tokens[2].trim()});
            }
          } else {
            gaps++;
          }

          // Put a spacer between plugins if specified
          if ((tokens.length == 2 && tokens[0].startsWith("Plugins")
              && tokens[1].trim().equals("\"-\"")) || line.length() == 0) {
            configPlugins.add(new String[] {"spacer", ""});
          }
        }
      }
    } catch (final IOException ex) {
      // Ignore
    }

    if (configPlugins.isEmpty()) {
      return false;
    }

    // Put a spacer on the menu
    ij.Menus.installPlugin("", ij.Menus.PLUGINS_MENU, "-", "", IJ.getInstance());

    // Arrange on a grid.
    final Panel mainPanel = new Panel();
    final GridBagLayout grid = new GridBagLayout();

    mainPanel.setLayout(grid);

    addSpacer = false;
    int col = 0;
    int row = 0;
    for (final String[] plugin : configPlugins) {
      if (plugin[0].equals("next")) {
        col++;
        row = 0;
      } else if (plugin[0].equals("spacer")) {
        addSpacer = true;
      } else {
        row = addPlugin(mainPanel, grid, plugin[0], plugin[1], col, row);
      }
    }

    // Allow scrollbars to handle small screens.
    // Appropriately size the scrollpane from the default of 100x100.
    // The preferred size is only obtained if the panel is packed.
    add(mainPanel);
    pack();
    final Dimension d = mainPanel.getPreferredSize();
    remove(0); // Assume this is the only component

    final ScrollPane scroll = new ScrollPane();
    scroll.getHAdjustable().setUnitIncrement(16);
    scroll.getVAdjustable().setUnitIncrement(16);
    scroll.add(mainPanel);
    add(scroll, BorderLayout.CENTER);

    // Scale to the screen size
    d.width = Math.min(d.width, screenDimension.width - 100);
    d.height = Math.min(d.height, screenDimension.height - 150);

    final Insets insets = scroll.getInsets();
    d.width += insets.left + insets.right;
    d.height += insets.top + insets.bottom;

    if (IJ.isMacintosh()) {
      // This is needed as the OSX scroll pane adds scrollbars when the panel
      // is close in size to the scroll pane
      final int padding = 15;
      d.width += padding;
      d.height += padding;
    }

    scroll.setPreferredSize(d);
    scroll.setSize(d);

    return true;
  }

  /**
   * Selectively ignore certain plugins.
   *
   * @param tokens The tokens from the plugins.config file
   * @return true if the plugin should be ignored
   */
  private static boolean ignore(String[] tokens) {
    // Only copy the entries from the Plugins menu
    if (!tokens[0].startsWith("Plugins")) {
      return true;
    }

    // This plugin cannot be run unless in a macro
    return tokens[1].contains("SMLM Macro Extensions");
  }

  private static InputStream getToolsPluginsConfig() {
    // Look for smlm.config in the plugin directory
    final String pluginsDir = IJ.getDirectory("plugins");
    final String filename = pluginsDir + File.separator + "smlm.config";
    if (new File(filename).exists()) {
      try {
        return new FileInputStream(filename);
      } catch (final FileNotFoundException ex) {
        // Ignore and resort to default
      }
    }

    // Fall back to the embedded config in the jar file
    return getPluginsConfig();
  }

  /**
   * Gets the plugins.config embedded in the jar.
   *
   * @return the plugins config
   */
  public static InputStream getPluginsConfig() {
    // Get the embedded config in the jar file
    final Class<SmlmTools> resourceClass = SmlmTools.class;
    return resourceClass.getResourceAsStream("/uk/ac/sussex/gdsc/smlm/plugins.config");
  }

  @SuppressWarnings("unused")
  private int addPlugin(Panel mainPanel, GridBagLayout grid, String commandName,
      final String command, int col, int row) {
    // Disect the ImageJ plugins.config string, e.g.:
    // Plugins>GDSC SMLM, "Peak Fit", uk.ac.sussex.gdsc.smlm.ij.plugins.PeakFit

    commandName = commandName.replaceAll("\"", "");
    final Button button = new Button(commandName);
    String className = command;
    String arg = "";
    final int index = command.indexOf('(');
    if (index > 0) {
      className = command.substring(0, index);
      final int argStart = command.indexOf('"');
      if (argStart > 0) {
        final int argEnd = command.lastIndexOf('"');
        arg = command.substring(argStart + 1, argEnd);
      }
    }

    // Add to Plugins menu so that the macros/toolset will work
    if (!ij.Menus.commandInUse(commandName)) {
      if (addSpacer) {
        try {
          ij.Menus.getImageJMenu("Plugins").addSeparator();
        } catch (final NoSuchMethodError ex) {
          // Ignore. This ImageJ method is from IJ 1.48+
        }
      }
      ij.Menus.installPlugin(command, ij.Menus.PLUGINS_MENU, commandName, "", IJ.getInstance());
    }

    // Store the command to be invoked when the button is clicked
    plugins.put(commandName, new String[] {className, arg});
    button.addActionListener(event -> {
      // Get the plugin from the button label and run it
      final Button source = (Button) event.getSource();
      final String label = source.getLabel();
      // Use the IJ executer to run in a background thread
      new Executer(label, null);
    });

    if (addSpacer) {
      addSpacer = false;
      if (row != 0) {
        row = add(mainPanel, grid, new Panel(), col, row);
      }
    }

    row = add(mainPanel, grid, button, col, row);

    return row;
  }

  private static int add(Panel mainPanel, GridBagLayout grid, Component comp, int col, int row) {
    final GridBagConstraints c = new GridBagConstraints();
    c.gridx = col;
    c.gridy = row++;
    c.fill = GridBagConstraints.BOTH;
    if (col > 0) {
      c.insets.left = 10;
    }
    grid.setConstraints(comp, c);
    mainPanel.add(comp);
    return row;
  }
}
