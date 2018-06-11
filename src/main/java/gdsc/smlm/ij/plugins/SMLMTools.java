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

import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Panel;
import java.awt.Point;
import java.awt.ScrollPane;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;

import gdsc.core.utils.UnicodeReader;
import ij.Executer;
import ij.IJ;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GUI;
import ij.plugin.frame.PlugInFrame;

/**
 * Build a frame window to run all the GDSC SMLM ImageJ plugins defined in gdsc/smlm/plugins.config. Also add these
 * commands to the plugins menu.
 */
public class SMLMTools extends PlugInFrame implements ActionListener
{
	private static final long serialVersionUID = -5457127382849923056L;
	private static final String TITLE = "GDSC SMLM ImageJ Plugins";
	private static final String OPT_LOCATION = "SMLM_Plugins.location";

	private static PlugInFrame instance;

	// Store the screen dimension
	private static Dimension screenDimension;
	static
	{
		screenDimension = IJ.getScreenSize();
	}

	private HashMap<String, String[]> plugins = new HashMap<String, String[]>();
	private boolean addSpacer = false;

	/**
	 * Constructor.
	 * <p>
	 * Create a frame showing all the available plugins within the user [ImageJ]/plugins/smlm.config file or the default
	 * gdsc/smlm/plugins.config file.
	 */
	public SMLMTools()
	{
		super(TITLE);

		// Only allow one instance to run
		if (isFrameVisible())
		{
			if (!(instance.getTitle().equals(getTitle())))
			{
				closeFrame();
			}
			else
			{
				instance.toFront();
				return;
			}
		}

		if (!createFrame())
			return;

		instance = this;
		WindowManager.addWindow(this);

		pack();
		Point loc = Prefs.getLocation(OPT_LOCATION);
		if (loc != null)
			setLocation(loc);
		else
		{
			GUI.center(this);
		}
		if (IJ.isMacOSX())
			setResizable(false);
		setVisible(true);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.frame.PlugInFrame#windowClosing(java.awt.event.WindowEvent)
	 */
	public void windowClosing(WindowEvent e)
	{
		Prefs.saveLocation(OPT_LOCATION, getLocation());
		instance = null;
		close();
	}

	/**
	 * @return True if the instance of the SMLM Tools Frame is visible
	 */
	public static boolean isFrameVisible()
	{
		return (instance != null && instance.isVisible());
	}

	/**
	 * Close the instance of the SMLM Tools Frame
	 */
	public static void closeFrame()
	{
		if (instance != null)
		{
			Prefs.saveLocation(OPT_LOCATION, instance.getLocation());
			instance.close();
			instance = null;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.frame.PlugInFrame#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		// Do nothing. The frame has been created and the buttons run the plugins.
	}

	private boolean createFrame()
	{
		// Locate all the GDSC SMLM plugins using the plugins.config:
		InputStream readmeStream = getToolsPluginsConfig();

		ij.Menus.installPlugin("", ij.Menus.PLUGINS_MENU, "-", "", IJ.getInstance());

		// Read into memory
		ArrayList<String[]> plugins = new ArrayList<String[]>();
		int gaps = 0;
		BufferedReader input = null;
		try
		{
			input = new BufferedReader(new UnicodeReader(readmeStream, null));
			String line;
			while ((line = input.readLine()) != null)
			{
				if (line.startsWith("#"))
					continue;
				String[] tokens = line.split(",");
				if (tokens.length == 3)
				{
					// Only copy the entries from the Plugins menu
					if (!ignore(tokens))
					{
						if (!plugins.isEmpty())
						{
							// Multiple gaps indicates a new column
							if (gaps > 1)
							{
								plugins.add(new String[] { "next", "" });
							}
						}
						gaps = 0;
						plugins.add(new String[] { tokens[1].trim(), tokens[2].trim() });
					}
				}
				else
					gaps++;

				// Put a spacer between plugins if specified
				if ((tokens.length == 2 && tokens[0].startsWith("Plugins") && tokens[1].trim().equals("\"-\"")) ||
						line.length() == 0)
				{
					plugins.add(new String[] { "spacer", "" });
				}
			}
		}
		catch (IOException e)
		{
			// Ignore 
		}
		finally
		{
			if (input != null)
			{
				try
				{
					input.close();
				}
				catch (IOException e)
				{
					// Ignore
				}
			}
		}

		if (plugins.isEmpty())
			return false;

		// Arrange on a grid.
		Panel mainPanel = new Panel();
		GridBagLayout grid = new GridBagLayout();

		mainPanel.setLayout(grid);

		addSpacer = false;
		int col = 0, row = 0;
		for (String[] plugin : plugins)
		{
			if (plugin[0].equals("next"))
			{
				col++;
				row = 0;
			}
			else if (plugin[0].equals("spacer"))
				addSpacer = true;
			else
				row = addPlugin(mainPanel, grid, plugin[0], plugin[1], col, row);
		}

		// Allow scrollbars to handle small screens.
		// Appropriately size the scrollpane from the default of 100x100.
		// The preferred size is only obtained if the panel is packed.
		add(mainPanel);
		pack();
		Dimension d = mainPanel.getPreferredSize();
		remove(0); // Assume this is the only component
		
		ScrollPane scroll = new ScrollPane();
		scroll.getHAdjustable().setUnitIncrement(16);
		scroll.getVAdjustable().setUnitIncrement(16);
		scroll.add(mainPanel);
		add(scroll, BorderLayout.CENTER);

		// Scale to the screen size
		d.width = Math.min(d.width, screenDimension.width - 100);
		d.height = Math.min(d.height, screenDimension.height - 150);

		Insets insets = scroll.getInsets();
		d.width += insets.left + insets.right;
		d.height += insets.top + insets.bottom;

		if (IJ.isMacintosh())
		{
			// This is needed as the OSX scroll pane adds scrollbars when the panel 
			// is close in size to the scroll pane
			int padding = 15;
			d.width += padding;
			d.height += padding;
		}

		scroll.setPreferredSize(d);
		scroll.setSize(d);
		
		return true;
	}

	/**
	 * Selectively ignore certain plugins
	 * 
	 * @param tokens
	 *            The tokens from the plugins.config file
	 * @return true if the plugin should be ignored
	 */
	private boolean ignore(String[] tokens)
	{
		// Only copy the entries from the Plugins menu
		if (!tokens[0].startsWith("Plugins"))
			return true;

		// This plugin cannot be run unless in a macro
		if (tokens[1].contains("SMLM Macro Extensions"))
			return true;

		return false;
	}

	private static InputStream getToolsPluginsConfig()
	{
		// Look for smlm.config in the plugin directory
		String pluginsDir = IJ.getDirectory("plugins");
		String filename = pluginsDir + File.separator + "smlm.config";
		if (new File(filename).exists())
		{
			try
			{
				return new FileInputStream(filename);
			}
			catch (FileNotFoundException e)
			{
				// Ignore and resort to default
			}
		}

		// Fall back to the embedded config in the jar file
		return getPluginsConfig();
	}

	public static InputStream getPluginsConfig()
	{
		// Get the embedded config in the jar file
		Class<SMLMTools> resourceClass = SMLMTools.class;
		InputStream readmeStream = resourceClass.getResourceAsStream("/gdsc/smlm/plugins.config");
		return readmeStream;
	}

	private int addPlugin(Panel mainPanel, GridBagLayout grid, String commandName, final String command, int col,
			int row)
	{
		// Disect the ImageJ plugins.config string, e.g.:
		// Plugins>GDSC SMLM, "Peak Fit", gdsc.smlm.ij.plugins.PeakFit

		commandName = commandName.replaceAll("\"", "");
		Button button = new Button(commandName);
		String className = command;
		String arg = "";
		int index = command.indexOf('(');
		if (index > 0)
		{
			className = command.substring(0, index);
			int argStart = command.indexOf('"');
			if (argStart > 0)
			{
				int argEnd = command.lastIndexOf('"');
				arg = command.substring(argStart + 1, argEnd);
			}
		}

		// Add to Plugins menu so that the macros/toolset will work
		if (!ij.Menus.commandInUse(commandName))
		{
			if (addSpacer)
			{
				try
				{
					ij.Menus.getImageJMenu("Plugins").addSeparator();
				}
				catch (NoSuchMethodError e)
				{
					// Ignore. This ImageJ method is from IJ 1.48+
				}
			}
			ij.Menus.installPlugin(command, ij.Menus.PLUGINS_MENU, commandName, "", IJ.getInstance());
		}

		// Store the command to be invoked when the button is clicked
		plugins.put(commandName, new String[] { className, arg });
		button.addActionListener(this);

		if (addSpacer)
		{
			addSpacer = false;
			if (row != 0)
				row = add(mainPanel, grid, new Panel(), col, row);
		}

		row = add(mainPanel, grid, button, col, row);

		return row;
	}

	private int add(Panel mainPanel, GridBagLayout grid, Component comp, int col, int row)
	{
		GridBagConstraints c = new GridBagConstraints();
		c.gridx = col;
		c.gridy = row++;
		c.fill = GridBagConstraints.BOTH;
		if (col > 0)
			c.insets.left = 10;
		grid.setConstraints(comp, c);
		mainPanel.add(comp);
		return row;
	}

	public void actionPerformed(ActionEvent e)
	{
		// Get the plugin from the button label and run it
		Button button = (Button) e.getSource();
		String commandName = button.getLabel();

		//String[] args = plugins.get(commandName);
		//IJ.runPlugIn(commandName, args[0], args[1]); // Only in IJ 1.47+
		//IJ.runPlugIn(args[0], args[1]);

		// Use the IJ executer to run in a background thread
		new Executer(commandName, null);
	}
}
