package gdsc.smlm.ij.plugins;

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

import gdsc.smlm.utils.UnicodeReader;
import gdsc.smlm.utils.XmlUtils;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.util.LinkedList;

/**
 * Contains help dialogs for the GDSC ImageJ plugins
 */
public class About implements PlugIn
{
	private static String TITLE = "GDSC SMLM ImageJ Plugins";
	public static String HELP_URL = "http://www.sussex.ac.uk/gdsc/intranet/microscopy/imagej/smlm_plugins";
	private static String YEAR = "2014";

	public void run(String arg)
	{
		if (arg.equals("about"))
		{
			showAbout();
			return;
		}

		if (arg.equals("uninstall"))
		{
			showUnintallDialog();
			return;
		}

		if (arg.equals("toolset"))
		{
			installResource("/macros/toolsets/SMLM Tools.txt", "macros",
					"toolsets" + File.separator + "SMLM Tools.txt", "SMLM toolset",
					"Select the toolset from the ImageJ 'More Tools' menu to load buttons on to the ImageJ menu bar.");
			return;
		}

		if (arg.equals("config"))
		{
			installResource("/gdsc/smlm/plugins.config", "plugins", "smlm.config", "SMLM Tools Configuration",
					"The configuration file is used to specify which plugins to display on the SMLM Tools window.");
			return;
		}

		showAbout();
	}

	public static void showUnintallDialog()
	{
		IJ.showMessage(TITLE, "To uninstall this plugin, move the SMLM jar out\n"
				+ "of the plugins folder and restart ImageJ.");
	}

	public static void showAbout()
	{
		// Locate the README.txt file and load that into the dialog. Include revision
		Class<About> resourceClass = About.class;
		InputStream readmeStream = resourceClass.getResourceAsStream("/gdsc/smlm/README.txt");

		StringBuilder msg = new StringBuilder();
		String helpURL = HELP_URL;
		String version = Version.getVersion();
		String buildDate = Version.getBuildDate();

		try
		{
			// Read the contents of the README file
			BufferedReader input = new BufferedReader(new UnicodeReader(readmeStream, null));
			String line;
			while ((line = input.readLine()) != null)
			{
				if (line.contains("http:"))
				{
					helpURL = line;
				}
				else
				{
					if (line.equals(""))
						line = " "; // Required to insert a line in the GenericDialog
					msg.append(line).append("\n");
				}
			}
		}
		catch (IOException e)
		{
			// Default message
			msg.append("GDSC SMLM Plugins for ImageJ\n");
			msg.append(" \n");
			msg.append("Copyright (C) ").append(YEAR).append(" Alex Herbert\n");
			msg.append("MRC Genome Damage and Stability Centre\n");
			msg.append("University of Sussex, UK\n");
		}

		// Build final message
		msg = new StringBuilder(msg.toString().trim());
		if (version != Version.UNKNOWN || buildDate != Version.UNKNOWN)
			msg.append("\n \n");
		if (version != Version.UNKNOWN)
			msg.append("Version : ").append(version).append("\n");
		if (buildDate != Version.UNKNOWN)
			msg.append("Build Date : ").append(buildDate).append("\n");
		if (helpURL != null)
			msg.append("\n \n(Click help for more information)");

		GenericDialog gd = new GenericDialog(TITLE);
		gd.addMessage(msg.toString());
		gd.addHelp(helpURL);
		gd.hideCancelButton();
		gd.showDialog();
	}

	private static void installResource(String resource, String ijDirectory, String destinationName,
			String resourceTitle, String notes)
	{
		Class<About> resourceClass = About.class;
		InputStream toolsetStream = resourceClass.getResourceAsStream(resource);
		if (toolsetStream == null)
			return;

		String dir = IJ.getDirectory(ijDirectory);
		if (dir == null)
		{
			IJ.error("Unable to locate " + ijDirectory + " directory");
		}
		GenericDialog gd = new GenericDialog(TITLE);
		String filename = dir + destinationName;
		boolean fileExists = new File(filename).exists();
		StringBuilder sb = new StringBuilder();
		sb.append("Install ").append(resourceTitle).append(" to:\n \n").append(filename);
		if (notes != null)
			sb.append("\n \n").append(XmlUtils.lineWrap(notes, 80, 0, null));
		if (fileExists)
		{
			sb.append("\n \n* Note this will over-write the existing file.");
		}
		gd.addMessage(sb.toString());
		if (fileExists)
			gd.addCheckbox("Remove existing file", false);
		gd.showDialog();

		if (gd.wasCanceled())
			return;

		if (fileExists && gd.getNextBoolean())
		{
			try
			{
				new File(filename).delete();
			}
			catch (SecurityException e)
			{
				IJ.error("Unable to remove existing file");
			}
			return;
		}

		BufferedReader input = null;
		BufferedWriter output = null;
		try
		{
			// Copy the contents of the file

			// Read
			LinkedList<String> contents = new LinkedList<String>();
			input = new BufferedReader(new UnicodeReader(toolsetStream, null));
			String line;
			while ((line = input.readLine()) != null)
			{
				contents.add(line);
			}

			// Write
			FileOutputStream fos = new FileOutputStream(filename);
			output = new BufferedWriter(new OutputStreamWriter(fos, "UTF-8"));
			for (String content : contents)
			{
				output.write(content);
				output.newLine();
			}
		}
		catch (IOException e)
		{
			IJ.error("Unable to install " + resourceTitle + ".\n \n" + e.getMessage());
		}
		finally
		{
			close(input);
			close(output);
		}
	}

	private static void close(BufferedReader input)
	{
		if (input != null)
		{
			try
			{
				input.close();
			}
			catch (IOException e)
			{
			}
		}
	}

	private static void close(BufferedWriter output)
	{
		if (output != null)
		{
			try
			{
				output.close();
			}
			catch (IOException e)
			{
			}
		}
	}
}
