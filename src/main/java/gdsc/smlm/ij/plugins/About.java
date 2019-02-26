package gdsc.smlm.ij.plugins;

import uk.ac.sussex.gdsc.core.utils.UnicodeReader;
import uk.ac.sussex.gdsc.core.utils.XmlUtils;

import gdsc.smlm.utils.XStreamXmlUtils;
import gdsc.smlm.Version;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.macro.ExtensionDescriptor;
import ij.macro.Functions;
import ij.macro.MacroExtension;
import ij.plugin.PlugIn;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.LinkedList;

/**
 * Contains help dialogs for the GDSC ImageJ plugins
 */
public class About implements PlugIn, MacroExtension
{
	private static String TITLE = "GDSC SMLM ImageJ Plugins";
	public static String HELP_URL = "http://www.sussex.ac.uk/gdsc/intranet/microscopy/imagej/smlm_plugins";
	private static String YEAR = "2016";

	enum ConfigureOption
	{
		INSTALL("Install"), REMOVE("Remove"), EDIT("Edit & Install");

		private String name;

		private ConfigureOption(String name)
		{
			this.name = name;
		}

		@Override
		public String toString()
		{
			return name;
		}
	}

	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);
		
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
					"Select the toolset from the ImageJ 'More Tools' menu to load buttons on to the ImageJ menu bar.",
					ConfigureOption.INSTALL, ConfigureOption.REMOVE);
			return;
		}

		if (arg.equals("config"))
		{
			int result = installResource(
					"/gdsc/smlm/plugins.config",
					"plugins",
					"smlm.config",
					"SMLM Tools Configuration",
					"The configuration file is used to specify which plugins to display on the SMLM Tools window. Creating a custom file will need to be repeated when the available plugins change.",
					ConfigureOption.INSTALL, ConfigureOption.EDIT, ConfigureOption.REMOVE);
			// If install/remove was successful then reload the GDSC SMLM Panel if it is showing.
			if (result != -1 && SMLMTools.isFrameVisible())
			{
				SMLMTools.closeFrame();
				new SMLMTools();
			}
			return;
		}

		if (arg.equals("ext"))
		{
			setupExtensions();
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

		BufferedReader input = null;
		try
		{
			// Read the contents of the README file
			input = new BufferedReader(new UnicodeReader(readmeStream, null));
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
		finally
		{
			try
			{
				input.close();
			}
			catch (IOException e)
			{
			}
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

	/**
	 * @param resource
	 * @param ijDirectory
	 * @param destinationName
	 * @param resourceTitle
	 * @param notes
	 * @param options
	 * @return -1 on error, 0 if installed, 1 if removed
	 */
	private static int installResource(String resource, String ijDirectory, String destinationName,
			String resourceTitle, String notes, ConfigureOption... options)
	{
		Class<About> resourceClass = About.class;
		InputStream toolsetStream = resourceClass.getResourceAsStream(resource);
		if (toolsetStream == null)
			return -1;

		String dir = IJ.getDirectory(ijDirectory);
		if (dir == null)
		{
			IJ.error("Unable to locate " + ijDirectory + " directory");
			return -1;
		}

		EnumSet<ConfigureOption> opt = EnumSet.of(options[0], options);

		GenericDialog gd = new GenericDialog(TITLE);
		String filename = dir + destinationName;
		boolean fileExists = new File(filename).exists();
		StringBuilder sb = new StringBuilder();
		sb.append("Configure resource '").append(resourceTitle).append("' at:\n \n").append(filename);
		if (notes != null)
			sb.append("\n \n").append(XmlUtils.lineWrap(notes, 80, 0, null));

		gd.addMessage(sb.toString());

		// Configure the options
		String[] choices = new String[3];
		ConfigureOption[] optChoices = new ConfigureOption[choices.length];
		int count = 0;
		if (opt.contains(ConfigureOption.INSTALL))
		{
			choices[count] = ConfigureOption.INSTALL.toString();
			if (fileExists)
				choices[count] += " (overwrite)";
			optChoices[count] = ConfigureOption.INSTALL;
			count++;
		}
		if (opt.contains(ConfigureOption.EDIT))
		{
			choices[count] = ConfigureOption.EDIT.toString();
			if (fileExists)
				choices[count] += " (overwrite)";
			optChoices[count] = ConfigureOption.EDIT;
			count++;
		}
		if (opt.contains(ConfigureOption.REMOVE) && fileExists)
		{
			choices[count] = ConfigureOption.REMOVE.toString();
			optChoices[count] = ConfigureOption.REMOVE;
			count++;
		}

		if (count == 0)
			return -1;
		choices = Arrays.copyOf(choices, count);
		gd.addChoice("Option", choices, choices[0]);

		gd.showDialog();

		if (gd.wasCanceled())
			return -1;

		ConfigureOption choice = optChoices[gd.getNextChoiceIndex()];

		if (choice == ConfigureOption.REMOVE)
		{
			try
			{
				new File(filename).delete();
				return 1;
			}
			catch (SecurityException e)
			{
				IJ.error("Unable to remove existing file");
			}
			return -1;
		}

		// Read the file
		LinkedList<String> contents = new LinkedList<String>();
		BufferedReader input = null;
		try
		{
			// Read
			input = new BufferedReader(new UnicodeReader(toolsetStream, null));
			String line;
			while ((line = input.readLine()) != null)
			{
				contents.add(line);
			}
		}
		catch (IOException e)
		{
			IJ.error("Unable to install " + resourceTitle + ".\n \n" + e.getMessage());
			return -1;
		}
		finally
		{
			close(input);
		}

		if (choice == ConfigureOption.EDIT)
		{
			// Allow the user to edit the file contents
			gd = new GenericDialog(TITLE);
			gd.addMessage("Edit the file contents before install:");
			sb.setLength(0);
			for (String line : contents)
				sb.append(line).append("\n");
			gd.addTextAreas(sb.toString(), null, 20, 80);
			gd.showDialog();
			if (gd.wasOKed())
			{
				contents.clear();
				String text = gd.getNextText();
				for (String line : text.split("\n"))
					contents.add(line);
			}
		}

		// Install the file
		BufferedWriter output = null;
		try
		{
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
			close(output);
		}
		return 0;
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

	private void setupExtensions()
	{
		Functions.registerExtensions(this);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.macro.MacroExtension#handleExtension(java.lang.String, java.lang.Object[])
	 */
	public String handleExtension(String name, Object[] args)
	{
		if (name == null)
			return "";
		if (name.equals("getNumberOfSpecies"))
		{
			return TraceDiffusion.getNumberOfSpecies(args);
		}
		if (name.equals("getD"))
		{
			return TraceDiffusion.getD(args);
		}
		if (name.equals("getF"))
		{
			return TraceDiffusion.getF(args);
		}
		if (name.equals("getSpecies"))
		{
			return TraceDiffusion.getSpecies(args);
		}
		return "";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.macro.MacroExtension#getExtensionFunctions()
	 */
	public ExtensionDescriptor[] getExtensionFunctions()
	{
		ArrayList<ExtensionDescriptor> list = new ArrayList<ExtensionDescriptor>(3);
		list.add(ExtensionDescriptor.newDescriptor("getNumberOfSpecies", this, MacroExtension.ARG_NUMBER +
				MacroExtension.ARG_OUTPUT));
		list.add(ExtensionDescriptor.newDescriptor("getD", this, MacroExtension.ARG_NUMBER, MacroExtension.ARG_NUMBER +
				MacroExtension.ARG_OUTPUT));
		list.add(ExtensionDescriptor.newDescriptor("getF", this, MacroExtension.ARG_NUMBER, MacroExtension.ARG_NUMBER +
				MacroExtension.ARG_OUTPUT));
		list.add(ExtensionDescriptor.newDescriptor("getSpecies", this, MacroExtension.ARG_NUMBER, MacroExtension.ARG_NUMBER +
				MacroExtension.ARG_OUTPUT, MacroExtension.ARG_NUMBER +
				MacroExtension.ARG_OUTPUT));
		return list.toArray(new ExtensionDescriptor[list.size()]);
	}
}
