/*----------------------------------------------------------------------------- 
 * GDSC Plugins for ImageJ
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/
package gdsc.common.ij;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;

import gdsc.common.analytics.ClientParameters;
import gdsc.common.analytics.ClientParametersManager;
import gdsc.common.analytics.ConsoleLogger;
import gdsc.common.analytics.JGoogleAnalyticsTracker;
import gdsc.common.analytics.JGoogleAnalyticsTracker.DispatchMode;
import gdsc.common.analytics.JGoogleAnalyticsTracker.MeasurementProtocolVersion;
import ij.IJ;
import ij.ImageJ;
import ij.Prefs;
import ij.gui.GenericDialog;

/**
 * Provide a global reference to JGoogleAnalyticsTracker
 * <p>
 * This class exists in both the codebase for the GDSC ImageJ Plugins and GDSC SMLM ImageJ Plugins simply to allow
 * tracking of both applications together with the same client Id and session in one instance of ImageJ.
 */
public class ImageJAnalyticsTracker
{
	private static final String PROPERTY_GA_CLIENT_ID = "gdsc.ga.clientId";
	private static final String PROPERTY_GA_STATE = "gdsc.ga.state";
	private static final String PROPERTY_GA_LAST_VERSION = "gdsc.ga.lastVersion";
	private static final String PROPERTY_GA_ANONYMIZE = "gdsc.ga.anonymize";
	private static final int DISABLED = -1;
	private static final int UNKNOWN = 0;
	private static final int ENABLED = 1;

	private static ClientParameters clientParameters = null;
	private static JGoogleAnalyticsTracker tracker = null;

	/**
	 * Flag indicating that the user has opted out of analytics
	 */
	private static int state = (int) Prefs.get(PROPERTY_GA_STATE, UNKNOWN);
	/**
	 * Flag indicating that the IP address of the sender will be anonymized
	 */
	private static int anonymized = (int) Prefs.get(PROPERTY_GA_ANONYMIZE, UNKNOWN);

	/**
	 * Track a page view
	 * 
	 * @param documentPath
	 *            The document path (must not be null)
	 * @param documentTitle
	 *            The document title
	 */
	public static void pageview(String pageUrl, String pageTitle)
	{
		initialiseTracker();
		if (isDisabled())
			return;
		tracker.pageview(pageUrl, pageTitle);
	}

	/**
	 * Create the client parameters for the tracker and populate with information about ImageJ and the system
	 */
	public static void initialise()
	{
		if (clientParameters == null)
		{
			synchronized (ImageJAnalyticsTracker.class)
			{
				// In case another thread was waiting to do this
				if (clientParameters != null)
					return;

				// Get the client parameters
				final String clientId = Prefs.get(PROPERTY_GA_CLIENT_ID, null);
				clientParameters = new ClientParameters("UA-74107394-1", clientId, "GDSC ImageJ Plugins");

				ClientParametersManager.populate(clientParameters);

				// Record for next time
				Prefs.set(PROPERTY_GA_CLIENT_ID, clientParameters.getClientId());

				// Record the version of analytics we are using
				clientParameters.setApplicationVersion(gdsc.common.analytics.Version.VERSION_X_X_X);

				// Use custom dimensions to record client data. These should be registered 
				// in the analytics account for the given tracking ID

				// Record the ImageJ information.
				clientParameters.addCustomDimension(1, getImageJInfo());

				// Java version
				clientParameters.addCustomDimension(2, System.getProperty("java.version"));

				// OS
				clientParameters.addCustomDimension(3, System.getProperty("os.name"));
				clientParameters.addCustomDimension(4, System.getProperty("os.version"));
				clientParameters.addCustomDimension(5, System.getProperty("os.arch"));
			}
		}
	}

	/**
	 * @return The ImageJ information (Application name and version)
	 */
	private static String getImageJInfo()
	{
		ImageJ ij = IJ.getInstance();
		if (ij == null)
		{
			// Run embedded without showing
			ij = new ImageJ(ImageJ.NO_SHOW);
		}

		// ImageJ version
		// This call should return a string like:
		//   ImageJ 1.48a; Java 1.7.0_11 [64-bit]; Windows 7 6.1; 29MB of 5376MB (<1%)
		// (This should also be different if we are running within Fiji)
		String info = ij.getInfo();
		if (info.indexOf(';') != -1)
			info = info.substring(0, info.indexOf(';'));
		return info;
	}

	/**
	 * Add a custom dimension
	 * <p>
	 * Note that custom dimensions have to be created for your site before they can be used in analytics reports.
	 * 
	 * @see https://support.google.com/analytics/answer/2709829
	 * 
	 * @param index
	 *            The dimension index (1-20 or 1-200 for premium accounts)
	 * @param value
	 *            The dimension value (must not be null)
	 */
	public static void addCustomDimension(int index, String value)
	{
		initialise();
		clientParameters.addCustomDimension(index, value);
	}

	/**
	 * Create the tracker
	 */
	private static void initialiseTracker()
	{
		if (tracker == null)
		{
			synchronized (ImageJAnalyticsTracker.class)
			{
				// Check again since this may be a second thread that was waiting  
				if (tracker != null)
					return;

				initialise();

				tracker = new JGoogleAnalyticsTracker(clientParameters, MeasurementProtocolVersion.V_1,
						DispatchMode.SINGLE_THREAD);

				tracker.setAnonymised(isAnonymized());

				// XXX - Disable in production code 
				// DEBUG: Enable logging
				tracker.setLogger(new ConsoleLogger());
			}
		}
	}

	/**
	 * Provide a method to read an ImageJ properties file and create the map between the ImageJ plugin class and
	 * argument and the ImageJ menu path and plugin title.
	 * 
	 * @param map
	 *            The map object to populate
	 * @param pluginsStream
	 *            The ImageJ properties file
	 */
	public static void buildPluginMap(HashMap<String, String[]> map, InputStream pluginsStream)
	{
		BufferedReader input = null;
		try
		{
			input = new BufferedReader(new InputStreamReader(pluginsStream));
			String line;
			while ((line = input.readLine()) != null)
			{
				if (line.startsWith("#"))
					continue;
				String[] tokens = line.split(",");
				if (tokens.length == 3)
				{
					// Plugins have [Menu path, Name, class(argument)], e.g.
					// Plugins>GDSC>Colocalisation, "CDA (macro)", gdsc.colocalisation.cda.CDA_Plugin("macro")

					String documentTitle = tokens[1].replaceAll("[\"']", "").trim();
					String documentPath = getDocumentPath(tokens[0], documentTitle);
					String key = getKey(tokens[2]);
					map.put(key, new String[] { documentPath, documentTitle });
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

	}

	/**
	 * Split the menu path string and create a document path
	 * 
	 * @param menuPath
	 *            The ImageJ menu path string
	 * @param documentTitle
	 * @return The document path
	 */
	private static String getDocumentPath(String menuPath, String documentTitle)
	{
		StringBuilder sb = new StringBuilder();
		for (String field : menuPath.split(">"))
		{
			sb.append('/').append(field.trim());
		}
		sb.append('/').append(documentTitle);
		return sb.toString();
	}

	/**
	 * Get the raw class name and string argument from the ImageJ 'org.package.Class("argument")' field
	 * 
	 * @param string
	 *            The field contents
	 * @return The hash key
	 */
	private static String getKey(String string)
	{
		String name = string.trim(), argument = null;
		final int index = name.indexOf('(');
		if (index != -1)
		{
			// Get the remaining text and remove the quotes " and brackets ()
			argument = name.substring(index).replaceAll("[\"'()]", "").trim();
			// Get the class name
			name = name.substring(0, index);
		}
		return getKey(name, argument);
	}

	/**
	 * Get the key used for the given name and argument in the plugin map
	 * 
	 * @param name
	 * @param argument
	 * @return The key
	 */
	public static String getKey(String name, String argument)
	{
		return (argument != null && argument.length() > 0) ? name + '.' + argument : name;
	}

	/**
	 * @return True if analytics is disabled
	 */
	public static boolean isDisabled()
	{
		return (state == DISABLED);
	}

	/**
	 * Set the state of the analytics tracker
	 * 
	 * @param disabled
	 *            True to disable analytics
	 */
	public static void setDisabled(boolean disabled)
	{
		final int oldState = ImageJAnalyticsTracker.state;
		ImageJAnalyticsTracker.state = (disabled) ? DISABLED : ENABLED;

		Prefs.set(PROPERTY_GA_LAST_VERSION, getVersion());

		if (oldState != state)
		{
			Prefs.set(PROPERTY_GA_STATE, state);

			// Record this opt in/out status change

			// If the state was previously unknown and they opt out then do nothing.
			// This is a user who never wants to be tracked.
			if (oldState == UNKNOWN && disabled)
				return;

			// Otherwise record the in/out status change, either:
			// - The user was previously opt in but now opts out; or
			// - The user was previously opt out but now opts in

			final boolean enabled = !disabled;

			initialiseTracker();
			// Reset the session if tracking has been enabled.
			if (enabled)
				tracker.resetSession();
			// Track the opt status change with an event
			tracker.event("Tracking", Boolean.toString(enabled), getVersion(), null);
		}
	}

	/**
	 * @return True if the IP address of the sender will be anonymized
	 */
	public static boolean isAnonymized()
	{
		return (anonymized == ENABLED);
	}

	/**
	 * Set the state of IP anonymization
	 * 
	 * @param anonymize
	 *            True if the IP address of the sender will be anonymized
	 */
	public static void setAnonymized(boolean anonymize)
	{
		final int oldAnonymized = anonymized;
		ImageJAnalyticsTracker.anonymized = (anonymize) ? ENABLED : DISABLED;

		Prefs.set(PROPERTY_GA_LAST_VERSION, getVersion());

		if (oldAnonymized != anonymized)
		{
			Prefs.set(PROPERTY_GA_ANONYMIZE, anonymized);

			// Make sure the tracker is informed
			if (tracker != null)
				tracker.setAnonymised(isAnonymized());
		}
	}

	/**
	 * @return The version of the code
	 */
	private static String getVersion()
	{
		return gdsc.common.analytics.Version.VERSION_X_X;
	}

	/**
	 * Check the opt-in/out status. If it is not known for this version of the analytics code then return true.
	 */
	public static boolean unknownStatus()
	{
		String lastVersion = Prefs.get(PROPERTY_GA_LAST_VERSION, "");
		String thisVersion = getVersion();
		return (state == UNKNOWN || anonymized == UNKNOWN || !lastVersion.equals(thisVersion));
	}

	/**
	 * Show a dialog allowing users to opt in/out of Google Analytics
	 * 
	 * @param title
	 *            The dialog title
	 * @param autoMessage
	 *            Set to true to display the message about automatically showing when the status is unknown
	 */
	public static void showDialog(String title, boolean autoMessage)
	{
		GenericDialog gd = new GenericDialog(title);
		// @formatter:off
		gd.addMessage("GDSC ImageJ Plugins\n \n" +
				"The use of these plugins is free and unconstrained.\n" +
				"The code uses Google Analytics to help us understand\n" +
				"how users are using the plugins.\n \n" +
				"Privacy Policy\n \n" +
				"No personal information or data within ImageJ is recorded.\n \n" +
				"We record the plugin name and the software running ImageJ.\n" +
				"This happens in the background when nothing else is active so will\n" +
				"not slow down ImageJ. IP anonymization will truncate your IP\n" +
				"address to a region, usually a country. For more details click\n" +
				"the Help button.\n \n" + 
				"Click here to opt-out from Google Analytics");
		gd.addHelp("http://www.sussex.ac.uk/gdsc/intranet/microscopy/imagej/tracking");
		gd.addCheckbox("Disabled", isDisabled());
		gd.addCheckbox("Anonymise IP", isAnonymized());
		if (autoMessage)
		{
		gd.addMessage(
				"Note: This dialog is automatically shown if your preferences\n" +
				"are not known, or after a release that changes the tracking data.");
		}
		// @formatter:on
		gd.hideCancelButton();
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		final boolean disabled = gd.getNextBoolean();
		final boolean anonymize = gd.getNextBoolean();
		
		// Anonymize first to respect the user choice if they still have tracking on
		setAnonymized(anonymize);
		setDisabled(disabled);
	}
}
