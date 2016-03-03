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
package gdsc.smlm.ij.plugins;

import java.util.HashMap;

import gdsc.core.ij.ImageJAnalyticsTracker;
import gdsc.smlm.Version;
import ij.plugin.PlugIn;

/**
 * Provide methods to track code usage within ImageJ
 */
public class SMLMUsageTracker implements PlugIn
{
	private static final String TITLE = "SMLM Usage Tracker";
	private static HashMap<String, String[]> map = new HashMap<String, String[]>();
	private static boolean trackerInitialised = false;
	private static boolean mapInitialised = false;

	/**
	 * Record the use of the ImageJ plugin using the raw class and argument passed by ImageJ. The plugins config
	 * file will be used to identify the correct ImageJ plugin path and title.
	 * 
	 * @param clazz
	 *            The class
	 * @param argument
	 *            The plugin argument
	 */
	public static void recordPlugin(@SuppressWarnings("rawtypes") Class clazz, String argument)
	{
		initialiseTracker();
		if (ImageJAnalyticsTracker.isDisabled())
			return;

		initialiseMap();

		final String[] pair = map.get(ImageJAnalyticsTracker.getKey(clazz.getName(), argument));
		if (pair == null)
		{
			recordPlugin(clazz.getName().replace('.', '/'), argument);
		}
		else
		{
			trackPageView(pair[0], pair[1]);
		}
	}

	/**
	 * Record the use of the ImageJ plugin
	 * 
	 * @param name
	 *            The plugin name
	 * @param argument
	 *            The plugin argument
	 */
	private static void recordPlugin(String name, String argument)
	{
		String url = '/' + name;
		if (argument != null && argument.length() > 0)
			url += "?arg=" + argument;

		trackPageView(url, name);
	}

	private static void trackPageView(String pageUrl, String pageTitle)
	{
		ImageJAnalyticsTracker.pageview(pageUrl, pageTitle);
	}

	/**
	 * Create the tracker and then verify the opt in/out status of the user if it unknown or a new major.minor version.
	 */
	private static void initialiseTracker()
	{
		if (trackerInitialised)
			return;

		synchronized (SMLMUsageTracker.class)
		{
			// Check again since this may be a second thread that was waiting  
			if (trackerInitialised)
				return;

			trackerInitialised = true;
			
			// Initialise analytics
			ImageJAnalyticsTracker.initialise();
			// Record the version of the GDSC SMLM plugins
			ImageJAnalyticsTracker.addCustomDimension(7, Version.getVersion());
			// Prompt the user to opt-in/out of analytics if the status is unknown
			if (ImageJAnalyticsTracker.unknownStatus())
				ImageJAnalyticsTracker.showDialog(TITLE, true);
		}
	}

	/**
	 * Create the map between the ImageJ plugin class and argument and the ImageJ menu path and plugin title
	 */
	private static void initialiseMap()
	{
		if (mapInitialised)
			return;

		synchronized (map)
		{
			// Check again since this may be a second thread that was waiting  
			if (mapInitialised)
				return;

			mapInitialised = true;
			ImageJAnalyticsTracker.buildPluginMap(map, SMLMTools.getPluginsConfig());
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		recordPlugin(this.getClass(), arg);
		ImageJAnalyticsTracker.showDialog(TITLE, false);
	}
}
