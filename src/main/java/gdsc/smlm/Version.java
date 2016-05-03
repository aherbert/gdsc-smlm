package gdsc.smlm;

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

import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Show the version information contained in the gdsc/smlm/Version.txt file.
 * <p>
 * Uses Semantic Versioning.
 * 
 * @see http://semver.org/
 */
public class Version
{
	public static final String UNKNOWN = "unknown";
	private static String version = null;
	private static String buildDate = null;

	static
	{
		// Locate the version file
		Class<Version> resourceClass = Version.class;
		InputStream propertiesStream = resourceClass.getResourceAsStream("/gdsc/smlm/Version.txt");

		try
		{
			// Read the version properties
			Properties props = new Properties();
			props.load(propertiesStream);
			version = props.getProperty("version");
			buildDate = props.getProperty("build.date");
		}
		catch (IOException e)
		{
			// Ignore
		}

		if (version == null || version.length() == 0)
			version = UNKNOWN;
		if (buildDate == null || buildDate.length() == 0)
			buildDate = UNKNOWN;
	}

	public static void main(String[] args)
	{
		StringBuilder msg = new StringBuilder();
		String newLine = System.getProperty("line.separator");
		msg.append("Version : ").append(version).append(newLine);
		msg.append("Build Date : ").append(buildDate).append(newLine);
		System.out.print(msg.toString());
	}

	/**
	 * Get the GDSC SMLM version
	 * 
	 * @return The gdsc.smlm package version
	 */
	public static String getVersion()
	{
		return version;
	}

	/**
	 * Get the GDSC SMLM package build date
	 * 
	 * @return The gdsc.smlm package build date
	 */
	public static String getBuildDate()
	{
		return buildDate;
	}

	/**
	 * Get the major version
	 * 
	 * @return The major version (or 0 if unknown)
	 */
	public static int getMajorVersion()
	{
		Pattern p = Pattern.compile("^\\d+");
		Matcher m = p.matcher(version);
		if (m.find())
			return Integer.parseInt(m.group());
		return 0;
	}

	/**
	 * Get the minor version
	 * 
	 * @return The minor version (or 0 if unknown)
	 */
	public static int getMinorVersion()
	{
		Pattern p = Pattern.compile("^\\d+\\.(\\d+)");
		Matcher m = p.matcher(version);
		if (m.find())
			return Integer.parseInt(m.group(1));
		return 0;
	}

	/**
	 * Get the patch version
	 * 
	 * @return The patch version (or 0 if unknown)
	 */
	public static int getPatchVersion()
	{
		Pattern p = Pattern.compile("^\\d+\\.\\d+\\.(\\d+)");
		Matcher m = p.matcher(version);
		if (m.find())
			return Integer.parseInt(m.group(1));
		return 0;
	}

	/**
	 * Get a string with the major, minor and patch versions
	 * 
	 * @return Major.Minor.Patch
	 */
	public static String getMajorMinorPatch()
	{
		Pattern p = Pattern.compile("^\\d+\\.\\d+\\.\\d+");
		Matcher m = p.matcher(version);
		if (m.find())
			return m.group();
		return "";
	}
}
