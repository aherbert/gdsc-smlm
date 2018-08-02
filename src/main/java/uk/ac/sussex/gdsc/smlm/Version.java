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
package uk.ac.sussex.gdsc.smlm;

import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Show the version information contained in the /uk/ac/sussex/gdsc/smlm/Version.txt file.
 * <p>
 * Uses Semantic Versioning.
 *
 * @see <a href="http://semver.org/">http://semver.org/</a>
 */
public class Version
{
    /** Constant for the string "unknown" */
    public static final String UNKNOWN = "unknown";
    private static String version = null;
    private static String buildDate = null;

    static
    {
        // Locate the version file
        final Class<Version> resourceClass = Version.class;
        try (final InputStream propertiesStream = resourceClass
                .getResourceAsStream("/uk/ac/sussex/gdsc/smlm/Version.txt"))
        {
            // Read the version properties
            final Properties props = new Properties();
            props.load(propertiesStream);
            version = props.getProperty("version");
            buildDate = props.getProperty("build.date");
        }
        catch (final IOException e)
        {
            // Ignore
        }

        if (version == null || version.length() == 0)
            version = UNKNOWN;
        if (buildDate == null || buildDate.length() == 0)
            buildDate = UNKNOWN;
    }

    /**
     * The main method. Output the version and build date.
     *
     * @param args
     *            the arguments
     */
    public static void main(String[] args)
    {
        final StringBuilder msg = new StringBuilder();
        final String newLine = System.getProperty("line.separator");
        msg.append("Version : ").append(version).append(newLine);
        msg.append("Build Date : ").append(buildDate).append(newLine);
        System.out.print(msg.toString());
    }

    /**
     * Get the GDSC SMLM version
     *
     * @return The uk.ac.sussex.gdsc.smlm package version
     */
    public static String getVersion()
    {
        return version;
    }

    /**
     * Get the GDSC SMLM package build date
     *
     * @return The uk.ac.sussex.gdsc.smlm package build date
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
        final Pattern p = Pattern.compile("^\\d+");
        final Matcher m = p.matcher(version);
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
        final Pattern p = Pattern.compile("^\\d+\\.(\\d+)");
        final Matcher m = p.matcher(version);
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
        final Pattern p = Pattern.compile("^\\d+\\.\\d+\\.(\\d+)");
        final Matcher m = p.matcher(version);
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
        final Pattern p = Pattern.compile("^\\d+\\.\\d+\\.\\d+");
        final Matcher m = p.matcher(version);
        if (m.find())
            return m.group();
        return "";
    }
}
