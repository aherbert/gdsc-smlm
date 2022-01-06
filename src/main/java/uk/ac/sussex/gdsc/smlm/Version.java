/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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
import java.net.URL;
import java.util.jar.Attributes;
import java.util.jar.Manifest;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Show the version information contained in the source jar manifest.
 */
public class Version {
  /** Constant for the string "unknown". */
  public static final String UNKNOWN = "";
  private static String version;
  private static String buildDate;
  private static String buildNumber;

  static {
    final Manifest manifest = loadManifest(Version.class);
    if (manifest != null) {
      final Attributes attributes = manifest.getMainAttributes();
      version = attributes.getValue("Specification-Version");
      buildDate = attributes.getValue("Implementation-Date");
      buildNumber = attributes.getValue("Implementation-Build");
    }

    if (version == null) {
      version = UNKNOWN;
    }
    if (buildDate == null) {
      buildDate = UNKNOWN;
    }
    if (buildNumber == null) {
      buildNumber = UNKNOWN;
    }
  }

  /**
   * The main method. Output the version and build date.
   *
   * @param args the arguments
   */
  public static void main(String[] args) {
    final StringBuilder msg = new StringBuilder();
    final String newLine = System.lineSeparator();
    msg.append("Version : ").append(version).append(newLine);
    msg.append("Build Date : ").append(buildDate).append(newLine);
    msg.append("Build Number : ").append(buildNumber).append(newLine);
    System.out.print(msg);
  }

  /**
   * Get the GDSC SMLM version.
   *
   * @return The uk.ac.sussex.gdsc.smlm package version
   */
  public static String getVersion() {
    return version;
  }

  /**
   * Get the GDSC SMLM package build date.
   *
   * @return The uk.ac.sussex.gdsc.smlm package build date
   */
  public static String getBuildDate() {
    return buildDate;
  }

  /**
   * Get the GDSC SMLM package build number.
   *
   * @return The uk.ac.sussex.gdsc.smlm package build number
   */
  public static String getBuildNumber() {
    return buildNumber;
  }

  /**
   * Get the major version.
   *
   * @return The major version (or 0 if unknown)
   */
  public static int getMajorVersion() {
    final Pattern p = Pattern.compile("^\\d+");
    final Matcher m = p.matcher(version);
    if (m.find()) {
      return Integer.parseInt(m.group());
    }
    return 0;
  }

  /**
   * Get the minor version.
   *
   * @return The minor version (or 0 if unknown)
   */
  public static int getMinorVersion() {
    final Pattern p = Pattern.compile("^\\d+\\.(\\d+)");
    final Matcher m = p.matcher(version);
    if (m.find()) {
      return Integer.parseInt(m.group(1));
    }
    return 0;
  }

  /**
   * Get the patch version.
   *
   * @return The patch version (or 0 if unknown)
   */
  public static int getPatchVersion() {
    final Pattern p = Pattern.compile("^\\d+\\.\\d+\\.(\\d+)");
    final Matcher m = p.matcher(version);
    if (m.find()) {
      return Integer.parseInt(m.group(1));
    }
    return 0;
  }

  /**
   * Get a string with the major, minor and patch versions.
   *
   * @return Major.Minor.Patch
   */
  public static String getMajorMinorPatch() {
    final Pattern p = Pattern.compile("^\\d+\\.\\d+\\.\\d+");
    final Matcher m = p.matcher(version);
    if (m.find()) {
      return m.group();
    }
    return "";
  }

  /**
   * Load the jar manifest for the given class.
   *
   * <p>If not from a jar or an IO exception occurs return null.
   *
   * @param clazz the class
   * @return the manifest (or null)
   */
  public static Manifest loadManifest(Class<?> clazz) {
    final String resource = "/" + clazz.getName().replace(".", "/") + ".class";
    final String classPath = clazz.getResource(resource).toString();
    if (!classPath.startsWith("jar")) {
      // Class not from JAR
      return null;
    }
    final String manifestPath =
        classPath.substring(0, classPath.lastIndexOf('!') + 1) + "/META-INF/MANIFEST.MF";
    try {
      try (InputStream in = new URL(manifestPath).openStream()) {
        return new Manifest(in);
      }
    } catch (final IOException ex) {
      // Ignore this
    }
    return null;
  }
}
