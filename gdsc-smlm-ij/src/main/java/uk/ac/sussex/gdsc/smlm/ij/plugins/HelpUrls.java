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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;
import java.util.function.BiConsumer;
import java.util.logging.Level;
import org.apache.commons.lang3.StringUtils;
import uk.ac.sussex.gdsc.core.data.VisibleForTesting;
import uk.ac.sussex.gdsc.core.ij.ImageJPluginLoggerHelper;
import uk.ac.sussex.gdsc.smlm.Version;

/**
 * Contains help URLs for the GDSC ImageJ plugins.
 */
public final class HelpUrls {
  /** The system property to define the base help url. */
  private static final String HELP_URL_BASE_PROPERTY = "gdsc.smlm.help.url.base";
  /** The system property to define the verion of the help url. */
  private static final String HELP_URL_VERSION_PROPERTY = "gdsc.smlm.help.url.version";
  /** The system property to define the help url properties configuration file. */
  private static final String HELP_URL_CONF_PROPERTY = "gdsc.smlm.help.url.conf";
  /** The help url for the SMLM plugins. */
  private static final String DEFAULT_HELP_URL = "https://gdsc-smlm.readthedocs.io/en/";
  /** The version for the SMLM plugins. */
  private static final String DEFAULT_HELP_VERSION = "latest";
  /** The url configuration for the SMLM plugins. */
  private static final String DEFAULT_HELP_CONF = "/uk/ac/sussex/gdsc/smlm/ij/help.config";

  /** The url base. */
  private static final String URL_BASE;

  /** The map storing sub-URLs to append to the base URL. */
  private static final Properties URL_PAGE_MAP;

  static {
    // url = base + version

    // Runtime specified.
    // No checks are made that the input URL is valid
    final String base = System.getProperty(HELP_URL_BASE_PROPERTY, DEFAULT_HELP_URL);
    String version = System.getProperty(HELP_URL_VERSION_PROPERTY);

    if (version == null) {
      version = getTagVersionOrDefault(Version.getVersion());
    }

    URL_BASE = createUrl(base, version);
    URL_PAGE_MAP = loadUrls(System.getProperty(HELP_URL_CONF_PROPERTY, DEFAULT_HELP_CONF));
  }

  /** No construction. */
  private HelpUrls() {}

  /**
   * Load the URLs from the given resource.
   *
   * @param resource the resource
   * @return the properties
   */
  @VisibleForTesting
  static Properties loadUrls(String resource) {
    final Properties props = new Properties();
    // Ignore blank input. This will load a file listing of the root resource class path.
    if (StringUtils.isNotBlank(resource)) {
      try (InputStream in = HelpUrls.class.getResourceAsStream(resource)) {
        if (in != null) {
          props.load(in);
        }
      } catch (IOException ex) {
        ImageJPluginLoggerHelper.getDefaultLogger().log(Level.WARNING, "Failed to load Help URLs",
            ex);
      }
    }
    return props;
  }

  /**
   * Gets the base help URL.
   *
   * @return the URL
   */
  public static String getUrl() {
    return URL_BASE;
  }

  /**
   * Gets the URL for the given key.
   *
   * @param key the key
   * @return the URL
   */
  public static String getUrl(String key) {
    return getUrl(key, URL_PAGE_MAP);
  }

  /**
   * Gets the URL for the given key. Returns the base url with the suffix specified by the key
   * appended. If there is no entry for the key then only the base url is returned.
   *
   * @param key the key
   * @param map the map containing the URL suffix
   * @return the URL
   */
  @VisibleForTesting
  static String getUrl(String key, Properties map) {
    final String path = map.getProperty(key);
    return path == null ? URL_BASE : URL_BASE + path;
  }

  /**
   * Enumerate all the (key, value) URLs.
   *
   * @param action the action
   */
  @VisibleForTesting
  static void forEach(BiConsumer<String, String> action) {
    URL_PAGE_MAP.forEach((k, v) -> action.accept((String) k, URL_BASE + v));
  }

  /**
   * Gets the help url version for a release based on knowledge of the release tag pattern, or else
   * the default help version ("latest").
   *
   * <p>The method assumes semver format using numbers for major.minor[.patch]. The tag for a
   * release will have a 'v' prefix.
   *
   * <pre>
   * 1.0-SNAPSHOT    latest
   * 1.0             v1.0
   * 1.0.0           v1.0.0
   * 1.0.1           v1.0.1
   * 1.1-SNAPSHOT    latest
   * 1.1             v1.1
   * </pre>
   *
   * <p>Note: The patch number is optional if zero. No logic is used to detect a zero patch number.
   * This respects the choice chosen for the semver.
   *
   * @param version the version
   * @return the tag version or default
   */
  @VisibleForTesting
  static String getTagVersionOrDefault(String version) {
    // Tags have a 'v' prefixed to the semver format of the version:
    // vX.Y[.Z]
    // Z is added if non-zero.
    // The following just assumes a tag for a release version will be the same with
    // a v prefix. The regex looks for 'X.Y' and optionally '.Z'.
    // No logic for a patch version of zero is made.
    // If 1.0.0 is the release version then the tag should be 'v1.0.0' not 'v1.0'
    if (version.matches("^(\\d+\\.\\d+)(\\.\\d+)?$")) {
      // Add the 'v' prefix
      return "v" + version;
    }
    ImageJPluginLoggerHelper.getLogger(HelpUrls.class)
        .fine(() -> String.format("GDSC SMLM help url version is not a release: <%s>, using: <%s>",
            version, DEFAULT_HELP_VERSION));
    // Fall back to the default
    return DEFAULT_HELP_VERSION;
  }

  /**
   * Create the URL by joining the base to the remaining path, adding the separator '/' if required.
   * A '/' is added at the end of the url.
   *
   * @param base the base
   * @param extra the extra
   * @return the url string
   */
  @VisibleForTesting
  static String createUrl(String base, String... extra) {
    final StringBuilder sb = new StringBuilder(base.length());
    append(sb, base);
    for (final String s : extra) {
      append(sb, s);
    }
    return sb.toString();
  }

  /**
   * Append the string if non-empty. Appends the separator '/' if missing from the string.
   *
   * @param sb the StringBuilder
   * @param s the string
   */
  private static void append(StringBuilder sb, String s) {
    if (!s.isEmpty()) {
      sb.append(s);
      if (sb.charAt(sb.length() - 1) != '/') {
        sb.append('/');
      }
    }
  }
}
