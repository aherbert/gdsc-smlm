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

/**
 * Contains help URLs for the GDSC ImageJ plugins.
 */
public final class HelpUrls {
  /** The system property to define the base help url. */
  private static final String HELP_URL_BASE_PROPERTY = "gdsc.smlm.help.url.base";
  /** The system property to define the help url properties configuration file. */
  private static final String HELP_URL_CONF_PROPERTY = "gdsc.smlm.help.url.conf";
  /** The help url for the SMLM plugins. */
  private static final String DEFAULT_HELP_URL = "https://gdsc-smlm.readthedocs.io/en/latest/";
  /** The url configuration for the SMLM plugins. */
  private static final String DEFAULT_HELP_CONF = "/uk/ac/sussex/gdsc/smlm/help.config";

  /** The url base. */
  private static String urlBase;

  /** The map storing sub-URLs to append to the base URL. */
  private static Properties map;

  static {
    urlBase = System.getProperty(HELP_URL_BASE_PROPERTY, DEFAULT_HELP_URL);
    map = loadUrls(System.getProperty(HELP_URL_CONF_PROPERTY, DEFAULT_HELP_CONF));
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
    return urlBase;
  }

  /**
   * Gets the URL for the given key.
   *
   * @param key the key
   * @return the URL
   */
  public static String getUrl(String key) {
    return getUrl(key, map);
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
    return path == null ? urlBase : urlBase + path;
  }

  /**
   * Enumerate all the (key, value) URLs.
   *
   * @param action the action
   */
  @VisibleForTesting
  static void forEach(BiConsumer<String, String> action) {
    map.forEach((k, v) -> action.accept((String) k, urlBase + v));
  }
}
