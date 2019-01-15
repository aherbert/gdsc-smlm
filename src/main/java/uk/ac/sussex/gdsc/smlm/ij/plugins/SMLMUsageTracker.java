/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

import uk.ac.sussex.gdsc.core.ij.ImageJAnalyticsUtils;
import uk.ac.sussex.gdsc.core.ij.ImageJPluginLoggerHelper;
import uk.ac.sussex.gdsc.smlm.Version;

import ij.plugin.PlugIn;

import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicBoolean;

/**
 * Provide methods to track code usage within ImageJ.
 */
public class SMLMUsageTracker implements PlugIn {
  private static final String TITLE = "SMLM Usage Tracker";

  /** A flag used when the dialog is shown. */
  private static final AtomicBoolean dialogShown = new AtomicBoolean();

  static {
    // This ensures all GDSC loggers redirect from the console to the ImageJ log window.
    // This is here to ensure all GDSC plugins that use this tracking class set-up logging
    // redirection.
    ImageJPluginLoggerHelper.getLogger(SMLMUsageTracker.class.getName());
  }

  /**
   * Initialise on demand the analytics code.
   *
   * <p>This is used to avoid synchronisation during initialisation.
   *
   * <a href="https://en.wikipedia.org/wiki/Initialization-on-demand_holder_idiom">Initialisation on
   * demand</a>
   */
  private static class LazyAnalyticsHolder {
    static {
      // Record the version of the GDSC SMLM plugins
      ImageJAnalyticsUtils.addCustomDimension(7, Version.getVersion());
      // Prompt the user to opt-in/out of analytics if the status is unknown
      if (ImageJAnalyticsUtils.unknownStatus()) {
        showDialog(true);
      }
    }

    /**
     * Checks if is disabled.
     *
     * @return true, if is disabled
     */
    static boolean isDisabled() {
      return ImageJAnalyticsUtils.isDisabled();
    }
  }

  /**
   * Initialise on demand the plugin map.
   *
   * <p>This is used to avoid synchronisation during initialisation.
   *
   * <a href="https://en.wikipedia.org/wiki/Initialization-on-demand_holder_idiom">Initialisation on
   * demand</a>
   */
  private static class LazyMapHolder {
    private static final HashMap<String, String[]> map;

    static {
      final HashMap<String, String[]> localMap = new HashMap<>();
      ImageJAnalyticsUtils.buildPluginMap(localMap, SMLMTools.getPluginsConfig(),
          StandardCharsets.UTF_8);
      map = localMap;
    }
  }

  /**
   * Record the use of the ImageJ plugin using the raw class and argument passed by ImageJ. The
   * plugins config file will be used to identify the correct ImageJ plugin path and title.
   *
   * @param clazz The class
   * @param argument The plugin argument
   */
  public static void recordPlugin(@SuppressWarnings("rawtypes") Class clazz, String argument) {
    if (LazyAnalyticsHolder.isDisabled()) {
      return;
    }

    final String[] pair =
        LazyMapHolder.map.get(ImageJAnalyticsUtils.getKey(clazz.getName(), argument));
    if (pair == null) {
      recordPlugin(clazz.getName().replace('.', '/'), argument);
    } else {
      trackPageView(pair[0], pair[1]);
    }
  }

  /**
   * Record the use of the ImageJ plugin.
   *
   * @param name The plugin name
   * @param argument The plugin argument
   */
  private static void recordPlugin(String name, String argument) {
    final StringBuilder url = new StringBuilder(name.length() + 16);
    // Assume plugin name has no '/' prefix
    url.append('/').append(name);
    if (argument != null && argument.length() > 0) {
      url.append("?arg=").append(argument);
    }
    trackPageView(url.toString(), name);
  }

  private static void trackPageView(String pageUrl, String pageTitle) {
    ImageJAnalyticsUtils.pageview(pageUrl, pageTitle);
  }

  @Override
  public void run(String arg) {
    // If this is the first plugin to call recordPlugin(...) then the dialog may be shown.
    dialogShown.set(false);
    recordPlugin(this.getClass(), arg);
    if (!dialogShown.get()) {
      showDialog(false);
    }
  }

  /**
   * Show a dialog allowing users to opt in/out of Google Analytics.
   *
   * @param autoMessage Set to true to display the message about automatically showing when the
   *        status is unknown
   */
  static void showDialog(boolean autoMessage) {
    dialogShown.set(true);
    ImageJAnalyticsUtils.showDialog(TITLE, autoMessage);
  }
}
