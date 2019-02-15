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

package uk.ac.sussex.gdsc.smlm.utils;

import com.thoughtworks.xstream.XStream;
import com.thoughtworks.xstream.XStreamException;
import com.thoughtworks.xstream.io.xml.DomDriver;

import java.util.regex.Pattern;

/**
 * Provide XML utilities using XStream.
 */
public class XStreamUtils {
  private static XStream xs;

  private static class PatternLoader {
    static final Pattern PACKAGE_PATTERN = Pattern.compile("(</?)gdsc.smlm");
  }

  /**
   * Convert an object to XML.
   *
   * @param obj the object
   * @return XML string representation
   */
  public static String toXml(Object obj) {
    init();
    if (xs != null) {
      try {
        return xs.toXML(obj);
      } catch (final XStreamException ex) {
        // ex.printStackTrace();
      }
    }
    return "";
  }

  /**
   * Load an object from the XML string representation.
   *
   * @param xml the xml
   * @return the object
   */
  public static Object fromXml(String xml) {
    init();
    if (xs != null) {
      try {
        return xs.fromXML(xml);
      } catch (final XStreamException ex) {
        // ex.printStackTrace();
      }
    }
    return null;
  }

  private static void init() {
    if (xs == null) {
      xs = new XStream(new DomDriver());
      XStream.setupDefaultSecurity(xs); // to be removed after 1.5
      xs.allowTypesByWildcard(new String[] {"uk.ac.sussex.gdsc.smlm.**"});
    }
  }

  /**
   * Update any XML elements using the old {@code <gdsc.smlm.*>} package name to the new
   * {@code <uk.ac.sussex.gdsc.smlm.*>} package name.
   *
   * @param xml the xml
   * @return the updated xml
   */
  public static String updateGdscPackageName(String xml) {
    // Fix for reading old versions:
    // Support package gdsc.smlm renamed to uk.ac.sussex.gdsc.smlm
    if (xml.contains("<gdsc.smlm")) {
      return PatternLoader.PACKAGE_PATTERN.matcher(xml).replaceAll("$1uk.ac.sussex.gdsc.smlm");
    }
    return xml;
  }
}
