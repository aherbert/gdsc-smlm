/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.results.filter;

import com.thoughtworks.xstream.XStream;
import com.thoughtworks.xstream.annotations.XStreamOmitField;
import com.thoughtworks.xstream.io.xml.DomDriver;
import java.util.logging.Level;
import java.util.logging.Logger;
import uk.ac.sussex.gdsc.core.annotation.Nullable;

/**
 * Uses {@link XStream} functionality for reading/writing {@link Filter} members as XML.
 */
public final class FilterXStreamUtils {
  private static final Logger logger = Logger.getLogger(FilterXStreamUtils.class.getName());

  @XStreamOmitField
  private static final XStream xs;

  static {
    xs = new XStream(new DomDriver());
    if (xs != null) {
      try {
        XStream.setupDefaultSecurity(xs); // to be removed after 1.5
        xs.allowTypesByWildcard(new String[] {"uk.ac.sussex.gdsc.smlm.**"});

        xs.autodetectAnnotations(true);

        addAlias(FilterSet.class);

        // Add aliases for all Filter classes
        addAlias(AndFilter.class);
        addAlias(AnrFilter.class);
        addAlias(CombinedFilter.class);
        addAlias(CoordinateFilter.class);
        addAlias(DirectFilter.class);
        addAlias(EShiftFilter.class);
        addAlias(Filter.class);
        addAlias(HysteresisFilter.class);
        addAlias(MultiFilter.class);
        addAlias(MultiFilter2.class);
        addAlias(MultiFilterCrlb.class);
        addAlias(MultiHysteresisFilter.class);
        addAlias(MultiHysteresisFilter2.class);
        addAlias(MultiPathFilter.class);
        addAlias(OrFilter.class);
        addAlias(PrecisionFilter.class);
        addAlias(PrecisionFilter2.class);
        addAlias(PrecisionCrlbFilter.class);
        addAlias(PrecisionHysteresisFilter.class);
        addAlias(PrecisionHysteresisFilter2.class);
        addAlias(SbrFilter.class);
        addAlias(ShiftFilter.class);
        addAlias(SignalFilter.class);
        addAlias(SnrFilter.class);
        addAlias(SnrHysteresisFilter.class);
        addAlias(TraceFilter.class);
        addAlias(WidthFilter.class);
        addAlias(WidthFilter2.class);
        addAlias(XyWidthFilter.class);
        addAlias(XyWidthFilter2.class);
        addAlias(ZCoordinateFilter.class);
      } catch (final Exception ex) {
        logger.log(Level.SEVERE, "Failed to initialise XStream", ex);
      }
    }
  }

  /**
   * No public construction.
   */
  private FilterXStreamUtils() {}

  /**
   * Add a class name alias to the global XStream object used for serialisation.
   *
   * <p>Should be called to produce neater XML output for new sub-class types prior to using
   * {@link #toXml(Object)} or {@link #fromXml(String)}.
   *
   * @param type The class
   */
  public static void addAlias(Class<?> type) {
    if (xs != null) {
      xs.alias(type.getSimpleName(), type);
    }
  }

  /**
   * Create an XML representation of this object.
   *
   * <p>If conversion fails then an empty string is returned.
   *
   * @param object the object
   * @return An XML representation of this object
   */
  public static String toXml(Object object) {
    if (xs != null) {
      try {
        return xs.toXML(object);
      } catch (final Exception ex) {
        logger.log(Level.WARNING, "Failed to serialise to XML", ex);
      }
    }
    return "";
  }

  /**
   * Create an object from the XML representation.
   *
   * @param xml the xml
   * @return the object
   */
  public static @Nullable Object fromXml(String xml) {
    if (xs != null) {
      try {
        return xs.fromXML(xml);
      } catch (final Exception ex) {
        logger.log(Level.WARNING, "Failed to deserialise from XML", ex);
      }
    }
    return null;
  }

  /**
   * Gets the single instance of {@link XStream}.
   *
   * <p>The instance has been initialised with the appropriate aliases for reading and writing
   * {@link Filter} members as XML.
   *
   * @return An XStream object for reading/writing Filter objects
   */
  public static XStream getXStreamInstance() {
    return xs;
  }
}
