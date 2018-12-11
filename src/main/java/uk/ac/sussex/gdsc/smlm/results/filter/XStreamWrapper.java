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

package uk.ac.sussex.gdsc.smlm.results.filter;

import com.thoughtworks.xstream.XStream;
import com.thoughtworks.xstream.XStreamException;
import com.thoughtworks.xstream.annotations.XStreamOmitField;
import com.thoughtworks.xstream.io.xml.DomDriver;

/**
 * Wraps the XStream functionality for reading/writing package members as XML. Initialises XStream
 * for tidy XML.
 */
public abstract class XStreamWrapper {
  @XStreamOmitField
  private static XStream xs = null;

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
        addAlias(ANRFilter.class);
        addAlias(CombinedFilter.class);
        addAlias(CoordinateFilter.class);
        addAlias(DirectFilter.class);
        addAlias(EShiftFilter.class);
        addAlias(Filter.class);
        addAlias(HysteresisFilter.class);
        addAlias(MultiFilter.class);
        addAlias(MultiFilter2.class);
        addAlias(MultiFilterCRLB.class);
        addAlias(MultiHysteresisFilter.class);
        addAlias(MultiHysteresisFilter2.class);
        addAlias(MultiPathFilter.class);
        addAlias(OrFilter.class);
        addAlias(PrecisionFilter.class);
        addAlias(PrecisionFilter2.class);
        addAlias(PrecisionCRLBFilter.class);
        addAlias(PrecisionHysteresisFilter.class);
        addAlias(PrecisionHysteresisFilter2.class);
        addAlias(SBRFilter.class);
        addAlias(ShiftFilter.class);
        addAlias(SignalFilter.class);
        addAlias(SNRFilter.class);
        addAlias(SNRHysteresisFilter.class);
        addAlias(TraceFilter.class);
        addAlias(WidthFilter.class);
        addAlias(WidthFilter2.class);
        addAlias(XYWidthFilter.class);
        addAlias(XYWidthFilter2.class);
        addAlias(ZCoordinateFilter.class);

        // Removed dependency on reflections since this has other jar dependencies
        // Reflections reflections = new Reflections("uk.ac.sussex.gdsc.smlm.results.filter");
        // Set<Class<? extends DirectFilter>> subTypes = reflections.getSubTypesOf(Filter.class);
        // for (Class<? extends DirectFilter> type : subTypes)
        // addAlias(type);
      } catch (final XStreamException ex) {
        ex.printStackTrace();
      } catch (final Exception ex) {
        ex.printStackTrace();
      }
    }
  }

  /**
   * Add a class name alias to the global XStream object used for serialisation.
   *
   * <p>Should be called to produce neater XML output for new sub-class types prior to using
   * {@link #toXML(Object)} or {@link #fromXML(String)}.
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
   * @param object the object
   * @return An XML representation of this object
   */
  public static String toXML(Object object) {
    if (xs != null) {
      try {
        return xs.toXML(object);
      } catch (final XStreamException ex) {
        ex.printStackTrace();
      } catch (final Exception ex) {
        ex.printStackTrace();
      }
    }
    return "";
  }

  /**
   * Create the filter from the XML representation.
   *
   * @param xml the xml
   * @return the filter
   */
  public static Object fromXML(String xml) {
    if (xs != null) {
      try {
        return xs.fromXML(xml);
      } catch (final XStreamException ex) {
        ex.printStackTrace();
      } catch (final Exception ex) {
        ex.printStackTrace();
      }
    }
    return null;
  }

  /**
   * @return An XStream object for reading/writing package members
   */
  public static XStream getInstance() {
    return xs;
  }
}
