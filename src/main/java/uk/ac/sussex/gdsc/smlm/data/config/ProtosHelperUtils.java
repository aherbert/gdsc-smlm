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

package uk.ac.sussex.gdsc.smlm.data.config;

/**
 * Contains helper functions for the ResultsProtos class.
 */
final class ProtosHelperUtils {

  private static final String UNKNOWN_NAME = "Unknown name: ";
  private static final String UNKNOWN_METHOD = "Unknown method: ";

  /** The constant "Unknown". */
  public static final String UNKNOWN = "Unknown";

  /** No public constructor. */
  private ProtosHelperUtils() {}

  /**
   * Create a message when the name of the object cannot be determined.
   *
   * @param object the object
   * @return the message
   */
  public static String unknownNameMessage(Object object) {
    return UNKNOWN_NAME + object;
  }

  /**
   * Create a message when the method of the object cannot be determined.
   *
   * @param object the object
   * @return the message
   */
  public static String unknownMethodMessage(Object object) {
    return UNKNOWN_METHOD + object;
  }
}
