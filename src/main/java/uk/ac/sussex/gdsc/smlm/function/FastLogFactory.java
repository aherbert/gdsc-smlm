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
package uk.ac.sussex.gdsc.smlm.function;

/**
 * Factory class for creating a fast log instance.
 */
public class FastLogFactory {
  private static FastLog fastLog = null;

  /**
   * Gets the global fast log instance.
   *
   * @return the fast log instance
   */
  public static FastLog getFastLog() {
    if (fastLog == null) {
      fastLog = new TurboLog();
    }
    return fastLog;
  }

  /**
   * Gets an instance of FastLog that uses Math.log, i.e. this is not fast.
   *
   * @return the log instance
   */
  public static FastLog getLog() {
    return NonFastLog.INSTANCE;
  }
}
