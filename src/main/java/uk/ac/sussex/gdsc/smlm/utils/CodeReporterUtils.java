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

package uk.ac.sussex.gdsc.smlm.utils;

import java.io.PrintStream;

/**
 * Allow the position of executing code to be reported.
 */
public final class CodeReporterUtils {
  /**
   * No public constructor.
   */
  private CodeReporterUtils() {}

  /**
   * Print the class name, method name and line number from the first element of the throwable's
   * stack trace to the output.
   *
   * @param throwable the throwable
   * @param out the out
   * @see Throwable#getStackTrace()
   * @see StackTraceElement#getClassName()
   * @see StackTraceElement#getMethodName()
   * @see StackTraceElement#getLineNumber()
   */
  public static void report(Throwable throwable, PrintStream out) {
    final StackTraceElement e = throwable.getStackTrace()[0];
    out.printf("%s:%s:%d%n", e.getClassName(), e.getMethodName(), e.getLineNumber());
  }
}
