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

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear;

import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;

/**
 * Exception to throw if a function solver failed.
 */
public class FunctionSolverException extends RuntimeException {
  /** The Constant serialVersionUID. */
  private static final long serialVersionUID = -5131234527135746186L;

  /** The fit status. This indicates why the solver failed. */
  public final FitStatus fitStatus;

  /**
   * Instantiates a new function solver exception.
   *
   * @param fitStatus the fit status
   */
  public FunctionSolverException(FitStatus fitStatus) {
    super();
    this.fitStatus = fitStatus;
  }

  /**
   * Instantiates a new function solver exception.
   *
   * @param fitStatus the fit status
   * @param message the message
   */
  public FunctionSolverException(FitStatus fitStatus, String message) {
    super(message);
    this.fitStatus = fitStatus;
  }

  /**
   * Instantiates a new function solver exception.
   *
   * @param fitStatus the fit status
   * @param message the message
   * @param cause the cause
   */
  public FunctionSolverException(FitStatus fitStatus, String message, Throwable cause) {
    super(message, cause);
    this.fitStatus = fitStatus;
  }

  /**
   * Instantiates a new function solver exception.
   *
   * @param fitStatus the fit status
   * @param cause the cause
   */
  public FunctionSolverException(FitStatus fitStatus, Throwable cause) {
    super(cause);
    this.fitStatus = fitStatus;
  }
}
