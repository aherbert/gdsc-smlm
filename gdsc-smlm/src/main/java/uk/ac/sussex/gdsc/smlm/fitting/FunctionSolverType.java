/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.fitting;

/**
 * Define the type of function solver.
 */
public enum FunctionSolverType {
  /** Least Squares Estimator. */
  LSE("Least Squares Estimator"),
  /** Weighted Least Squares Estimator. */
  WLSE("Weighted Least Squares Estimator"),
  /** Maximum Likelihood Estimator. */
  MLE("Maximum Likelihood Estimator");

  private String niceName;

  FunctionSolverType(String name) {
    niceName = name;
  }

  @Override
  public String toString() {
    return getName();
  }

  /**
   * Gets the name.
   *
   * @return the name
   */
  public String getName() {
    return niceName;
  }
}
