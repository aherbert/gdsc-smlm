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

package uk.ac.sussex.gdsc.smlm.ga;

/**
 * Exception to throw if the population size is invalid ({@code <2} individuals).
 */
public class InvalidPopulationSize extends RuntimeException {
  private static final long serialVersionUID = -7425196200200953611L;

  /**
   * Instantiates a new invalid population size.
   *
   * @param size the size
   * @param minSize the min size
   */
  public InvalidPopulationSize(int size, int minSize) {
    super("Population size (" + size + ") must be at least " + minSize);
  }
}
