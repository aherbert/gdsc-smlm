/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.math3.optim;

import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimplePointChecker;

/**
 * Check if the position has converged.
 */
public class PositionChecker extends SimplePointChecker<PointValuePair>
    implements OptimizationData {
  /**
   * Build an instance with specified thresholds.
   *
   * <p>In order to perform only relative checks, the absolute tolerance must be set to a negative
   * value. In order to perform only absolute checks, the relative tolerance must be set to a
   * negative value.
   *
   * @param relative relative tolerance threshold
   * @param absolute absolute tolerance threshold
   */
  public PositionChecker(double relative, double absolute) {
    super(relative, absolute);
  }
}
