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

package uk.ac.sussex.gdsc.smlm.results.procedures;

/**
 * Interface for accessing the localisation precision of Gaussian 2D fitting computed using the
 * Mortensen formula for Maximum Likelihood Estimation using local background photons to estimate
 * noise.
 *
 * <p>See Mortensen, et al (2010) Nature Methods 7, 377-383, equation 6.
 */
public interface MlePrecisionBProcedure {
  /**
   * Executes this procedure.
   *
   * @param precision the precision
   */
  void executeMlePrecisionB(double precision);
}
