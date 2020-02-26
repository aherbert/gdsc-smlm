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

/**
 * Specifies a peak fitting result for use in filtering. Any result implementing this interface can
 * be directly filtered without requiring the filter to be initialised with calibration data. This
 * result can be assigned matches to actual data.
 */
public interface AssignablePreprocessedPeakResult extends PreprocessedPeakResult {
  /**
   * Set the assignments between this result and the true data.
   *
   * @param assignments The assignments
   */
  void setAssignments(ResultAssignment[] assignments);

  /**
   * Sets the ignore flag. If true then the result should be ignored from the total counts when
   * scoring.
   *
   * @param ignore the new ignore flag
   */
  void setIgnore(boolean ignore);
}
