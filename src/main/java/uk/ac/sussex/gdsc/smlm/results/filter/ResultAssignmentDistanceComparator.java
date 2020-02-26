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

import java.io.Serializable;
import java.util.Comparator;

/**
 * Compares the {@link ResultAssignment} using the distance, lowest first.
 */
public class ResultAssignmentDistanceComparator
    implements Comparator<ResultAssignment>, Serializable {
  private static final long serialVersionUID = 1L;

  /** The instance. */
  public static final ResultAssignmentDistanceComparator INSTANCE =
      new ResultAssignmentDistanceComparator();

  @Override
  public int compare(ResultAssignment o1, ResultAssignment o2) {
    return Double.compare(o1.distance, o2.distance);
  }
}
