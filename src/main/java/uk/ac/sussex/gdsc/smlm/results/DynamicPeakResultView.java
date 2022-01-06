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

package uk.ac.sussex.gdsc.smlm.results;

import uk.ac.sussex.gdsc.smlm.results.predicates.FramePeakResultPredicate;
import uk.ac.sussex.gdsc.smlm.results.predicates.IdPeakResultPredicate;

/**
 * Provides a dynamic view of the results. Changes to the underlying results store will be reflected
 * in the view.
 */
public class DynamicPeakResultView implements PeakResultView {
  private final PeakResultStore store;

  /**
   * Instantiates a new cached peak result view.
   *
   * @param store the store
   */
  public DynamicPeakResultView(PeakResultStore store) {
    this.store = store;
  }

  /**
   * Gets the results by frame.
   *
   * @param frame the frame
   * @return the results
   */
  @Override
  public PeakResult[] getResultsByFrame(int frame) {
    return store.subset(new FramePeakResultPredicate(frame));
  }

  /**
   * Gets the results by id.
   *
   * @param id the id
   * @return the results
   */
  @Override
  public PeakResult[] getResultsById(int id) {
    return store.subset(new IdPeakResultPredicate(id));
  }
}
