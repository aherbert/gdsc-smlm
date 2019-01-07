/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

/**
 * Store the score from analysis of the non-filter parameters during direct filter analysis.
 */
public class ParameterScoreResult {
  /** The score. */
  final double score;
  /** The criteria. */
  final double criteria;
  /** The parameters. */
  final double[] parameters;
  /** The text. */
  final String text;

  /**
   * Instantiates a new parameter score result.
   *
   * @param score the score
   * @param criteria the criteria
   * @param parameters the parameters
   * @param text the text
   */
  public ParameterScoreResult(double score, double criteria, double[] parameters, String text) {
    this.score = score;
    this.criteria = criteria;
    this.parameters = parameters;
    this.text = text;
  }
}
