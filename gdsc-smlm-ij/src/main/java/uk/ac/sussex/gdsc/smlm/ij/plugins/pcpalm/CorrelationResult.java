/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Package
 *
 * Software for single molecule localisation microscopy (SMLM) in ImageJ
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

package uk.ac.sussex.gdsc.smlm.ij.plugins.pcpalm;

/**
 * Used to store the correlation (g(r)) result for the PC-PALM analysis.
 */
class CorrelationResult {
  /** The id. */
  int id;

  /** The source for the result. */
  String source;

  /** The minx. */
  double minx;

  /** The miny. */
  double miny;

  /** The maxx. */
  double maxx;

  /** The maxy. */
  double maxy;

  /** The nm per pixel. */
  double nmPerPixel;

  /** The peak density. */
  double peakDensity;

  /** The number of unique points. */
  double uniquePoints;

  /**
   * Set to true if pixels in the image have a 1/0 value. Otherwise it is assumed the image has
   * continuous real values, i.e. pixels can have 1 or more localisations.
   */
  boolean binaryImage;

  /** The correlation curve. */
  double[][] gr;

  /** The spatial domain. */
  boolean spatialDomain;

  /**
   * Instantiates a new correlation result.
   *
   * @param id the id
   * @param source the source
   * @param minx the minx
   * @param miny the miny
   * @param maxx the maxx
   * @param maxy the maxy
   * @param uniquePoints the unique points
   * @param nmPerPixel the nm per pixel
   * @param peakDensity the peak density
   * @param binaryImage Set to true if pixels in the image have a 1/0 value. Otherwise it is assumed
   *        the image has continuous real values, i.e. pixels can have 1 or more localisations.
   * @param gr the correlation curve
   * @param spatialDomain the spatial domain
   */
  CorrelationResult(int id, String source, double minx, double miny, double maxx, double maxy,
      double uniquePoints, double nmPerPixel, double peakDensity, boolean binaryImage,
      double[][] gr, boolean spatialDomain) {
    this.id = id;
    this.source = source;
    this.minx = minx;
    this.miny = miny;
    this.maxx = maxx;
    this.maxy = maxy;
    this.uniquePoints = uniquePoints;
    this.nmPerPixel = nmPerPixel;
    this.peakDensity = peakDensity;
    this.binaryImage = binaryImage;
    this.gr = gr;
    this.spatialDomain = spatialDomain;
  }

  /**
   * Compare the results using the Id.
   *
   * @param r1 the first result
   * @param r2 the second result
   * @return the comparison (-1, 0, or 1)
   */
  static int compare(CorrelationResult r1, CorrelationResult r2) {
    return Integer.compare(r1.id, r2.id);
  }
}
