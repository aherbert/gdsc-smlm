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

package uk.ac.sussex.gdsc.smlm.results;

import org.apache.commons.math3.util.FastMath;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;

/**
 * Define a cluster of localisations.
 */
public class Cluster {
  /**
   * The method for computing the centroid.
   */
  public enum CentroidMethod {
    /** Weight each localisation evenly. */
    STANDARD,
    /** Weight each localisation using the signal. */
    SIGNAL_WEIGHTED
  }

  /** The results. */
  protected PeakResultStoreList results = new ArrayPeakResultStore(2);
  private float[] centroid;
  private int id;

  /**
   * Instantiates a new cluster.
   */
  public Cluster() {}

  /**
   * Instantiates a new cluster.
   *
   * @param result the result
   */
  public Cluster(PeakResult result) {
    add(result);
  }

  /**
   * Get the size.
   *
   * @return the size
   */
  public int size() {
    return results.size();
  }

  /**
   * Gets the points.
   *
   * @return the points
   */
  public PeakResultStoreList getPoints() {
    return results;
  }

  /**
   * Adds the result.
   *
   * @param result the result
   */
  public void add(PeakResult result) {
    results.add(result);
    centroid = null;
  }

  /**
   * Checks if is empty.
   *
   * @return true, if is empty
   */
  public boolean isEmpty() {
    return size() == 0;
  }

  /**
   * Gets the centroid.
   *
   * @param method the method
   * @return the centroid
   */
  public float[] getCentroid(CentroidMethod method) {
    if (centroid == null && results.size() != 0) {
      switch (method) {
        case SIGNAL_WEIGHTED:
          final float[] weights = new float[results.size()];
          for (int i = 0; i < results.size(); i++) {
            weights[i] = Math.abs(results.get(i).getIntensity());
          }
          // Normalise weights?
          return getCentroid(results, weights);

        case STANDARD:
        default:
          return getCentroid();
      }
    }
    return centroid;
  }

  private float[] getCentroid(PeakResultStoreList results, float[] weights) {
    centroid = new float[2];
    double sum = 0;
    for (int i = 0; i < results.size(); i++) {
      final PeakResult result = results.get(i);
      final float w = weights[i++];
      sum += w;
      centroid[0] += result.getXPosition() * w;
      centroid[1] += result.getYPosition() * w;
    }
    centroid[0] /= sum;
    centroid[1] /= sum;
    return centroid;
  }

  /**
   * Gets the centroid.
   *
   * @return the centroid
   */
  public float[] getCentroid() {
    if (centroid == null && results.size() != 0) {
      centroid = new float[2];
      for (int i = 0; i < results.size(); i++) {
        final PeakResult result = results.get(i);
        centroid[0] += result.getXPosition();
        centroid[1] += result.getYPosition();
      }
      centroid[0] /= results.size();
      centroid[1] /= results.size();
    }
    return centroid;
  }

  /**
   * Remove the calculated centroid from memory (forces a refresh of the centroid).
   */
  public void resetCentroid() {
    centroid = null;
  }

  /**
   * Gets the standard deviation of the distances from the centroid.
   *
   * @return The standard deviation of the distances from the centroid
   */
  public float getStandardDeviation() {
    if (results.size() < 2) {
      return 0;
    }
    getCentroid();
    double ssx = 0;
    for (int i = 0; i < results.size(); i++) {
      final PeakResult result = results.get(i);
      final double dx = result.getXPosition() - centroid[0];
      final double dy = result.getYPosition() - centroid[1];
      final double d2 = dx * dx + dy * dy;
      ssx += d2;
    }
    return (float) Math.sqrt(ssx / (results.size() - 1));
  }

  /**
   * Calculate the weighted localisation precision using the PC-PALM formula of Sengupta, et al
   * (2013) Nature Protocols 8, 345.
   *
   * <p>Also sets the centroid if it has not been calculated using the signal weighted
   * centre-of-mass
   *
   * <p>Note that the PeakResult must have valid values in the precision field, otherwise a value of
   * 1 is used.
   *
   * @param converter the converter to convert the distances to nm
   * @return The weighted localisation precision of the group peak (in nm)
   */
  public double getLocalisationPrecision(TypeConverter<DistanceUnit> converter) {
    if (converter == null || converter.to() != DistanceUnit.NM) {
      return 0;
    }

    final int n = size();
    if (n == 0) {
      centroid = null;
      return 0;
    }

    if (n == 1) {
      final PeakResult result = results.get(0);
      if (centroid == null) {
        centroid = new float[] {result.getXPosition(), result.getYPosition()};
      }
      return checkPrecision(result.getPrecision());
    }

    final float[] photons = new float[results.size()];
    for (int i = 0; i < results.size(); i++) {
      final PeakResult result = results.get(i);
      photons[i++] = Math.abs(result.getIntensity());
    }

    double sumNi = 0;
    double xm = 0;
    double ym = 0;
    for (int i = 0; i < results.size(); i++) {
      final PeakResult result = results.get(i);
      final float Ni = photons[i++];
      sumNi += Ni;
      xm += result.getXPosition() * Ni;
      ym += result.getYPosition() * Ni;
    }
    xm /= sumNi;
    ym /= sumNi;

    if (centroid == null) {
      centroid = new float[] {(float) xm, (float) ym};
    }

    double sumXi2Ni = 0;
    double sumYi2Ni = 0;
    double sumS2 = 0;
    for (int i = 0; i < results.size(); i++) {
      final PeakResult result = results.get(i);
      final float Ni = photons[i++];

      final double dx = converter.convert(result.getXPosition() - xm);
      final double dy = converter.convert(result.getYPosition() - ym);

      sumXi2Ni += dx * dx * Ni;
      sumYi2Ni += dy * dy * Ni;
      sumS2 += MathUtils.pow2(checkPrecision(result.getPrecision())) * Ni;
    }

    final double sumNin = sumNi * n;
    final double sumS2_sumNin = sumS2 / sumNin;
    final double sxm = Math.sqrt(sumXi2Ni / sumNin + sumS2_sumNin) / 1.414213562;
    final double sym = Math.sqrt(sumYi2Ni / sumNin + sumS2_sumNin) / 1.414213562;

    final double sPeak = FastMath.max(sxm, sym);

    return sPeak;
  }

  private static double checkPrecision(double precision) {
    return (precision > 0 && precision < Double.MAX_VALUE) ? precision : 1;

  }

  /**
   * Gets the first PeakResult in the cluster (or null).
   *
   * @return The first PeakResult in the cluster (or null)
   */
  public PeakResult getHead() {
    if (isEmpty()) {
      return null;
    }
    return results.get(0);
  }

  /**
   * Gets the last PeakResult in the cluster (or null).
   *
   * @return The last PeakResult in the cluster (or null)
   */
  public PeakResult getTail() {
    if (isEmpty()) {
      return null;
    }
    return results.get(results.size() - 1);
  }

  /**
   * Gets the result from the set.
   *
   * @param index the index
   * @return the peak result
   */
  public PeakResult get(int index) {
    return results.get(index);
  }

  /**
   * Gets the signal.
   *
   * @return The total signal
   */
  public double getSignal() {
    double sum = 0;
    for (int i = 0; i < results.size(); i++) {
      final PeakResult result = results.get(i);
      sum += result.getIntensity();
    }
    return sum;
  }

  /**
   * Sort in time order.
   */
  public void sort() {
    results.sort();
  }

  /**
   * Gets the id.
   *
   * @return the id
   */
  public int getId() {
    return id;
  }

  /**
   * Sets the id.
   *
   * @param id the new id
   */
  public void setId(int id) {
    this.id = id;
  }

  /**
   * Expand any localisations that have a different start and end frame into a series. Note that
   * this will increase the size of the cluster.
   *
   * <p>The results are copies save for the end frame. This makes analysis of the signal invalid as
   * it will have been increased n-fold for each localisation that spans n frames. The original
   * multi-frame result is removed.
   */
  public void expandToSingles() {
    // Check for expansion
    int singles = 0;
    while (singles < results.size()) {
      final PeakResult result = results.get(singles);
      if (result.getFrame() != result.getEndFrame()) {
        break;
      }
      singles++;
    }

    if (singles == size()) {
      return;
    }

    final PeakResultStoreList newResults = new ArrayPeakResultStore(size());
    for (int i = 0; i < singles; i++) {
      newResults.add(results.get(i));
    }

    for (int i = singles; i < results.size(); i++) {
      final PeakResult result = results.get(i);
      if (result.getFrame() != result.getEndFrame()) {
        for (int peak = result.getFrame(); peak <= result.getEndFrame(); peak++) {
          newResults.add(new ExtendedPeakResult(peak, result.getOrigX(), result.getOrigY(),
              result.getOrigValue(), result.getError(), result.getNoise(),
              result.getMeanIntensity(), result.getParameters(), result.getParameterDeviations(),
              peak, result.getId()));
        }
      } else {
        newResults.add(result);
      }
    }

    results = newResults;
    resetCentroid();
  }

  /**
   * Remove the first and last result. If the size is 2 or less then the new size will be zero.
   */
  public void removeEnds() {
    if (size() <= 2) {
      results.clear();
    } else {
      final PeakResultStoreList newResults = new ArrayPeakResultStore(size() - 2);
      for (int i = 1, size = size() - 1; i < size; i++) {
        newResults.add(results.get(i));
      }
      results = newResults;
    }
    resetCentroid();
  }

  /**
   * Gets the mean-squared distance (MSD) between adjacent localisations.
   *
   * @return The mean squared-distance between adjacent localisations.
   */
  public double getMsd() {
    if (size() < 2) {
      return 0;
    }
    double sum = 0;
    PeakResult last = results.get(0);
    for (int i = 1; i < results.size(); i++) {
      final PeakResult result = results.get(i);
      sum += last.distance2(result);
      last = result;
    }
    return sum / (size() - 1);
  }

  /**
   * Gets the mean distance between adjacent localisations.
   *
   * @return The mean distance between adjacent localisations.
   */
  public double getMeanDistance() {
    if (size() < 2) {
      return 0;
    }
    double sum = 0;
    PeakResult last = results.get(0);
    for (int i = 1; i < results.size(); i++) {
      final PeakResult result = results.get(i);
      sum += last.distance(result);
      last = result;
    }
    return sum / (size() - 1);
  }
}
