/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Contains a collection of localisations associated with a discrete-time. This can be used to model
 * a moving localisation simulated on a more refined time-scale to a single localisation.
 */
public class LocalisationModelSet implements Comparable<LocalisationModelSet> {
  private int id;
  private final int time;
  private final List<LocalisationModel> localisations = new ArrayList<>();
  private double[] data;
  private LocalisationModelSet previous;
  private LocalisationModelSet next;

  /**
   * Create a new localisation.
   *
   * @param id the id
   * @param time the time
   */
  public LocalisationModelSet(int id, int time) {
    this.id = id;
    this.time = time;
  }

  /**
   * Adds the localisation.
   *
   * @param localisation the localisation
   */
  public void add(LocalisationModel localisation) {
    localisations.add(localisation);
  }

  /**
   * Get the size.
   *
   * @return the size
   */
  public int size() {
    return localisations.size();
  }

  /**
   * Gets the localisations.
   *
   * @return the localisations
   */
  public List<LocalisationModel> getLocalisations() {
    return localisations;
  }

  /**
   * Gets the id.
   *
   * @return The Id.
   */
  public int getId() {
    return id;
  }

  /**
   * Allow the package to set the id.
   *
   * @param id The Id
   */
  void setId(int id) {
    this.id = id;
  }

  /**
   * Gets the time.
   *
   * @return The time.
   */
  public int getTime() {
    return time;
  }

  /** {@inheritDoc} */
  @Override
  public int compareTo(LocalisationModelSet other) {
    return Integer.compare(time, other.time);
  }

  /**
   * Checks if this localisation is on for the entire duration of the time interval.
   *
   * @return True if this localisation is on for the entire duration of the time interval.
   */
  public boolean isContinuous() {
    if (localisations.isEmpty()) {
      return false;
    }

    // All localisations must be continuous and consecutive in time
    final int[] t = new int[localisations.size()];
    int count = 0;
    for (final LocalisationModel l : localisations) {
      t[count++] = l.getTime();
      if (!l.isContinuous()) {
        return false;
      }
    }

    // Check consecutive in time
    Arrays.sort(t);

    return ((t[t.length - 1] - t[0] + 1) <= count);
  }

  /**
   * Gets the data.
   *
   * @return the data.
   */
  public double[] getData() {
    return data;
  }

  /**
   * Sets the data.
   *
   * @param data the data to set
   */
  public void setData(double[] data) {
    this.data = data;
  }

  /**
   * Convert the set of localisations to a single localisation with the combined signal and the
   * centroid location (centre-of-mass weighted by intensity).
   *
   * @return the localisation model
   */
  public LocalisationModel toLocalisation() {
    double intensity = 0;
    final double[] xyz = new double[3];
    for (final LocalisationModel l : localisations) {
      final double s = l.getIntensity();
      intensity += s;
      final double[] xyz2 = l.getCoordinates();
      for (int i = 0; i < 3; i++) {
        xyz[i] += xyz2[i] * s;
      }
    }
    if (!localisations.isEmpty()) {
      for (int i = 0; i < 3; i++) {
        xyz[i] /= intensity;
      }
    }

    final LocalisationModel l = new LocalisationModel(id, time, xyz, intensity,
        isContinuous() ? LocalisationModel.CONTINUOUS : LocalisationModel.SINGLE);
    l.setData(data);
    return l;
  }

  /**
   * Gets the total intensity.
   *
   * @return The total intensity.
   */
  public double getIntensity() {
    double intensity = 0;
    for (final LocalisationModel l : localisations) {
      intensity += l.getIntensity();
    }
    return intensity;
  }

  /**
   * Gets the previous model set.
   *
   * @return the previous model set
   */
  public LocalisationModelSet getPrevious() {
    return previous;
  }

  /**
   * Sets the previous model set.
   *
   * @param previous the previous to set
   */
  public void setPrevious(LocalisationModelSet previous) {
    this.previous = previous;
    if (previous != null) {
      previous.next = this;
    }
  }

  /**
   * Gets the next model set.
   *
   * @return the next model set
   */
  public LocalisationModelSet getNext() {
    return next;
  }

  /**
   * Sets the next model set.
   *
   * @param next the next to set
   */
  public void setNext(LocalisationModelSet next) {
    this.next = next;
    if (next != null) {
      next.previous = this;
    }
  }

  /**
   * Checks if either of the previous/next pointers are not null.
   *
   * @return True if either of the previous/next pointers are not null
   */
  public boolean hasNeighbour() {
    return next != null || previous != null;
  }
}
