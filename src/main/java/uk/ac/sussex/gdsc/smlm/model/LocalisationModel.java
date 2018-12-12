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

import java.util.Arrays;

/**
 * Contains a discrete-time model for the position and intensity of a localisation.
 */
public class LocalisationModel implements Comparable<LocalisationModel> {
  private double[] xyz;
  private int id;
  private int time;
  private double intensity;
  private int state;
  private LocalisationModel previous;
  private LocalisationModel next;
  private double[] data;
  private int label;

  /** A single localisation. */
  public static final int SINGLE = 0;
  /** A localisation with a localisation in the previous frame. */
  public static final int PREVIOUS = 1;
  /** A localisation with a localisation in the next frame. */
  public static final int NEXT = 2;
  /** A localisation with a localisation in the previous or next frame. */
  public static final int NEIGHBOUR = PREVIOUS | NEXT;
  /** A localisation with a localisation in both the previous and next frame. */
  public static final int CONTINUOUS = 4 | PREVIOUS | NEXT;

  /**
   * Create a new localisation.
   *
   * @param id the id
   * @param time the time
   * @param x the x
   * @param y the y
   * @param z the z
   * @param intensity the intensity
   * @param state the state
   */
  public LocalisationModel(int id, int time, double x, double y, double z, double intensity,
      int state) {
    init(id, time, new double[] {x, y, z}, intensity, state, null, null);
  }

  /**
   * Create a new localisation.
   *
   * @param id the id
   * @param time the time
   * @param xyz Coordinates (will be deep copied)
   * @param intensity the intensity
   * @param state the state
   */
  public LocalisationModel(int id, int time, double[] xyz, double intensity, int state) {
    init(id, time, Arrays.copyOf(xyz, xyz.length), intensity, state, null, null);
  }

  /**
   * Create a new localisation.
   *
   * @param id the id
   * @param time the time
   * @param x the x
   * @param y the y
   * @param z the z
   * @param intensity the intensity
   * @param state the state
   * @param previous The previous localisation in this pulse (should be continuous time points)
   * @param next The next localisation in this pulse (should be continuous time points)
   */
  public LocalisationModel(int id, int time, double x, double y, double z, double intensity,
      int state, LocalisationModel previous, LocalisationModel next) {
    init(id, time, new double[] {x, y, z}, intensity, state, previous, next);
  }

  /**
   * Create a new localisation.
   *
   * @param id the id
   * @param time the time
   * @param xyz Coordinates (will be deep copied)
   * @param intensity the intensity
   * @param state the state
   * @param previous The previous localisation in this pulse (should be continuous time points)
   * @param next The next localisation in this pulse (should be continuous time points)
   */
  public LocalisationModel(int id, int time, double[] xyz, double intensity, int state,
      LocalisationModel previous, LocalisationModel next) {
    init(id, time, Arrays.copyOf(xyz, xyz.length), intensity, state, previous, next);
  }

  private void init(int id, int time, double[] xyz, double intensity, int state,
      LocalisationModel previous, LocalisationModel next) {
    this.xyz = xyz;
    this.id = id;
    this.time = time;
    this.intensity = intensity;
    this.state = state;
    this.previous = previous;
    this.next = next;
  }

  /**
   * Gets the x.
   *
   * @return the x
   */
  public double getX() {
    return xyz[0];
  }

  /**
   * Gets the y.
   *
   * @return the y
   */
  public double getY() {
    return xyz[1];
  }

  /**
   * Gets the z.
   *
   * @return the z
   */
  public double getZ() {
    return xyz[2];
  }

  /**
   * Gets the coordinates.
   *
   * @return The coordinates (x,y,z)
   */
  public double[] getCoordinates() {
    return xyz;
  }

  /**
   * Gets the id.
   *
   * @return The Id
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
   * @return The time
   */
  public int getTime() {
    return time;
  }

  /**
   * Gets the intensity.
   *
   * @return the intensity
   */
  public double getIntensity() {
    return intensity;
  }

  /**
   * Sets the intensity.
   *
   * @param intensity The new intensity
   */
  public void setIntensity(double intensity) {
    this.intensity = intensity;
  }

  /** {@inheritDoc} */
  @Override
  public int compareTo(LocalisationModel o) {
    if (time == o.time) {
      return Double.compare(o.intensity, intensity);
    }
    return (time < o.time) ? -1 : 1;
  }

  /**
   * Checks for previous.
   *
   * @return True if this localisation is on in the previous time interval
   */
  public boolean hasPrevious() {
    return (state & PREVIOUS) == PREVIOUS;
  }

  /**
   * Checks for next.
   *
   * @return True if this localisation is on in the next time interval
   */
  public boolean hasNext() {
    return (state & NEXT) == NEXT;
  }

  /**
   * Checks for neighbour.
   *
   * @return True if this localisation is on in the previous or next time interval
   */
  public boolean hasNeighbour() {
    return (state & (NEIGHBOUR)) != 0;
  }

  /**
   * Checks if is continuous.
   *
   * @return True if this localisation is on for the entire duration of the time interval
   */
  public boolean isContinuous() {
    return (state & CONTINUOUS) == CONTINUOUS;
  }

  /**
   * Gets the previous.
   *
   * @return the previous
   */
  public LocalisationModel getPrevious() {
    return previous;
  }

  /**
   * Sets the previous.
   *
   * @param previous the previous to set
   */
  public void setPrevious(LocalisationModel previous) {
    this.previous = previous;
    if (previous != null) {
      previous.next = this;
    }
  }

  /**
   * Gets the next.
   *
   * @return the next
   */
  public LocalisationModel getNext() {
    return next;
  }

  /**
   * Sets the next.
   *
   * @param next the next to set
   */
  public void setNext(LocalisationModel next) {
    this.next = next;
    if (next != null) {
      next.previous = this;
    }
  }

  /**
   * Gets the data.
   *
   * @return the data
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
   * Gets the label.
   *
   * @return the label
   */
  public int getLabel() {
    return label;
  }

  /**
   * Sets the label. This can be used to identify subsets of molecules.
   *
   * @param label the new label
   */
  public void setLabel(int label) {
    this.label = label;
  }
}
