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

package uk.ac.sussex.gdsc.smlm.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Contains a model for a blinking fluorophore.
 *
 * <p>Based on the work of Coltharp et al (2012) Accurate Construction of photoactivated
 * localization microscopy images for quantitative measurements. PLOS One 7, Issue 12, pp 1-15
 */
public abstract class FluorophoreSequenceModel extends MoleculeModel {
  /**
   * Instantiates a new fluorophore sequence model.
   *
   * @param id the id
   * @param x the x
   * @param y the y
   * @param z the z
   */
  public FluorophoreSequenceModel(int id, double x, double y, double z) {
    super(id, x, y, z);
  }

  /**
   * Instantiates a new fluorophore sequence model.
   *
   * @param id the id
   * @param xyz the xyz
   */
  public FluorophoreSequenceModel(int id, double[] xyz) {
    super(id, xyz);
  }

  /**
   * The number of times the molecule went into the dark state.
   */
  private int blinks;
  /**
   * A sequence of fluorescent bursts in pairs of {on,off} times. The burst sequence will be length
   * = 2 * (blinks+1)
   */
  private double[] burstSequence = new double[] {0, 0};

  /**
   * Sets the burst sequence.
   *
   * @param sequence the new burst sequence
   */
  protected void setBurstSequence(double[] sequence) {
    if (sequence != null && sequence.length > 1) {
      blinks = (sequence.length / 2) - 1;

      // Ensure the sequence array is an even number in length
      final int length = 2 * (blinks + 1);
      if (sequence.length == length) {
        burstSequence = sequence;
      } else {
        burstSequence = Arrays.copyOf(sequence, length);
      }
    }
  }

  /**
   * Gets the number of blinks.
   *
   * @return The number of times the fluorophore blinked
   */
  public int getNumberOfBlinks() {
    return blinks;
  }

  /**
   * Get the start time, i.e. when the molecule activated.
   *
   * <p>Note that a molecule will always have a start time even if it has no blinks. This models a
   * molecule that turns on and then bleaches immediately.
   *
   * @return The start time
   */
  public double getStartTime() {
    return burstSequence[0];
  }

  /**
   * Get the end time, i.e. when the molecule bleached.
   *
   * @return The end time
   */
  public double getEndTime() {
    return burstSequence[burstSequence.length - 1];
  }

  /**
   * Gets the burst sequence.
   *
   * @return Fluorescent bursts arranged as list of on/off times: {onT,offT}
   */
  public List<double[]> getBurstSequence() {
    final ArrayList<double[]> data = new ArrayList<>(blinks + 1);
    for (int i = 0; i <= blinks; i++) {
      data.add(new double[] {burstSequence[i * 2], burstSequence[i * 2 + 1]});
    }
    return data;
  }

  /**
   * Gets the sampled burst sequence.
   *
   * @return Fluorescent bursts arranged as list of on/off times in integer sampling intervals:
   *         {onT,offT}
   */
  public List<int[]> getSampledBurstSequence() {
    final ArrayList<int[]> data = new ArrayList<>(blinks + 1);
    for (int i = 0; i <= blinks; i++) {
      data.add(new int[] {(int) (burstSequence[i * 2]), (int) (burstSequence[i * 2 + 1])});
    }
    return data;
  }

  /**
   * Order by time ascending.
   *
   * @param o1 the first model
   * @param o2 The second model
   * @return The comparison result -1, 0, or 1
   */
  public static int compare(FluorophoreSequenceModel o1, FluorophoreSequenceModel o2) {
    return Double.compare(o1.getStartTime(), o2.getStartTime());
  }

  /**
   * Gets the on times.
   *
   * @return The duration of the on times
   */
  public double[] getOnTimes() {
    final double[] onTimes = new double[blinks + 1];
    for (int i = 0; i <= blinks; i++) {
      onTimes[i] = burstSequence[i * 2 + 1] - burstSequence[i * 2];
    }
    return onTimes;
  }

  /**
   * Gets the off times.
   *
   * @return The duration of the off times
   */
  public double[] getOffTimes() {
    if (blinks < 1) {
      return new double[0];
    }

    final double[] offTimes = new double[blinks];
    for (int i = 1; i <= blinks; i++) {
      offTimes[i - 1] = burstSequence[i * 2] - burstSequence[i * 2 - 1];
    }
    return offTimes;
  }

  /**
   * Gets the sampled on times.
   *
   * @return The duration of the on times if sampled at integer time intervals
   */
  public int[] getSampledOnTimes() {
    if (blinks == 0) {
      return new int[] {end(burstSequence[1]) - start(burstSequence[0])};
    }

    // Process all blinks. Join together blinks with an off-time that would not be noticed,
    // i.e. where the molecule was on in consecutive frames.
    final int[] onTimes = new int[blinks + 1];
    int count = 0;
    int tstart = (int) burstSequence[0];
    for (int i = 0; i < blinks; i++) {
      final int end1 = end(burstSequence[i * 2 + 1]);
      final int start2 = start(burstSequence[(i + 1) * 2]);

      if (start2 - end1 > 0) {
        onTimes[count++] = end1 - tstart;
        tstart = start2;
      }
    }
    onTimes[count++] = end(getEndTime()) - tstart;

    return Arrays.copyOf(onTimes, count);
  }

  /**
   * Convert the start time to an integer.
   *
   * @param time the time
   * @return the integer start time
   */
  private static int start(double time) {
    return (int) time;
  }

  /**
   * Convert the end time to an integer.
   *
   * @param time the time
   * @return the integer end time
   */
  private static int end(double time) {
    return (int) (Math.ceil(time));
  }

  /**
   * Gets the sampled off times.
   *
   * @return The duration of the off times if sampled at integer time intervals
   */
  public int[] getSampledOffTimes() {
    if (blinks == 0) {
      return new int[0];
    }

    // Process all blinks. Join together blinks with an off-time that would not be noticed,
    // i.e. where the molecule was on in consecutive frames.
    final int[] offTimes = new int[blinks];
    int count = 0;
    for (int i = 0; i < blinks; i++) {
      final int end1 = end(burstSequence[i * 2 + 1]);
      final int start2 = start(burstSequence[(i + 1) * 2]);

      if (start2 - end1 > 0) {
        offTimes[count++] = start2 - end1;
      }
    }

    return Arrays.copyOf(offTimes, count);
  }

  /**
   * Gets the on frames.
   *
   * @return An array of frames when the molecule was on
   */
  public int[] getOnFrames() {
    final int sequenceStartT = (int) getStartTime();
    final int sequenceEndT = (int) getEndTime();

    int count = 0;
    final int[] onFrames = new int[sequenceEndT - sequenceStartT + 1];
    for (int i = 0; i <= blinks; i++) {
      final int on = (int) (burstSequence[i * 2]);
      final int off = (int) (burstSequence[i * 2 + 1]);

      for (int t = on; t <= off; t++) {
        onFrames[count++] = t;
      }
    }

    return Arrays.copyOf(onFrames, count);
  }

  /**
   * Scale the times using the specified factor. Allows adjusting the relative time of the sequence.
   *
   * @param scale the scale
   */
  public void adjustTime(double scale) {
    if (scale < 0) {
      throw new IllegalArgumentException("Scale factor must be above zero");
    }
    for (int i = 0; i < burstSequence.length; i++) {
      burstSequence[i] *= scale;
    }
  }
}
