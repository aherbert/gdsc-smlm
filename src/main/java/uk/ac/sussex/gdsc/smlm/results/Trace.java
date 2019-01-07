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

import gnu.trove.list.linked.TIntLinkedList;

/**
 * Define a cluster of localisations from different frames that represent a single molecule trace.
 */
public class Trace extends Cluster {
  private int blinks = -1;
  private int[] onTimes;
  private int[] offTimes;

  /**
   * Instantiates a new trace.
   */
  public Trace() {
    super();
  }

  /**
   * Instantiates a new trace.
   *
   * @param peakResult the peak result
   */
  public Trace(PeakResult peakResult) {
    super(peakResult);
  }

  @Override
  public void add(PeakResult result) {
    super.add(result);
    blinks = -1; // Invalidate the analysis
  }

  private void analyse() {
    if (blinks == -1) {
      if (isEmpty()) {
        blinks = 0;
        onTimes = offTimes = null;
        return;
      }
      if (results.size() == 1) {
        blinks = 1;
        onTimes = new int[] {1};
        offTimes = null;
        return;
      }

      // Ensure in the correct time-order
      sort();
      final TIntLinkedList on = new TIntLinkedList();
      final TIntLinkedList off = new TIntLinkedList();

      blinks = 1;
      int t1 = results.get(0).getFrame();
      int onStart = t1;
      for (int i = 0; i < results.size() - 1; i++) {
        final int t2 = results.get(i + 1).getFrame();
        final int diff = t2 - t1;
        if (diff > 1) {
          off.add(diff - 1);
          on.add(t1 - onStart + 1);
          blinks++;
          onStart = t2;
        }
        t1 = t2;
      }
      on.add(t1 - onStart + 1);

      onTimes = on.toArray();
      offTimes = off.toArray();
    }
  }

  /**
   * Gets the number of times the molecule blinked.
   *
   * @return The number of times the molecule blinked.
   */
  public int getBlinks() {
    analyse();
    return blinks;
  }

  /**
   * Gets the average on time for the molecule.
   *
   * @return The average on time for the molecule.
   */
  public double getOnTime() {
    analyse();
    return getAverage(onTimes);
  }

  /**
   * Gets the average off time for the molecule.
   *
   * @return The average off time for the molecule.
   */
  public double getOffTime() {
    analyse();
    return getAverage(offTimes);
  }

  private static double getAverage(int[] times) {
    if (times != null) {
      double av = 0;
      for (final int t : times) {
        av += t;
      }
      return av / times.length;
    }
    return 0;
  }

  /**
   * Gets the on times.
   *
   * @return the on-times.
   */
  public int[] getOnTimes() {
    analyse();
    return onTimes;
  }

  /**
   * Gets the off times.
   *
   * @return the off-times.
   */
  public int[] getOffTimes() {
    analyse();
    return offTimes;
  }

  @Override
  public void removeEnds() {
    super.removeEnds();
    blinks = -1; // Invalidate the analysis
  }
}
