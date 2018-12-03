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
package uk.ac.sussex.gdsc.smlm.ij.plugins;

import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.TestUtils;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import org.opentest4j.AssertionFailedError;

import java.util.Iterator;

/**
 * Test the PulseActivationAnalysis unmixing functions.
 */
@SuppressWarnings({"javadoc"})
public class PulseActivationAnalysisTest {
  @SeededTest
  public void canLinearlyUnmix2Channels(RandomSeed seed) {
    final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
    for (int n = 0; n <= 2; n++) {
      for (int m = 0; m <= 2; m++) {
        canLinearlyUnmix2Channels(r, n, m);
      }
    }
  }

  private static void canLinearlyUnmix2Channels(UniformRandomProvider r, int n, int m) {
    try {
      for (int loop = 0; loop < 10; loop++) {
        // A rough mix of each channel
        final double[] d = create(2, r, 1, 100);

        // Crosstalk should be below 50%
        final double[] c = create(2, r, 0, 0.5);

        // Enumerate
        final Iterator<int[]> it = CombinatoricsUtils.combinationsIterator(2, n);
        while (it.hasNext()) {
          final int[] channels = it.next();
          final double[] dd = new double[2];
          for (final int i : channels) {
            dd[i] = d[i];
          }

          final Iterator<int[]> it2 = CombinatoricsUtils.combinationsIterator(2, m);
          while (it2.hasNext()) {
            final int[] crosstalk = it2.next();
            final double[] cc = new double[2];
            for (final int i : crosstalk) {
              cc[i] = c[i];
            }

            canLinearlyUnmix2Channels(dd[0], dd[1], cc[0], cc[1]);
          }
        }
      }
    } catch (final AssertionFailedError ex) {
      TestUtils.wrapAssertionFailedError(ex, () -> String.format("channels=%d, crosstalk=%d", n, m));
    }
  }

  private static double[] create(int size, UniformRandomProvider r, double min, double range) {
    final double[] d = new double[size];
    for (int i = 0; i < size; i++) {
      d[i] = min + r.nextDouble() * range;
    }
    return d;
  }

  private static void canLinearlyUnmix2Channels(double d1, double d2, double C21, double C12) {
    // Solving:
    // D1 = d1 + C21 * d2
    // D2 = d2 + C12 * d1

    final double D1 = d1 + C21 * d2;
    final double D2 = d2 + C12 * d1;

    final double[] d = PulseActivationAnalysis.unmix(D1, D2, C21, C12);

    Assertions.assertEquals(d1, d[0], delta(d1), "d1");
    Assertions.assertEquals(d2, d[1], delta(d2), "d2");
  }

  private static double delta(double d) {
    d *= 1e-10;
    // Set a limit where we are close enough. This mainly applies when d is zero.
    if (d < 1e-8) {
      return 1e-8;
    }
    return d;
  }

  @SeededTest
  public void canLinearlyUnmix3Channels(RandomSeed seed) {
    final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
    for (int n = 0; n <= 3; n++) {
      for (int m = 0; m <= 6; m++) {
        canLinearlyUnmix3Channels(r, n, m);
      }
    }
  }

  private void canLinearlyUnmix3Channels(UniformRandomProvider r, int n, int m) {
    try {
      for (int loop = 0; loop < 10; loop++) {
        // A rough mix of each channel
        final double[] d = create(6, r, 100, 100);

        // Total crosstalk per channel should be below 50%
        final double[] c = create(6, r, 0, 0.25);

        // Enumerate
        final Iterator<int[]> it = CombinatoricsUtils.combinationsIterator(3, n);
        while (it.hasNext()) {
          final int[] channels = it.next();
          final double[] dd = new double[3];
          for (final int i : channels) {
            dd[i] = d[i];
          }

          final Iterator<int[]> it2 = CombinatoricsUtils.combinationsIterator(6, m);
          while (it2.hasNext()) {
            final int[] crosstalk = it2.next();
            final double[] cc = new double[6];
            for (final int i : crosstalk) {
              cc[i] = c[i];
            }

            canLinearlyUnmix3Channels(dd[0], dd[1], dd[2], cc[0], cc[1], cc[2], cc[3], cc[4],
                cc[5]);
          }
        }
      }
    } catch (final AssertionFailedError ex) {
      TestUtils.wrapAssertionFailedError(ex, () -> String.format("channels=%d, crosstalk=%d", n, m));
    }
  }

  public void canLinearlyUnmix3Channels(double d1, double d2, double d3, double C21, double C31,
      double C12, double C32, double C13, double C23) {
    // Solving:
    // D1 = d1 + C21 * d2 + C31 * d3
    // D2 = d2 + C12 * d1 + C32 * d3
    // D3 = d3 + C13 * d1 + C23 * d2

    final double D1 = d1 + C21 * d2 + C31 * d3;
    final double D2 = d2 + C12 * d1 + C32 * d3;
    final double D3 = d3 + C13 * d1 + C23 * d2;

    final double[] d = PulseActivationAnalysis.unmix(D1, D2, D3, C21, C31, C12, C32, C13, C23);

    Assertions.assertEquals(d1, d[0], delta(d1), "d1");
    Assertions.assertEquals(d2, d[1], delta(d2), "d2");
    Assertions.assertEquals(d3, d[2], delta(d3), "d3");
  }
}
