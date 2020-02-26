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

import java.util.Iterator;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import org.opentest4j.AssertionFailedError;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.TestUtils;

/**
 * Test the PulseActivationAnalysis unmixing functions.
 */
@SuppressWarnings({"javadoc"})
public class PulseActivationAnalysisTest {
  // CHECKSTYLE.OFF: ParameterName

  @SeededTest
  public void canLinearlyUnmix2Channels(RandomSeed seed) {
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    for (int n = 1; n <= 2; n++) {
      for (int m = 1; m <= 2; m++) {
        canLinearlyUnmix2Channels(r, n, m);
      }
    }
  }

  /**
   * Check if can linearly unmix 2 channels.
   *
   * <p>Randomly creates density data for each channel and then crosstalk between them. Only some of
   * the density and cross talk are sampled to test unmixing random combinations.
   *
   * @param rng the random generator
   * @param n the number of channels that are selected to have density
   * @param m the number of crosstalk values that are selected from the total
   */
  private static void canLinearlyUnmix2Channels(UniformRandomProvider rng, int n, int m) {
    try {
      for (int loop = 0; loop < 10; loop++) {
        // A rough mix of each channel
        final double[] d = createData(2, rng, 1, 100);

        // Crosstalk should be below 50%
        final double[] c = createData(2, rng, 0, 0.5);

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
      TestUtils.wrapAssertionFailedError(ex,
          () -> String.format("channels=%d, crosstalk=%d", n, m));
    }
  }

  private static void canLinearlyUnmix2Channels(double d1, double d2, double c21, double c12) {
    // Solving:
    // D1 = d1 + c21 * d2
    // D2 = d2 + c12 * d1

    final double D1 = d1 + c21 * d2;
    final double D2 = d2 + c12 * d1;

    final double[] d = PulseActivationAnalysis.unmix(D1, D2, c21, c12);

    Assertions.assertEquals(d1, d[0], delta(d1), "d1");
    Assertions.assertEquals(d2, d[1], delta(d2), "d2");
  }

  private static double delta(double value) {
    value *= 1e-10;
    // Set a limit where we are close enough. This mainly applies when d is zero.
    if (value < 1e-8) {
      return 1e-8;
    }
    return value;
  }

  @SeededTest
  public void canLinearlyUnmix3Channels(RandomSeed seed) {
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    for (int n = 1; n <= 3; n++) {
      for (int m = 1; m <= 6; m++) {
        canLinearlyUnmix3Channels(r, n, m);
      }
    }
  }

  /**
   * Check if can linearly unmix 3 channels.
   *
   * <p>Randomly creates density data for each channel and then crosstalk between them. Only some of
   * the density and cross talk are sampled to test unmixing random combinations.
   *
   * @param rng the random generator
   * @param n the number of channels that are selected to have density
   * @param m the number of crosstalk values that are selected from the total
   */
  private static void canLinearlyUnmix3Channels(UniformRandomProvider rng, int n, int m) {
    try {
      for (int loop = 0; loop < 10; loop++) {
        // A rough mix of each channel
        final double[] d = createData(6, rng, 100, 100);

        // Total crosstalk per channel should be below 50%
        final double[] c = createData(6, rng, 0, 0.25);

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
      TestUtils.wrapAssertionFailedError(ex,
          () -> String.format("channels=%d, crosstalk=%d", n, m));
    }
  }

  private static void canLinearlyUnmix3Channels(double d1, double d2, double d3, double c21,
      double c31, double c12, double c32, double c13, double c23) {
    // Solving:
    // D1 = d1 + c21 * d2 + c31 * d3
    // D2 = d2 + c12 * d1 + c32 * d3
    // D3 = d3 + c13 * d1 + c23 * d2

    final double D1 = d1 + c21 * d2 + c31 * d3;
    final double D2 = d2 + c12 * d1 + c32 * d3;
    final double D3 = d3 + c13 * d1 + c23 * d2;

    final double[] d = PulseActivationAnalysis.unmix(D1, D2, D3, c21, c31, c12, c32, c13, c23);

    Assertions.assertEquals(d1, d[0], delta(d1), "d1");
    Assertions.assertEquals(d2, d[1], delta(d2), "d2");
    Assertions.assertEquals(d3, d[2], delta(d3), "d3");
  }

  private static double[] createData(int size, UniformRandomProvider rng, double min,
      double range) {
    final double[] data = new double[size];
    for (int i = 0; i < size; i++) {
      data[i] = min + rng.nextDouble() * range;
    }
    return data;
  }
}
