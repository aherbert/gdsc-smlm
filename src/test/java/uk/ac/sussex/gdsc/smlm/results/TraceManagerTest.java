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

import uk.ac.sussex.gdsc.core.utils.rng.RandomUtils;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;

import gnu.trove.list.array.TIntArrayList;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;

import java.util.Arrays;
import java.util.Comparator;

@SuppressWarnings({"javadoc"})
public class TraceManagerTest {
  @SeededTest
  public void canTraceSinglePulseWithFixedCoords(RandomSeed seed) {
    final UniformRandomProvider rand = RngUtils.create(seed.getSeedAsLong());
    final float[] params = createParams(rand);
    final Trace trace = new Trace();
    for (int i = 0; i < 5; i++) {
      trace.add(new PeakResult(i, 0, 0, 0, 0, 0, 0, params, null));
    }

    runTracing(rand, 0, 1, trace);
  }

  @SeededTest
  public void canTraceSinglePulseWithMovingCoords(RandomSeed seed) {
    final UniformRandomProvider rand = RngUtils.create(seed.getSeedAsLong());
    final float distance = 0.5f;

    final float[] params = createParams(rand);
    final Trace trace = new Trace();
    for (int i = 0; i < 5; i++) {
      move(rand, params, distance);
      trace.add(new PeakResult(i, 0, 0, 0, 0, 0, 0, params, null));
    }

    runTracing(rand, distance, 1, trace);
  }

  @SeededTest
  public void canTraceMultiplePulseWithFixedCoords(RandomSeed seed) {
    final UniformRandomProvider rand = RngUtils.create(seed.getSeedAsLong());

    final float[] params = createParams(rand);
    final Trace trace = new Trace();
    int time = 0;
    for (int i = 0; i < 5; i++) {
      trace.add(new PeakResult(time++, 0, 0, 0, 0, 0, 0, params, null));
    }
    time++;
    for (int i = 0; i < 5; i++) {
      trace.add(new PeakResult(time++, 0, 0, 0, 0, 0, 0, params, null));
    }

    runTracing(rand, 0, 2, trace);
  }

  @SeededTest
  public void canTraceMultiplePulseWithMovingCoords(RandomSeed seed) {
    final UniformRandomProvider rand = RngUtils.create(seed.getSeedAsLong());
    final float distance = 0.5f;

    final float[] params = createParams(rand);
    final Trace trace = new Trace();
    int time = 0;
    for (int i = 0; i < 5; i++) {
      move(rand, params, distance);
      trace.add(new PeakResult(time++, 0, 0, 0, 0, 0, 0, params, null));
    }
    time++;
    for (int i = 0; i < 5; i++) {
      move(rand, params, distance);
      trace.add(new PeakResult(time++, 0, 0, 0, 0, 0, 0, params, null));
    }

    runTracing(rand, distance, 2, trace);
  }

  private static void runTracing(final UniformRandomProvider rnd, double distance, int time,
      Trace... expected) {
    final MemoryPeakResults results = toPeakResults(rnd, expected);
    final TraceManager tm = new TraceManager(results);
    final int n = tm.traceMolecules(distance, time);
    Assertions.assertEquals(expected.length, n, "Incorrect number of traces");

    final Trace[] actual = tm.getTraces();
    areEqual(expected, actual);
  }

  @SeededTest
  public void canTraceMultipleFluorophoresWithFixedCoords(RandomSeed seed) {
    simulate(seed, 1000, 1, 5, 2, 0);
  }

  @SeededTest
  public void canTraceMultiplePulsingFluorophoresWithFixedCoords(RandomSeed seed) {
    simulate(seed, 1000, 5, 5, 10, 0);
  }

  @SeededTest
  public void canTraceMultipleFluorophoresWithMovingCoords(RandomSeed seed) {
    simulateMoving(seed, 1000, 1, 5, 2);
  }

  @SeededTest
  public void canTraceMultiplePulsingFluorophoresWithMovingCoords(RandomSeed seed) {
    simulateMoving(seed, 100, 5, 5, 10);
  }

  private static void simulate(RandomSeed seed, int molecules, int maxPulses, int maxOnTime,
      int maxOffTime, float distance) {
    final UniformRandomProvider rand = RngUtils.create(seed.getSeedAsLong());
    final Trace[] expected = new Trace[molecules];
    for (int j = 0; j < expected.length; j++) {
      final float[] params = createParams(rand);
      int time = 1 + rand.nextInt(200);
      final Trace trace = new Trace();
      final int pulses = 1 + rand.nextInt(maxPulses);
      for (int p = 0; p < pulses; p++) {
        final int length = 1 + rand.nextInt(maxOnTime);
        for (int i = 0; i < length; i++) {
          move(rand, params, distance);
          trace.add(new PeakResult(time++, 0, 0, 0, 0, 0, 0, params, null));
        }
        time += 1 + rand.nextInt(maxOffTime);
      }
      expected[j] = trace;
    }

    final double d = (distance > 0) ? Math.sqrt(2.05 * distance * distance) : 0;
    runTracing(rand, d, maxOffTime + 1, expected);
  }

  private static void simulateMoving(RandomSeed seed, int molecules, int maxPulses, int maxOnTime,
      int maxOffTime) {
    final UniformRandomProvider rand = RngUtils.create(seed.getSeedAsLong());

    // When the molecules are moving their paths may intersect.
    // Thus each molecule is allocated a 2x2 square to move within
    // with a 2 pixel border (i.e. a 4x4 square per molecule) ensuring
    // no overlap at the clustering distance of sqrt(2)
    final int n = (int) Math.sqrt(molecules);

    final TIntArrayList list = new TIntArrayList(molecules);
    for (int y = 0; list.size() < molecules; y += 4) {
      for (int x = 0; x < n && list.size() < molecules; x += 4) {
        list.add(y * n + x);
      }
    }
    final int[] positions = list.toArray();
    RandomUtils.shuffle(positions, rand);

    // Offsets for movement around the 2x2 region
    final int[][] offsets = new int[][] {{0, 0}, {0, 1}, {1, 0}, {1, 1},};

    final Trace[] expected = new Trace[molecules];
    for (int j = 0; j < expected.length; j++) {
      final int x = positions[j] % n;
      final int y = positions[j] / n;
      final float[] params = Gaussian2DPeakResultHelper.createOneAxisParams(0, 1, x, y, 0, 1);
      int time = 1 + rand.nextInt(200);
      final Trace trace = new Trace();
      final int pulses = 1 + rand.nextInt(maxPulses);
      for (int p = 0; p < pulses; p++) {
        final int length = 1 + rand.nextInt(maxOnTime);
        for (int i = 0; i < length; i++) {
          final int[] offset = offsets[rand.nextInt(4)];
          params[Gaussian2DFunction.X_POSITION] = x + offset[0];
          params[Gaussian2DFunction.Y_POSITION] = y + offset[1];
          trace.add(new PeakResult(time++, 0, 0, 0, 0, 0, 0, params, null));
        }
        time += 1 + rand.nextInt(maxOffTime);
      }
      expected[j] = trace;
    }

    // distance should be bigger than sqrt(2) but less than 2
    runTracing(rand, 1.5, maxOffTime + 1, expected);
  }

  private static float[] createParams(UniformRandomProvider rand) {
    return Gaussian2DPeakResultHelper.createOneAxisParams(0, 1, rand.nextFloat() * 256f,
        rand.nextFloat() * 256f, 0, 1);
  }

  private static MemoryPeakResults toPeakResults(final UniformRandomProvider rnd, Trace... traces) {
    final PeakResultStoreList results = new ArrayPeakResultStore(traces.length);
    for (final Trace t : traces) {
      results.addStore(t.getPoints());
    }
    // Shuffle
    results.shuffle(rnd);
    return new MemoryPeakResults(results);
  }

  private static void areEqual(Trace[] expected, Trace[] actual) {
    Assertions.assertNotNull(expected);
    Assertions.assertNotNull(actual);
    Assertions.assertEquals(expected.length, actual.length, "Traces are different lengths");

    sort(expected);
    sort(actual);

    for (int i = 0; i < expected.length; i++) {
      final int ii = i;
      final PeakResultStoreList e = expected[i].getPoints();
      final PeakResultStoreList a = actual[i].getPoints();
      Assertions.assertEquals(e.size(), a.size(),
          () -> "Points are different lengths [" + ii + "]");
      for (int j = 0; j < e.size(); j++) {
        final int jj = j;
        final PeakResult p1 = e.get(j);
        final PeakResult p2 = a.get(j);
        Assertions.assertEquals(p1.getFrame(), p2.getFrame(),
            () -> "Frames different [" + ii + "][" + jj + "]");
        Assertions.assertEquals(p1.getXPosition(), p2.getXPosition(), 1e-3f,
            () -> "X different [" + ii + "][" + jj + "]");
        Assertions.assertEquals(p1.getYPosition(), p2.getYPosition(), 1e-3f,
            () -> "Y different [" + ii + "][" + jj + "]");
      }
    }
  }

  private static void sort(Trace[] traces) {
    Arrays.sort(traces, new Comparator<Trace>() {
      @Override
      public int compare(Trace o1, Trace o2) {
        final PeakResult p1 = o1.getHead();
        final PeakResult p2 = o2.getHead();
        final int result = p1.getFrame() - p2.getFrame();
        if (result != 0) {
          return result;
        }
        if (p1.getXPosition() < p2.getXPosition()) {
          return -1;
        }
        if (p1.getXPosition() > p2.getXPosition()) {
          return 1;
        }
        if (p1.getYPosition() < p2.getYPosition()) {
          return -1;
        }
        if (p1.getYPosition() > p2.getYPosition()) {
          return 1;
        }
        return 0;
      }
    });
  }

  private static void move(UniformRandomProvider rand, float[] params, float distance) {
    if (distance > 0) {
      params[PeakResult.X] += -distance + rand.nextDouble() * 2 * distance;
      params[PeakResult.Y] += -distance + rand.nextDouble() * 2 * distance;
    }
  }
}
