/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
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

package uk.ac.sussex.gdsc.smlm.results;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Locale;
import java.util.function.Consumer;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.logging.TrackProgress;
import uk.ac.sussex.gdsc.core.logging.TrackProgressAdapter;
import uk.ac.sussex.gdsc.core.utils.rng.RandomUtils;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.results.TraceManager.TraceMode;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngFactory;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;

@SuppressWarnings({"javadoc"})
class TraceManagerTest {

  @Test
  void testTraceModeName() {
    for (final TraceMode mode : TraceMode.values()) {
      Assertions.assertEquals(mode.name().toLowerCase(Locale.ROOT).replace('_', ' '),
          mode.toString().toLowerCase(Locale.ROOT).replace('_', ' '));
    }
  }

  @Test
  void testEmptyResults() {
    Assertions.assertThrows(IllegalArgumentException.class, () -> new TraceManager(null));
    final MemoryPeakResults results = new MemoryPeakResults();
    Assertions.assertThrows(IllegalArgumentException.class, () -> new TraceManager(results));
  }

  @Test
  void testInvalidTimeRange() {
    Assertions.assertThrows(IllegalArgumentException.class, () -> new TraceManager(null));
    final MemoryPeakResults results = new MemoryPeakResults();
    final PeakResult peak = new PeakResult(0, 0, 100) {
      @Override
      public int getFrame() {
        return 2;
      }

      @Override
      public int getEndFrame() {
        return 1;
      }
    };
    results.add(peak);
    Assertions.assertThrows(IllegalArgumentException.class, () -> new TraceManager(results));
  }

  @Test
  void testInvalidStartTime() {
    Assertions.assertThrows(IllegalArgumentException.class, () -> new TraceManager(null));
    final MemoryPeakResults results = new MemoryPeakResults();
    final PeakResult peak = new PeakResult(0, 0, 100) {
      @Override
      public int getFrame() {
        return -1;
      }
    };
    results.add(peak);
    Assertions.assertThrows(IllegalArgumentException.class, () -> new TraceManager(results));
  }

  @Test
  void testProperties() {
    final MemoryPeakResults results = new MemoryPeakResults();
    final PeakResult peak = new PeakResult(0, 0, 100);
    results.add(peak);
    final TraceManager tm = new TraceManager(results);
    Assertions.assertEquals(0, tm.getActivationFrameInterval());
    Assertions.assertEquals(0, tm.getActivationFrameWindow());
    Assertions.assertEquals(0, tm.getDistanceExclusion());
    Assertions.assertEquals(0, tm.getPulseInterval());
    Assertions.assertEquals(TraceMode.LATEST_FORERUNNER, tm.getTraceMode());
    tm.setActivationFrameInterval(10);
    tm.setActivationFrameWindow(2);
    tm.setDistanceExclusion(44);
    tm.setPulseInterval(5);
    tm.setTraceMode(TraceMode.SINGLE_LINKAGE);
    Assertions.assertEquals(10, tm.getActivationFrameInterval());
    Assertions.assertEquals(2, tm.getActivationFrameWindow());
    Assertions.assertEquals(44, tm.getDistanceExclusion());
    Assertions.assertEquals(5, tm.getPulseInterval());
    Assertions.assertEquals(TraceMode.SINGLE_LINKAGE, tm.getTraceMode());
    // Hit edge case of deactivating the filtering of activation frames
    tm.setActivationFrameInterval(0);
    tm.setActivationFrameWindow(0);
    Assertions.assertEquals(0, tm.getActivationFrameInterval());
    Assertions.assertEquals(0, tm.getActivationFrameWindow());
  }

  @SeededTest
  void canTraceSinglePulseWithFixedCoords(RandomSeed seed) {
    final UniformRandomProvider rand = RngFactory.create(seed.get());
    final float[] params = createParams(rand);
    final Trace trace = new Trace();
    for (int i = 0; i < 5; i++) {
      trace.add(new PeakResult(i, 0, 0, 0, 0, 0, 0, params, null));
    }

    runTracing(rand, 0, 1, trace);
  }

  @SeededTest
  void canTraceSinglePulseWithMovingCoords(RandomSeed seed) {
    final UniformRandomProvider rand = RngFactory.create(seed.get());
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
  void canTraceMultiplePulseWithFixedCoords(RandomSeed seed) {
    final UniformRandomProvider rand = RngFactory.create(seed.get());

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
  void canTraceMultiplePulseWithMovingCoords(RandomSeed seed) {
    final UniformRandomProvider rand = RngFactory.create(seed.get());
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

  @SeededTest
  void canTraceUsingNoThresholds(RandomSeed seed) {
    final UniformRandomProvider rand = RngFactory.create(seed.get());
    final float intensity = 100;
    final float distance = 0.5f;

    // p1 and p2 cannot be linked
    // p3 can link to either, but p2 is closer
    final PeakResult p1 = new PeakResult(1, -distance, 0, intensity);
    final PeakResult p2 = new PeakResult(2, distance * 0.5f, 0, intensity);
    final PeakResult p3 = new PeakResult(3, 0, 0, intensity);

    // no tracing
    final int time = 1;
    runTracing(rand, 0, time, createTrace(p1), createTrace(p2), createTrace(p3));
    runTracing(rand, -1, time, createTrace(p1), createTrace(p2), createTrace(p3));
    runTracing(rand, distance, 0, createTrace(p1), createTrace(p2), createTrace(p3));
    runTracing(rand, distance, -1, createTrace(p1), createTrace(p2), createTrace(p3));
  }

  @SeededTest
  void canTraceUsingTraceMode(RandomSeed seed) {
    final UniformRandomProvider rand = RngFactory.create(seed.get());
    final float intensity = 100;
    final float distance = 0.5f;

    // p1 and p2 cannot be linked
    // p3 can link to either, but p2 is closer
    final PeakResult p1 = new PeakResult(1, -distance, 0, intensity);
    final PeakResult p2 = new PeakResult(2, distance * 0.5f, 0, intensity);
    final PeakResult p3 = new PeakResult(3, 0, 0, intensity);

    final int time = 2;
    runTracing(rand, distance, time, tm -> tm.setTraceMode(TraceMode.LATEST_FORERUNNER),
        createTrace(p1), createTrace(p2, p3));
    runTracing(rand, distance, time, tm -> tm.setTraceMode(TraceMode.EARLIEST_FORERUNNER),
        createTrace(p1, p3), createTrace(p2));
    runTracing(rand, distance, time - 1, tm -> tm.setTraceMode(TraceMode.EARLIEST_FORERUNNER),
        createTrace(p1), createTrace(p2, p3));
    runTracing(rand, distance, time, tm -> tm.setTraceMode(TraceMode.SINGLE_LINKAGE),
        createTrace(p2, p3), createTrace(p1));
  }

  @SeededTest
  void canTraceUsingTraceMode2(RandomSeed seed) {
    final UniformRandomProvider rand = RngFactory.create(seed.get());
    final float intensity = 100;
    final float distance = 0.5f;

    // p1 and p2 cannot be linked
    // p3 can link to either, but p1 is closer
    final PeakResult p1 = new PeakResult(1, distance * 0.5f, 0, intensity);
    final PeakResult p2 = new PeakResult(2, -distance, 0, intensity);
    final PeakResult p3 = new PeakResult(3, 0, 0, intensity);

    final int time = 2;
    runTracing(rand, distance, time, tm -> tm.setTraceMode(TraceMode.LATEST_FORERUNNER),
        createTrace(p1), createTrace(p2, p3));
    runTracing(rand, distance, time, tm -> tm.setTraceMode(TraceMode.EARLIEST_FORERUNNER),
        createTrace(p1, p3), createTrace(p2));
    runTracing(rand, distance, time - 1, tm -> tm.setTraceMode(TraceMode.EARLIEST_FORERUNNER),
        createTrace(p1), createTrace(p2, p3));
    runTracing(rand, distance, time, tm -> tm.setTraceMode(TraceMode.SINGLE_LINKAGE),
        createTrace(p1, p3), createTrace(p2));
  }

  @SeededTest
  void canTraceUsingTraceMode3(RandomSeed seed) {
    final UniformRandomProvider rand = RngFactory.create(seed.get());
    final float intensity = 100;
    final float distance = 0.5f;

    // Edge case where no points can be linked.
    // This tests the iteration over the existing tracks can terminate without error.
    final PeakResult p1a = new PeakResult(1, 0, 0, intensity);
    final PeakResult p1b = new PeakResult(1, 0, 0, intensity);
    final PeakResult p2a = new PeakResult(2, 2 * distance, 0, intensity);
    final PeakResult p2b = new PeakResult(2, 0, 2 * distance, intensity);
    final PeakResult p3a = new PeakResult(3, 4 * distance, 0, intensity);
    final PeakResult p3b = new PeakResult(3, 0, 4 * distance, intensity);

    final int time = 2;
    final Trace[] traces = new Trace[] {createTrace(p1a), createTrace(p1b), createTrace(p2a),
        createTrace(p2b), createTrace(p3a), createTrace(p3b)};
    runTracing(rand, distance, time, tm -> tm.setTraceMode(TraceMode.LATEST_FORERUNNER), traces);
    runTracing(rand, distance, time, tm -> tm.setTraceMode(TraceMode.EARLIEST_FORERUNNER), traces);
  }

  @SeededTest
  void canTraceUsingExclusionDistance(RandomSeed seed) {
    final UniformRandomProvider rand = RngFactory.create(seed.get());
    final float intensity = 100;
    final float distance = 0.5f;

    // p1 and p2 cannot be linked
    // p3 can link to either, but p2 is closer
    final PeakResult p1 = new PeakResult(1, 0, 0, intensity);
    final PeakResult p2 = new PeakResult(1, distance * 2, 0, intensity);
    final PeakResult p3 = new PeakResult(5, 0, 0, intensity);

    final int time = 10;
    runTracing(rand, distance, time, createTrace(p1, p3), createTrace(p2));
    runTracing(rand, distance, time, tm -> tm.setDistanceExclusion(distance), createTrace(p1, p3),
        createTrace(p2));
    runTracing(rand, distance, time, tm -> tm.setDistanceExclusion(distance * 3), createTrace(p1),
        createTrace(p2), createTrace(p3));
    runTracing(rand, distance, time, tm -> {
      tm.setDistanceExclusion(distance * 3);
      tm.setActivationFrameInterval(3);
      tm.setActivationFrameWindow(1);
    }, createTrace(p1), createTrace(p2), createIgnoredTrace(p3));
  }

  @SeededTest
  void canSplitTracesUsingPulseWindow(RandomSeed seed) {
    final UniformRandomProvider rand = RngFactory.create(seed.get());
    final float intensity = 100;
    final float distance = 0.5f;

    final PeakResult p1 = new PeakResult(1, 0, 0, intensity);
    final PeakResult p2 = new PeakResult(2, 0, 0, intensity);
    final PeakResult p3 = new PeakResult(8, 0, 0, intensity);

    int time = 10;
    runTracing(rand, distance, time, createTrace(p1, p2, p3));
    runTracing(rand, distance, time, tm -> tm.setPulseInterval(5), createTrace(p1, p2),
        createTrace(p3));
    // activation interval but no window
    runTracing(rand, distance, time, tm -> tm.setActivationFrameInterval(5),
        createTrace(p1, p2, p3));

    // Two traces
    time = 5;
    runTracing(rand, distance, time, createTrace(p1, p2), createTrace(p3));
    // Second trace start is not in the activation window
    runTracing(rand, distance, time, tm -> {
      tm.setActivationFrameInterval(5);
      tm.setActivationFrameWindow(2);
    }, createTrace(p1, p2), createIgnoredTrace(p3));

    // No tracing but filtering is active
    final MemoryPeakResults results = new MemoryPeakResults();
    results.addAll(Arrays.asList(p1, p2, p3));
    final TraceManager tm = new TraceManager(results);
    tm.setActivationFrameInterval(5);
    tm.setActivationFrameWindow(2);
    areEqual(new Trace[] {createTrace(p1), createTrace(p2)}, tm.getTraces());
  }

  @SeededTest
  void canTraceToExistingTrajectories(RandomSeed seed) {
    final UniformRandomProvider rand = RngFactory.create(seed.get());
    final float intensity = 100;
    final float distance = 0.5f;

    final Trace trace1 = new Trace();
    trace1.add(new PeakResult(1, 0, 0, intensity));
    trace1.add(new PeakResult(2, distance, 0, intensity));

    // Spot 3 is within the trace distance of spot 1. But spot 1 should
    // be joined in a trajectory to spot 2. Spot 3 is not within the distance
    // to spot 2.
    final Trace trace2 = new Trace();
    trace2.add(new PeakResult(3, 0, distance, intensity));

    final int time = 2;
    runTracing(rand, distance, time, trace1, trace2);
  }

  @SeededTest
  void canTraceUsingTracker(RandomSeed seed) {
    final UniformRandomProvider rand = RngFactory.create(seed.get());
    final float intensity = 100;
    final float distance = 0.5f;

    // p1 and p2 can be linked
    final PeakResult p1 = new PeakResult(1, distance, 0, intensity);
    final PeakResult p2 = new PeakResult(2, 0, 0, intensity);

    final int[] calls = {0, 0};
    final TrackProgress tracker = new TrackProgressAdapter() {
      @Override
      public void progress(double fraction) {
        calls[0]++;
      }

      @Override
      public void progress(long position, long total) {
        calls[1]++;
      }
    };

    final int time = 1;
    runTracing(rand, distance, time, tm -> {
      Assertions.assertNull(tm.getTracker());
      tm.setTracker(tracker);
      Assertions.assertSame(tracker, tm.getTracker());
    }, createTrace(p1, p2));

    Assertions.assertNotEquals(0, calls[0]);
    Assertions.assertNotEquals(0, calls[1]);
  }

  @SeededTest
  void canTraceUsingEndFrame(RandomSeed seed) {
    final UniformRandomProvider rand = RngFactory.create(seed.get());
    final float intensity = 100;
    final float distance = 0.5f;

    // p1 and p2 cannot be linked as p1 ends after p2 starts.
    // p3 can link to p2, not p1 as p1 ends after p3 starts.
    // p4 can link to p3 or p1, p1 is closer
    final PeakResult p1 = createWithEndFrame(1, 3, 0, 0);
    final PeakResult p2 = new PeakResult(2, distance, 0, intensity);
    final PeakResult p3 = new PeakResult(3, distance * 0.5f, 0, intensity);
    final PeakResult p4 = new PeakResult(4, 0, 0, intensity);

    final int time = 10;
    runTracing(rand, distance, time, createTrace(p1, p4), createTrace(p2, p3));
  }

  private static PeakResult createWithEndFrame(int frame, int endFrame, float x, float y) {
    final ExtendedPeakResult p = new ExtendedPeakResult(frame, x, y, 100, 0);
    p.setEndFrame(endFrame);
    return p;
  }

  @SeededTest
  void canTraceUsingSingleFrame(RandomSeed seed) {
    final UniformRandomProvider rand = RngFactory.create(seed.get());
    final float intensity = 100;
    final float distance = 0.5f;

    // p1 and p2 cannot be linked as p1 ends after p2 starts.
    // p3 can link to p2, not p1 as p1 ends after p3 starts.
    // p4 can link to p3 or p1, p1 is closer
    final PeakResult p1 = new PeakResult(2, 0, 0, intensity);
    final PeakResult p2 = new PeakResult(2, distance, 0, intensity);
    final PeakResult p3 = new PeakResult(2, distance * 0.5f, 0, intensity);

    final int time = 10;
    runTracing(rand, distance, time, createTrace(p1));
    runTracing(rand, distance, time, createTrace(p1), createTrace(p2));
    runTracing(rand, distance, time, createTrace(p1), createTrace(p2), createTrace(p3));
  }

  private static void runTracing(final UniformRandomProvider rnd, double distance, int time,
      Trace... expected) {
    runTracing(rnd, distance, time, null, expected);
  }

  private static void runTracing(final UniformRandomProvider rnd, double distance, int time,
      Consumer<TraceManager> setupAction, Trace... traces) {
    final MemoryPeakResults results = toPeakResults(rnd, traces);
    final TraceManager tm = new TraceManager(results);
    if (setupAction != null) {
      setupAction.accept(tm);
    }
    final int[] filtered = {0};
    final Trace[] expected = Arrays.stream(traces).filter(t -> {
      if (t.getId() >= 0) {
        return true;
      }
      filtered[0] += t.size();
      return false;
    }).toArray(Trace[]::new);
    final int n = tm.traceMolecules(distance, time);
    Assertions.assertEquals(expected.length, n, "Incorrect number of traces");

    final Trace[] actual = tm.getTraces();
    areEqual(expected, actual);
    Assertions.assertEquals(actual.length, tm.getTotalTraces());
    Assertions.assertEquals(filtered[0], tm.getTotalFiltered());
  }

  @SeededTest
  void canTraceMultipleFluorophoresWithFixedCoords(RandomSeed seed) {
    simulate(seed, 1000, 1, 5, 2, 0);
  }

  @SeededTest
  void canTraceMultiplePulsingFluorophoresWithFixedCoords(RandomSeed seed) {
    simulate(seed, 1000, 5, 5, 10, 0);
  }

  @SeededTest
  void canTraceMultipleFluorophoresWithMovingCoords(RandomSeed seed) {
    simulateMoving(seed, 1000, 1, 5, 2);
  }

  @SeededTest
  void canTraceMultiplePulsingFluorophoresWithMovingCoords(RandomSeed seed) {
    simulateMoving(seed, 100, 5, 5, 10);
  }

  private static void simulate(RandomSeed seed, int molecules, int maxPulses, int maxOnTime,
      int maxOffTime, float distance) {
    final UniformRandomProvider rand = RngFactory.create(seed.get());
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
    final UniformRandomProvider rand = RngFactory.create(seed.get());

    // When the molecules are moving their paths may intersect.
    // Thus each molecule is allocated a 2x2 square to move within
    // with a 2 pixel border (i.e. a 4x4 square per molecule) ensuring
    // no overlap at the clustering distance of sqrt(2)
    final int n = (int) Math.sqrt(molecules);

    final IntArrayList list = new IntArrayList(molecules);
    for (int y = 0; list.size() < molecules; y += 4) {
      for (int x = 0; x < n && list.size() < molecules; x += 4) {
        list.add(y * n + x);
      }
    }
    final int[] positions = list.toIntArray();
    RandomUtils.shuffle(positions, rand);

    // Offsets for movement around the 2x2 region
    final int[][] offsets = new int[][] {{0, 0}, {0, 1}, {1, 0}, {1, 1}};

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

  private static Trace createTrace(PeakResult... peakResults) {
    final Trace trace = new Trace();
    for (final PeakResult r : peakResults) {
      trace.add(r);
    }
    return trace;
  }

  private static Trace createIgnoredTrace(PeakResult... peakResults) {
    final Trace trace = createTrace(peakResults);
    trace.setId(-1);
    return trace;
  }

  @Test
  void canFilterTraces() {
    final PeakResult p1 = new PeakResult(1, 0, 0, 100);
    final PeakResult p2 = new PeakResult(2, 0, 0, 100);
    final PeakResult p3 = new PeakResult(3, 0, 0, 100);
    final Trace[] traces =
        new Trace[] {createTrace(p1), createTrace(p2), createTrace(p3), createTrace()};

    final MemoryPeakResults results = new MemoryPeakResults();
    results.add(new PeakResult(44, 0, 0, 100));
    final TraceManager tm = new TraceManager(results);
    areEqual(new Trace[] {createTrace(p1)}, tm.filterTraces(traces, 10));
    areEqual(new Trace[] {createTrace(p1), createTrace(p3)}, tm.filterTraces(traces, 2));
  }
}
