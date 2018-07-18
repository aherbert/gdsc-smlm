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
package uk.ac.sussex.gdsc.smlm.results;

import java.util.Arrays;
import java.util.Comparator;

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import gnu.trove.list.array.TIntArrayList;
import uk.ac.sussex.gdsc.core.utils.Random;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.test.TestSettings;

@SuppressWarnings({ "javadoc" })
public class TraceManagerTest
{
	@Test
	public void canTraceSinglePulseWithFixedCoords()
	{
		final RandomGenerator rand = TestSettings.getRandomGenerator();
		final float[] params = createParams(rand);
		final Trace trace = new Trace();
		for (int i = 0; i < 5; i++)
			trace.add(new PeakResult(i, 0, 0, 0, 0, 0, 0, params, null));

		runTracing(0, 1, trace);
	}

	@Test
	public void canTraceSinglePulseWithMovingCoords()
	{
		final RandomGenerator rand = TestSettings.getRandomGenerator();
		final float distance = 0.5f;

		final float[] params = createParams(rand);
		final Trace trace = new Trace();
		for (int i = 0; i < 5; i++)
		{
			move(rand, params, distance);
			trace.add(new PeakResult(i, 0, 0, 0, 0, 0, 0, params, null));
		}

		runTracing(distance, 1, trace);
	}

	@Test
	public void canTraceMultiplePulseWithFixedCoords()
	{
		final RandomGenerator rand = TestSettings.getRandomGenerator();

		final float[] params = createParams(rand);
		final Trace trace = new Trace();
		int t = 0;
		for (int i = 0; i < 5; i++)
			trace.add(new PeakResult(t++, 0, 0, 0, 0, 0, 0, params, null));
		t++;
		for (int i = 0; i < 5; i++)
			trace.add(new PeakResult(t++, 0, 0, 0, 0, 0, 0, params, null));

		runTracing(0, 2, trace);
	}

	@Test
	public void canTraceMultiplePulseWithMovingCoords()
	{
		final RandomGenerator rand = TestSettings.getRandomGenerator();
		final float distance = 0.5f;

		final float[] params = createParams(rand);
		final Trace trace = new Trace();
		int t = 0;
		for (int i = 0; i < 5; i++)
		{
			move(rand, params, distance);
			trace.add(new PeakResult(t++, 0, 0, 0, 0, 0, 0, params, null));
		}
		t++;
		for (int i = 0; i < 5; i++)
		{
			move(rand, params, distance);
			trace.add(new PeakResult(t++, 0, 0, 0, 0, 0, 0, params, null));
		}

		runTracing(distance, 2, trace);
	}

	private static void runTracing(double d, int t, Trace... expected)
	{
		final MemoryPeakResults results = toPeakResults(expected);
		final TraceManager tm = new TraceManager(results);
		final int n = tm.traceMolecules(d, t);
		Assert.assertEquals("Incorrect number of traces", expected.length, n);

		final Trace[] actual = tm.getTraces();
		areEqual(expected, actual);
	}

	@Test
	public void canTraceMultipleFluorophoresWithFixedCoords()
	{
		simulate(TestSettings.getRandomGenerator(), 1000, 1, 5, 2, 0);
	}

	@Test
	public void canTraceMultiplePulsingFluorophoresWithFixedCoords()
	{
		simulate(TestSettings.getRandomGenerator(), 1000, 5, 5, 10, 0);
	}

	@Test
	public void canTraceMultipleFluorophoresWithMovingCoords()
	{
		simulateMoving(TestSettings.getRandomGenerator(), 1000, 1, 5, 2);
	}

	@Test
	public void canTraceMultiplePulsingFluorophoresWithMovingCoords()
	{
		simulateMoving(TestSettings.getRandomGenerator(), 100, 5, 5, 10);
	}

	private static void simulate(RandomGenerator rand, int molecules, int maxPulses, int maxOnTime, int maxOffTime,
			float distance)
	{
		final Trace[] expected = new Trace[molecules];
		for (int j = 0; j < expected.length; j++)
		{
			final float[] params = createParams(rand);
			int t = 1 + rand.nextInt(200);
			final Trace trace = new Trace();
			final int pulses = 1 + rand.nextInt(maxPulses);
			for (int p = 0; p < pulses; p++)
			{
				final int length = 1 + rand.nextInt(maxOnTime);
				for (int i = 0; i < length; i++)
				{
					move(rand, params, distance);
					trace.add(new PeakResult(t++, 0, 0, 0, 0, 0, 0, params, null));
				}
				t += 1 + rand.nextInt(maxOffTime);
			}
			expected[j] = trace;
		}

		final double d = (distance > 0) ? Math.sqrt(2.05 * distance * distance) : 0;
		runTracing(d, maxOffTime + 1, expected);
	}

	private static void simulateMoving(RandomGenerator rand, int molecules, int maxPulses, int maxOnTime, int maxOffTime)
	{
		// When the molecules are moving their paths may intersect.
		// Thus each molecule is allocated a 2x2 square to move within
		// with a 2 pixel border (i.e. a 4x4 square per molecule) ensuring
		// no overlap at the clustering distance of sqrt(2)
		final int n = (int) Math.sqrt(molecules);

		final TIntArrayList list = new TIntArrayList(molecules);
		for (int y = 0; list.size() < molecules; y += 4)
			for (int x = 0; x < n && list.size() < molecules; x += 4)
				list.add(y * n + x);
		final int[] positions = list.toArray();
		Random.shuffle(positions, rand);

		// Offsets for movement around the 2x2 region
		final int[][] offsets = new int[][] { { 0, 0 }, { 0, 1 }, { 1, 0 }, { 1, 1 }, };

		final Trace[] expected = new Trace[molecules];
		for (int j = 0; j < expected.length; j++)
		{
			final int x = positions[j] % n;
			final int y = positions[j] / n;
			final float[] params = Gaussian2DPeakResultHelper.createOneAxisParams(0, 1, x, y, 0, 1);
			int t = 1 + rand.nextInt(200);
			final Trace trace = new Trace();
			final int pulses = 1 + rand.nextInt(maxPulses);
			for (int p = 0; p < pulses; p++)
			{
				final int length = 1 + rand.nextInt(maxOnTime);
				for (int i = 0; i < length; i++)
				{
					final int[] offset = offsets[rand.nextInt(4)];
					params[Gaussian2DFunction.X_POSITION] = x + offset[0];
					params[Gaussian2DFunction.Y_POSITION] = y + offset[1];
					trace.add(new PeakResult(t++, 0, 0, 0, 0, 0, 0, params, null));
				}
				t += 1 + rand.nextInt(maxOffTime);
			}
			expected[j] = trace;
		}

		// distance should be bigger than sqrt(2) but less than 2
		runTracing(1.5, maxOffTime + 1, expected);
	}

	private static float[] createParams(RandomGenerator rand)
	{
		return Gaussian2DPeakResultHelper.createOneAxisParams(0, 1, rand.nextFloat() * 256f, rand.nextFloat() * 256f, 0,
				1);
	}

	private static MemoryPeakResults toPeakResults(Trace... traces)
	{
		final PeakResultStoreList results = new ArrayPeakResultStore(traces.length);
		for (final Trace t : traces)
			results.addStore(t.getPoints());
		// Shuffle
		final RandomGenerator rnd = TestSettings.getRandomGenerator();
		results.shuffle(rnd);
		return new MemoryPeakResults(results);
	}

	private static void areEqual(Trace[] expected, Trace[] actual)
	{
		Assert.assertNotNull(expected);
		Assert.assertNotNull(actual);
		Assert.assertEquals("Traces are different lengths", expected.length, actual.length);

		sort(expected);
		sort(actual);

		for (int i = 0; i < expected.length; i++)
		{
			final PeakResultStoreList e = expected[i].getPoints();
			final PeakResultStoreList a = actual[i].getPoints();
			Assert.assertEquals("Points are different lengths [" + i + "]", e.size(), a.size());
			for (int j = 0; j < e.size(); j++)
			{
				final PeakResult p1 = e.get(j);
				final PeakResult p2 = a.get(j);
				Assert.assertEquals("Frames different", p1.getFrame(), p2.getFrame());
				Assert.assertEquals("X different", p1.getXPosition(), p2.getXPosition(), 1e-3f);
				Assert.assertEquals("Y different", p1.getYPosition(), p2.getYPosition(), 1e-3f);
			}
		}
	}

	private static void sort(Trace[] traces)
	{
		Arrays.sort(traces, new Comparator<Trace>()
		{
			@Override
			public int compare(Trace o1, Trace o2)
			{
				final PeakResult p1 = o1.getHead();
				final PeakResult p2 = o2.getHead();
				final int result = p1.getFrame() - p2.getFrame();
				if (result != 0)
					return result;
				if (p1.getXPosition() < p2.getXPosition())
					return -1;
				if (p1.getXPosition() > p2.getXPosition())
					return 1;
				if (p1.getYPosition() < p2.getYPosition())
					return -1;
				if (p1.getYPosition() > p2.getYPosition())
					return 1;
				return 0;
			}
		});
	}

	private static void move(RandomGenerator rand, float[] params, float distance)
	{
		if (distance > 0)
		{
			params[PeakResult.X] += -distance + rand.nextDouble() * 2 * distance;
			params[PeakResult.Y] += -distance + rand.nextDouble() * 2 * distance;
		}
	}
}