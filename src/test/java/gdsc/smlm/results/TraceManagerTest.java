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
package gdsc.smlm.results;

import java.util.Arrays;
import java.util.Comparator;

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import gdsc.test.TestSettings;

public class TraceManagerTest
{
	@Test
	public void canTraceSinglePulseWithFixedCoords()
	{
		RandomGenerator rand = TestSettings.getRandomGenerator();
		float[] params = createParams(rand);
		Trace trace = new Trace();
		for (int i = 0; i < 5; i++)
		{
			trace.add(new PeakResult(i, 0, 0, 0, 0, 0, 0, params, null));
		}

		runTracing(0, 1, trace);
	}

	@Test
	public void canTraceSinglePulseWithMovingCoords()
	{
		RandomGenerator rand = TestSettings.getRandomGenerator();
		float distance = 0.5f;

		float[] params = createParams(rand);
		Trace trace = new Trace();
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
		RandomGenerator rand = TestSettings.getRandomGenerator();
		int maxOffTime = 5;

		float[] params = createParams(rand);
		Trace trace = new Trace();
		for (int i = 0; i < 5; i++)
		{
			trace.add(new PeakResult(i, 0, 0, 0, 0, 0, 0, params, null));
		}
		for (int i = 0; i < 5; i++)
		{
			trace.add(new PeakResult(i + maxOffTime, 0, 0, 0, 0, 0, 0, params, null));
		}

		runTracing(0, maxOffTime + 1, trace);
	}

	@Test
	public void canTraceMultiplePulseWithMovingCoords()
	{
		RandomGenerator rand = TestSettings.getRandomGenerator();
		float distance = 0.5f;
		int maxOffTime = 5;

		float[] params = createParams(rand);
		Trace trace = new Trace();
		for (int i = 0; i < 5; i++)
		{
			move(rand, params, distance);
			trace.add(new PeakResult(i, 0, 0, 0, 0, 0, 0, params, null));
		}
		for (int i = 0; i < 5; i++)
		{
			move(rand, params, distance);
			trace.add(new PeakResult(i + maxOffTime, 0, 0, 0, 0, 0, 0, params, null));
		}

		runTracing(distance, maxOffTime + 1, trace);
	}

	private void runTracing(double d, int t, Trace... expected)
	{
		MemoryPeakResults results = toPeakResults(expected);
		TraceManager tm = new TraceManager(results);
		int n = tm.traceMolecules(d, t);
		Assert.assertEquals("Incorrect number of traces", expected.length, n);

		Trace[] actual = tm.getTraces();
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
		// This test can fail if the moving fluorophores paths intersect
		// so we use a known seed that is ok
		simulate(TestSettings.getRandomGenerator(3), 1000, 1, 5, 2, 0.5f);
	}

	@Test
	public void canTraceMultiplePulsingFluorophoresWithMovingCoords()
	{
		// This test can fail if the moving fluorophores paths intersect
		// so we use a known seed that is ok
		simulate(TestSettings.getRandomGenerator(3005), 100, 5, 5, 10, 0.5f);
	}

	private void simulate(RandomGenerator rand, int molecules, int maxPulses, int maxOnTime, int maxOffTime,
			float distance)
	{
		Trace[] expected = new Trace[molecules];
		for (int j = 0; j < expected.length; j++)
		{
			float[] params = createParams(rand);
			int t = 1 + rand.nextInt(200);
			Trace trace = new Trace();
			int pulses = 1 + rand.nextInt(maxPulses);
			for (int p = 0; p < pulses; p++)
			{
				int length = 1 + rand.nextInt(maxOnTime);
				for (int i = 0; i < length; i++)
				{
					move(rand, params, distance);
					trace.add(new PeakResult(t++, 0, 0, 0, 0, 0, 0, params, null));
				}
				t += 1 + rand.nextInt(maxOffTime);
			}
			expected[j] = trace;
		}

		double d = (distance > 0) ? Math.sqrt(2.05 * distance * distance) : 0;
		runTracing(d, maxOffTime + 1, expected);
	}

	private static float[] createParams(RandomGenerator rand)
	{
		return Gaussian2DPeakResultHelper.createOneAxisParams(0, 1, rand.nextFloat() * 256f, rand.nextFloat() * 256f, 0,
				1);
	}

	private MemoryPeakResults toPeakResults(Trace... traces)
	{
		PeakResultStoreList results = new ArrayPeakResultStore(traces.length);
		for (Trace t : traces)
		{
			results.addStore(t.getPoints());
		}
		// Shuffle
		RandomGenerator rnd = TestSettings.getRandomGenerator();
		results.shuffle(rnd);
		return new MemoryPeakResults(results);
	}

	private void areEqual(Trace[] expected, Trace[] actual)
	{
		Assert.assertNotNull(expected);
		Assert.assertNotNull(actual);
		Assert.assertEquals("Traces are different lengths", expected.length, actual.length);

		sort(expected);
		sort(actual);

		for (int i = 0; i < expected.length; i++)
		{
			PeakResultStoreList e = expected[i].getPoints();
			PeakResultStoreList a = actual[i].getPoints();
			Assert.assertEquals("Points are different lengths [" + i + "]", e.size(), a.size());
			for (int j = 0; j < e.size(); j++)
			{
				PeakResult p1 = e.get(j);
				PeakResult p2 = a.get(j);
				Assert.assertEquals("Frames different", p1.getFrame(), p2.getFrame());
				Assert.assertEquals("X different", p1.getXPosition(), p2.getXPosition(), 1e-3f);
				Assert.assertEquals("Y different", p1.getYPosition(), p2.getYPosition(), 1e-3f);
			}
		}
	}

	private void sort(Trace[] traces)
	{
		Arrays.sort(traces, new Comparator<Trace>()
		{
			@Override
			public int compare(Trace o1, Trace o2)
			{
				PeakResult p1 = o1.getHead();
				PeakResult p2 = o2.getHead();
				int result = p1.getFrame() - p2.getFrame();
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
