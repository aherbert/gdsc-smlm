package gdsc.smlm.results;

import gdsc.core.utils.Random;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;

import org.junit.Assert;
import org.junit.Test;

public class TraceManagerTest
{
	@Test
	public void canTraceSinglePulseWithFixedCoords()
	{
		Random rand = new Random(30051977);
		float[] params = createParams(rand);
		Trace trace = new Trace();
		for (int i = 0; i < 5; i++)
		{
			trace.add(new PeakResult(i, 0, 0, 0, 0, 0, params, null));
		}

		runTracing(0, 1, trace);
	}

	@Test
	public void canTraceSinglePulseWithMovingCoords()
	{
		Random rand = new Random(30051977);
		float distance = 0.5f;

		float[] params = createParams(rand);
		Trace trace = new Trace();
		for (int i = 0; i < 5; i++)
		{
			move(rand, params, distance);
			trace.add(new PeakResult(i, 0, 0, 0, 0, 0, params, null));
		}

		runTracing(distance, 1, trace);
	}

	@Test
	public void canTraceMultiplePulseWithFixedCoords()
	{
		Random rand = new Random(30051977);
		int maxOffTime = 5;

		float[] params = createParams(rand);
		Trace trace = new Trace();
		for (int i = 0; i < 5; i++)
		{
			trace.add(new PeakResult(i, 0, 0, 0, 0, 0, params, null));
		}
		for (int i = 0; i < 5; i++)
		{
			trace.add(new PeakResult(i + maxOffTime, 0, 0, 0, 0, 0, params, null));
		}

		runTracing(0, maxOffTime + 1, trace);
	}

	@Test
	public void canTraceMultiplePulseWithMovingCoords()
	{
		Random rand = new Random(30051977);
		float distance = 0.5f;
		int maxOffTime = 5;

		float[] params = createParams(rand);
		Trace trace = new Trace();
		for (int i = 0; i < 5; i++)
		{
			move(rand, params, distance);
			trace.add(new PeakResult(i, 0, 0, 0, 0, 0, params, null));
		}
		for (int i = 0; i < 5; i++)
		{
			move(rand, params, distance);
			trace.add(new PeakResult(i + maxOffTime, 0, 0, 0, 0, 0, params, null));
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
		simulate(new Random(), 1000, 1, 5, 2, 0);
	}

	@Test
	public void canTraceMultiplePulsingFluorophoresWithFixedCoords()
	{
		simulate(new Random(), 1000, 5, 5, 10, 0);
	}

	@Test
	public void canTraceMultipleFluorophoresWithMovingCoords()
	{
		// This test can fail if the moving fluorophores paths intersect
		// so we used a known seed that is ok
		simulate(new Random(3), 1000, 1, 5, 2, 0.5f);
	}

	@Test
	public void canTraceMultiplePulsingFluorophoresWithMovingCoords()
	{
		// This test can fail if the moving fluorophores paths intersect
		// so we used a known seed that is ok
		simulate(new Random(3005), 100, 5, 5, 10, 0.5f);
	}

	private void simulate(Random rand, int molecules, int maxPulses, int maxOnTime, int maxOffTime, float distance)
	{
		Trace[] expected = new Trace[molecules];
		for (int j = 0; j < expected.length; j++)
		{
			float[] params = createParams(rand);
			int t = rand.nextInt(1, 200);
			Trace trace = new Trace();
			int pulses = rand.nextInt(1, maxPulses);
			for (int p = 0; p < pulses; p++)
			{
				int length = rand.nextInt(1, maxOnTime);
				for (int i = 0; i < length; i++)
				{
					move(rand, params, distance);
					trace.add(new PeakResult(t++, 0, 0, 0, 0, 0, params, null));
				}
				t += rand.nextInt(1, maxOffTime);
			}
			expected[j] = trace;
		}

		double d = (distance > 0) ? Math.sqrt(2.05 * distance * distance) : 0;
		runTracing(d, maxOffTime + 1, expected);
	}

	private static float[] createParams(Random rand)
	{
		float[] params = new float[7];
		params[Gaussian2DFunction.X_POSITION] = rand.next() * 256f;
		params[Gaussian2DFunction.Y_POSITION] = rand.next() * 256f;
		return params;
	}

	private MemoryPeakResults toPeakResults(Trace... traces)
	{
		ArrayList<PeakResult> results= new ArrayList<PeakResult>();
		for (Trace t : traces)
		{
			results.addAll(t.getPoints());
		}
		Collections.shuffle(results);
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
			ArrayList<PeakResult> e = expected[i].getPoints();
			ArrayList<PeakResult> a = actual[i].getPoints();
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

	private static void move(Random rand, float[] params, float distance)
	{
		if (distance > 0)
		{
			params[Gaussian2DFunction.X_POSITION] += -distance + rand.next() * 2 * distance;
			params[Gaussian2DFunction.Y_POSITION] += -distance + rand.next() * 2 * distance;
		}
	}
}
