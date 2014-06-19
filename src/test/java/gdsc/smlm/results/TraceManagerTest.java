package gdsc.smlm.results;

import gdsc.smlm.fitting.function.Gaussian2DFunction;
import gdsc.smlm.utils.Random;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;

import org.junit.Assert;
import org.junit.Test;

public class TraceManagerTest
{
	private gdsc.smlm.utils.Random rand = new Random();

	@Test
	public void canTraceSinglePulseWithFixedCoords()
	{
		float[] params = createParams();
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
		float distance = 0.5f;

		float[] params = createParams();
		Trace trace = new Trace();
		for (int i = 0; i < 5; i++)
		{
			move(params, distance);
			trace.add(new PeakResult(i, 0, 0, 0, 0, 0, params, null));
		}

		runTracing(distance, 1, trace);
	}

	@Test
	public void canTraceMultiplePulseWithFixedCoords()
	{
		int maxOffTime = 5;

		float[] params = createParams();
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
		float distance = 0.5f;
		int maxOffTime = 5;

		float[] params = createParams();
		Trace trace = new Trace();
		for (int i = 0; i < 5; i++)
		{
			move(params, distance);
			trace.add(new PeakResult(i, 0, 0, 0, 0, 0, params, null));
		}
		for (int i = 0; i < 5; i++)
		{
			move(params, distance);
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
		simulate(1000, 1, 5, 2, 0);
	}

	@Test
	public void canTraceMultiplePulsingFluorophoresWithFixedCoords()
	{
		simulate(1000, 5, 5, 10, 0);
	}

	@Test
	public void canTraceMultipleFluorophoresWithMovingCoords()
	{
		// This test can fail if the moving fluorophores paths intersect
		simulate(1000, 1, 5, 2, 0.5f);
	}

	@Test
	public void canTraceMultiplePulsingFluorophoresWithMovingCoords()
	{
		// This test can fail if the moving fluorophores paths intersect
		simulate(100, 5, 5, 10, 0.5f);
	}

	private void simulate(int molecules, int maxPulses, int maxOnTime, int maxOffTime, float distance)
	{
		Trace[] expected = new Trace[molecules];
		for (int j = 0; j < expected.length; j++)
		{
			float[] params = createParams();
			int t = nextInt(200);
			Trace trace = new Trace();
			int pulses = nextInt(maxPulses);
			for (int p = 0; p < pulses; p++)
			{
				int length = nextInt(maxOnTime);
				for (int i = 0; i < length; i++)
				{
					move(params, distance);
					trace.add(new PeakResult(t++, 0, 0, 0, 0, 0, params, null));
				}
				t += nextInt(maxOffTime);
			}
			expected[j] = trace;
		}

		double d = (distance > 0) ? Math.sqrt(2.05 * distance * distance) : 0;
		runTracing(d, maxOffTime + 1, expected);
	}

	private float[] createParams()
	{
		float[] params = new float[7];
		params[Gaussian2DFunction.X_POSITION] = rand.next() * 256f;
		params[Gaussian2DFunction.Y_POSITION] = rand.next() * 256f;
		return params;
	}

	private MemoryPeakResults toPeakResults(Trace... traces)
	{
		MemoryPeakResults results = new MemoryPeakResults();
		for (Trace t : traces)
		{
			results.addAll(t.getPoints());
		}
		Collections.shuffle(results.getResults());
		return results;
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
				Assert.assertEquals("Frames different", p1.peak, p2.peak);
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
				int result = p1.peak - p2.peak;
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

	private void move(float[] params, float distance)
	{
		if (distance > 0)
		{
			params[Gaussian2DFunction.X_POSITION] += -distance + rand.next() * 2 * distance;
			params[Gaussian2DFunction.Y_POSITION] += -distance + rand.next() * 2 * distance;
		}
	}

	/**
	 * Return an integer in the range 1 - max
	 * 
	 * @param max
	 * @return
	 */
	private int nextInt(int max)
	{
		return 1 + (int) (rand.next() * max);
	}
}
