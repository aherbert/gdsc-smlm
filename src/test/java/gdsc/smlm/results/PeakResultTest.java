package gdsc.smlm.results;

import org.junit.Assert;
import org.junit.Test;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;

public class PeakResultTest
{
	double[] test_a = { 100, 130, 160 };
	double[] test_s = { 80, 100, 140 };
	double[] test_N = { 1, 10, 30, 100, 1000 };
	double[] test_b2 = { 0, 1, 2, 4, 8 };
	int minPoints = 3, maxPoints = 20;

	@Test
	public void canCalculateMaximumLikelihoodVariance()
	{
		for (double a : test_a)
			for (double s : test_s)
				for (double N : test_N)
					for (double b2 : test_b2)
						for (int points = 3; points <= 20; points++)
						{
							PeakResult.getMLVarianceX(a, s, N, b2, true, points);
						}
	}

	@Test
	public void lowerIntegrationPointsApproximateMaximumLikelihoodVariance()
	{
		double[] sum = new double[maxPoints + 1];
		int count = 0;
		for (double a : test_a)
			for (double s : test_s)
				for (double N : test_N)
					for (double b2 : test_b2)
					{
						count++;
						double e = PeakResult.getMLVarianceX(a, s, N, b2, true, 30);
						for (int points = minPoints; points <= maxPoints; points++)
						{
							double o = PeakResult.getMLVarianceX(a, s, N, b2, true, points);
							double error = DoubleEquality.relativeError(e, o);
							sum[points] += error;
							if (error > 1e-2)
							{
								String msg = String.format("a=%f, s=%f, N=%f, b2=%f, points=%d : %f != %f : %f\n", a,
										s, N, b2, points, e, o, error);
								Assert.assertTrue(msg, false);
							}
						}
					}

		for (int points = minPoints; points <= maxPoints; points++)
		{
			System.out.printf("Points = %d, Av error = %f\n", points, sum[points] / count);
		}
	}

	@Test
	public void runSpeedTest()
	{
		// Test with realistic parameters

		// Warm-up
		for (double a : new double[] { 108 })
			for (double s : new double[] { 120 })
				for (double N : new double[] { 50, 100, 300 })
					for (double b2 : new double[] { 0.5, 1, 2 })
						for (int points = 3; points <= 20; points++)
						{
							PeakResult.getMLVarianceX(a, s, N, b2, true, points);
						}

		// Get average performance
		double[] sum = new double[maxPoints + 1];
		double[] sum2 = new double[sum.length];
		long[] time = new long[sum.length];
		long count = 0, count2 = 0;

		for (double a : new double[] { 108 })
			for (double s : new double[] { 120 })
				for (double N : new double[] { 50, 100, 300 })
					for (double b2 : new double[] { 0.5, 1, 2 })
					{
						long min = Long.MAX_VALUE;
						for (int points = 3; points <= 20; points++)
						{
							long t = System.nanoTime();
							for (int i = 0; i < 1000; i++)
								PeakResult.getMLVarianceX(a, s, N, b2, true, points);
							t = time[points] = System.nanoTime() - t;
							if (min > t)
								min = t;
						}
						// Proportional weighting to the calculation that takes the longest
						count++;
						count2 += min;

						// Store relative performance
						double factor = 1.0 / min;
						for (int points = 3; points <= 20; points++)
						{
							sum[points] += time[points] * factor;
							sum2[points] += time[points];
						}
					}

		for (int points = minPoints; points <= maxPoints; points++)
		{
			System.out.printf("Points = %d, Av relative time = %f, Slow down factor = %f\n", points, sum[points] /
					count, sum2[points] / count2);
		}
	}
}
