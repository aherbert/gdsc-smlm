package gdsc.smlm.results;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.junit.Assert;
import org.junit.Assume;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.PSFProtos.PSFType;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

public class Gaussian2DPeakResultHelperTest
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
							Gaussian2DPeakResultHelper.getMLVarianceX(a, s, N, b2, true, points);
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
						double e = Gaussian2DPeakResultHelper.getMLVarianceX(a, s, N, b2, true, 30);
						for (int points = minPoints; points <= maxPoints; points++)
						{
							double o = Gaussian2DPeakResultHelper.getMLVarianceX(a, s, N, b2, true, points);
							double error = DoubleEquality.relativeError(e, o);
							sum[points] += error;
							if (error > 1e-2)
							{
								String msg = String.format("a=%f, s=%f, N=%f, b2=%f, points=%d : %f != %f : %f\n", a, s,
										N, b2, points, e, o, error);
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
							Gaussian2DPeakResultHelper.getMLVarianceX(a, s, N, b2, true, points);
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
								Gaussian2DPeakResultHelper.getMLVarianceX(a, s, N, b2, true, points);
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
			System.out.printf("Points = %d, Av relative time = %f, Slow down factor = %f\n", points,
					sum[points] / count, sum2[points] / count2);
		}
	}

	@Test
	public void canComputePixelAmplitude()
	{
		float[] x = new float[] { 0f, 0.1f, 0.3f, 0.5f, 0.7f, 1f };
		float[] s = new float[] { 0.8f, 1f, 1.5f, 2.2f };

		float[] paramsf = new float[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		paramsf[Gaussian2DFunction.BACKGROUND] = 0;
		paramsf[Gaussian2DFunction.SIGNAL] = 105;

		Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, 1, 1, GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE,
				null);

		SimpleRegression r = new SimpleRegression(false);

		for (float tx : x)
			for (float ty : x)
				for (float sx : s)
					for (float sy : s)
					{
						paramsf[Gaussian2DFunction.X_POSITION] = tx;
						paramsf[Gaussian2DFunction.Y_POSITION] = ty;
						paramsf[Gaussian2DFunction.X_SD] = sx;
						paramsf[Gaussian2DFunction.Y_SD] = sy;

						// Get the answer using a single pixel image
						// Note the Gaussian2D functions set the centre of the pixel as 0,0 so offset
						double[] params = SimpleArrayUtils.toDouble(paramsf);
						params[Gaussian2DFunction.X_POSITION] -= 0.5;
						params[Gaussian2DFunction.Y_POSITION] -= 0.5;
						f.initialise0(params);
						double e = f.eval(0);

						PSF psf = PSFHelper.create(PSFType.TWO_AXIS_GAUSSIAN_2D);
						CalibrationWriter calibration = new CalibrationWriter();
						calibration.setCountPerPhoton(1);
						calibration.setIntensityUnit(IntensityUnit.PHOTON);
						calibration.setNmPerPixel(1);
						calibration.setDistanceUnit(DistanceUnit.PIXEL);
						Gaussian2DPeakResultCalculator calc = Gaussian2DPeakResultHelper.create(psf, calibration,
								Gaussian2DPeakResultHelper.AMPLITUDE | Gaussian2DPeakResultHelper.PIXEL_AMPLITUDE);
						double o1 = calc.getAmplitude(paramsf);
						double o2 = calc.getPixelAmplitude(paramsf);

						//System.out.printf("e=%f, o1=%f, o2=%f\n", e, o1, o2);
						Assert.assertEquals(e, o2, 1e-3);
						r.addData(e, o1);
					}

		//System.out.printf("Regression: pixel amplitude vs amplitude = %f, slope=%f, n=%d\n", r.getR(), r.getSlope(),
		//		r.getN());
		// The simple amplitude over estimates the actual pixel amplitude
		Assert.assertTrue(r.getSlope() > 1);
	}

	@Test
	public void canComputeCumulative()
	{
		Assert.assertEquals(0.6827, Gaussian2DPeakResultHelper.cumulative(1), 1e-3);
		Assert.assertEquals(0.9545, Gaussian2DPeakResultHelper.cumulative(2), 1e-3);
		Assert.assertEquals(0.9974, Gaussian2DPeakResultHelper.cumulative(3), 1e-3);
	}

	@Test
	public void canComputeSNR()
	{
		//@formatter:off
		Assert.assertEquals(0.6827 / (Math.PI), Gaussian2DPeakResultHelper.getSNR(1, 1, 1, 1, 1), 1e-3);
		Assert.assertEquals(15 * 0.6827 / (Math.PI * 2 * 1.5), Gaussian2DPeakResultHelper.getSNR(15, 2, 1.5, 1, 1), 1e-3);
		Assert.assertEquals(15 * 0.6827 / (Math.PI * 2 * 1.5 * 1.2), Gaussian2DPeakResultHelper.getSNR(15, 2, 1.5, 1, 1.2), 1e-3);
		Assert.assertEquals(15 * 0.9545 / (Math.PI * 2 * 2 * 1.5 * 2 * 1.2), Gaussian2DPeakResultHelper.getSNR(15, 2, 1.5, 2, 1.2), 1e-3);
		//@formatter:on

		// Test fixed versions verse dynamic
		RandomGenerator r = new Well19937c();
		for (int i = 0; i < 10; i++)
		{
			double intensity = r.nextDouble() * 100;
			double sx = r.nextDouble() * 2;
			double sy = r.nextDouble() * 2;
			double noise = r.nextDouble() * 3;
			Assert.assertEquals(Gaussian2DPeakResultHelper.getSNR(intensity, sx, sy, 1, noise),
					Gaussian2DPeakResultHelper.getSNR1(intensity, sx, sy, noise), 1e-3);
			Assert.assertEquals(Gaussian2DPeakResultHelper.getSNR(intensity, sx, sy, 2, noise),
					Gaussian2DPeakResultHelper.getSNR2(intensity, sx, sy, noise), 1e-3);
		}
	}

	@Test
	public void canComputeSNRVersesTotalSNR()
	{
		Assume.assumeTrue(false);

		double intensity = 100;
		double noise = 3;

		for (int i = 0; i <= 10; i++)
		{
			double sx = 1 + i / 10.0;
			for (int j = 0; j <= 10; j++)
			{
				double sy = 1 + j / 10.0;
				System.out.printf("%g,%g  %g  : %g  %g  %g  %g\n", sx, sy, intensity / noise,
						Gaussian2DPeakResultHelper.getSNR(intensity, sx, sy, 1, noise),
						Gaussian2DPeakResultHelper.getSNR(intensity, sx, sy, 1.5, noise),
						Gaussian2DPeakResultHelper.getSNR(intensity, sx, sy, 2, noise),
						Gaussian2DPeakResultHelper.getSNR(intensity, sx, sy, 3, noise));
				;
			}
		}
	}
}
