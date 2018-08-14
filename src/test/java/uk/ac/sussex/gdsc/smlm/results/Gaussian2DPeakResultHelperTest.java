package uk.ac.sussex.gdsc.smlm.results;

import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.Maths;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.PSFHelper;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.test.TestComplexity;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.functions.FunctionUtils;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssumptions;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;

@SuppressWarnings({ "javadoc" })
public class Gaussian2DPeakResultHelperTest
{
    private static Logger logger;

    @BeforeAll
    public static void beforeAll()
    {
        logger = Logger.getLogger(Gaussian2DPeakResultHelperTest.class.getName());
    }

    @AfterAll
    public static void afterAll()
    {
        logger = null;
    }

    double[] test_a = { 100, 130, 160 };
    double[] test_s = { 80, 100, 140 };
    double[] test_N = { 1, 10, 30, 100, 1000 };
    double[] test_b2 = { 0, 1, 2, 4, 8 };
    int minPoints = 3, maxPoints = 20;

    @Test
    public void canCalculateMaximumLikelihoodVariance()
    {
        int min = Gaussian2DPeakResultHelper.POINTS;
        int max = min;
        if (TestSettings.allow(TestComplexity.HIGH))
        {
            min = 3;
            max = 20;
        }
        for (final double a : test_a)
            for (final double s : test_s)
                for (final double N : test_N)
                    for (final double b2 : test_b2)
                        for (int points = min; points <= max; points++)
                            Gaussian2DPeakResultHelper.getMLVarianceX(a, s, N, b2, true, points);
    }

    @Test
    public void lowerIntegrationPointsApproximateMaximumLikelihoodVariance()
    {
        ExtraAssumptions.assume(logger, Level.INFO);
        ExtraAssumptions.assume(TestComplexity.HIGH);

        final double[] sum = new double[maxPoints + 1];
        int count = 0;
        for (final double a : test_a)
            for (final double s : test_s)
                for (final double N : test_N)
                    for (final double b2 : test_b2)
                    {
                        count++;
                        final double e = Gaussian2DPeakResultHelper.getMLVarianceX(a, s, N, b2, true, 30);
                        for (int points = minPoints; points <= maxPoints; points++)
                        {
                            final double o = Gaussian2DPeakResultHelper.getMLVarianceX(a, s, N, b2, true, points);
                            final double error = DoubleEquality.relativeError(e, o);
                            sum[points] += error;
                            if (error > 1e-2)
                                Assertions.fail(String.format("a=%f, s=%f, N=%f, b2=%f, points=%d : %f != %f : %f", a, s, N, b2,
                                        points, e, o, error));
                        }
                    }

        for (int points = minPoints; points <= maxPoints; points++)
            logger.info(FunctionUtils.getSupplier("Points = %d, Av error = %f", points, sum[points] / count));
    }

    @Test
    public void runSpeedTest()
    {
        ExtraAssumptions.assume(logger, Level.INFO);
        ExtraAssumptions.assume(TestComplexity.HIGH);

        // Test with realistic parameters

        // Warm-up
        for (final double a : new double[] { 108 })
            for (final double s : new double[] { 120 })
                for (final double N : new double[] { 50, 100, 300 })
                    for (final double b2 : new double[] { 0.5, 1, 2 })
                        for (int points = 3; points <= 20; points++)
                            Gaussian2DPeakResultHelper.getMLVarianceX(a, s, N, b2, true, points);

        // Get average performance
        final double[] sum = new double[maxPoints + 1];
        final double[] sum2 = new double[sum.length];
        final long[] time = new long[sum.length];
        long count = 0, count2 = 0;

        for (final double a : new double[] { 108 })
            for (final double s : new double[] { 120 })
                for (final double N : new double[] { 50, 100, 300 })
                    for (final double b2 : new double[] { 0.5, 1, 2 })
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
                        final double factor = 1.0 / min;
                        for (int points = 3; points <= 20; points++)
                        {
                            sum[points] += time[points] * factor;
                            sum2[points] += time[points];
                        }
                    }

        for (int points = minPoints; points <= maxPoints; points++)
            logger.info(FunctionUtils.getSupplier("Points = %d, Av relative time = %f, Slow down factor = %f", points,
                    sum[points] / count, sum2[points] / count2));
    }

    @Test
    public void canComputePixelAmplitude()
    {
        final float[] x = new float[] { 0f, 0.1f, 0.3f, 0.5f, 0.7f, 1f };
        final float[] s = new float[] { 0.8f, 1f, 1.5f, 2.2f };

        final float[] paramsf = new float[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
        paramsf[Gaussian2DFunction.BACKGROUND] = 0;
        paramsf[Gaussian2DFunction.SIGNAL] = 105;

        final Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, 1, 1,
                GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);

        final SimpleRegression r = new SimpleRegression(false);

        for (final float tx : x)
            for (final float ty : x)
                for (final float sx : s)
                    for (final float sy : s)
                    {
                        paramsf[Gaussian2DFunction.X_POSITION] = tx;
                        paramsf[Gaussian2DFunction.Y_POSITION] = ty;
                        paramsf[Gaussian2DFunction.X_SD] = sx;
                        paramsf[Gaussian2DFunction.Y_SD] = sy;

                        // Get the answer using a single pixel image
                        // Note the Gaussian2D functions set the centre of the pixel as 0,0 so offset
                        final double[] params = SimpleArrayUtils.toDouble(paramsf);
                        params[Gaussian2DFunction.X_POSITION] -= 0.5;
                        params[Gaussian2DFunction.Y_POSITION] -= 0.5;
                        f.initialise0(params);
                        final double e = f.eval(0);

                        final PSF psf = PSFHelper.create(PSFType.TWO_AXIS_GAUSSIAN_2D);
                        final CalibrationWriter calibration = new CalibrationWriter();
                        calibration.setCountPerPhoton(1);
                        calibration.setIntensityUnit(IntensityUnit.PHOTON);
                        calibration.setNmPerPixel(1);
                        calibration.setDistanceUnit(DistanceUnit.PIXEL);
                        final Gaussian2DPeakResultCalculator calc = Gaussian2DPeakResultHelper.create(psf, calibration,
                                Gaussian2DPeakResultHelper.AMPLITUDE | Gaussian2DPeakResultHelper.PIXEL_AMPLITUDE);
                        final double o1 = calc.getAmplitude(paramsf);
                        final double o2 = calc.getPixelAmplitude(paramsf);

                        //logger.fine(FunctionUtils.getSupplier("e=%f, o1=%f, o2=%f", e, o1, o2));
                        Assertions.assertEquals(e, o2, 1e-3);
                        r.addData(e, o1);
                    }

        //logger.fine(FunctionUtils.getSupplier("Regression: pixel amplitude vs amplitude = %f, slope=%f, n=%d", r.getR(), r.getSlope(),
        //		r.getN()));
        // The simple amplitude over estimates the actual pixel amplitude
        Assertions.assertTrue(r.getSlope() > 1);
    }

    @Test
    public void canComputeCumulative()
    {
        Assertions.assertEquals(0, Gaussian2DPeakResultHelper.cumulative(0));
        Assertions.assertEquals(0.6827, Gaussian2DPeakResultHelper.cumulative(1), 1e-3);
        Assertions.assertEquals(0.9545, Gaussian2DPeakResultHelper.cumulative(2), 1e-3);
        Assertions.assertEquals(0.9974, Gaussian2DPeakResultHelper.cumulative(3), 1e-3);
        Assertions.assertTrue(1 == Gaussian2DPeakResultHelper.cumulative(Double.POSITIVE_INFINITY));
    }

    @Test
    public void canComputeCumulative2DAndInverse()
    {
        Assertions.assertEquals(0, Gaussian2DPeakResultHelper.cumulative2D(0));
        Assertions.assertTrue(1 == Gaussian2DPeakResultHelper.cumulative2D(Double.POSITIVE_INFINITY));
        Assertions.assertEquals(0, Gaussian2DPeakResultHelper.inverseCumulative2D(0));
        Assertions.assertTrue(Double.POSITIVE_INFINITY == Gaussian2DPeakResultHelper.inverseCumulative2D(1));
        for (int i = 1; i <= 10; i++)
        {
            final double r = i / 10.0;
            final double p = Gaussian2DPeakResultHelper.cumulative2D(r);
            final double r2 = Gaussian2DPeakResultHelper.inverseCumulative2D(p);
            ExtraAssertions.assertEqualsRelative(r, r2, 1e-8);
        }
    }

    @SeededTest
    public void canComputeMeanSignalUsingR(RandomSeed seed)
    {
        final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());

        for (int i = 0; i < 10; i++)
        {
            final double intensity = rg.nextDouble() * 100;
            final double sx = rg.nextDouble() * 2;
            final double sy = rg.nextDouble() * 2;
            final double r = rg.nextDouble() * 5;
            assertEquals(intensity * Gaussian2DPeakResultHelper.cumulative2D(r) / (Math.PI * r * r * sx * sy),
                    Gaussian2DPeakResultHelper.getMeanSignalUsingR(intensity, sx, sy, r));

            // Test fixed versions verse dynamic
            assertEquals(Gaussian2DPeakResultHelper.getMeanSignalUsingR(intensity, sx, sy, 1),
                    Gaussian2DPeakResultHelper.getMeanSignalUsingR1(intensity, sx, sy));
            assertEquals(Gaussian2DPeakResultHelper.getMeanSignalUsingR(intensity, sx, sy, 2),
                    Gaussian2DPeakResultHelper.getMeanSignalUsingR2(intensity, sx, sy));
        }
    }

    private static void assertEquals(double e, double o)
    {
        ExtraAssertions.assertEqualsRelative(e, o, 1e-10);
    }

    @SeededTest
    public void canComputeMeanSignalUsingP(RandomSeed seed)
    {
        final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());

        for (int i = 0; i < 10; i++)
        {
            final double intensity = rg.nextDouble() * 100;
            final double sx = rg.nextDouble() * 2;
            final double sy = rg.nextDouble() * 2;
            final double p = rg.nextDouble();
            double e = intensity * p /
                    (Math.PI * Maths.pow2(Gaussian2DPeakResultHelper.inverseCumulative2D(p)) * sx * sy);
            double o = Gaussian2DPeakResultHelper.getMeanSignalUsingP(intensity, sx, sy, p);
            assertEquals(e, o);

            // Test fixed versions verse dynamic
            e = Gaussian2DPeakResultHelper.getMeanSignalUsingP(intensity, sx, sy, 0.5);
            o = Gaussian2DPeakResultHelper.getMeanSignalUsingP05(intensity, sx, sy);
            assertEquals(e, o);
        }
    }
}
