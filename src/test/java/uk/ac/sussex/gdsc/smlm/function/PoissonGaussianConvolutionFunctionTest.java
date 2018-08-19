package uk.ac.sussex.gdsc.smlm.function;

import java.util.function.Supplier;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssumptions;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RNGFactory;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLog;

@SuppressWarnings({ "javadoc" })
public class PoissonGaussianConvolutionFunctionTest
{
    private static Logger logger;

    @BeforeAll
    public static void beforeAll()
    {
        logger = Logger.getLogger(PoissonGaussianConvolutionFunctionTest.class.getName());
    }

    @AfterAll
    public static void afterAll()
    {
        logger = null;
    }

    double[] gain = PoissonGaussianFunctionTest.gain;
    double[] photons = PoissonGaussianFunctionTest.photons;
    double[] noise = PoissonGaussianFunctionTest.noise;

    @Test
    public void cumulativeProbabilityIsOneWithPDF()
    {
        for (final double g : gain)
            for (final double p : photons)
                for (final double s : noise)
                    cumulativeProbabilityIsOne(g, p, s, false);
    }

    @Test
    public void cumulativeProbabilityIsOneWithPMF()
    {
        for (final double g : gain)
            for (final double p : photons)
                for (final double s : noise)
                    cumulativeProbabilityIsOne(g, p, s, true);
    }

    @Test
    public void probabilityMatchesLogProbabilityWithPDF()
    {
        for (final double g : gain)
            for (final double p : photons)
                for (final double s : noise)
                    probabilityMatchesLogProbability(g, p, s, false);
    }

    @Test
    public void probabilityMatchesLogProbabilityWithPMF()
    {
        for (final double g : gain)
            for (final double p : photons)
                for (final double s : noise)
                    probabilityMatchesLogProbability(g, p, s, true);
    }

    private static void cumulativeProbabilityIsOne(final double gain, final double mu, final double s,
            boolean computePMF)
    {
        final double p2 = cumulativeProbability(gain, mu, s, computePMF);
        Assertions.assertEquals(1, p2, 0.02, () -> String.format("g=%f, mu=%f, s=%f, erf=%b", gain, mu, s, computePMF));
    }

    private static double cumulativeProbability(final double gain, final double mu, final double s, boolean computePMF)
    {
        // Note: The input s parameter is pre-gain.
        final PoissonGaussianConvolutionFunction f = PoissonGaussianConvolutionFunction
                .createWithStandardDeviation(1.0 / gain, s * gain);
        f.setComputePMF(computePMF);

        //final PoissonGaussianConvolutionFunction f2 = PoissonGaussianConvolutionFunction.createWithStandardDeviation(1.0 / gain, mu*gain, s * gain);
        //f2.setUsePicardApproximation(usePicard);

        double p = 0;
        int min = 1;
        int max = 0;

        // Note: The input mu parameter is pre-gain.
        final double e = mu;

        // Evaluate an initial range.
        // Gaussian should have >99% within +/- s
        // Poisson will have mean mu with a variance mu.
        // At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
        if (mu > 0)
        {
            final int[] range = PoissonGaussianFunctionTest.getRange(gain, mu, s);
            min = range[0];
            max = range[1];
            for (int x = min; x <= max; x++)
            {
                final double pp = f.likelihood(x, e);
                //logger.fine(FunctionUtils.getSupplier("x=%d, p=%g", x, pp);
                //logger.fine(FunctionUtils.getSupplier("x=%d, p=%f   %f", x, pp, f2.probability(x));
                p += pp;
            }
            //if (p > 1.01)
            //	Assertions.fail("P > 1: " + p);
        }

        // We have most of the likelihood density.
        // Now keep evaluating up and down until no difference
        final double changeTolerance = 1e-6;
        for (int x = min - 1;; x--)
        {
            min = x;
            final double pp = f.likelihood(x, e);
            //logger.fine(FunctionUtils.getSupplier("x=%d, p=%g", x, pp);
            p += pp;
            if (pp == 0 || pp / p < changeTolerance)
                break;
        }
        for (int x = max + 1;; x++)
        {
            max = x;
            final double pp = f.likelihood(x, e);
            //logger.fine(FunctionUtils.getSupplier("x=%d, p=%g", x, pp);
            p += pp;
            if (pp == 0 || pp / p < changeTolerance)
                break;
        }

        double p2 = p;
        if (!computePMF)
        {
            // Do a formal integration if the PDF
            //if (p < 0.98 || p > 1.02)
            //	logger.fine(FunctionUtils.getSupplier("g=%f, mu=%f, s=%f p=%f", gain, mu, s, p);
            final UnivariateIntegrator in = new SimpsonIntegrator(1e-4, 1e-6, 4,
                    SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
            p2 = in.integrate(Integer.MAX_VALUE, new UnivariateFunction()
            {
                @Override
                public double value(double x)
                {
                    return f.likelihood(x, e);
                }
            }, min, max);
        }

        if (p2 < 0.98 || p2 > 1.02)
            logger.log(TestLog.getRecord(Level.INFO, "g=%f, mu=%f, s=%f p=%f  %f", gain, mu, s, p, p2));

        return p2;
    }

    private static void probabilityMatchesLogProbability(final double gain, double mu, final double s,
            boolean computePMF)
    {
        // Note: The input s parameter is pre-gain.
        final PoissonGaussianConvolutionFunction f = PoissonGaussianConvolutionFunction
                .createWithStandardDeviation(1.0 / gain, s * gain);
        f.setComputePMF(computePMF);

        // Evaluate an initial range.
        // Gaussian should have >99% within +/- s
        // Poisson will have mean mu with a variance mu.
        // At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
        final int[] range = PoissonGaussianFunctionTest.getRange(gain, mu, s);
        final int min = range[0];
        final int max = range[1];
        // Note: The input mu parameter is pre-gain.
        final double e = mu;
        final Supplier<String> msg = () -> String.format("g=%f, mu=%f, s=%f, erf=%b", gain, mu, s, computePMF);
        for (int x = min; x <= max; x++)
        {
            final double p = f.likelihood(x, e);
            if (p == 0)
                continue;
            final double logP = f.logLikelihood(x, e);
            ExtraAssertions.assertEqualsRelative(Math.log(p), logP, 1e-3, msg);
        }
    }

    @SpeedTag
    @SeededTest
    public void pdfFasterThanPMF(RandomSeed seed)
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);

        // Realistic CCD parameters for speed test
        final double s = 7.16;
        final double g = 3.1;

        final PoissonGaussianConvolutionFunction f1 = PoissonGaussianConvolutionFunction
                .createWithStandardDeviation(1 / g, s);
        f1.setComputePMF(true);

        final PoissonGaussianConvolutionFunction f2 = PoissonGaussianConvolutionFunction
                .createWithStandardDeviation(1 / g, s);
        f2.setComputePMF(false);

        final UniformRandomProvider rg = RNGFactory.create(seed.getSeed());

        // Generate realistic data from the probability mass function
        final double[][] samples = new double[photons.length][];
        for (int j = 0; j < photons.length; j++)
        {
            final int start = (int) (4 * -s);
            int u = start;
            final StoredDataStatistics stats = new StoredDataStatistics();
            while (stats.getSum() < 0.995)
            {
                final double p = f1.likelihood(u, photons[j]);
                stats.add(p);
                if (u > 10 && p / stats.getSum() < 1e-6)
                    break;
                u++;
            }

            // Generate cumulative probability
            final double[] data = stats.getValues();
            for (int i = 1; i < data.length; i++)
                data[i] += data[i - 1];
            // Normalise
            for (int i = 0, end = data.length - 1; i < data.length; i++)
                data[i] /= data[end];

            // Sample
            final double[] sample = new double[1000];
            for (int i = 0; i < sample.length; i++)
            {
                final double p = rg.nextDouble();
                int x = 0;
                while (x < data.length && data[x] < p)
                    x++;
                sample[i] = start + x;
            }
            samples[j] = sample;
        }

        // Warm-up
        run(f1, samples, photons);
        run(f2, samples, photons);

        long t1 = 0;
        for (int i = 0; i < 5; i++)
            t1 += run(f1, samples, photons);

        long t2 = 0;
        for (int i = 0; i < 5; i++)
            t2 += run(f2, samples, photons);

        logger.log(TestLog.getTimingRecord("cdf", t1, "pdf", t2));
    }

    private static long run(PoissonGaussianConvolutionFunction f, double[][] samples, double[] photons)
    {
        final long start = System.nanoTime();
        for (int j = 0; j < photons.length; j++)
        {
            final double p = photons[j];
            for (final double x : samples[j])
                f.likelihood(x, p);
        }
        return System.nanoTime() - start;
    }

}
