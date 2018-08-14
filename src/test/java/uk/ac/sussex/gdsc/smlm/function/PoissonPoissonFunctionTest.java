package uk.ac.sussex.gdsc.smlm.function;

import java.util.Arrays;
import java.util.function.Supplier;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;

@SuppressWarnings({ "javadoc" })
public class PoissonPoissonFunctionTest
{
    private static Logger logger;

    @BeforeAll
    public static void beforeAll()
    {
        logger = Logger.getLogger(PoissonPoissonFunctionTest.class.getName());
    }

    @AfterAll
    public static void afterAll()
    {
        logger = null;
    }

    static double[] gain = PoissonGaussianFunctionTest.gain;
    static double[] photons = PoissonGaussianFunctionTest.photons;
    static
    {
        int c = 0;

        // No gain below 1
        c = 0;
        for (int i = 0; i < gain.length; i++)
            if (gain[i] >= 1)
                gain[c++] = gain[i];
        gain = Arrays.copyOf(gain, c);

        // No negative photons
        c = 0;
        for (int i = 0; i < photons.length; i++)
            if (photons[i] >= 0)
                photons[c++] = photons[i];
        photons = Arrays.copyOf(photons, c);
    }
    static double[] noise = PoissonGaussianFunctionTest.noise;

    @Test
    public void cumulativeProbabilityIsOne()
    {
        for (final double g : gain)
            for (final double p : photons)
                for (final double s : noise)
                    cumulativeProbabilityIsOne(g, p, s);
    }

    @Test
    public void probabilityMatchesLogProbability()
    {
        for (final double g : gain)
            for (final double p : photons)
                for (final double s : noise)
                    probabilityMatchesLogProbability(g, p, s);
    }

    private static void cumulativeProbabilityIsOne(final double gain, final double mu, final double s)
    {
        final double p2 = cumulativeProbability(gain, mu, s);
        // Only true with continuous distribution if the combined Poisson mean is above 4
        if (mu + s / gain > 4)
            Assertions.assertEquals(1, p2, 0.02, () -> String.format("g=%f, mu=%f, s=%f", gain, mu, s));
    }

    private static double cumulativeProbability(final double gain, final double mu, final double s)
    {
        // Note: The input s parameter is pre-gain.
        final PoissonPoissonFunction f = PoissonPoissonFunction.createWithStandardDeviation(1.0 / gain, s * gain);

        //final PoissonGaussianFunction f2 = PoissonGaussianFunction.createWithStandardDeviation(1.0 / gain, mu*gain, s * gain);
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
                //logger.fine(FunctionUtils.getSupplier("x=%d, p=%f", x, pp);
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
            //logger.fine(FunctionUtils.getSupplier("x=%d, p=%f", x, pp);
            p += pp;
            if (pp == 0 || pp / p < changeTolerance)
                break;
        }
        for (int x = max + 1;; x++)
        {
            max = x;
            final double pp = f.likelihood(x, e);
            //logger.fine(FunctionUtils.getSupplier("x=%d, p=%f", x, pp);
            p += pp;
            if (pp == 0 || pp / p < changeTolerance)
                break;
        }

        // Do a formal integration
        double p2 = 0;
        final UnivariateIntegrator in = new SimpsonIntegrator(1e-4, 1e-6, 3,
                SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
        p2 = in.integrate(Integer.MAX_VALUE, new UnivariateFunction()
        {
            @Override
            public double value(double x)
            {
                return f.likelihood(x, e);
            }
        }, min, max);

        logger.log(TestLog.getRecord(Level.INFO, "g=%f, mu=%f, s=%f p=%f  %f", gain, mu, s, p, p2));

        return p2;
    }

    private static void probabilityMatchesLogProbability(final double gain, double mu, final double s)
    {
        // Note: The input s parameter is pre-gain.
        final PoissonPoissonFunction f = PoissonPoissonFunction.createWithStandardDeviation(1.0 / gain, s * gain);

        // Evaluate an initial range.
        // Gaussian should have >99% within +/- s
        // Poisson will have mean mu with a variance mu.
        // At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
        final int[] range = PoissonGaussianFunctionTest.getRange(gain, mu, s);
        final int min = range[0];
        final int max = range[1];
        // Note: The input mu parameter is pre-gain.
        final double e = mu;
        final Supplier<String> msg = () -> String.format("g=%f, mu=%f, s=%f", gain, mu, s);
        for (int x = min; x <= max; x++)
        {
            final double p = f.likelihood(x, e);
            if (p == 0)
                continue;
            final double logP = f.logLikelihood(x, e);
            ExtraAssertions.assertEqualsRelative(Math.log(p), logP, 1e-3, msg);
        }
    }
}
