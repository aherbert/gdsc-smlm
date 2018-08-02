package uk.ac.sussex.gdsc.smlm.function;

import java.util.logging.Level;
import java.util.logging.Logger;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import uk.ac.sussex.gdsc.test.TestComplexity;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssumptions;

@SuppressWarnings({ "javadoc" })
public class PoissonGammaGaussianFisherInformationTest
{
    private static Logger logger;

    @BeforeAll
    public static void beforeAll()
    {
        logger = Logger.getLogger(PoissonGammaGaussianFisherInformationTest.class.getName());
    }

    @AfterAll
    public static void afterAll()
    {
        logger = null;
    }

    @Test
    public void canFindMaximumAndUpperLimit()
    {
        ExtraAssumptions.assume(logger, Level.INFO);
        ExtraAssumptions.assume(TestComplexity.MEDIUM);

        final double[] M = { 20, 500 };

        for (final double m : M)
            canFindMaximumAndUpperLimit(m);
    }

    private static void canFindMaximumAndUpperLimit(double m)
    {
        final PoissonGammaGaussianFisherInformation f = new PoissonGammaGaussianFisherInformation(m, 1);
        f.setMeanThreshold(Double.POSITIVE_INFINITY);

        // Due to a limited convolution (for an infinite range)
        // the class works up to mean of about 300. Above that the approximation using
        // half the Poisson Fisher information should be used instead.
        // 10^2 == 100 (OK), 10^2.5 == 316 (Fail)
        for (int exp = -12; exp <= 4; exp++)
            canFindMaximumAndUpperLimit(f, Math.pow(10, exp * 0.5));
        for (int i = 0; i <= 10; i++)
            canFindMaximumAndUpperLimit(f, 0.001 + i * 0.1);
    }

    private static void canFindMaximumAndUpperLimit(PoissonGammaGaussianFisherInformation f, double u)
    {
        final double[] max = f.findMaximum(u, 1e-6);
        final double[] upper = f.findUpperLimit(u, max, 1e-6);
        logger.fine(TestLog.getSupplier("m=%g u=%g max=%s %s (%s)  upper=%s %s (%s)", f.m, u, max[0], max[1], max[2],
                upper[0], upper[1], upper[2]));
    }

    @Test
    public void canComputeFisherInformation()
    {
        //canComputeFisherInformation(250, 13);
        //double m1 = 500;
        //double u = 350;
        //double s1 = 13;
        //PoissonGammaGaussianFisherInformation f = new PoissonGammaGaussianFisherInformation(m1, s1);
        //f.setMeanThreshold(Double.POSITIVE_INFINITY);
        //double I = f.getPoissonGammaGaussianI(u);
        //double upper = PoissonFisherInformation.getPoissonI(u);
        //logger.fine(TestLog.getSupplier("s=%g u=%g I=%s PoissonI=%s alpha=%s", s1, u, I, upper, I / upper);
        //if (true)
        //	return;

        final double[] M = { 20, 500 };
        final double[] S = { 3, 13 };

        for (final double m : M)
            for (final double s : S)
                canComputeFisherInformation(m, s);
    }

    private static void canComputeFisherInformation(double m, double s)
    {
        canComputeFisherInformation(new PoissonGammaGaussianFisherInformation(m, s));
    }

    private static void canComputeFisherInformation(PoissonGammaGaussianFisherInformation f)
    {
        f.setMeanThreshold(Double.POSITIVE_INFINITY);

        // Due to a limited convolution (for an infinite range)
        // the class works up to mean of about 300. Above that the approximation using
        // half the Poisson Fisher information should be used instead.
        // 10^2 == 100 (OK), 10^2.5 == 316 (Fail)

        // exp == -12 => p = 1e-6
        // exp ==  -8 => p = 1e-4
        // exp ==   0 => p = 1
        // exp ==   4 => p = 100
        for (int exp = -8; exp <= 4; exp++)
            //System.out.println(Math.pow(10, exp * 0.5));
            canComputeFisherInformation(f, Math.pow(10, exp * 0.5));
    }

    private static void canComputeFisherInformation(PoissonGammaGaussianFisherInformation f, double u)
    {
        final double I = f.getPoissonGammaGaussianI(u);
        final double upper = PoissonFisherInformation.getPoissonI(u);
        //logger.fine(TestLog.getSupplier("m=%g s=%g u=%g I=%s PoissonI=%s alpha=%s", f.m, f.s, u, I, upper, I / upper);
        Assertions.assertTrue(I <= upper, "Not less than Poisson information");
        // This is true at higher mean
        if (u > 10)
            Assertions.assertTrue(I >= 0.4999 * upper, "Not above half the Poisson information");
    }

    @Test
    public void canComputeAlpha()
    {
        // This is a report as nothing is asserted
        ExtraAssumptions.assume(logger, Level.INFO);
        ExtraAssumptions.assume(TestComplexity.VERY_HIGH);

        // Compute the alpha using a range of gain and standard deviation

        final int minm = 100;
        final int maxm = 200;

        // When Poisson mean is high s does not matter as the Dirac is insignificant
        // and the convolution is mute. This may not be true when m is low but at higher
        // m the function is smooth and convolution has little effect.
        final int mins = 5;
        final int maxs = 5;

        // The relative Fisher information plateaus at low photons.
        // Output where this is.
        final double p = 1e-100;
        final double p2 = 1e-200; //Double.MIN_NORMAL; //Double.longBitsToDouble(0x4000000000001L);

        final double upper = PoissonFisherInformation.getPoissonI(p);
        final double upper2 = PoissonFisherInformation.getPoissonI(p2);

        for (int s = mins; s <= maxs; s++)
        {
            double lastAlpha = 1;
            for (int m = minm; m <= maxm; m++)
            {
                final PoissonGammaGaussianFisherInformation f = new PoissonGammaGaussianFisherInformation(m, s);
                final double I = f.getPoissonGammaGaussianI(p);
                final double I2 = f.getPoissonGammaGaussianI(p2);
                final double alpha = I / upper;
                final double alpha2 = I2 / upper2;
                final double change = lastAlpha / alpha;
                logger.fine(TestLog.getSupplier("p=%g  p2=%s   m=%s  s=%s   I=%s (%s)  alpha=%s (%s)  (delta=%s)", p,
                        p2, m, s, I, I2, alpha, alpha2, change));
                lastAlpha = alpha;
            }
        }
    }

    // Note: In practice the low mean is not used during fitting as the background photons
    // will always contribute to each pixel. Values below 0.01 are rare.

    @Test
    public void canComputeFisherInformationWithLowMean()
    {
        ExtraAssumptions.assumeLowComplexity();
        computeFisherInformationWithMean(1e-100);
    }

    @Test
    public void canComputeFisherInformationWithVeryLowMean()
    {
        ExtraAssumptions.assumeMediumComplexity();
        computeFisherInformationWithMean(1e-300);
    }

    @Test
    public void canComputeFisherInformationWithLowestPossibleMean()
    {
        ExtraAssumptions.assumeHighComplexity();

        // Lowest value where the reciprocal is not infinity.
        double u = Double.longBitsToDouble(0x4000000000001L);

        // Binary search for the min value
        final boolean doSearch = false;
        if (doSearch)
        {
            // This is the full 52-bit mantissa of a double with zero for the unbiased exponent,
            // i.e. the largest sub-normal number.
            long upper = (1L << 52) - 1;
            Assertions.assertTrue(1.0 / Double.longBitsToDouble(upper) != Double.POSITIVE_INFINITY);
            long lower = 1;
            while (lower + 1 < upper)
            {
                // 1/Upper is not infinity
                // Test mid-point
                final long mid = (upper + lower) / 2;
                u = Double.longBitsToDouble(mid);
                if (1 / u == Double.POSITIVE_INFINITY)
                    lower = mid;
                else
                    // Mid point
                    upper = mid;
            }

            u = Double.longBitsToDouble(upper);
            logger.info(TestLog.getSupplier("(upper = 0x%s = %s", Long.toHexString(upper), u));
        }

        Assertions.assertTrue(1.0 / u != Double.POSITIVE_INFINITY);
        Assertions.assertTrue(1.0 / Math.nextDown(u) == Double.POSITIVE_INFINITY);

        computeFisherInformationWithMean(u);
    }

    private static void computeFisherInformationWithMean(double u)
    {
        final double[] M = { 20, 500 };
        final double[] S = { 3, 13 };

        for (final double m : M)
            for (final double s : S)
            {
                final PoissonGammaGaussianFisherInformation f = new PoissonGammaGaussianFisherInformation(m, s);
                final double I = f.getPoissonGammaGaussianI(u);
                final double upper = PoissonFisherInformation.getPoissonI(u);
                final double alpha = I / upper;
                TestLog.fine(logger, "m=%g s=%g u=%g I=%s PoissonI=%s alpha=%s", f.m, f.s, u, I, upper, alpha);
                ExtraAssertions.assertTrue(I < upper, "Fisher information (%s) is not below upper limit: %s", I, upper);
                Assertions.assertTrue(alpha > 0, "Alpha is not above zero");
            }
    }
}
