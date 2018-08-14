package uk.ac.sussex.gdsc.smlm.function;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
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
import org.opentest4j.AssertionFailedError;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.smlm.function.PoissonGammaGaussianFunction.ConvolutionMode;
import uk.ac.sussex.gdsc.test.TestComplexity;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssumptions;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;

@SuppressWarnings({ "javadoc" })
public class PoissonGammaGaussianFunctionTest
{
    private static Logger logger;

    @BeforeAll
    public static void beforeAll()
    {
        logger = Logger.getLogger(PoissonGammaGaussianFunctionTest.class.getName());
    }

    @AfterAll
    public static void afterAll()
    {
        logger = null;
    }

    // Noise is in Counts and gain is total gain.
    // This makes more sense when testing as the
    // PoissonGammaGaussianFunction accepts 1/gain and noise as parameters.

    // Poisson-Gamma convolution sums to above 1 at lower gain.
    // due to the Dirac delta function, i.e. the Poisson-Gamma convolution
    // is a PDF and the Dirac delta is a PMF from the Poisson PMF at c=0.
    // This summing on the integer intervals (for a PMF) is invalid.
    // Store the expected sum at different gain below 10 for testing.
    static double[] pgSum = new double[11];
    static
    {
        // Compute the sum at expected photons around 1. This produces
        // the highest sum as the contribution from the Poisson-Gamma to c=0
        // will be the greatest.
        // These are rounded up to 3 d.p. provide a safer upper bound.
        final boolean compute = false;
        if (compute)
        {
            final MathContext mc = new MathContext(4, RoundingMode.UP);

            final StringBuilder sb = new StringBuilder();
            final String newLine = System.getProperty("line.separator");
            for (int g = 1; g <= 10; g++)
            {
                double max = 0;
                final int steps = 10;
                for (int i = 0; i <= steps; i++)
                {
                    final double e = 0.5 * i / steps;
                    // Compute half for the first interval
                    double sum = PoissonGammaFunction.poissonGammaN(0, e, g) * 0.5 + PoissonGammaFunction.dirac(e);
                    for (int c = 1;; c++)
                    {
                        final double p = PoissonGammaFunction.poissonGamma(c, e, g);
                        sum += p;
                        if (p / sum < 1e-6)
                            break;
                    }
                    max = Math.max(max, sum);
                }
                pgSum[g] = max;
                final BigDecimal bd = new BigDecimal(max);
                sb.append(newLine);
                sb.append(String.format("pgSum[%d] = %.3f;", g, bd.round(mc).doubleValue()));
            }
            logger.info(sb.toString());
        }

        pgSum[1] = 1.019;
        pgSum[2] = 1.005;
        pgSum[3] = 1.003;
        pgSum[4] = 1.002;
        pgSum[5] = 1.001;
        pgSum[6] = 1.001;
        pgSum[7] = 1.001;
        pgSum[8] = 1.001;
        pgSum[9] = 1.001;
        pgSum[10] = 1.001;
    }

    double[] photons = { 0, 0.25, 0.5, 1, 2, 4, 10, 100 };
    double[] highPhotons = { 5000 };
    double[] lowPhotons = { 1e-2, /* 1e-4, */ 1e-6 };
    double[] noise = { 3, 10 }; // in counts
    double[] lowNoise = { 0.3, 1 }; // in counts
    double[] totalGain = { 6.5, 45 };

    @Test
    public void cumulativeGaussianProbabilityIsCorrect()
    {
        for (final double s : noise)
            for (final double g : totalGain)
                cumulativeGaussianProbabilityIsCorrect(s, g);
    }

    private static void cumulativeGaussianProbabilityIsCorrect(double s, double g)
    {
        // Read noise should be in proportion to the camera gain
        final PoissonGammaGaussianFunction f = new PoissonGammaGaussianFunction(1 / g, s);
        final double range = 5 * s;
        final int upper = (int) Math.ceil(range);
        final int lower = (int) Math.floor(-range);
        final SimpsonIntegrator in = new SimpsonIntegrator(1e-4, 1e-8, 3, 32);
        final UnivariateFunction uf = new UnivariateFunction()
        {
            @Override
            public double value(double x)
            {
                return f.gaussianPDF(x);
            }
        };
        for (int u = lower; u <= upper; u++)
        {
            final double ux = u + 0.5;
            final double lx = u - 0.5;
            final double e = in.integrate(20000, uf, lx, ux);
            final double o = f.gaussianCDF(ux) - f.gaussianCDF(lx);
            final double o2 = f.gaussianCDF(lx, ux);
            ExtraAssertions.assertEqualsRelative(e, o, 0.1);
            ExtraAssertions.assertEqualsRelative(o, o2, 1e-6);
        }
    }

    // The Poisson-Gamma has a delta function at c=0. This causes problems
    // if not correctly integrated.
    // Some modes create a PMF, others a PDF so handle appropriately.

    @Test
    public void cumulativeProbabilityIsOneWithDiscretePMFIntegrationAsPMF()
    {
        for (final double p : photons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PMF, true);
    }

    @Test
    public void cumulativeProbabilityIsOneWithDiscretePMFIntegrationAsPDF()
    {
        for (final double p : photons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PMF, false);
    }

    @Test
    public void cumulativeProbabilityIsOneWithDiscretePDFIntegrationAsPMF()
    {
        for (final double p : photons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PDF, true);
    }

    @Test
    public void cumulativeProbabilityIsOneWithDiscretePDFIntegrationAsPDF()
    {
        for (final double p : photons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PDF, false);
    }

    @Test
    public void cumulativeProbabilityIsOneWithApproximationAtGainAbove10AsPDF()
    {
        for (final double p : photons)
            for (final double s : noise)
                for (final double g : totalGain)
                {
                    if (g < 10)
                        continue;
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.APPROXIMATION, false);
                }
    }

    @Test
    public void cumulativeProbabilityIsOneWithApproximationAtGainAbove10AsPMF()
    {
        for (final double p : photons)
            for (final double s : noise)
                for (final double g : totalGain)
                {
                    if (g < 10)
                        continue;
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.APPROXIMATION, true);
                }
    }

    @Test
    public void cumulativeProbabilityIsOneWithSimpsonIntegrationAsPDF()
    {
        for (final double p : photons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.SIMPSON_PDF, false);
    }

    @Test
    public void cumulativeProbabilityIsOneWithSimpsonIntegrationAsPMF()
    {
        for (final double p : photons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.SIMPSON_PDF, true);
    }

    @Test
    public void cumulativeProbabilityIsOneWithLegendreGaussIntegrationAsPDF()
    {
        ExtraAssumptions.assume(TestComplexity.HIGH);
        for (final double p : photons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.LEGENDRE_GAUSS_PDF, false);
    }

    @Test
    public void cumulativeProbabilityIsOneWithLegendreGaussIntegrationAsPMF()
    {
        ExtraAssumptions.assume(TestComplexity.HIGH);
        for (final double p : photons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.LEGENDRE_GAUSS_PDF, true);
    }

    @Test
    public void cumulativeProbabilityIsOneWithDiscretePMFIntegrationAsPMFAtLowNoise()
    {
        ExtraAssumptions.assume(TestComplexity.HIGH);
        for (final double p : photons)
            for (final double s : lowNoise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PMF, true);
    }

    @Test
    public void cumulativeProbabilityIsOneWithDiscretePMFIntegrationAsPDFAtLowNoise()
    {
        ExtraAssumptions.assume(TestComplexity.HIGH);
        for (final double p : photons)
            for (final double s : lowNoise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PMF, false);
    }

    @Test
    public void cumulativeProbabilityIsOneWithDiscretePDFIntegrationAsPMFAtLowNoise()
    {
        ExtraAssumptions.assume(TestComplexity.HIGH);
        for (final double p : photons)
            for (final double s : lowNoise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PDF, true);
    }

    @Test
    public void cumulativeProbabilityIsOneWithDiscretePDFIntegrationAsPDFAtLowNoise()
    {
        ExtraAssumptions.assume(TestComplexity.HIGH);
        for (final double p : photons)
            for (final double s : lowNoise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PDF, false);
    }

    @Test
    public void cumulativeProbabilityIsOneWithApproximationAtGainAbove10AsPDFAtLowNoise()
    {
        ExtraAssumptions.assume(TestComplexity.HIGH);
        for (final double p : photons)
            for (final double s : lowNoise)
                for (final double g : totalGain)
                {
                    if (g < 10)
                        continue;
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.APPROXIMATION, false);
                }
    }

    // The approximation is not meant to be used as a PMF
    @Test
    public void cumulativeProbabilityIsNotOneWithApproximationAtGainAbove10AsPMFAtLowNoise()
    {
        ExtraAssumptions.assume(TestComplexity.HIGH);

        // This is the whole test. It could be separated into parameters that fail
        // and parameters that are OK.
        Assertions.assertThrows(AssertionFailedError.class, () -> {
            for (final double p : photons)
                for (final double s : lowNoise)
                    for (final double g : totalGain)
                    {
                        if (g < 10)
                            continue;
                        cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.APPROXIMATION, true);
                    }
        });
    }

    @Test
    public void cumulativeProbabilityIsOneWithSimpsonIntegrationAsPDFAtLowNoise()
    {
        ExtraAssumptions.assume(TestComplexity.HIGH);
        for (final double p : photons)
            for (final double s : lowNoise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.SIMPSON_PDF, false);
    }

    @Test
    public void cumulativeProbabilityIsOneWithSimpsonIntegrationAsPMFAtLowNoise()
    {
        ExtraAssumptions.assume(TestComplexity.HIGH);
        for (final double p : photons)
            for (final double s : lowNoise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.SIMPSON_PDF, true);
    }

    @Test
    public void cumulativeProbabilityIsOneWithLegendreGaussIntegrationAsPDFAtLowNoise()
    {
        ExtraAssumptions.assume(TestComplexity.HIGH);
        for (final double p : photons)
            for (final double s : lowNoise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.LEGENDRE_GAUSS_PDF, false);
    }

    @Test
    public void cumulativeProbabilityIsOneWithLegendreGaussIntegrationAsPMFAtLowNoise()
    {
        ExtraAssumptions.assume(TestComplexity.HIGH);
        for (final double p : photons)
            for (final double s : lowNoise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.LEGENDRE_GAUSS_PDF, true);
    }

    @Test
    public void cumulativeProbabilityIsOneWithDiscretePMFIntegrationAsPMFAtHighPhotons()
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        for (final double p : highPhotons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PMF, true);
    }

    @Test
    public void cumulativeProbabilityIsOneWithDiscretePMFIntegrationAsPDFAtHighPhotons()
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        for (final double p : highPhotons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PMF, false);
    }

    @Test
    public void cumulativeProbabilityIsOneWithDiscretePDFIntegrationAsPMFAtHighPhotons()
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        for (final double p : highPhotons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PDF, true);
    }

    @Test
    public void cumulativeProbabilityIsOneWithDiscretePDFIntegrationAsPDFAtHighPhotons()
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        for (final double p : highPhotons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PDF, false);
    }

    @Test
    public void cumulativeProbabilityIsOneWithApproximationAsPDFAtHighPhotons()
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        for (final double p : highPhotons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.APPROXIMATION, false);
    }

    @Test
    public void cumulativeProbabilityIsOneWithApproximationAsPMFAtHighPhotons()
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        for (final double p : highPhotons)
            for (final double s : noise)
                for (final double g : totalGain)
                {
                    if (g < 10)
                        continue;
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.APPROXIMATION, true);
                }
    }

    @Test
    public void cumulativeProbabilityIsOneWithSimpsonIntegrationAsPDFAtHighPhotons()
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        for (final double p : highPhotons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.SIMPSON_PDF, false);
    }

    @Test
    public void cumulativeProbabilityIsOneWithSimpsonIntegrationAsPMFAtHighPhotons()
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        for (final double p : highPhotons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.SIMPSON_PDF, true);
    }

    @Test
    public void cumulativeProbabilityIsOneWithLegendreGaussIntegrationAsPDFAtHighPhotons()
    {
        ExtraAssumptions.assume(TestComplexity.HIGH);
        for (final double p : highPhotons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.LEGENDRE_GAUSS_PDF, false);
    }

    @Test
    public void cumulativeProbabilityIsOneWithLegendreGaussIntegrationAsPMFAtHighPhotons()
    {
        ExtraAssumptions.assume(TestComplexity.HIGH);
        for (final double p : highPhotons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.LEGENDRE_GAUSS_PDF, true);
    }

    @Test
    public void cumulativeProbabilityIsOneWithDiscretePMFIntegrationAsPMFAtLowPhotons()
    {
        for (final double p : lowPhotons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PMF, true);
    }

    @Test
    public void cumulativeProbabilityIsOneWithDiscretePMFIntegrationAsPDFAtLowPhotons()
    {
        for (final double p : lowPhotons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PMF, false);
    }

    @Test
    public void cumulativeProbabilityIsOneWithDiscretePDFIntegrationAsPMFAtLowPhotons()
    {
        for (final double p : lowPhotons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PDF, true);
    }

    @Test
    public void cumulativeProbabilityIsOneWithDiscretePDFIntegrationAsPDFAtLowPhotons()
    {
        for (final double p : lowPhotons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PDF, false);
    }

    @Test
    public void cumulativeProbabilityIsOneWithApproximationAtGainAbove10AsPDFAtLowPhotons()
    {
        for (final double p : lowPhotons)
            for (final double s : noise)
                for (final double g : totalGain)
                {
                    if (g < 10)
                        continue;
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.APPROXIMATION, false);
                }
    }

    @Test
    public void cumulativeProbabilityIsOneWithApproximationAtGainAbove10AsPMFAtLowPhotons()
    {
        for (final double p : lowPhotons)
            for (final double s : noise)
                for (final double g : totalGain)
                {
                    if (g < 10)
                        continue;
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.APPROXIMATION, true);
                }
    }

    @Test
    public void cumulativeProbabilityIsOneWithSimpsonIntegrationAsPDFAtLowPhotons()
    {
        for (final double p : lowPhotons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.SIMPSON_PDF, false);
    }

    @Test
    public void cumulativeProbabilityIsOneWithSimpsonIntegrationAsPMFAtLowPhotons()
    {
        for (final double p : lowPhotons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.SIMPSON_PDF, true);
    }

    @Test
    public void cumulativeProbabilityIsOneWithLegendreGaussIntegrationAsPDFAtLowPhotons()
    {
        for (final double p : lowPhotons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.LEGENDRE_GAUSS_PDF, false);
    }

    @Test
    public void cumulativeProbabilityIsOneWithLegendreGaussIntegrationAsPMFAtLowPhotons()
    {
        for (final double p : lowPhotons)
            for (final double s : noise)
                for (final double g : totalGain)
                    cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.LEGENDRE_GAUSS_PDF, true);
    }

    @Test
    public void discretePDFCloselyMatchesPMFIntegration()
    {
        final double[] e = closelyMatchesPMFIntegration(0.34, ConvolutionMode.DISCRETE_PDF);
        logger.log(TestLog.getRecord(Level.FINE, "Discrete integration max error : rel = %g : abs = %g", e[0], e[1]));
    }

    @Test
    public void discretePMFCloselyMatchesPMFIntegration()
    {
        final double[] e = closelyMatchesPMFIntegration(0.22, ConvolutionMode.DISCRETE_PMF);
        logger.log(TestLog.getRecord(Level.FINE, "Discrete integration max error : rel = %g : abs = %g", e[0], e[1]));
    }

    @Test
    public void approximationCloselyMatchesPMFIntegration()
    {
        final double[] e = closelyMatchesPMFIntegration(0.22, ConvolutionMode.APPROXIMATION);
        logger.log(TestLog.getRecord(Level.FINE, "Approximation max error : rel = %g : abs = %g", e[0], e[1]));
    }

    @Test
    public void legedreGaussPDFMatchesPMFIntegration()
    {
        final double[] e = closelyMatchesPMFIntegration(0.03, ConvolutionMode.LEGENDRE_GAUSS_PDF);
        logger.log(TestLog.getRecord(Level.FINE, "Simpson integration max error : rel = %g : abs = %g", e[0], e[1]));
    }

    // Speed order is roughly: Approx, Simpson, Discrete PDF, Legendre, Discrete PMF
    // The most accurate over most settings p<<1e-5, p=1, p>>10 is the Simpson.

    @SpeedTag
    @SeededTest
    public void approximationFasterThanSimpsonIntegration(RandomSeed seed)
    {
        fasterThan(seed, ConvolutionMode.SIMPSON_PDF, ConvolutionMode.APPROXIMATION);
    }

    @SpeedTag
    @SeededTest
    public void simpsonIntegrationFasterThanDiscretePDFIntegration(RandomSeed seed)
    {
        fasterThan(seed, ConvolutionMode.DISCRETE_PDF, ConvolutionMode.SIMPSON_PDF);
    }

    @SpeedTag
    @SeededTest
    public void simpsonIntegrationFasterThanLegendreGaussIntegration(RandomSeed seed)
    {
        fasterThan(seed, ConvolutionMode.LEGENDRE_GAUSS_PDF, ConvolutionMode.SIMPSON_PDF);
    }

    @SpeedTag
    @SeededTest
    public void discretePDFIntegrationFasterThanDiscretePMFIntegration(RandomSeed seed)
    {
        fasterThan(seed, ConvolutionMode.DISCRETE_PMF, ConvolutionMode.DISCRETE_PDF);
    }

    private static void cumulativeProbabilityIsOne(final double mu, final double s, final double g,
            ConvolutionMode convolutionMode, boolean pmfMode)
    {
        final double p = cumulativeProbability(mu, s, g, convolutionMode, pmfMode);
        logger.log(TestLog.getRecord(Level.INFO, "%s : mu=%f, s=%f, g=%f, p=%f", getName(convolutionMode), mu, s, g, p));

        // Poisson-Gamma convolution approximation does not sum to 1 at lower gain
        // so account for this during the test.
        final double delta = 0.02;
        double upper = 1 + delta;
        final double lower = 1 - delta;
        if (g < 10)
            upper = pgSum[(int) g] + delta;

        if (p < lower || p > upper)
            Assertions.fail(String.format("mu=%f, s=%f, g=%f, p=%g", mu, s, g, p));
    }

    private static double cumulativeProbability(final double mu, final double s, final double g,
            ConvolutionMode convolutionMode, boolean pmfMode)
    {
        final PoissonGammaGaussianFunction f = new PoissonGammaGaussianFunction(1 / g, s);
        f.setConvolutionMode(convolutionMode);
        f.setPmfMode(pmfMode);
        f.setMinimumProbability(0);
        double p = 0;
        int min = 1;
        int max = 0;

        // Evaluate an initial range.
        // Gaussian should have >99% within +/- 3s
        // Poisson will have mean mu with a variance mu.
        // At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
        if (mu > 0)
        {
            final int[] range = PoissonGaussianFunctionTest.getRange(g, mu, s);
            min = range[0];
            max = range[1];
            for (int x = min; x <= max; x++)
            {
                final double pp = f.likelihood(x, mu);
                //TestLog.fine(logger,"x=%d, p=%g", x, pp);
                p += pp;
            }
            //if (p > 1.01)
            //	Assertions.fail("P > 1: " + p);
        }

        // We have most of the probability density.
        // Now keep evaluating up and down until no difference
        final double changeTolerance = 1e-6;
        for (int x = min - 1;; x--)
        {
            min = x;
            final double pp = f.likelihood(x, mu);
            //TestLog.fine(logger,"x=%d, p=%g", x, pp);
            p += pp;
            if (pp / p < changeTolerance)
                break;
        }
        for (int x = max + 1;; x++)
        {
            max = x;
            final double pp = f.likelihood(x, mu);
            //TestLog.fine(logger,"x=%d, p=%g", x, pp);
            p += pp;
            if (pp / p < changeTolerance)
                break;
        }

        // This is a simple integral. Compute the full integral if necessary.
        if (!pmfMode && (p < 0.98 || p > 1.02))
        {
            // Do a formal integration
            final UnivariateIntegrator in = new SimpsonIntegrator(1e-6, 1e-6, 4,
                    SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
            final double pp = in.integrate(Integer.MAX_VALUE, new UnivariateFunction()
            {
                @Override
                public double value(double x)
                {
                    return f.likelihood(x, mu);
                }
            }, min, max);
            logger.log(TestLog.getRecord(Level.FINE, "%s : mu=%f, rn=%f, cg=%f, s=%f, g=%f, p=%g => %g", getName(convolutionMode), mu, s, g,
                    s, g, p, pp));
            p = pp;
        }

        return p;
    }

    private double[] closelyMatchesPMFIntegration(double error, ConvolutionMode convolutionMode)
    {
        //DoubleEquality eq = new DoubleEquality(error, 1e-7);
        final double[] maxError = new double[2];
        for (final double s : noise)
            for (final double g : totalGain)
            {
                if (g < 10)
                    continue;

                // This is the reference for a PMF-type result
                final PoissonGammaGaussianFunction f1 = new PoissonGammaGaussianFunction(1 / g, s);
                f1.setConvolutionMode(ConvolutionMode.SIMPSON_PDF);
                f1.setPmfMode(true);
                f1.setMinimumProbability(0);

                final PoissonGammaGaussianFunction f2 = new PoissonGammaGaussianFunction(1 / g, s);
                f2.setConvolutionMode(convolutionMode);
                f2.setPmfMode(true);
                f2.setMinimumProbability(0);

                for (final double p : photons)
                {
                    final double pg = p * g;
                    final double min = pg * 0.5 - 5 * s;
                    final double max = 2 * pg;
                    for (double x = min; x < max; x += 1)
                    {
                        final double p1 = f1.likelihood(x, p);
                        final double p2 = f2.likelihood(x, p);

                        final double relativeError = DoubleEquality.relativeError(p1, p2);
                        final double absError = Math.abs(p1 - p2);
                        final boolean equal = relativeError <= error; //eq.almostEqualRelativeOrAbsolute(p1, p2);
                        if (!equal)
                        {
                            // Ignore small probabilities
                            if (p1 < 1e-3)
                                continue;

                            Assertions.fail(String.format("s=%g, g=%g, p=%g, x=%g: %g != %g (%g)", s, g, p, x, p1, p2,
                                    relativeError));
                        }
                        if (maxError[0] < relativeError)
                            maxError[0] = relativeError;
                        if (maxError[1] < absError)
                            maxError[1] = absError;
                    }
                }
            }
        return maxError;
    }

    private void fasterThan(RandomSeed seed, ConvolutionMode slow, ConvolutionMode fast)
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);

        // Realistic EM-CCD parameters for speed test
        final double s = 7.16;
        final double g = 39.1;

        final PoissonGammaGaussianFunction f1 = new PoissonGammaGaussianFunction(1 / g, s);
        f1.setConvolutionMode(slow);

        final PoissonGammaGaussianFunction f2 = new PoissonGammaGaussianFunction(1 / g, s);
        f2.setConvolutionMode(fast);

        final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());

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

        logger.log(TestLog.getTimingRecord(getName(f1), t1, getName(f2), t2));
    }

    private static long run(PoissonGammaGaussianFunction f, double[][] samples, double[] photons)
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

    private static String getName(PoissonGammaGaussianFunction f)
    {
        return getName(f.getConvolutionMode());
    }

    private static String getName(ConvolutionMode convolutionMode)
    {
        return convolutionMode.toString();
    }
}
