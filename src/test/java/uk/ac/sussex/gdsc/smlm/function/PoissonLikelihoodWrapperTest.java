package uk.ac.sussex.gdsc.smlm.function;

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Precision;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssumptions;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLog;

@SuppressWarnings({ "javadoc" })
public class PoissonLikelihoodWrapperTest
{
    private static Logger logger;

    @BeforeAll
    public static void beforeAll()
    {
        logger = Logger.getLogger(PoissonLikelihoodWrapperTest.class.getName());
    }

    @AfterAll
    public static void afterAll()
    {
        logger = null;
    }

    double alpha = 1 / 40.0;
    double[] photons = { 0.25, 0.5, 1, 2, 4, 10, 100, 1000 };
    // Set this at the range output from cumulativeProbabilityIsOneWithIntegerData
    int[] maxRange = { 6, 7, 10, 13, 17, 29, 149, 1141 };

    DoubleEquality eqPerDatum = new DoubleEquality(5e-2, 0.01);
    DoubleEquality eq = new DoubleEquality(5e-3, 0.001);

    static String[] NAME;
    static
    {
        NAME = new String[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
        for (int i = 0; i < NAME.length; i++)
            NAME[i] = Gaussian2DFunction.getName(i);
    }

    // Compute as per Numerical Recipes 5.7.
    // Approximate error accuracy in single precision: Ef
    // Step size for derivatives:
    // h ~ (Ef)^(1/3) * xc
    // xc is the characteristic scale over which x changes, assumed to be 1 (not x as per NR since x is close to zero)
    final double h_ = 0.01; //(double) (Math.pow(1e-3f, 1.0 / 3));

    int[] testx = new int[] { 4, 5, 6 };
    int[] testy = new int[] { 4, 5, 6 };
    // Do not test zero background since this is an edge case for the likelihood function
    double[] testbackground_ = new double[] { 10, 400 };

    double[] testsignal1_ = new double[] { 15, 55, 105 };
    double[] testcx1_ = new double[] { 4.9, 5.3 };
    double[] testcy1_ = new double[] { 4.8, 5.2 };
    double[] testcz1_ = new double[] { -1.5, 1.0 };
    double[][] testw1_ = new double[][] { { 1.1, 1.4 }, { 1.1, 1.7 }, { 1.5, 1.2 }, { 1.3, 1.7 }, };
    double[] testangle1_ = new double[] { Math.PI / 5, Math.PI / 3 };

    double[] testbackground, testsignal1, testangle1, testcx1, testcy1, testcz1;
    double[][] testw1;

    int maxx = 10;
    double background = 50;
    double angle = 0;
    double width = 5;

    @Test
    public void fitFixedComputesGradientPerDatum()
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_FIXED);
    }

    @Test
    public void fitCircleComputesGradientPerDatum()
    {
        functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_CIRCLE);
    }

    @Test
    public void fitFreeCircleComputesGradientPerDatum()
    {
        functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_FREE_CIRCLE);
    }

    @Test
    public void fitEllipticalComputesGradientPerDatum()
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_ELLIPTICAL);
    }

    @Test
    public void fitNBFixedComputesGradientPerDatum()
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED);
    }

    @Test
    public void fitNBCircleComputesGradientPerDatum()
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_SIMPLE_NB_CIRCLE);
    }

    @Test
    public void fitNBFreeCircleComputesGradientPerDatum()
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
    }

    @Test
    public void fitNBEllipticalComputesGradientPerDatum()
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
    }

    private void functionComputesGradientPerDatum(int flags)
    {
        final Gaussian2DFunction f1 = GaussianFunctionFactory.create2D(1, maxx, maxx, flags, null);
        // Setup
        testbackground = testbackground_;
        testsignal1 = testsignal1_;
        testcx1 = testcx1_;
        testcy1 = testcy1_;
        testcz1 = testcz1_;
        testw1 = testw1_;
        testangle1 = testangle1_;
        if (!f1.evaluatesBackground())
            testbackground = new double[] { testbackground[0] };
        if (!f1.evaluatesSignal())
            testsignal1 = new double[] { testsignal1[0] };

        if (!f1.evaluatesZ())
            testcz1 = new double[] { 0 };
        boolean noSecondWidth = false;
        if (!f1.evaluatesSD0())
        {
            // Just use 1 width
            testw1 = new double[][] { testw1[0] };
            // If no width 0 then assume we have no width 1 as well
            noSecondWidth = true;
        }
        else if (!f1.evaluatesSD1())
        {
            // No evaluation of second width needs only variation in width 0 so truncate
            testw1 = Arrays.copyOf(testw1, 2);
            noSecondWidth = true;
        }
        if (noSecondWidth)
            for (int i = 0; i < testw1.length; i++)
                testw1[i][1] = testw1[i][0];
        if (!f1.evaluatesAngle())
            testangle1 = new double[] { 0 };

        if (f1.evaluatesBackground())
            functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.BACKGROUND);
        if (f1.evaluatesSignal())
            functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.SIGNAL);
        functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.X_POSITION);
        functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.Y_POSITION);
        if (f1.evaluatesZ())
            functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.Z_POSITION);
        if (f1.evaluatesSD0())
            functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.X_SD);
        if (f1.evaluatesSD1())
            functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.Y_SD);
        if (f1.evaluatesAngle())
            functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.ANGLE);
    }

    private void functionComputesTargetGradientPerDatum(Gaussian2DFunction f1, int targetParameter)
    {
        final int[] indices = f1.gradientIndices();
        final int gradientIndex = findGradientIndex(f1, targetParameter);
        final double[] dyda = new double[indices.length];
        double[] a;

        PoissonLikelihoodWrapper ff1;

        final int n = maxx * maxx;
        int count = 0, total = 0;

        for (final double background : testbackground)
            for (final double signal1 : testsignal1)
                for (final double cx1 : testcx1)
                    for (final double cy1 : testcy1)
                        for (final double cz1 : testcz1)
                            for (final double[] w1 : testw1)
                                for (final double angle1 : testangle1)
                                {
                                    a = createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1);

                                    // Create y as a function we would want to move towards
                                    final double[] a2 = a.clone();
                                    a2[targetParameter] *= 1.1;
                                    f1.initialise(a2);
                                    final double[] data = new double[maxx * maxx];
                                    for (int i = 0; i < n; i++)
                                        data[i] = f1.eval(i);

                                    ff1 = new PoissonLikelihoodWrapper(f1, a, data, n, alpha);

                                    // Numerically solve gradient.
                                    // Calculate the step size h to be an exact numerical representation
                                    final double xx = a[targetParameter];

                                    // Get h to minimise roundoff error
                                    final double h = Precision.representableDelta(xx, h_);

                                    for (final int x : testx)
                                        for (final int y : testy)
                                        {
                                            final int i = y * maxx + x;
                                            a[targetParameter] = xx;
                                            ff1.likelihood(getVariables(indices, a), dyda, i);

                                            // Evaluate at (x+h) and (x-h)
                                            a[targetParameter] = xx + h;
                                            final double value2 = ff1.likelihood(getVariables(indices, a), i);

                                            a[targetParameter] = xx - h;
                                            final double value3 = ff1.likelihood(getVariables(indices, a), i);

                                            final double gradient = (value2 - value3) / (2 * h);
                                            boolean ok = Math.signum(gradient) == Math.signum(dyda[gradientIndex]) ||
                                                    Math.abs(gradient - dyda[gradientIndex]) < 0.1;
                                            //logger.fine(FunctionUtils.getSupplier("[%s-%s]/2*%g : %g == %g", "" + value2, "" + value3, h, gradient,
                                            //		dyda[gradientIndex]));
                                            if (!ok)
                                                Assertions.fail(NAME[targetParameter] + ": " + gradient + " != " +
                                                        dyda[gradientIndex]);
                                            ok = eqPerDatum.almostEqualRelativeOrAbsolute(gradient,
                                                    dyda[gradientIndex]);
                                            if (ok)
                                                count++;
                                            total++;
                                        }
                                }
        final double p = (100.0 * count) / total;
        logger.log(TestLog.getRecord(Level.INFO, "Per Datum %s : %s = %d / %d (%.2f)", f1.getClass().getSimpleName(), NAME[targetParameter],
                count, total, p));
        Assertions.assertTrue(p > 90, () -> NAME[targetParameter] + " fraction too low per datum: " + p);
    }

    @Test
    public void fitFixedComputesGradient()
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        functionComputesGradient(GaussianFunctionFactory.FIT_FIXED);
    }

    @Test
    public void fitCircleComputesGradient()
    {
        functionComputesGradient(GaussianFunctionFactory.FIT_CIRCLE);
    }

    @Test
    public void fitFreeCircleComputesGradient()
    {
        functionComputesGradient(GaussianFunctionFactory.FIT_FREE_CIRCLE);
    }

    @Test
    public void fitEllipticalComputesGradient()
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        // The elliptical function gradient evaluation is worse
        final DoubleEquality tmp = eq;
        eq = eqPerDatum;
        functionComputesGradient(GaussianFunctionFactory.FIT_ELLIPTICAL);
        eq = tmp;
    }

    @Test
    public void fitNBFixedComputesGradient()
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        functionComputesGradient(GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED);
    }

    @Test
    public void fitNBCircleComputesGradient()
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        functionComputesGradient(GaussianFunctionFactory.FIT_SIMPLE_NB_CIRCLE);
    }

    @Test
    public void fitNBFreeCircleComputesGradient()
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        functionComputesGradient(GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
    }

    @Test
    public void fitNBEllipticalComputesGradient()
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        // The elliptical function gradient evaluation is worse
        final DoubleEquality tmp = eq;
        eq = eqPerDatum;
        functionComputesGradient(GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
        eq = tmp;
    }

    private void functionComputesGradient(int flags)
    {
        final Gaussian2DFunction f1 = GaussianFunctionFactory.create2D(1, maxx, maxx, flags, null);
        testbackground = testbackground_;
        testsignal1 = testsignal1_;
        testcx1 = testcx1_;
        testcy1 = testcy1_;
        testcz1 = testcz1_;
        testw1 = testw1_;
        testangle1 = testangle1_;
        if (!f1.evaluatesBackground())
            testbackground = new double[] { testbackground[0] };
        if (!f1.evaluatesSignal())
            testsignal1 = new double[] { testsignal1[0] };

        if (!f1.evaluatesZ())
            testcz1 = new double[] { 0 };
        boolean noSecondWidth = false;
        if (!f1.evaluatesSD0())
        {
            // Just use 1 width
            testw1 = new double[][] { testw1[0] };
            // If no width 0 then assume we have no width 1 as well
            noSecondWidth = true;
        }
        else if (!f1.evaluatesSD1())
        {
            // No evaluation of second width needs only variation in width 0 so truncate
            testw1 = Arrays.copyOf(testw1, 2);
            noSecondWidth = true;
        }
        if (noSecondWidth)
            for (int i = 0; i < testw1.length; i++)
                testw1[i][1] = testw1[i][0];
        if (!f1.evaluatesAngle())
            testangle1 = new double[] { 0 };

        final double fraction = 90;
        if (f1.evaluatesBackground())
            functionComputesTargetGradient(f1, Gaussian2DFunction.BACKGROUND, fraction);
        if (f1.evaluatesSignal())
            functionComputesTargetGradient(f1, Gaussian2DFunction.SIGNAL, fraction);
        functionComputesTargetGradient(f1, Gaussian2DFunction.X_POSITION, fraction);
        functionComputesTargetGradient(f1, Gaussian2DFunction.Y_POSITION, fraction);
        if (f1.evaluatesZ())
            functionComputesTargetGradient(f1, Gaussian2DFunction.Z_POSITION, fraction);
        if (f1.evaluatesSD0())
            functionComputesTargetGradient(f1, Gaussian2DFunction.X_SD, fraction);
        if (f1.evaluatesSD1())
            functionComputesTargetGradient(f1, Gaussian2DFunction.Y_SD, fraction);
        if (f1.evaluatesAngle())
            functionComputesTargetGradient(f1, Gaussian2DFunction.ANGLE, fraction);
    }

    private void functionComputesTargetGradient(Gaussian2DFunction f1, int targetParameter, double threshold)
    {
        final int[] indices = f1.gradientIndices();
        final int gradientIndex = findGradientIndex(f1, targetParameter);
        final double[] dyda = new double[indices.length];
        double[] a;

        PoissonLikelihoodWrapper ff1;

        final int n = maxx * maxx;
        int count = 0, total = 0;

        for (final double background : testbackground)
            for (final double signal1 : testsignal1)
                for (final double cx1 : testcx1)
                    for (final double cy1 : testcy1)
                        for (final double cz1 : testcz1)
                            for (final double[] w1 : testw1)
                                for (final double angle1 : testangle1)
                                {
                                    a = createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1);

                                    // Create y as a function we would want to move towards
                                    final double[] a2 = a.clone();
                                    a2[targetParameter] *= 1.3;
                                    f1.initialise(a2);
                                    final double[] data = new double[maxx * maxx];
                                    for (int i = 0; i < n; i++)
                                        data[i] = f1.eval(i);

                                    ff1 = new PoissonLikelihoodWrapper(f1, a, data, n, alpha);

                                    // Numerically solve gradient.
                                    // Calculate the step size h to be an exact numerical representation
                                    final double xx = a[targetParameter];

                                    // Get h to minimise roundoff error
                                    final double h = Precision.representableDelta(xx, h_);

                                    ff1.likelihood(getVariables(indices, a), dyda);

                                    // Evaluate at (x+h) and (x-h)
                                    a[targetParameter] = xx + h;
                                    final double value2 = ff1.likelihood(getVariables(indices, a));

                                    a[targetParameter] = xx - h;
                                    final double value3 = ff1.likelihood(getVariables(indices, a));

                                    final double gradient = (value2 - value3) / (2 * h);
                                    boolean ok = Math.signum(gradient) == Math.signum(dyda[gradientIndex]) ||
                                            Math.abs(gradient - dyda[gradientIndex]) < 0.1;
                                    //logger.fine(FunctionUtils.getSupplier("[%s-%s]/2*%g : %g == %g", "" + value2, "" + value3, h, gradient,
                                    //		dyda[gradientIndex]));
                                    if (!ok)
                                        Assertions.fail(
                                                NAME[targetParameter] + ": " + gradient + " != " + dyda[gradientIndex]);
                                    ok = eq.almostEqualRelativeOrAbsolute(gradient, dyda[gradientIndex]);
                                    if (ok)
                                        count++;
                                    total++;

                                }
        final double p = (100.0 * count) / total;
        logger.log(TestLog.getRecord(Level.INFO, "%s : %s = %d / %d (%.2f)", f1.getClass().getSimpleName(), NAME[targetParameter], count,
                total, p));
        Assertions.assertTrue(p > threshold, () -> NAME[targetParameter] + " fraction too low: " + p);
    }

    private static double[] getVariables(int[] indices, double[] a)
    {
        final double[] variables = new double[indices.length];
        for (int i = 0; i < indices.length; i++)
            variables[i] = a[indices[i]];
        return variables;
    }

    private static int findGradientIndex(Gaussian2DFunction f, int targetParameter)
    {
        final int i = f.findGradientIndex(targetParameter);
        Assertions.assertTrue(i >= 0, "Cannot find gradient index");
        return i;
    }

    double[] createParameters(double... args)
    {
        return args;
    }

    @Test
    public void cumulativeProbabilityIsOneWithIntegerData()
    {
        // Initialise for large observed count
        PoissonLikelihoodWrapper.likelihood(1, photons[photons.length - 1] * 2);

        for (final double p : photons)
            cumulativeProbabilityIsOneWithIntegerData(p);
    }

    private static void cumulativeProbabilityIsOneWithIntegerData(final double mu)
    {
        double p = 0;
        int x = 0;

        // Evaluate an initial range.
        // Poisson will have mean mu with a variance mu.
        // At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
        if (mu > 0)
        {
            final int max = (int) Math.ceil(mu + 3 * Math.sqrt(mu));
            for (; x <= max; x++)
            {
                final double pp = PoissonLikelihoodWrapper.likelihood(mu, x);
                //logger.fine(FunctionUtils.getSupplier("x=%d, p=%f", x, pp);
                p += pp;
            }
            if (p > 1.01)
                Assertions.fail("P > 1: " + p);
        }

        // We have most of the probability density.
        // Now keep evaluating up until no difference
        final double changeTolerance = 1e-6;
        for (;; x++)
        {
            final double pp = PoissonLikelihoodWrapper.likelihood(mu, x);
            //logger.fine(FunctionUtils.getSupplier("x=%d, p=%f", x, pp);
            p += pp;
            if (pp / p < changeTolerance)
                break;
        }
        logger.log(TestLog.getRecord(Level.INFO, "mu=%f, p=%f, max=%d", mu, p, x));
        Assertions.assertEquals(1, p, 0.02, () -> String.format("mu=%f", mu));
    }

    @Test
    public void cumulativeProbabilityIsOneWithRealDataForCountAbove4()
    {
        //for (int i = photons.length; i-- > 0;)
        for (int i = 0; i < photons.length; i++)
            if (photons[i] >= 4)
                cumulativeProbabilityIsOneWithRealData(photons[i], maxRange[i] + 1);
    }

    private static void cumulativeProbabilityIsOneWithRealData(final double mu, int max)
    {
        double p = 0;

        final SimpsonIntegrator in = new SimpsonIntegrator();

        p = in.integrate(20000, new UnivariateFunction()
        {
            @Override
            public double value(double x)
            {
                return PoissonLikelihoodWrapper.likelihood(mu, x);
            }
        }, 0, max);

        logger.log(TestLog.getRecord(Level.INFO, "mu=%f, p=%f", mu, p));
        Assertions.assertEquals(1, p, 0.02, () -> String.format("mu=%f", mu));
    }

    @Test
    public void cumulativeProbabilityIsOneFromLikelihoodForCountAbove4()
    {
        for (int i = 0; i < photons.length; i++)
            if (photons[i] >= 4)
                cumulativeProbabilityIsOneFromLikelihood(photons[i]);
    }

    private void cumulativeProbabilityIsOneFromLikelihood(final double mu)
    {
        // Determine upper limit for a Poisson
        final int limit = new PoissonDistribution(mu).inverseCumulativeProbability(0.999);

        // Expand to allow for the gain
        final int n = (int) Math.ceil(limit / alpha);

        // Evaluate all values from zero to large n
        final double[] k = SimpleArrayUtils.newArray(n, 0, 1.0);
        final double[] a = new double[0];
        final double[] g = new double[0];

        final NonLinearFunction nlf = new NonLinearFunction()
        {
            @Override
            public void initialise(double[] a)
            {
                // Do nothing
            }

            @Override
            public int[] gradientIndices()
            {
                return new int[0];
            }

            @Override
            public double eval(int x, double[] dyda, double[] w)
            {
                return 0;
            }

            @Override
            public double eval(int x)
            {
                return mu;
            }

            @Override
            public double eval(int x, double[] dyda)
            {
                return mu;
            }

            @Override
            public boolean canComputeWeights()
            {
                return false;
            }

            @Override
            public double evalw(int x, double[] w)
            {
                return 0;
            }

            @Override
            public int getNumberOfGradients()
            {
                return 0;
            }
        };
        PoissonLikelihoodWrapper f = new PoissonLikelihoodWrapper(nlf, a, Arrays.copyOf(k, n), n, alpha);

        // Keep evaluating up until no difference
        final double changeTolerance = 1e-6;
        double total = 0;
        double p = 0;
        int i = 0;
        while (i < n)
        {
            final double nll = f.computeLikelihood(i);
            final double nll2 = f.computeLikelihood(g, i);
            i++;
            Assertions.assertEquals(nll, nll2, 1e-10, "computeLikelihood(i)");
            total += nll;
            final double pp = FastMath.exp(-nll);
            //logger.fine(FunctionUtils.getSupplier("mu=%f, o=%f, i=%d, pp=%f", mu, mu / alpha, i, pp);
            p += pp;
            if (p > 0.5 && pp / p < changeTolerance)
                break;
        }

        logger.log(TestLog.getRecord(Level.INFO, "mu=%f, limit=%d, p=%f", mu, limit, p));
        Assertions.assertEquals(1, p, 0.02, () -> String.format("mu=%f", mu));

        // Check the function can compute the same total
        f = new PoissonLikelihoodWrapper(nlf, a, k, i, alpha);
        final double sum = f.computeLikelihood();
        final double sum2 = f.computeLikelihood(g);
        final double delta = 1e-10;
        ExtraAssertions.assertEqualsRelative(total, sum, delta, "computeLikelihood");
        ExtraAssertions.assertEqualsRelative(total, sum2, delta, "computeLikelihood with gradient");
    }
}
