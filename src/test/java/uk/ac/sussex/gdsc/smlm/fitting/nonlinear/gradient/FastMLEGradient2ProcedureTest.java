package uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.math3.util.Precision;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.RandomGeneratorAdapter;
import uk.ac.sussex.gdsc.smlm.function.DummyGradientFunction;
import uk.ac.sussex.gdsc.smlm.function.FakeGradientFunction;
import uk.ac.sussex.gdsc.smlm.function.Gradient2Function;
import uk.ac.sussex.gdsc.smlm.function.OffsetGradient2Function;
import uk.ac.sussex.gdsc.smlm.function.ValueProcedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.function.gaussian.HoltzerAstigmatismZModel;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.SingleAstigmatismErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.math3.distribution.CustomPoissonDistribution;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssumptions;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RNGFactory;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestCounter;
import uk.ac.sussex.gdsc.test.utils.TestLog;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;
import uk.ac.sussex.gdsc.test.utils.functions.IndexSupplier;
import uk.ac.sussex.gdsc.test.utils.functions.IntArrayFormatSupplier;

/**
 * Contains speed tests for the methods for calculating the Hessian and gradient vector
 * for use in the LVM algorithm.
 */
@SuppressWarnings({ "javadoc" })
public class FastMLEGradient2ProcedureTest
{
    private static Logger logger;

    @BeforeAll
    public static void beforeAll()
    {
        logger = Logger.getLogger(FastMLEGradient2ProcedureTest.class.getName());
    }

    @AfterAll
    public static void afterAll()
    {
        logger = null;
    }

    DoubleEquality eq = new DoubleEquality(1e-6, 1e-16);

    int MAX_ITER = 20000;
    int blockWidth = 10;
    double background = 0.5;
    double signal = 100;
    double angle = Math.PI;
    double xpos = 5;
    double ypos = 5;
    double xwidth = 1.2;
    double ywidth = 1.2;

    private static double nextUniform(UniformRandomProvider r, double min, double max)
    {
        return min + r.nextDouble() * (max - min);
    }

    private static double random(UniformRandomProvider r, double d)
    {
        return d - d * 0.1 + r.nextDouble() * 0.2;
    }

    @Test
    public void gradientProcedureFactoryCreatesOptimisedProcedures()
    {
        final double[] y = new double[0];
        Assertions.assertEquals(
                FastMLEGradient2ProcedureFactory.createUnrolled(y, new DummyGradientFunction(4)).getClass(),
                FastMLEGradient2Procedure4.class);
        Assertions.assertEquals(
                FastMLEGradient2ProcedureFactory.createUnrolled(y, new DummyGradientFunction(5)).getClass(),
                FastMLEGradient2Procedure5.class);
        Assertions.assertEquals(
                FastMLEGradient2ProcedureFactory.createUnrolled(y, new DummyGradientFunction(6)).getClass(),
                FastMLEGradient2Procedure6.class);
    }

    @SeededTest
    public void gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator(RandomSeed seed)
    {
        gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator(seed, 4);
        gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator(seed, 5);
        gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator(seed, 6);
        gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator(seed, 11);
        gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator(seed, 21);
    }

    @SpeedTag
    @SeededTest
    public void gradientProcedureIsNotSlowerThanGradientCalculator(RandomSeed seed)
    {
        // Note: The procedure does not have a lot of work within loops. It is only a single loop
        // so unrolling does not produce performance gains. The JVM can optimise this.

        gradientProcedureIsNotSlowerThanGradientCalculator(seed, 4);
        gradientProcedureIsNotSlowerThanGradientCalculator(seed, 5);
        gradientProcedureIsNotSlowerThanGradientCalculator(seed, 6);
        gradientProcedureIsNotSlowerThanGradientCalculator(seed, 11);
        gradientProcedureIsNotSlowerThanGradientCalculator(seed, 21);
    }

    private void gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator(RandomSeed seed, int nparams)
    {
        final int iter = 10;

        final ArrayList<double[]> paramsList = new ArrayList<>(iter);
        final ArrayList<double[]> yList = new ArrayList<>(iter);

        createFakeData(RNGFactory.create(seed.getSeed()), nparams, iter, paramsList, yList);
        final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

        final MLEGradientCalculator calc = (MLEGradientCalculator) GradientCalculatorFactory.newCalculator(nparams,
                true);

        final IndexSupplier msg = new IndexSupplier(1, "["+nparams+"] Result: not same @ ", null);
        for (int i = 0; i < paramsList.size(); i++)
        {
            final FastMLEGradient2Procedure p = FastMLEGradient2ProcedureFactory.createUnrolled(yList.get(i), func);
            final double s = p.computeLogLikelihood(paramsList.get(i));
            final double s2 = calc.logLikelihood(yList.get(i), paramsList.get(i), func);
            // Virtually the same ...
            ExtraAssertions.assertEqualsRelative(s, s2, 1e-5, msg.set(0, i));
        }
    }

    @SeededTest
    public void gradientProcedureComputesSameWithPrecomputed(RandomSeed seed)
    {
        final int iter = 10;
        final UniformRandomProvider r = RNGFactory.create(seed.getSeed());

        final ErfGaussian2DFunction f1 = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, 10, 10,
                GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
        final ErfGaussian2DFunction f2 = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(2, 10, 10,
                GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);

        final double[] a1 = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
        final double[] a2 = new double[1 + 2 * Gaussian2DFunction.PARAMETERS_PER_PEAK];

        final double[] x = new double[f1.size()];
        final double[] b = new double[f1.size()];

        final CustomPoissonDistribution pd = new CustomPoissonDistribution(new RandomGeneratorAdapter(r), 1);

        for (int i = 0; i < iter; i++)
        {
            final int ii = i;

            a2[Gaussian2DFunction.BACKGROUND] = nextUniform(r, 0.1, 0.3);
            a2[Gaussian2DFunction.SIGNAL] = nextUniform(r, 100, 300);
            a2[Gaussian2DFunction.X_POSITION] = nextUniform(r, 3, 5);
            a2[Gaussian2DFunction.Y_POSITION] = nextUniform(r, 3, 5);
            a2[Gaussian2DFunction.Z_POSITION] = nextUniform(r, -2, 2);
            a2[Gaussian2DFunction.X_SD] = nextUniform(r, 1, 1.3);
            a2[Gaussian2DFunction.Y_SD] = nextUniform(r, 1, 1.3);
            a2[Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.SIGNAL] = nextUniform(r, 100, 300);
            a2[Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_POSITION] = nextUniform(r, 5, 7);
            a2[Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Y_POSITION] = nextUniform(r, 5, 7);
            a2[Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Z_POSITION] = nextUniform(r, -3, 1);
            a2[Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_SD] = nextUniform(r, 1, 1.3);
            a2[Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Y_SD] = nextUniform(r, 1, 1.3);

            // Simulate Poisson data
            f2.initialise0(a2);
            f1.forEach(new ValueProcedure()
            {
                int k = 0;

                @Override
                public void execute(double value)
                {
                    if (value > 0)
                    {
                        pd.setMeanUnsafe(value);
                        x[k++] = pd.sample();
                    }
                    else
                        x[k++] = 0;
                }
            });

            // Precompute peak 2 (no background)
            a1[Gaussian2DFunction.BACKGROUND] = 0;
            for (int j = 1; j < 7; j++)
                a1[j] = a2[Gaussian2DFunction.PARAMETERS_PER_PEAK + j];
            f1.initialise0(a1);
            f1.forEach(new ValueProcedure()
            {
                int k = 0;

                @Override
                public void execute(double value)
                {
                    b[k++] = value;
                }
            });

            // Reset to peak 1
            for (int j = 0; j < 7; j++)
                a1[j] = a2[j];

            // Compute peak 1+2
            final FastMLEGradient2Procedure p12 = FastMLEGradient2ProcedureFactory.create(x, f2);
            p12.computeSecondDerivative(a2);
            final double[] d11 = Arrays.copyOf(p12.d1, f1.getNumberOfGradients());
            final double[] d21 = Arrays.copyOf(p12.d2, f1.getNumberOfGradients());

            // Compute peak 1+(precomputed 2)
            final FastMLEGradient2Procedure p1b2 = FastMLEGradient2ProcedureFactory.create(x,
                    OffsetGradient2Function.wrapGradient2Function(f1, b));
            p1b2.computeSecondDerivative(a1);
            final double[] d12 = p1b2.d1;
            final double[] d22 = p1b2.d2;

            Assertions.assertArrayEquals(p12.u, p1b2.u, 1e-10, () -> " Result: Not same @ " + ii);
            Assertions.assertArrayEquals(d11, d12, 1e-10, () -> " D1: Not same @ " + ii);
            Assertions.assertArrayEquals(d21, d22, 1e-10, () -> " D2: Not same @ " + ii);

            final double[] v1 = p12.computeValue(a2);
            final double[] v2 = p1b2.computeValue(a1);
            Assertions.assertArrayEquals(v1, v2, 1e-10, () -> " Value: Not same @ " + ii);

            final double[] d1 = Arrays.copyOf(p12.computeFirstDerivative(a2), f1.getNumberOfGradients());
            final double[] d2 = p1b2.computeFirstDerivative(a1);
            Assertions.assertArrayEquals(d1, d2, 1e-10, () -> " 1st derivative: Not same @ " + ii);
        }
    }

    private abstract class Timer
    {
        private int loops;
        int min = 5;

        Timer()
        {
        }

        Timer(int min)
        {
            this.min = min;
        }

        long getTime()
        {
            // Run till stable timing
            long t1 = time();
            for (int i = 0; i < 10; i++)
            {
                final long t2 = t1;
                t1 = time();
                if (loops >= min && DoubleEquality.relativeError(t1, t2) < 0.02) // 2% difference
                    break;
            }
            return t1;
        }

        long time()
        {
            loops++;
            long t = System.nanoTime();
            run();
            t = System.nanoTime() - t;
            return t;
        }

        abstract void run();
    }

    private void gradientProcedureIsNotSlowerThanGradientCalculator(RandomSeed seed, final int nparams)
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);

        final int iter = 1000;

        final ArrayList<double[]> paramsList = new ArrayList<>(iter);
        final ArrayList<double[]> yList = new ArrayList<>(iter);

        createFakeData(RNGFactory.create(seed.getSeed()), nparams, iter, paramsList, yList);
        final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

        final MLEGradientCalculator calc = (MLEGradientCalculator) GradientCalculatorFactory.newCalculator(nparams,
                true);

        for (int i = 0; i < paramsList.size(); i++)
            calc.logLikelihood(yList.get(i), paramsList.get(i), func);

        for (int i = 0; i < paramsList.size(); i++)
        {
            final FastMLEGradient2Procedure p = FastMLEGradient2ProcedureFactory.createUnrolled(yList.get(i), func);
            p.computeLogLikelihood(paramsList.get(i));
        }

        // Realistic loops for an optimisation
        final int loops = 15;

        // Run till stable timing
        final Timer t1 = new Timer()
        {
            @Override
            void run()
            {
                for (int i = 0, k = 0; i < iter; i++)
                {
                    final MLEGradientCalculator calc = (MLEGradientCalculator) GradientCalculatorFactory
                            .newCalculator(nparams, true);
                    for (int j = loops; j-- > 0;)
                        calc.logLikelihood(yList.get(i), paramsList.get(k++ % iter), func);
                }
            }
        };
        final long time1 = t1.getTime();

        final Timer t2 = new Timer(t1.loops)
        {
            @Override
            void run()
            {
                for (int i = 0, k = 0; i < iter; i++)
                {
                    final FastMLEGradient2Procedure p = FastMLEGradient2ProcedureFactory.createUnrolled(yList.get(i),
                            func);
                    for (int j = loops; j-- > 0;)
                        p.computeLogLikelihood(paramsList.get(k++ % iter));
                }
            }
        };
        final long time2 = t2.getTime();

        logger.log(TestLog.getTimingRecord("GradientCalculator " + nparams, time1, "FastMLEGradient2Procedure", time2));
    }

    @SeededTest
    public void gradientProcedureUnrolledComputesSameAsGradientProcedure(RandomSeed seed)
    {
        gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 4);
        gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 5);
        gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 6);
    }

    private void gradientProcedureUnrolledComputesSameAsGradientProcedure(RandomSeed seed, int nparams)
    {
        final int iter = 10;

        final ArrayList<double[]> paramsList = new ArrayList<>(iter);
        final ArrayList<double[]> yList = new ArrayList<>(iter);

        createFakeData(RNGFactory.create(seed.getSeed()), nparams, iter, paramsList, yList);
        final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

        // Create messages
        final IntArrayFormatSupplier msgLL = getMessage(nparams, "[%d] LL: Not same @ %d");
        final IntArrayFormatSupplier msg1Dv = getMessage(nparams, "[%d] first derivative: Not same value @ %d");
        final IntArrayFormatSupplier msg1Dd1 = getMessage(nparams, "[%d] first derivative: Not same d1 @ %d");
        final IntArrayFormatSupplier msg2Dv = getMessage(nparams, "[%d] second derivative: Not same value @ %d");
        final IntArrayFormatSupplier msg2Dd1 = getMessage(nparams, "[%d] second derivative: Not same d1 @ %d");
        final IntArrayFormatSupplier msg2Dd2 = getMessage(nparams, "[%d] second derivative: Not same d2 @ %d");
        
        FastMLEGradient2Procedure p1, p2;
        for (int i = 0; i < paramsList.size(); i++)
        {
            p1 = new FastMLEGradient2Procedure(yList.get(i), func);
            p2 = FastMLEGradient2ProcedureFactory.createUnrolled(yList.get(i), func);
            final double[] a = paramsList.get(i);

            final double ll1 = p1.computeLogLikelihood(a);
            final double ll2 = p2.computeLogLikelihood(a);
            Assertions.assertEquals(ll1, ll2, msgLL.set(1, i));

            p1 = new FastMLEGradient2Procedure(yList.get(i), func);
            p2 = FastMLEGradient2ProcedureFactory.createUnrolled(yList.get(i), func);
            p1.computeFirstDerivative(a);
            p2.computeFirstDerivative(a);
            Assertions.assertArrayEquals(p1.u, p2.u, msg1Dv.set(1, i));
            Assertions.assertArrayEquals(p1.d1, p2.d1, msg1Dd1.set(1, i));

            p1 = new FastMLEGradient2Procedure(yList.get(i), func);
            p2 = FastMLEGradient2ProcedureFactory.createUnrolled(yList.get(i), func);
            p1.computeSecondDerivative(a);
            p2.computeSecondDerivative(a);
            Assertions.assertArrayEquals(p1.u, p2.u, msg2Dv.set(1, i));
            Assertions.assertArrayEquals(p1.d1, p2.d1, msg2Dd1.set(1, i));
            Assertions.assertArrayEquals(p1.d2, p2.d2, msg2Dd2.set(1, i));
        }
    }
    
    private static IntArrayFormatSupplier getMessage(int nparams, String format) {
        final IntArrayFormatSupplier msg = new IntArrayFormatSupplier(format, 2);
        msg.set(0, nparams);
        return msg;
    }

    @SpeedTag
    @SeededTest
    public void gradientProcedureIsFasterUnrolledThanGradientProcedure(RandomSeed seed)
    {
        gradientProcedureLinearIsFasterThanGradientProcedure(seed, 4);
        gradientProcedureLinearIsFasterThanGradientProcedure(seed, 5);
        gradientProcedureLinearIsFasterThanGradientProcedure(seed, 6);
    }

    private void gradientProcedureLinearIsFasterThanGradientProcedure(RandomSeed seed, final int nparams)
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);

        final int iter = 100;

        final ArrayList<double[]> paramsList = new ArrayList<>(iter);
        final ArrayList<double[]> yList = new ArrayList<>(iter);

        createData(RNGFactory.create(seed.getSeed()), 1, iter, paramsList, yList);

        // Remove the timing of the function call by creating a dummy function
        final Gradient2Function func = new FakeGradientFunction(blockWidth, nparams);

        for (int i = 0; i < paramsList.size(); i++)
        {
            final FastMLEGradient2Procedure p1 = new FastMLEGradient2Procedure(yList.get(i), func);
            p1.computeSecondDerivative(paramsList.get(i));
            p1.computeSecondDerivative(paramsList.get(i));

            final FastMLEGradient2Procedure p2 = FastMLEGradient2ProcedureFactory.createUnrolled(yList.get(i), func);
            p2.computeSecondDerivative(paramsList.get(i));
            p2.computeSecondDerivative(paramsList.get(i));

            // Check they are the same
            final int ii = i;
            Assertions.assertArrayEquals(p1.d1, p2.d1, () -> "D1 " + ii);
            Assertions.assertArrayEquals(p1.d2, p2.d2, () -> "D2 " + ii);
        }

        // Realistic loops for an optimisation
        final int loops = 15;

        // Run till stable timing
        final Timer t1 = new Timer()
        {
            @Override
            void run()
            {
                for (int i = 0, k = 0; i < paramsList.size(); i++)
                {
                    final FastMLEGradient2Procedure p1 = new FastMLEGradient2Procedure(yList.get(i), func);
                    for (int j = loops; j-- > 0;)
                        p1.computeSecondDerivative(paramsList.get(k++ % iter));
                }
            }
        };
        final long time1 = t1.getTime();

        final Timer t2 = new Timer(t1.loops)
        {
            @Override
            void run()
            {
                for (int i = 0, k = 0; i < paramsList.size(); i++)
                {
                    final FastMLEGradient2Procedure p2 = FastMLEGradient2ProcedureFactory.createUnrolled(yList.get(i),
                            func);
                    for (int j = loops; j-- > 0;)
                        p2.computeSecondDerivative(paramsList.get(k++ % iter));
                }
            }
        };
        final long time2 = t2.getTime();

        logger.log(TestLog.getRecord(Level.INFO, "Standard = %d : Unrolled %d = %d : %fx", time1, nparams, time2,
                (1.0 * time1) / time2));
        Assertions.assertTrue(time2 < time1 * 1.5);
    }

    @SeededTest
    public void gradientCalculatorComputesGradient(RandomSeed seed)
    {
        gradientCalculatorComputesGradient(seed, new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth));

        // Use a reasonable z-depth function from the Smith, et al (2010) paper (page 377)
        final double sx = 1.08;
        final double sy = 1.01;
        final double gamma = 0.389;
        final double d = 0.531;
        final double Ax = -0.0708;
        final double Bx = -0.073;
        final double Ay = 0.164;
        final double By = 0.0417;
        final HoltzerAstigmatismZModel zModel = HoltzerAstigmatismZModel.create(sx, sy, gamma, d, Ax, Bx, Ay, By);
        gradientCalculatorComputesGradient(seed,
                new SingleAstigmatismErfGaussian2DFunction(blockWidth, blockWidth, zModel));
    }

    private void gradientCalculatorComputesGradient(RandomSeed seed, ErfGaussian2DFunction func)
    {
        // Check the first and second derivatives
        final int nparams = func.getNumberOfGradients();
        final int[] indices = func.gradientIndices();

        final int iter = 100;

        final ArrayList<double[]> paramsList = new ArrayList<>(iter);
        final ArrayList<double[]> yList = new ArrayList<>(iter);

        createData(RNGFactory.create(seed.getSeed()), 1, iter, paramsList, yList, true);

        // for the gradients
        final double delta = 1e-4;
        final DoubleEquality eq = new DoubleEquality(5e-2, 1e-16);

        // Must compute most of the time
        final int failureLimit = TestCounter.computeFailureLimit(iter, 0.1);
        final TestCounter failCounter = new TestCounter(failureLimit, nparams);

        for (int i = 0; i < paramsList.size(); i++)
        {
            final int ii = i;
            final double[] y = yList.get(i);
            final double[] a = paramsList.get(i);
            final double[] a2 = a.clone();
            final FastMLEGradient2Procedure p = FastMLEGradient2ProcedureFactory.create(y, func);
            //double ll = p.computeLogLikelihood(a);
            p.computeSecondDerivative(a);
            final double[] d1 = p.d1.clone();
            final double[] d2 = p.d2.clone();
            for (int j = 0; j < nparams; j++)
            {
                final int j_ = j;
                final int k = indices[j];
                final double d = Precision.representableDelta(a[k], (a[k] == 0) ? delta : a[k] * delta);
                a2[k] = a[k] + d;
                final double llh = p.computeLogLikelihood(a2);
                p.computeFirstDerivative(a2);
                final double[] d1h = p.d1.clone();
                a2[k] = a[k] - d;
                final double lll = p.computeLogLikelihood(a2);
                p.computeFirstDerivative(a2);
                final double[] d1l = p.d1.clone();
                a2[k] = a[k];

                final double gradient1 = (llh - lll) / (2 * d);
                final double gradient2 = (d1h[j] - d1l[j]) / (2 * d);
                //logger.fine("[%d,%d] ll - %f  (%s %f+/-%f) d1 %f ?= %f : d2 %f ?= %f", i, k, ll, func.getName(k), a[k], d,
                //		gradient1, d1[j], gradient2, d2[j]);
                failCounter.run(j, () -> {
                    return eq.almostEqualRelativeOrAbsolute(gradient1, d1[j_]);
                }, () -> {
                    Assertions.fail(FunctionUtils.getSupplier("Not same gradient1 @ %d,%d: %s != %s (error=%s)", ii, j_, gradient1, d1[j_],
                            DoubleEquality.relativeError(gradient1, d1[j_])));
                });
                failCounter.run(nparams + j, () -> {
                    return eq.almostEqualRelativeOrAbsolute(gradient2, d2[j_]);
                }, () -> {
                    Assertions.fail(FunctionUtils.getSupplier("Not same gradient2 @ %d,%d: %s != %s (error=%s)", ii, j_, gradient2, d2[j_],
                            DoubleEquality.relativeError(gradient2, d2[j_])));
                });
            }
        }
    }

    /**
     * Create random elliptical Gaussian data an returns the data plus an estimate of the parameters.
     * Only the chosen parameters are randomised and returned for a maximum of (background, amplitude, angle, xpos,
     * ypos, xwidth, ywidth }
     *
     * @param r
     *            the random
     * @param npeaks
     *            the npeaks
     * @param params
     *            set on output
     * @param randomiseParams
     *            Set to true to randomise the params
     * @return the double[]
     */
    private double[] doubleCreateGaussianData(UniformRandomProvider r, int npeaks, double[] params,
            boolean randomiseParams)
    {
        final int n = blockWidth * blockWidth;

        // Generate a 2D Gaussian
        final ErfGaussian2DFunction func = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(npeaks, blockWidth,
                blockWidth, GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
        params[0] = random(r, background);
        for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK)
        {
            params[j + Gaussian2DFunction.SIGNAL] = random(r, signal);
            params[j + Gaussian2DFunction.X_POSITION] = random(r, xpos);
            params[j + Gaussian2DFunction.Y_POSITION] = random(r, ypos);
            params[j + Gaussian2DFunction.X_SD] = random(r, xwidth);
            params[j + Gaussian2DFunction.Y_SD] = random(r, ywidth);
        }

        final double[] y = new double[n];
        func.initialise(params);
        final CustomPoissonDistribution pd = new CustomPoissonDistribution(new RandomGeneratorAdapter(r), 1);
        for (int i = 0; i < y.length; i++)
        {
            // Add random Poisson noise
            final double u = func.eval(i);
            pd.setMean(u);
            y[i] = pd.sample();
        }

        if (randomiseParams)
        {
            params[0] = random(r, params[0]);
            for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK)
            {
                params[j + Gaussian2DFunction.SIGNAL] = random(r, params[j + Gaussian2DFunction.SIGNAL]);
                params[j + Gaussian2DFunction.X_POSITION] = random(r, params[j + Gaussian2DFunction.X_POSITION]);
                params[j + Gaussian2DFunction.Y_POSITION] = random(r, params[j + Gaussian2DFunction.Y_POSITION]);
                params[j + Gaussian2DFunction.X_SD] = random(r, params[j + Gaussian2DFunction.X_SD]);
                params[j + Gaussian2DFunction.Y_SD] = random(r, params[j + Gaussian2DFunction.Y_SD]); //params[j + 4];
            }
        }

        return y;
    }

    protected int[] createData(UniformRandomProvider r, int npeaks, int iter, ArrayList<double[]> paramsList,
            ArrayList<double[]> yList)
    {
        return createData(r, npeaks, iter, paramsList, yList, true);
    }

    protected int[] createData(UniformRandomProvider r, int npeaks, int iter, ArrayList<double[]> paramsList,
            ArrayList<double[]> yList, boolean randomiseParams)
    {
        final int[] x = new int[blockWidth * blockWidth];
        for (int i = 0; i < x.length; i++)
            x[i] = i;
        for (int i = 0; i < iter; i++)
        {
            final double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK * npeaks];
            final double[] y = doubleCreateGaussianData(r, npeaks, params, randomiseParams);
            paramsList.add(params);
            yList.add(y);
        }
        return x;
    }

    protected int[] createFakeData(UniformRandomProvider r, int nparams, int iter, ArrayList<double[]> paramsList,
            ArrayList<double[]> yList)
    {
        final int[] x = new int[blockWidth * blockWidth];
        for (int i = 0; i < x.length; i++)
            x[i] = i;
        for (int i = 0; i < iter; i++)
        {
            final double[] params = new double[nparams];
            final double[] y = createFakeData(r, params);
            paramsList.add(params);
            yList.add(y);
        }
        return x;
    }

    private double[] createFakeData(UniformRandomProvider r, double[] params)
    {
        final int n = blockWidth * blockWidth;

        for (int i = 0; i < params.length; i++)
            params[i] = r.nextDouble();

        final double[] y = new double[n];
        for (int i = 0; i < y.length; i++)
            y[i] = r.nextDouble() * 10;

        return y;
    }

    protected ArrayList<double[]> copyList(ArrayList<double[]> paramsList)
    {
        final ArrayList<double[]> params2List = new ArrayList<>(paramsList.size());
        for (int i = 0; i < paramsList.size(); i++)
            params2List.add(paramsList.get(i).clone());
        return params2List;
    }
}
