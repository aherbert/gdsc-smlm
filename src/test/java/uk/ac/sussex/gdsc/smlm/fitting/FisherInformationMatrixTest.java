package uk.ac.sussex.gdsc.smlm.fitting;

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.rng.UniformRandomProvider;
import org.ejml.data.DenseMatrix64F;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.opentest4j.AssertionFailedError;

import uk.ac.sussex.gdsc.core.utils.Maths;
import uk.ac.sussex.gdsc.core.utils.Random;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.GradientCalculatorFactory;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RNGFactory;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;

@SuppressWarnings({ "javadoc" })
public class FisherInformationMatrixTest
{
    private static Logger logger;

    @BeforeAll
    public static void beforeAll()
    {
        logger = Logger.getLogger(FisherInformationMatrixTest.class.getName());
    }

    @AfterAll
    public static void afterAll()
    {
        logger = null;
    }

    /** The level for logging output. */
    private final Level level = Level.FINE;

    @SeededTest
    public void canComputeCRLB(RandomSeed seed)
    {
        final UniformRandomProvider rg = RNGFactory.create(seed.getSeed());
        for (int n = 1; n < 10; n++)
            canComputeCRLB(rg, n, 0, true);
    }

    @SeededTest
    public void canComputeCRLBWithZeros(RandomSeed seed)
    {
        final UniformRandomProvider rg = RNGFactory.create(seed.getSeed());
        for (int n = 2; n < 10; n++)
        {
            canComputeCRLB(rg, n, 1, true);
            canComputeCRLB(rg, n, n / 2, true);
        }
    }

    @SeededTest
    public void canComputeCRLBWithReciprocal(RandomSeed seed)
    {
        final UniformRandomProvider rg = RNGFactory.create(seed.getSeed());
        for (int n = 1; n < 10; n++)
            canComputeCRLB(rg, n, 0, false);
    }

    @SeededTest
    public void canComputeCRLBWithReciprocalWithZeros(RandomSeed seed)
    {
        final UniformRandomProvider rg = RNGFactory.create(seed.getSeed());
        for (int n = 2; n < 10; n++)
        {
            canComputeCRLB(rg, n, 1, false);
            canComputeCRLB(rg, n, n / 2, false);
        }
    }

    @SeededTest
    public void inversionDoesNotMatchReciprocal(RandomSeed seed)
    {
        final UniformRandomProvider rg = RNGFactory.create(seed.getSeed());
        for (int n = 1; n < 10; n++)
        {
            final FisherInformationMatrix m = createFisherInformationMatrix(rg, n, 0);
            final double[] crlb = m.crlb();
            final double[] crlb2 = m.crlbReciprocal();
            // These increasingly do not match with increasing number of parameters.
            if (logger.isLoggable(level))
                logger.log(level, FunctionUtils.getSupplier("%s =? %s", Arrays.toString(crlb), Arrays.toString(crlb2)));
            if (n > 1)
                // Just do a sum so we have a test
                Assertions.assertThrows(AssertionFailedError.class, () -> {
                    Assertions.assertEquals(Maths.sum(crlb), Maths.sum(crlb2));
                });
        }
    }

    private double[] canComputeCRLB(UniformRandomProvider rg, int n, int k, boolean invert)
    {
        final FisherInformationMatrix m = createFisherInformationMatrix(rg, n, k);

        // Invert for CRLB
        final double[] crlb = (invert) ? m.crlb() : m.crlbReciprocal();
        if (logger.isLoggable(level))
            logger.log(level, FunctionUtils.getSupplier("n=%d, k=%d : %s", n, k, Arrays.toString(crlb)));
        Assertions.assertNotNull(crlb, () -> String.format("CRLB failed: n=%d, k=%d", n, k));
        return crlb;
    }

    private static FisherInformationMatrix createFisherInformationMatrix(UniformRandomProvider rg, int n, int k)
    {
        final int maxx = 10;
        final int size = maxx * maxx;

        // Use a real Gaussian function here to compute the Fisher information.
        // The matrix may be sensitive to the type of equation used.
        int npeaks = 1;
        Gaussian2DFunction f = createFunction(maxx, npeaks);
        while (f.getNumberOfGradients() < n)
        {
            npeaks++;
            f = createFunction(maxx, npeaks);
        }

        final double[] a = new double[1 + npeaks * Gaussian2DFunction.PARAMETERS_PER_PEAK];
        a[Gaussian2DFunction.BACKGROUND] = nextUniform(rg, 1, 5);
        for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK)
        {
            a[j + Gaussian2DFunction.SIGNAL] = nextUniform(rg, 100, 300);
            // Non-overlapping peaks otherwise the CRLB are poor
            a[j + Gaussian2DFunction.X_POSITION] = nextUniform(rg, 2 + i * 2, 4 + i * 2);
            a[j + Gaussian2DFunction.Y_POSITION] = nextUniform(rg, 2 + i * 2, 4 + i * 2);
            a[j + Gaussian2DFunction.X_SD] = nextUniform(rg, 1.5, 2);
            a[j + Gaussian2DFunction.Y_SD] = nextUniform(rg, 1.5, 2);
        }
        f.initialise(a);

        final GradientCalculator c = GradientCalculatorFactory.newCalculator(f.getNumberOfGradients());
        double[][] I = c.fisherInformationMatrix(size, a, f);

        //TestLog.fine(logger,"n=%d, k=%d, I=", n, k);
        //for (int i = 0; i < I.length; i++)
        //	TestLog.fine(logger,Arrays.toString(I[i]));

        // Reduce to the desired size
        I = Arrays.copyOf(I, n);
        for (int i = 0; i < n; i++)
            I[i] = Arrays.copyOf(I[i], n);

        // Zero selected columns
        if (k > 0)
        {
            final int[] zero = Random.sample(k, n, rg); // new RandomDataGenerator(UniformRandomProvider).nextPermutation(n, k);
            for (final int i : zero)
                for (int j = 0; j < n; j++)
                    I[i][j] = I[j][i] = 0;
        }

        //TestLog.fine(logger,"n=%d, k=%d", n, k);
        //for (int i = 0; i < n; i++)
        //	TestLog.fine(logger,Arrays.toString(I[i]));

        // Create matrix
        return new FisherInformationMatrix(I, 1e-3);
    }

    private static double nextUniform(UniformRandomProvider r, double min, double max)
    {
        return min + r.nextDouble() * (max - min);
    }

    private static Gaussian2DFunction createFunction(int maxx, int npeaks)
    {
        final Gaussian2DFunction f = GaussianFunctionFactory.create2D(npeaks, maxx, maxx,
                GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
        return f;
    }

    @SeededTest
    public void canProduceSubset(RandomSeed seed)
    {
        final int k = 5;
        final int n = 10;

        final UniformRandomProvider UniformRandomProvider = RNGFactory.create(seed.getSeed());
        final FisherInformationMatrix m = createRandomMatrix(UniformRandomProvider, n);
        final DenseMatrix64F e = m.getMatrix();
        if (logger.isLoggable(level))
            logger.log(level, String.valueOf(e));

        for (int run = 1; run < 10; run++)
        {
            final int[] indices = Random.sample(k, n, UniformRandomProvider);
            Arrays.sort(indices);
            final DenseMatrix64F o = m.subset(indices).getMatrix();
            if (logger.isLoggable(level))
            {
                logger.log(level, FunctionUtils.getSupplier(Arrays.toString(indices)));
                logger.log(level, String.valueOf(o));
            }
            for (int i = 0; i < indices.length; i++)
                for (int j = 0; j < indices.length; j++)
                    Assertions.assertEquals(e.get(indices[i], indices[j]), o.get(i, j));
        }
    }

    private static FisherInformationMatrix createRandomMatrix(UniformRandomProvider UniformRandomProvider, int n)
    {
        final double[] data = new double[n * n];
        for (int i = 0; i < data.length; i++)
            data[i] = UniformRandomProvider.nextDouble();
        return new FisherInformationMatrix(data, n);
    }

    @SeededTest
    public void computeWithSubsetReducesTheCRLB(RandomSeed seed)
    {
        final UniformRandomProvider rg = RNGFactory.create(seed.getSeed());
        final Gaussian2DFunction f = createFunction(10, 1);
        final int perPeak = f.getGradientParametersPerPeak();
        // Create a matrix with 2 peaks + background
        final FisherInformationMatrix m = createFisherInformationMatrix(rg, 1 + 2 * perPeak, 0);
        // Subset each peak
        final int[] indices = SimpleArrayUtils.newArray(1 + perPeak, 0, 1);
        final FisherInformationMatrix m1 = m.subset(indices);
        for (int i = 1; i < indices.length; i++)
            indices[i] += perPeak;
        final FisherInformationMatrix m2 = m.subset(indices);

        //TestLog.fine(logger,m.getMatrix());
        //TestLog.fine(logger,m1.getMatrix());
        //TestLog.fine(logger,m2.getMatrix());

        final double[] crlb = m.crlb();
        final double[] crlb1 = m1.crlb();
        final double[] crlb2 = m2.crlb();
        final double[] crlbB = Arrays.copyOf(crlb1, crlb.length);
        System.arraycopy(crlb2, 1, crlbB, crlb1.length, perPeak);

        //TestLog.fine(logger,Arrays.toString(crlb));
        //TestLog.fine(logger,Arrays.toString(crlb1));
        //TestLog.fine(logger,Arrays.toString(crlb2));
        //TestLog.fine(logger,Arrays.toString(crlbB));

        // Removing the interaction between fit parameters lowers the bounds
        for (int i = 0; i < crlb.length; i++)
            Assertions.assertTrue(crlbB[i] < crlb[i]);
    }
}
