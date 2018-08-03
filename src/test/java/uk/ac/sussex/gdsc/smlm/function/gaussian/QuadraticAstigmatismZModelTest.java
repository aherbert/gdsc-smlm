package uk.ac.sussex.gdsc.smlm.function.gaussian;

import java.util.logging.Level;
import java.util.logging.Logger;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;

@SuppressWarnings({ "javadoc" })
public class QuadraticAstigmatismZModelTest
{
    private static Logger logger;

    @BeforeAll
    public static void beforeAll()
    {
        logger = Logger.getLogger(QuadraticAstigmatismZModelTest.class.getName());
    }

    @AfterAll
    public static void afterAll()
    {
        logger = null;
    }

    protected DoubleEquality eq = new DoubleEquality(1e-5, 1e-7);

    // Compute as per Numerical Recipes 5.7.
    // Approximate error accuracy in single precision: Ef
    // Step size for derivatives:
    // h ~ (Ef)^(1/3) * xc
    // xc is the characteristic scale over which x changes, assumed to be 1 (not x as per NR since x is close to zero)
    protected double h_ = 0.0001; //(double) (Math.pow(1e-3f, 1.0 / 3));

    @Test
    public void canStaticComputeGradient()
    {
        canStaticComputeGradient(1.2);
        canStaticComputeGradient(0.3);
    }

    private void canStaticComputeGradient(double zDepth)
    {
        final double[] ds_dz = new double[1];
        final double[] ds_dz2 = new double[2];
        final double[] ds_duz = new double[1];
        final double[] ds_dlz = new double[1];
        final boolean record = logger.isLoggable(Level.INFO);
        for (double z = -0.5; z < 0.5; z += 0.01)
        {
            final double s0 = QuadraticAstigmatismZModel.getS(z, zDepth);
            final double s1 = QuadraticAstigmatismZModel.getS1(z, zDepth, ds_dz);
            final double s2 = QuadraticAstigmatismZModel.getS2(z, zDepth, ds_dz2);

            Assertions.assertEquals(s0, s1);
            Assertions.assertEquals(s0, s2);
            Assertions.assertEquals(ds_dz[0], ds_dz2[0]);

            final double uz = z + h_;
            final double lz = z - h_;
            final double upper = QuadraticAstigmatismZModel.getS1(uz, zDepth, ds_duz);
            final double lower = QuadraticAstigmatismZModel.getS1(lz, zDepth, ds_dlz);

            final double e1 = (upper - lower) / (uz - lz);
            final double o1 = ds_dz[0];

            // Second gradient
            final double e2 = (ds_duz[0] - ds_dlz[0]) / (uz - lz);
            final double o2 = ds_dz2[1];

            if (record)
                logger.log(TestLog.getRecord(Level.INFO, "z=%f s=%f : ds_dz=%g  %g  (%g): d2s_dz2=%g   %g  (%g)", z, s0, e1, o1,
                        DoubleEquality.relativeError(o1, e1), e2, o2, DoubleEquality.relativeError(o2, e2)));

            //double error = DoubleEquality.relativeError(o, e);
            if (Math.abs(z) > 0.02)
                ExtraAssertions.assertTrue((e1 * o1) >= 0, "%s sign != %s", e1, o1);
            ExtraAssertions.assertTrue(eq.almostEqualRelativeOrAbsolute(e1, o1), "%s != %s", e1, o1);

            if (Math.abs(z) > 0.02)
                ExtraAssertions.assertTrue((e2 * o2) >= 0, "%s sign != %s", e2, o2);
            ExtraAssertions.assertTrue(eq.almostEqualRelativeOrAbsolute(e2, o2), "%s != %s", e2, o2);
        }
    }
}
