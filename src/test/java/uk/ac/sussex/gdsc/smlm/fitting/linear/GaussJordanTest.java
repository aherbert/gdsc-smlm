package uk.ac.sussex.gdsc.smlm.fitting.linear;

import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;

import java.util.Arrays;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import uk.ac.sussex.gdsc.test.LogLevel;
import uk.ac.sussex.gdsc.test.TestSettings;

@SuppressWarnings({ "javadoc" })
public class GaussJordanTest
{
    private static Logger logger;

    @BeforeAll
    public static void beforeAll()
    {
        logger = Logger.getLogger(GaussJordanTest.class.getName());
    }

    @AfterAll
    public static void afterAll()
    {
        logger = null;
    }

	@Test
	public void canSolveLinearEquation()
	{
		final GaussJordan solver = new GaussJordan();

		// Solves (one) linear equation, a x = b, for x[n]

		// Example taken from http://eigen.tuxfamily.org/dox/TutorialLinearAlgebra.html
		final float[][] a = new float[][] { new float[] { 1, 2, 3 }, new float[] { 4, 5, 6 }, new float[] { 7, 8, 10 } };
		final float[] b = new float[] { 3, 3, 4 };
		final float[] expecteds = new float[] { -2, 1, 1 };

		final boolean result = solver.solve(a, b);

		Assertions.assertTrue(result);
		Assertions.assertArrayEquals(expecteds, b, 1e-4f);

		if (logger.isLoggable(Level.INFO))
		{
			log("x = %s\n", Arrays.toString(b));
			for (int i = 0; i < b.length; i++)
				log("a[%d] = %s\n", i, Arrays.toString(a[i]));
		}
	}

	void log(String format, Object... args)
	{
		System.out.printf(format, args);
	}
}
