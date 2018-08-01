package uk.ac.sussex.gdsc.smlm.utils;

import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import uk.ac.sussex.gdsc.test.LogLevel;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;

@SuppressWarnings({ "javadoc" })
public class TensorTest
{
    private static Logger logger;

    @BeforeAll
    public static void beforeAll()
    {
        logger = Logger.getLogger(TensorTest.class.getName());
    }

    @AfterAll
    public static void afterAll()
    {
        logger = null;
    }

	@Test
	public void canComputeTensor3D()
	{
		//@formatter:off
		final float[][] data = new float[][] {
			{ 2, 1, 0, 1, 2, 1, 0, 1, 2 },
			//{ 1, 0, 0, 0, 1, 0, 0, 0, 1 },
			//{ 1, 0, 0, 0, 1, 0, 0, 0, 1 },
		};
		//@formatter:on
		final Tensor3D t = new Tensor3D(data, 3, 3);
		Assertions.assertTrue(t.hasDecomposition());
		final double[] com = t.getCentreOfMass();
		Assertions.assertArrayEquals(new double[] { 1, 1, 0 }, com);
		final double[] v = t.getEigenValues();
		final double[][] vv = t.getEigenVectors();
		print(com, v, vv);
		for (int i = 1; i < v.length; i++)
			Assertions.assertTrue(v[i - 1] >= v[i]);
	}

	private static void print(double[] com, double[] v, double[][] vv)
	{
		if (logger.isLoggable(Level.INFO))
		{
			System.out.printf("com = %s\n", java.util.Arrays.toString(com));
			for (int i = 0; i < v.length; i++)
				System.out.printf("[%d] %f = %s  %.2f\n", i, v[i], java.util.Arrays.toString(vv[i]),
						180.0 * Math.atan2(vv[i][1], vv[i][0]) / Math.PI);
		}
	}

	@Test
	public void canComputeTensor2D()
	{
		//@formatter:off
		// Line through [0][0], [1][1], [2][2]
		// longest axis of object is -45 degrees
		final float[] data = new float[] {
				//1, 0, 0, 0, 1, 0, 0, 0, 1
				2, 1, 0, 1, 2, 1, 0, 1, 2
				//2, 0, 0, 0, 0, 0, 0, 0, 2
				};
		//@formatter:on
		final Tensor2D t = new Tensor2D(data, 3, 3);
		Assertions.assertTrue(t.hasDecomposition());
		final double[] com = t.getCentreOfMass();
		Assertions.assertArrayEquals(new double[] { 1, 1 }, com);
		final double[] v = t.getEigenValues();
		final double[][] vv = t.getEigenVectors();
		print(com, v, vv);
		for (int i = 1; i < v.length; i++)
			Assertions.assertTrue(v[i - 1] >= v[i]);
	}

	@SeededTest
	public void canComputeSameTensor(RandomSeed seed)
	{
		final UniformRandomProvider random = TestSettings.getRandomGenerator(seed.getSeed());
		final int w = 3, h = 4;
		final float[] data = new float[w * h];
		for (int i = 0; i < 10; i++)
		{
			for (int j = data.length; j-- > 0;)
				data[j] = random.nextFloat();

			final Tensor2D t2 = new Tensor2D(data, w, h);
			final Tensor3D t3 = new Tensor3D(new float[][] { data }, w, h);

			final double[] com2 = t2.getCentreOfMass();
			final double[] v2 = t2.getEigenValues();
			final double[][] vv2 = t2.getEigenVectors();
			final double[] com3 = t3.getCentreOfMass();
			final double[] v3 = t3.getEigenValues();
			final double[][] vv3 = t3.getEigenVectors();
			for (int k = 0; k < 2; k++)
			{
				Assertions.assertEquals(com2[k], com3[k]);
				ExtraAssertions.assertEqualsRelative(v2[k], v3[k + 1], 1e-6);
				for (int kk = 0; kk < 2; kk++)
				{
					// Swap vector direction
					if (Math.signum(vv2[k][kk]) != Math.signum(vv3[k + 1][kk]))
					{
						vv2[k][0] = -vv2[k][0];
						vv2[k][1] = -vv2[k][1];
					}
					ExtraAssertions.assertEqualsRelative(vv2[k][kk], vv3[k + 1][kk], 1e-6);
				}
			}
		}
	}
}
