package gdsc.smlm.utils;

import gdsc.smlm.TestSettings;
import gdsc.smlm.utils.DoubleEquality;
import gdsc.smlm.utils.FloatEquality;

import java.util.Random;

import org.junit.Assert;
import org.junit.Test;

public class EqualityTest
{
	int MAX_ITER = 2000000;

	@Test
	public void canComputeEquality()
	{
		float maxRelativeError = 1e-2f;
		float maxAbsoluteError = 1e-16f;
		FloatEquality equality = new FloatEquality(maxRelativeError, maxAbsoluteError, 3);

		for (int i = 0; i < 100; i++)
		{
			float f = (float) (i / 10000.0);
			Assert.assertTrue("not equal " + f, equality.almostEqualRelativeOrAbsolute(f, f));
			Assert.assertTrue("not equal " + f,
					equality.almostEqualRelativeOrAbsolute(f, f * (1.00f + maxRelativeError - 1e-3f)));
			if (i > 0)
				Assert.assertFalse("equal " + f,
						equality.almostEqualRelativeOrAbsolute(f, f * (1.0f + 2.0f * maxRelativeError)));
			Assert.assertEquals("compare == " + f, 0, equality.compareComplement(f, f));
			Assert.assertEquals("compare < " + f, -1, equality.compareComplement(f, f + 0.01f));
			Assert.assertEquals("compare > " + f, 1, equality.compareComplement(f + 0.01f, f));
			Assert.assertEquals("compare ~= " + f, 0, equality.compareComplement(f, f * 1.000001f));
			Assert.assertEquals("compare =~ " + f, 0, equality.compareComplement(f * 1.000001f, f));
		}

		intBits(100f);
		intBits(10f);
		intBits(1f);
		intBits(1e-1f);
		intBits(1e-2f);
		intBits(1e-3f);
		intBits(1e-4f);
		intBits(1e-5f);
		intBits(1e-6f);
		intBits(1e-7f);
		intBits(1e-8f);
		intBits(1e-9f);
		intBits(1e-10f);
		intBits(1e-11f);
		intBits(1e-12f);
		intBits(1e-13f);
		intBits(1e-14f);
		intBits(1e-15f);
		intBits(1e-16f);
		intBits(1e-26f);
		intBits(1e-36f);

		for (int i = 0; i < 18; i++)
			log("sig = %d -> %d : %d\n", i, FloatEquality.getUlps(i), DoubleEquality.getUlps(i));

		log("%d\n", DoubleEquality.complement(0e-16, -1.11022e-16));
		log("%d\n", DoubleEquality.complement(-1.11022e-15, -1.11022e-16));
	}

	/**
	 * Used to check what the int difference between float actually is
	 * 
	 * @param f
	 * @param f2
	 */
	private void intBits(float f)
	{
		float f3 = f + f * 1e-2f;
		float f4 = f - f * 1e-2f;
		System.out.printf("%g -> %g = %d : %d (%g : %g)\n", f, f3, FloatEquality.complement(f3, f),
				DoubleEquality.complement(f3, f), FloatEquality.relativeError(f, f3),
				DoubleEquality.relativeError(f, f3));
		System.out.printf("%g -> %g = %d : %d (%g : %g)\n", f, f4, FloatEquality.complement(f4, f),
				DoubleEquality.complement(f4, f), FloatEquality.relativeError(f, f4),
				DoubleEquality.relativeError(f, f4));
	}

	@Test
	public void floatRelativeIsSlowerThanFloatComplement()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		float maxRelativeError = 1e-2f;
		float maxAbsoluteError = 1e-16f;
		int significantDigits = 3;
		FloatEquality equality = new FloatEquality(maxRelativeError, maxAbsoluteError, significantDigits);

		// Create data
		Random rand = new Random(30051977);
		float[] data = new float[MAX_ITER];
		float[] data2 = new float[data.length];

		for (int i = 0; i < data.length; i++)
		{
			data[i] = 1e-10f * ((rand.nextFloat() > 0.5) ? 1.001f : 1.1f);
			data2[i] *= (rand.nextFloat() > 0.5) ? 1.001f : 1.1f;
		}

		relative(equality, data, data2);
		complement(equality, data, data2);

		long start1 = System.nanoTime();
		relative(equality, data, data2);
		start1 = System.nanoTime() - start1;

		long start2 = System.nanoTime();
		complement(equality, data, data2);
		start2 = System.nanoTime() - start2;

		log("floatRelative = %d : floatComplement = %d : %fx\n", start1, start2, (1.0 * start1) / start2);
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(start2 < start1);
	}

	@Test
	public void doubleRelativeIsSlowerThanFloatComplement()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		double maxRelativeError = 1e-2;
		double maxAbsoluteError = 1e-16;
		int significantDigits = 3;
		DoubleEquality equality = new DoubleEquality(maxRelativeError, maxAbsoluteError, significantDigits);

		// Create data
		Random rand = new Random(30051977);
		double[] data = new double[MAX_ITER];
		double[] data2 = new double[data.length];

		for (int i = 0; i < data.length; i++)
		{
			data[i] = 1e-10f * ((rand.nextFloat() > 0.5) ? 1.001f : 1.1f);
			data2[i] *= (rand.nextFloat() > 0.5) ? 1.001f : 1.1f;
		}

		relative(equality, data, data2);
		complement(equality, data, data2);

		long start1 = System.nanoTime();
		relative(equality, data, data2);
		start1 = System.nanoTime() - start1;

		long start2 = System.nanoTime();
		complement(equality, data, data2);
		start2 = System.nanoTime() - start2;

		log("doubleRelative = %d : doubleComplement = %d : %fx\n", start1, start2, (1.0 * start1) / start2);
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(start2 < start1);
	}

	@Test
	public void floatRelativeIsSlowerThanDoubleRelative()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		float maxRelativeError = 1e-2f;
		float maxAbsoluteError = 1e-16f;
		int significantDigits = 3;
		FloatEquality equality = new FloatEquality(maxRelativeError, maxAbsoluteError, significantDigits);
		DoubleEquality equality2 = new DoubleEquality(maxRelativeError, maxAbsoluteError, significantDigits);

		// Create data
		Random rand = new Random(30051977);
		float[] data = new float[MAX_ITER];
		float[] data2 = new float[data.length];
		double[] data3 = new double[data.length];
		double[] data4 = new double[data.length];

		for (int i = 0; i < data.length; i++)
		{
			data[i] = 1e-10f * ((rand.nextFloat() > 0.5) ? 1.001f : 1.1f);
			data2[i] *= (rand.nextFloat() > 0.5) ? 1.001f : 1.1f;
			data3[i] = data[i];
			data4[i] = data2[i];
		}

		relative(equality, data, data2);
		relative(equality2, data3, data4);

		long start1 = System.nanoTime();
		relative(equality, data, data2);
		start1 = System.nanoTime() - start1;

		long start2 = System.nanoTime();
		relative(equality2, data3, data4);
		start2 = System.nanoTime() - start2;

		log("floatRelative = %d : doubleRelative = %d : %fx\n", start1, start2, (1.0 * start1) / start2);
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(start2 < start1);
	}

	@Test
	public void floatComplementIsSlowerThanDoubleComplement()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		float maxComplementError = 1e-2f;
		float maxAbsoluteError = 1e-16f;
		int significantDigits = 3;
		FloatEquality equality = new FloatEquality(maxComplementError, maxAbsoluteError, significantDigits);
		DoubleEquality equality2 = new DoubleEquality(maxComplementError, maxAbsoluteError, significantDigits);

		// Create data
		Random rand = new Random(30051977);
		float[] data = new float[MAX_ITER];
		float[] data2 = new float[data.length];
		double[] data3 = new double[data.length];
		double[] data4 = new double[data.length];

		for (int i = 0; i < data.length; i++)
		{
			data[i] = 1e-10f * ((rand.nextFloat() > 0.5) ? 1.001f : 1.1f);
			data2[i] *= (rand.nextFloat() > 0.5) ? 1.001f : 1.1f;
			data3[i] = data[i];
			data4[i] = data2[i];
		}

		complement(equality, data, data2);
		complement(equality2, data3, data4);

		long start1 = System.nanoTime();
		complement(equality, data, data2);
		start1 = System.nanoTime() - start1;

		long start2 = System.nanoTime();
		complement(equality2, data3, data4);
		start2 = System.nanoTime() - start2;

		log("floatComplement = %d : doubleComplement = %d : %fx\n", start1, start2, (1.0 * start1) / start2);
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(start2 < start1);
	}

	private void relative(FloatEquality equality, float[] data, float[] data2)
	{
		for (int i = 0; i < data.length; i++)
			equality.almostEqualRelativeOrAbsolute(data[i], data2[i]);
	}

	private void relative(DoubleEquality equality, double[] data, double[] data2)
	{
		for (int i = 0; i < data.length; i++)
			equality.almostEqualRelativeOrAbsolute(data[i], data2[i]);
	}

	private void complement(FloatEquality equality, float[] data, float[] data2)
	{
		for (int i = 0; i < data.length; i++)
			equality.almostEqualComplement(data[i], data2[i]);
	}

	private void complement(DoubleEquality equality, double[] data, double[] data2)
	{
		for (int i = 0; i < data.length; i++)
			equality.almostEqualComplement(data[i], data2[i]);
	}

	void log(String format, Object... args)
	{
		System.out.printf(format, args);
	}
}
