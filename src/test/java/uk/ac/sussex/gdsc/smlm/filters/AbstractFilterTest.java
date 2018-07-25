/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package uk.ac.sussex.gdsc.smlm.filters;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.random.RandomGenerator;

import uk.ac.sussex.gdsc.core.utils.FloatEquality;
import uk.ac.sussex.gdsc.core.utils.Random;
import uk.ac.sussex.gdsc.test.LogLevel;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;

@SuppressWarnings({ "javadoc" })
public class AbstractFilterTest
{
	static final boolean debug = TestSettings.allow(LogLevel.DEBUG);

	// TODO - The test data should be representative of the final use case

	/** The primes used for the width/height of images during filter testing. */
	static int[] primes = new int[] { 113, /* 97, 53, */ 29 };

	/** The primes used for the width/height of images during speed testing. */
	static int[] speedPrimes = new int[] { 113 };

	/**
	 * The box sizes used during filter testing.
	 * 15 is required to make the box larger than the smallest image.
	 */
	static int[] boxSizes = new int[] { 15, 5, 3, 2, 1 };

	/** The box sizes used during filter testing for filters that can use non-integer sizes. */
	static float[] fBoxSizes;
	static
	{
		fBoxSizes = new float[boxSizes.length];
		for (int i = 0; i < boxSizes.length; i++)
			fBoxSizes[i] = boxSizes[i] - 0.5f;
	}

	/** The check internal flags [true,false]. */
	static boolean[] checkInternal = new boolean[] { true, false };

	/**
	 * Create random float data
	 *
	 * @param rg
	 *            the random generator
	 * @param width
	 *            the width
	 * @param height
	 *            the height
	 * @return the float data
	 */
	static float[] createData(RandomGenerator rg, int width, int height)
	{
		final float[] data = new float[width * height];
		for (int i = data.length; i-- > 0;)
			data[i] = rg.nextFloat();
		return data;
	}

	/**
	 * Create random int data
	 *
	 * @param rg
	 *            the random generator
	 * @param width
	 *            the width
	 * @param height
	 *            the height
	 * @return the int data
	 */
	static int[] createIntData(RandomGenerator rg, int width, int height)
	{
		final int[] data = new int[width * height];
		for (int i = data.length; i-- > 0;)
			data[i] = i;
		Random.shuffle(data, rg);
		return data;
	}

	static ArrayList<float[]> dataSet = new ArrayList<>();
	static RandomGenerator rg = null;

	/**
	 * Create random float data.
	 *
	 * @param size
	 *            the number of datasets
	 * @return the array list of random data
	 */
	static ArrayList<float[]> getSpeedData(int size)
	{
		synchronized (dataSet)
		{
			if (rg == null)
				rg = TestSettings.getRandomGenerator();
			while (dataSet.size() < size)
				dataSet.add(createData(rg, primes[0], primes[0]));
		}

		final ArrayList<float[]> dataSet2 = new ArrayList<>(size);
		for (int i = 0; i < size; i++)
			dataSet2.add(dataSet.get(i).clone());
		return dataSet2;
	}

	/**
	 * Create random int data.
	 *
	 * @param size
	 *            the number of datasets
	 * @return the array list of random data
	 */
	static ArrayList<int[]> getIntSpeedData(int size)
	{
		synchronized (dataSet)
		{
			if (rg == null)
				rg = TestSettings.getRandomGenerator();
			while (dataSet.size() < size)
				dataSet.add(createData(rg, primes[0], primes[0]));
		}

		final ArrayList<int[]> dataSet2 = new ArrayList<>(size);
		for (int i = 0; i < size; i++)
		{
			final float[] f = dataSet.get(i);
			final int[] d = new int[f.length];
			for (int j = 0; j < d.length; j++)
				d[j] = (int) (4096 * f[j]);
			dataSet2.add(d);
		}
		return dataSet2;
	}

	static double speedUpFactor(long slowTotal, long fastTotal)
	{
		return (1.0 * slowTotal) / fastTotal;
	}

	static float[] floatClone(float[] data1)
	{
		return Arrays.copyOf(data1, data1.length);
	}

	static float[] floatClone(int[] data1)
	{
		final float[] data2 = new float[data1.length];
		for (int i = data2.length; i-- > 0;)
			data2[i] = data1[i];
		return data2;
	}

	static int[] intClone(int[] data1)
	{
		return Arrays.copyOf(data1, data1.length);
	}

	static void floatArrayEquals(FloatEquality eq, float[] data1, float[] data2, int maxx, int maxy, float boxSize,
			String format, Object... args)
	{
		// Ignore the border
		final int border = (int) Math.ceil(boxSize);
		for (int y = border; y < maxy - border - 1; y++)
		{
			int index = y * maxx + border;
			for (int x = border; x < maxx - border - 1; x++, index++)
				if (!eq.almostEqualRelativeOrAbsolute(data1[index], data2[index]))
				{
					final String message = String.format(format, args);
					ExtraAssertions.fail("%s [%d,%d] %f != %f  (%g)", message, x, y, data1[index], data2[index],
							FloatEquality.relativeError(data1[index], data2[index]));
				}
		}
	}

	static void intArrayEquals(int[] data1, int[] data2, int maxx, int maxy, float boxSize, String format,
			Object... args)
	{
		// Ignore the border
		final int border = (int) Math.ceil(boxSize);
		for (int y = border; y < maxy - border - 1; y++)
		{
			int index = y * maxx + border;
			for (int x = border; x < maxx - border - 1; x++, index++)
				if (data1[index] != data2[index])
				{
					final String message = String.format(format, args);
					ExtraAssertions.fail("%s [%d,%d] %f != %f  (%g)", message, x, y, data1[index], data2[index],
							FloatEquality.relativeError(data1[index], data2[index]));
				}
		}
	}
}
