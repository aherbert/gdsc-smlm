package gdsc.smlm.filters;

import java.util.Arrays;

//import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Maths;
import gnu.trove.list.array.TDoubleArrayList;

public abstract class WeightedFilterTest
{
	gdsc.core.utils.Random rand;

	int[] primes = new int[] { 113, 97, 53 };
	int[] boxSizes = new int[] { 9, 5, 3, 2, 1 };
	float[] offsets = new float[] { 0, 0.3f, 0.6f };
	boolean[] checkInternal = new boolean[] { true, false };

	float[] createData(int width, int height)
	{
		float[] data = new float[width * height];
		for (int i = data.length; i-- > 0;)
			data[i] = i;

		rand.shuffle(data);

		return data;
	}

	abstract DataFilter createDataFilter();

	@Test
	public void evenWeightsDoesNotAlterFiltering()
	{
		rand = new gdsc.core.utils.Random(-30051976);

		DataFilter filter1 = createDataFilter();
		DataFilter filter2 = createDataFilter();
		float[] offsets = (filter1.isInterpolated) ? this.offsets : new float[1];

		int[] primes = Arrays.copyOf(this.primes, this.primes.length - 1);

		for (int width : primes)
			for (int height : primes)
			{
				float[] data = createData(width, height);

				// Uniform weights
				float[] w = new float[width * height];
				Arrays.fill(w, 0.5f);
				filter2.setWeights(w, width, height);

				for (float boxSize : boxSizes)
					for (float offset : offsets)
						for (boolean internal : checkInternal)
						{
							float[] e = filter(data, width, height, boxSize - offset, internal, filter1);
							float[] o = filter(data, width, height, boxSize - offset, internal, filter2);
							try
							{
								Assert.assertArrayEquals(e, o, 0.1f);
							}
							catch (AssertionError ex)
							{
								String msg = String.format("%s : [%dx%d] @ %.1f [internal=%b]", filter2.name, width,
										height, boxSize, internal);
								throw new AssertionError(msg, ex);
							}
						}
			}
	}

	@Test
	public void filterDoesNotAlterImageMean()
	{
		rand = new gdsc.core.utils.Random(-30051976);
		//ExponentialDistribution ed = new ExponentialDistribution(rand, 57,
		//		ExponentialDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);

		DataFilter filter = createDataFilter();
		float[] offsets = (filter.isInterpolated) ? this.offsets : new float[1];

		TDoubleArrayList l1 = new TDoubleArrayList();

		for (int width : primes)
			for (int height : primes)
			{
				float[] data = createData(width, height);
				l1.reset();

				filter.setWeights(null, width, height);
				for (float boxSize : boxSizes)
					for (float offset : offsets)
						for (boolean internal : checkInternal)
						{
							l1.add(getMean(data, width, height, boxSize - offset, internal, filter));
						}

				double[] e = l1.toArray();
				int ei = 0;

				// Uniform weights
				float[] w = new float[width * height];

				Arrays.fill(w, 0.5f);
				filter.setWeights(w, width, height);
				for (float boxSize : boxSizes)
					for (float offset : offsets)
						for (boolean internal : checkInternal)
						{
							testMean(data, width, height, boxSize - offset, internal, filter, "w=0.5", e[ei++], 1e-5);
						}

				// Weights simulating the variance of sCMOS pixels
				for (int i = 0; i < w.length; i++)
				{
					//w[i] = (float) (1.0 / Math.max(0.01, ed.sample()));
					w[i] = (float) (1.0 / Math.max(0.01, rand.nextGaussian() * 0.2 + 2));
					//w[i] = 0.5f;
				}

				ei = 0;
				filter.setWeights(w, width, height);
				for (float boxSize : boxSizes)
					for (float offset : offsets)
						for (boolean internal : checkInternal)
						{
							testMean(data, width, height, boxSize - offset, internal, filter, "w=?", e[ei++], 1e-2);
						}
			}
	}

	private static float[] filter(float[] data, int width, int height, float boxSize, boolean internal,
			DataFilter filter)
	{
		float[] data1 = data.clone();
		if (internal)
		{
			filter.filterInternal(data1, width, height, boxSize);
		}
		else
		{
			filter.filter(data1, width, height, boxSize);
		}
		return data1;
	}

	private static double getMean(float[] data, int width, int height, float boxSize, boolean internal,
			DataFilter filter)
	{
		return Maths.sum(filter(data, width, height, boxSize, internal, filter)) / data.length;
	}

	private static double testMean(float[] data, int width, int height, float boxSize, boolean internal,
			DataFilter filter, String title, double u1, double tol)
	{
		double u2 = getMean(data, width, height, boxSize, internal, filter);
		double error = DoubleEquality.relativeError(u1, u2);
		//System.out.printf("%s : %s [%dx%d] @ %.1f [internal=%b] : %g => %g  (%g)\n", filter.name, title,
		//		width, height, boxSize, internal, u1, u2, error);
		if (error > tol)
		{
			String msg = String.format("%s : %s [%dx%d] @ %.1f [internal=%b] : %g => %g  (%g)", filter.name, title,
					width, height, boxSize, internal, u1, u2, error);
			Assert.fail(msg);
		}
		return u2;
	}
}
