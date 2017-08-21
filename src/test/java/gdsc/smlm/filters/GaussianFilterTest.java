package gdsc.smlm.filters;

import java.util.Arrays;

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.junit.Test;
import org.junit.internal.ArrayComparisonFailure;

import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Maths;

public class GaussianFilterTest
{
	private gdsc.core.utils.Random rand;

	// TODO - The test data should be representative of the final use case
	int[] primes = new int[] { 113, 97, 53, 29 };
	//int[] primes = new int[] { 1024 };
	int[] boxSizes = new int[] { 15, 9, 5, 3, 2, 1 };
	boolean[] checkInternal = new boolean[] { true, false };

	@Test
	public void filterDoesNotAlterImageMean()
	{
		rand = new gdsc.core.utils.Random(-30051976);
		ExponentialDistribution ed = new ExponentialDistribution(rand, 57,
				ExponentialDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
		GaussianFilter filter = new GaussianFilter();

		for (int width : primes)
			for (int height : primes)
			{
				float[] data = createData(width, height);
				double u1 = Maths.sum(data) / data.length;

				// TODO - rearrange so that the w=null, w=1 and weights are run sequentially
				// for the same box size.
				
				float[] w = null;
				filter.setWeights(w, width, height);
				for (float boxSize : boxSizes)
					for (boolean internal : checkInternal)
					{
						filterDoesNotAlterImageMean("w=null", data, w, width, height, boxSize, internal, u1, filter);
						filterDoesNotAlterImageMean("w=null", data, w, width, height, boxSize - 0.3f, internal, u1,
								filter);
						filterDoesNotAlterImageMean("w=null", data, w, width, height, boxSize - 0.6f, internal, u1,
								filter);
					}

				// Uniform weights
				w = new float[width * height];
				Arrays.fill(w, 1);
				filter.setWeights(w, width, height);
				for (float boxSize : boxSizes)
					for (boolean internal : checkInternal)
					{
						filterDoesNotAlterImageMean("w=1", data, w, width, height, boxSize, internal, u1, filter);
						filterDoesNotAlterImageMean("w=1", data, w, width, height, boxSize - 0.3f, internal, u1,
								filter);
						filterDoesNotAlterImageMean("w=1", data, w, width, height, boxSize - 0.6f, internal, u1,
								filter);
					}

				// Weights simulating the variance of sCMOS pixels
				for (int i = 0; i < w.length; i++)
				{
					w[i] = (float) (1.0 / Math.max(0.01, ed.sample()));
				}

				filter.setWeights(w, width, height);
				for (float boxSize : boxSizes)
					for (boolean internal : checkInternal)
					{
						filterDoesNotAlterImageMean("Weighted", data, w, width, height, boxSize, internal, u1, filter);
						filterDoesNotAlterImageMean("Weighted", data, w, width, height, boxSize - 0.3f, internal, u1,
								filter);
						filterDoesNotAlterImageMean("Weighted", data, w, width, height, boxSize - 0.6f, internal, u1,
								filter);
					}
			}
	}

	private void filterDoesNotAlterImageMean(String title, float[] data, float[] w, int width, int height,
			float boxSize, boolean internal, double u1, GaussianFilter filter) throws ArrayComparisonFailure
	{
		float[] data1 = data.clone();
		if (internal)
		{
			filter.convolveInternal(data1, width, height, boxSize);
		}
		else
		{
			filter.convolve(data1, width, height, boxSize);
		}

		// Check the weights do not alter the image mean
		double u2 = Maths.sum(data1) / data.length;
		System.out.printf("[%dx%d] @ %.1f [%s,internal=%b] : %g => %g  (%g)\n", width, height, boxSize, title, internal, u1, u2,
				DoubleEquality.relativeError(u1, u2));
	}

	private float[] createData(int width, int height)
	{
		float[] data = new float[width * height];
		for (int i = data.length; i-- > 0;)
			data[i] = i;

		rand.shuffle(data);

		return data;
	}
}
