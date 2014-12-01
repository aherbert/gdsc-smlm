package gdsc.smlm.model;

import java.awt.Rectangle;

import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.ij.results.IJImagePeakResults;
import gdsc.smlm.results.MemoryPeakResults;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

public class SphericalDistributionTest
{
	private RandomGenerator rand = new Well19937c(System.currentTimeMillis() + System.identityHashCode(this));

	@Test
	public void canSampleUsingTransformationMethod()
	{
		double radius = 10 + rand.nextDouble() * 10;
		SphericalDistribution dist = new SphericalDistribution(radius, rand);
		dist.setUseRejectionMethod(false);
		for (int i = 100; i-- > 0;)
			dist.next();
	}

	@Test
	public void canSampleUsingRejectionMethod()
	{
		double radius = 10 + rand.nextDouble() * 10;
		SphericalDistribution dist = new SphericalDistribution(radius, rand);
		dist.setUseRejectionMethod(true);
		for (int i = 100; i-- > 0;)
			dist.next();
	}

	@Test
	public void rejectionMethodIsFasterThanTransformationMethod()
	{
		double radius = 10 + rand.nextDouble() * 10;
		SphericalDistribution dist = new SphericalDistribution(radius, rand);
		dist.setUseRejectionMethod(false);
		for (int i = 100; i-- > 0;)
			dist.next();
		dist.setUseRejectionMethod(true);
		for (int i = 100; i-- > 0;)
			dist.next();

		dist.setUseRejectionMethod(false);
		long time1 = getRunTime(dist);
		dist.setUseRejectionMethod(true);
		long time2 = getRunTime(dist);
		String msg = String.format("Rejection = %d, Transformation = %d\n", time2, time1);
		Assert.assertTrue("Rejection method is slower: " + msg, time1 > time2);
		System.out.printf(msg);
	}

	private long getRunTime(SphericalDistribution dist)
	{
		long start = System.nanoTime();
		for (int i = 1000000; i-- > 0;)
			dist.next();
		return System.nanoTime() - start;
	}

	// These are not tests. They draw an image and use classes outside the package.
	// Comment out for production code.

	//@Test
	public void rejectionMethodSamplesEvenly()
	{
		drawImage(true);
	}

	//@Test
	public void transformationMethodSamplesEvenly()
	{
		drawImage(false);
	}

	private void drawImage(boolean useRejctionMethod)
	{
		MemoryPeakResults results = new MemoryPeakResults();
		results.setSortAfterEnd(true);
		int radius = 10;
		Rectangle bounds = new Rectangle(0, 0, radius * 2, radius * 2);
		SphericalDistribution dist = new SphericalDistribution(radius, rand);
		dist.setUseRejectionMethod(useRejctionMethod);
		float scale = 10;
		results.begin();
		for (int i = 100000; i-- > 0;)
		{
			double[] xyz = dist.next();
			int peak = (int) (1 + scale*radius + Math.round(scale * xyz[2]));
			float[] params = new float[7];
			params[Gaussian2DFunction.X_POSITION] = radius + (float) xyz[0];
			params[Gaussian2DFunction.Y_POSITION] = radius + (float) xyz[1];
			results.add(peak, 0, 0, 0, 0, 0, params, null);
		}
		results.end();
		IJImagePeakResults image = new IJImagePeakResults((useRejctionMethod) ? "Rejection Method"
				: "Transformation Method", bounds, scale);
		image.setRollingWindowSize(1);
		image.begin();
		image.addAll(results.getResults());
		// Place breakpoint here in debug mode to view the image. 
		// It should have an even colour through the stack.
		image.end();
	}
}
