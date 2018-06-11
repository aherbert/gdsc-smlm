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
package gdsc.smlm.model;

import java.awt.Rectangle;
import java.util.Arrays;

import gdsc.smlm.ij.results.IJImagePeakResults;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;

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
		float intensity = 1;
		for (int i = 100000; i-- > 0;)
		{
			double[] xyz = dist.next();
			int frame = (int) (1 + scale * radius + Math.round(scale * xyz[2]));
			float x = radius + (float) xyz[0];
			float y = radius + (float) xyz[1];
			results.add(new PeakResult(frame, x, y, intensity));
		}
		results.end();
		IJImagePeakResults image = new IJImagePeakResults(
				(useRejctionMethod) ? "Rejection Method" : "Transformation Method", bounds, scale);
		image.setRollingWindowSize(1);
		image.begin();
		image.addAll(Arrays.asList(results.toArray()));
		// Place breakpoint here in debug mode to view the image. 
		// It should have an even colour through the stack.
		image.end();
	}
}
