package gdsc.smlm.utils;

import org.junit.Assert;
import org.junit.Test;

public class GaussianKernelTest
{
	@Test
	public void canComputeGaussianKernelIncScaleIncRange()
	{
		for (int i = 0; i < 5; i++)
		{
			double s = 0.33 * (1 << i);

			GaussianKernel k = new GaussianKernel(s);

			for (int scale = 1; scale < 16; scale *= 2)
				for (int range = 3; range < 5; range++)
				{
					//System.out.printf("s=%g scale=%d range=%d\n", s, scale, range);

					for (boolean edgeCorrection : new boolean[] { true, false })
					{
						double[] e = GaussianKernel.makeGaussianKernel(s * scale, range, edgeCorrection);
						double[] o = k.getGaussianKernel(scale, range, edgeCorrection);

						Assert.assertArrayEquals(e, o, 0.0);
					}
				}
		}
	}

	@Test
	public void canComputeGaussianKernelDecScaleIncRange()
	{
		for (int i = 0; i < 5; i++)
		{
			double s = 0.33 * (1 << i);

			GaussianKernel k = new GaussianKernel(s);

			for (int scale = 16; scale > 1; scale /= 2)
				for (int range = 3; range < 5; range++)
				{
					//System.out.printf("s=%g scale=%d range=%d\n", s, scale, range);

					for (boolean edgeCorrection : new boolean[] { true, false })
					{
						double[] e = GaussianKernel.makeGaussianKernel(s * scale, range, edgeCorrection);
						double[] o = k.getGaussianKernel(scale, range, edgeCorrection);

						Assert.assertArrayEquals(e, o, 0.0);
					}
				}
		}
	}

	@Test
	public void canComputeGaussianKernelIncRangeIncScale()
	{
		for (int i = 0; i < 5; i++)
		{
			double s = 0.33 * (1 << i);

			GaussianKernel k = new GaussianKernel(s);

			for (int range = 3; range < 5; range++)
				for (int scale = 1; scale < 16; scale *= 2)
				{
					//System.out.printf("s=%g scale=%d range=%d\n", s, scale, range);

					for (boolean edgeCorrection : new boolean[] { true, false })
					{
						double[] e = GaussianKernel.makeGaussianKernel(s * scale, range, edgeCorrection);
						double[] o = k.getGaussianKernel(scale, range, edgeCorrection);

						Assert.assertArrayEquals(e, o, 0.0);
					}
				}
		}
	}

	@Test
	public void canComputeGaussianKernelIncRangeDecScale()
	{
		for (int i = 0; i < 5; i++)
		{
			double s = 0.33 * (1 << i);

			GaussianKernel k = new GaussianKernel(s);

			for (int range = 3; range < 5; range++)
				for (int scale = 16; scale > 1; scale /= 2)
				{
					//System.out.printf("s=%g scale=%d range=%d\n", s, scale, range);

					for (boolean edgeCorrection : new boolean[] { true, false })
					{
						double[] e = GaussianKernel.makeGaussianKernel(s * scale, range, edgeCorrection);
						double[] o = k.getGaussianKernel(scale, range, edgeCorrection);

						Assert.assertArrayEquals(e, o, 0.0);
					}
				}
		}
	}
}
