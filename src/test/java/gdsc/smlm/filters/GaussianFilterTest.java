package gdsc.smlm.filters;

import org.junit.Assert;
import org.junit.Test;

import gdsc.core.test.BaseTimingTask;
import gdsc.core.test.TimingService;

public class GaussianFilterTest
{
	private gdsc.core.utils.Random rand;

	int[] primes = new int[] { 113, 97, 53 };
	double[] sigmas = new double[] { 9.3, 5, 3.2, 2.1, 0.5 };
	int size = 256;

	@Test
	public void floatFilterIsSameAsDoubleFilter()
	{
		rand = new gdsc.core.utils.Random(-30051976);

		GaussianFilter ff = new GaussianFilter();
		DoubleGaussianFilter df = new DoubleGaussianFilter();
		float[] e, o;

		float tolerance = 1e-2f;

		for (int width : primes)
			for (int height : primes)
			{
				float[] data = createData(width, height);

				for (double sigma : sigmas)
				{
					e = data.clone();
					df.convolve(e, width, height, sigma);
					o = data.clone();
					ff.convolve(o, width, height, sigma);
					Assert.assertArrayEquals("Full", e, o, tolerance);

					e = data.clone();
					df.convolveInternal(e, width, height, sigma);
					o = data.clone();
					ff.convolveInternal(o, width, height, sigma);
					Assert.assertArrayEquals("Internal", e, o, tolerance);
				}
			}
	}

	private abstract class MyTimingTask extends BaseTimingTask
	{
		float[][] data;
		double sigma;

		public MyTimingTask(String name, float[][] data, double sigma)
		{
			super(name);
			this.data = data;
			this.sigma = sigma;
		}

		public int getSize()
		{
			return data.length;
		}

		public Object getData(int i)
		{
			return data[i].clone();
		}

		public Object run(Object data)
		{
			float[] d = (float[]) data;
			return filter(d);
		}

		abstract float[] filter(float[] d);
	}

	private class FloatTimingTask extends MyTimingTask
	{
		GaussianFilter gf = new GaussianFilter();

		public FloatTimingTask(float[][] data, double sigma)
		{
			super("float " + sigma, data, sigma);
		}

		@Override
		float[] filter(float[] d)
		{
			gf.convolve(d, size, size, sigma);
			return d;
		}
	}

	private class DoubleTimingTask extends MyTimingTask
	{
		DoubleGaussianFilter gf = new DoubleGaussianFilter();

		public DoubleTimingTask(float[][] data, double sigma)
		{
			super("double " + sigma, data, sigma);
		}

		@Override
		float[] filter(float[] d)
		{
			gf.convolve(d, size, size, sigma);
			return d;
		}
	}

	@Test
	public void floatFilterIsFasterThanDoubleFilter()
	{
		rand = new gdsc.core.utils.Random(-30051976);

		float[][] data = new float[10][];
		for (int i = 0; i < data.length; i++)
			data[i] = createData(size, size);

		// TODO - Test with a float[] pixels array but double sums ...
		
		TimingService ts = new TimingService();
		for (double sigma : sigmas)
		{
			ts.execute(new FloatTimingTask(data, sigma));
			ts.execute(new DoubleTimingTask(data, sigma));
		}
		int size = ts.getSize();
		ts.repeat();
		ts.report(size);
		for (int i = 0, j = size; i < sigmas.length; i++, j += 2)
		{
			Assert.assertTrue(ts.get(j).getMean() < ts.get(j + 1).getMean() * 1.1);
		}
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
