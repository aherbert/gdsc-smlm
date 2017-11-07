package gdsc.smlm.filters;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.jtransforms.dht.FloatDHT_2D;
import org.jtransforms.fft.FloatFFT_2D;
import org.jtransforms.utils.CommonUtils;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.test.BaseTimingTask;
import gdsc.core.test.TimingService;
import ij.plugin.filter.EDM;
import ij.process.ByteProcessor;
import ij.process.FHT2;
import ij.process.FloatProcessor;

public class JTransformsTest
{
	private FloatProcessor createProcessor(int size, int x, int y, int w, int h, RandomGenerator r)
	{
		ByteProcessor bp = new ByteProcessor(size, size);
		bp.setColor(255);
		bp.fillOval(x, y, w, h);
		EDM e = new EDM();
		FloatProcessor fp = e.makeFloatEDM(bp, 0, true);
		if (r != null)
		{
			float[] d = (float[]) fp.getPixels();
			for (int i = 0; i < d.length; i++)
				d[i] += r.nextFloat() * 0.01;
		}
		return fp;
	}

	@Test
	public void canCorrelateUsingFFT()
	{
		canComputeUsingFFT(false);
	}

	@Test
	public void canConvolveUsingFFT()
	{
		canComputeUsingFFT(true);
	}

	private void canComputeUsingFFT(boolean convolution)
	{
		int size = 16;
		int ex = 5, ey = 7;
		int ox = 1, oy = 2;
		FloatProcessor fp1 = createProcessor(size, ex, ey, 4, 4, null);
		FloatProcessor fp2 = createProcessor(size, size / 2 + ox, size / 2 + oy, 4, 4, null);

		float[] input1 = (float[]) fp1.getPixels();
		float[] input2 = (float[]) fp2.getPixels();

		FHTFilter ff = new FHTFilter(input2.clone(), size, size);
		ff.setConvolution(convolution);
		float[] e = input1.clone();
		ff.filter(e, size, size);

		// Do the same with JTransforms
		float[] data1 = new float[input1.length * 2];
		FloatFFT_2D fft = new FloatFFT_2D(size, size);
		System.arraycopy(input1, 0, data1, 0, input1.length);
		float[] data2 = new float[data1.length];
		System.arraycopy(input2, 0, data2, 0, input2.length);

		fft.realForwardFull(data1);
		fft.realForwardFull(data2);

		// Multiply
		// https://en.wikipedia.org/wiki/Complex_number#Multiplication_and_division
		for (int i = 0; i < data2.length; i += 2)
		{
			float a = data1[i];
			float b = data1[i + 1];
			float c = data2[i];
			float d = (convolution) ? data2[i + 1] : -data2[i + 1]; // Get the conjugate for correlation
			data1[i] = a * c - b * d;
			data1[i + 1] = b * c + a * d;
		}

		fft.complexInverse(data1, true);

		float[] o = new float[e.length];
		for (int i = 0, j = 0; i < o.length; i++, j += 2)
			o[i] = data1[j];
		FHT2.swapQuadrants(new FloatProcessor(size, size, o));

		Assert.assertArrayEquals(e, o, 1e-3f);
	}

	@Test
	public void canComputeFHTUsingJTransforms()
	{
		// Note: no need to test the correlation as the transformed data 
		// is the same format as FHT so we just test that.
		
		int size = 16;
		int ex = 5, ey = 7;
		int ox = 1, oy = 2;
		FloatProcessor fp1 = createProcessor(size, ex, ey, 4, 4, null);
		FloatProcessor fp2 = createProcessor(size, size / 2 + ox, size / 2 + oy, 4, 4, null);

		float[] input1 = (float[]) fp1.getPixels();
		float[] input2 = (float[]) fp2.getPixels();

		FHT2 fht1 = new FHT2(input1.clone(), size, false);
		FHT2 fht2 = new FHT2(input2.clone(), size, false);

		fht1.transform();
		fht2.transform();

		// Do the same with JTransforms
		FloatDHT_2D dht = new FloatDHT_2D(size, size);

		dht.forward(input1);
		dht.forward(input2);

		Assert.assertArrayEquals(fht1.getData(), input1, 1e-5f);
		Assert.assertArrayEquals(fht2.getData(), input2, 1e-5f);
	}

	private abstract class DHTSpeedTask extends BaseTimingTask
	{
		int maxN;
		float[][] data;

		public DHTSpeedTask(String name, int maxN, float[][] data)
		{
			super(name);
			this.maxN = maxN;
			this.data = data;
		}

		public int getSize()
		{
			return 1;
		}

		public Object getData(int i)
		{
			return clone(data);
		}

		private float[][] clone(float[][] data)
		{
			float[][] data2 = new float[data.length][];
			for (int i = 0; i < data.length; i++)
				data2[i] = data[i].clone();
			return data2;
		}

		public Object run(Object data)
		{
			return run((float[][]) data);
		}

		abstract Object run(float[][] data);
	}

	private class IJDHTSpeedTask extends DHTSpeedTask
	{
		FHT2 fht2;

		public IJDHTSpeedTask(int maxN, float[][] data)
		{
			super(FHT2.class.getSimpleName(), maxN, data);
			// Create one so we have the pre-computed tables
			fht2 = new FHT2(data[0].clone(), maxN, false);
			fht2.transform();
		}

		@Override
		Object run(float[][] data)
		{
			for (int i = 0; i < data.length; i += 2)
			{
				// Forward
				FHT2 fht = new FHT2(data[i], maxN, false);
				fht.copyTables(fht2);
				fht.transform();
				// Reverse
				fht = new FHT2(data[i + 1], maxN, true);
				fht.copyTables(fht2);
				fht.transform();
			}
			return null;
		}
	}

	private class JTransformsDHTSpeedTask extends DHTSpeedTask
	{
		FloatDHT_2D dht;

		public JTransformsDHTSpeedTask(int maxN, float[][] data)
		{
			super(FloatDHT_2D.class.getSimpleName(), maxN, data);
			dht = new FloatDHT_2D(maxN, maxN);
		}

		@Override
		Object run(float[][] data)
		{
			for (int i = 0; i < data.length; i += 2)
			{
				// Forward
				dht.forward(data[i]);
				// Reverse
				dht.inverse(data[i + 1], true);
			}
			return null;
		}
	}

	@Test
	public void jTransforms2DDHTIsFasterThanFHT2()
	{
		// Test the forward DHT of data. and reverse transform or the pre-computed correlation.

		int size = 256;
		int w = size / 4;
		RandomGenerator r = new Well19937c(30051977);
		RandomDataGenerator rdg = new RandomDataGenerator(r);

		// Blob in the centre
		FloatProcessor fp = createProcessor(size, size / 2 - w / 2, size / 2 - w / 2, w, w, null);
		FHT2 fht2 = new FHT2((float[]) fp.getPixels(), size, false);
		fht2.transform();
		fht2.initialiseFastMultiply();

		// Random blobs, original and correlated
		int N = 40;
		float[][] data = new float[N * 2][];
		int lower = w;
		int upper = size - w;
		for (int i = 0, j = 0; i < N; i++)
		{
			int x = rdg.nextInt(lower, upper);
			int y = rdg.nextInt(lower, upper);
			fp = createProcessor(size, x, y, w, w, r);
			float[] pixels = (float[]) fp.getPixels();
			data[j++] = pixels.clone();
			FHT2 fht1 = new FHT2(pixels, size, false);
			fht1.copyTables(fht2);
			fht2.transform();
			float[] pixels2 = new float[pixels.length];
			fht2.conjugateMultiply(fht2, pixels2);
			data[j++] = pixels2;
		}
		
		//CommonUtils.setThreadsBeginN_1D_FFT_2Threads(Long.MAX_VALUE);
		//CommonUtils.setThreadsBeginN_1D_FFT_4Threads(Long.MAX_VALUE);
		CommonUtils.setThreadsBeginN_2D(Long.MAX_VALUE);
		
		TimingService ts = new TimingService();
		ts.execute(new IJDHTSpeedTask(size, data));
		ts.execute(new JTransformsDHTSpeedTask(size, data));
		ts.repeat();
		ts.report();
		
		Assert.assertTrue(ts.get(-1).getMean() < ts.get(-2).getMean());
	}
}
