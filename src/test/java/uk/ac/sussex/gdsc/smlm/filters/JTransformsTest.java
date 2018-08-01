package uk.ac.sussex.gdsc.smlm.filters;

import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;

import org.apache.commons.rng.UniformRandomProvider;
import org.jtransforms.dht.FloatDHT_2D;
import org.jtransforms.fft.FloatFFT_2D;
import org.jtransforms.utils.CommonUtils;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import ij.plugin.filter.EDM;
import ij.process.ByteProcessor;
import ij.process.FHT;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import uk.ac.sussex.gdsc.core.ij.process.FHT2;
import uk.ac.sussex.gdsc.smlm.filters.FHTFilter.Operation;
import uk.ac.sussex.gdsc.test.BaseTimingTask;
import uk.ac.sussex.gdsc.test.LogLevel;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.TimingService;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssumptions;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;

@SuppressWarnings({ "javadoc" })
public class JTransformsTest
{
    private static Logger logger;

    @BeforeAll
    public static void beforeAll()
    {
        logger = Logger.getLogger(JTransformsTest.class.getName());
    }

    @AfterAll
    public static void afterAll()
    {
        logger = null;
    }

	private static FloatProcessor createProcessor(int size, int x, int y, int w, int h, UniformRandomProvider r)
	{
		final ByteProcessor bp = new ByteProcessor(size, size);
		bp.setColor(255);
		bp.fillOval(x, y, w, h);
		final EDM e = new EDM();
		final FloatProcessor fp = e.makeFloatEDM(bp, 0, true);
		if (r != null)
		{
			final float[] d = (float[]) fp.getPixels();
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

	private static void canComputeUsingFFT(boolean convolution)
	{
		final int size = 16;
		final int ex = 5, ey = 7;
		final int ox = 1, oy = 2;
		final FloatProcessor fp1 = createProcessor(size, ex, ey, 4, 4, null);
		final FloatProcessor fp2 = createProcessor(size, size / 2 + ox, size / 2 + oy, 4, 4, null);

		final float[] input1 = (float[]) fp1.getPixels();
		final float[] input2 = (float[]) fp2.getPixels();

		final FHTFilter ff = new FHTFilter(input2.clone(), size, size);
		ff.setOperation((convolution) ? Operation.CONVOLUTION : Operation.CORRELATION);
		final float[] e = input1.clone();
		ff.filter(e, size, size);

		// Do the same with JTransforms
		final float[] data1 = new float[input1.length * 2];
		final FloatFFT_2D fft = new FloatFFT_2D(size, size);
		System.arraycopy(input1, 0, data1, 0, input1.length);
		final float[] data2 = new float[data1.length];
		System.arraycopy(input2, 0, data2, 0, input2.length);

		fft.realForwardFull(data1);
		fft.realForwardFull(data2);

		// Multiply
		// https://en.wikipedia.org/wiki/Complex_number#Multiplication_and_division
		for (int i = 0; i < data2.length; i += 2)
		{
			final float a = data1[i];
			final float b = data1[i + 1];
			final float c = data2[i];
			final float d = (convolution) ? data2[i + 1] : -data2[i + 1]; // Get the conjugate for correlation
			data1[i] = a * c - b * d;
			data1[i + 1] = b * c + a * d;
		}

		fft.complexInverse(data1, true);

		final float[] o = new float[e.length];
		for (int i = 0, j = 0; i < o.length; i++, j += 2)
			o[i] = data1[j];
		FHT2.swapQuadrants(new FloatProcessor(size, size, o));

		Assertions.assertArrayEquals(e, o, 1e-3f);
	}

	@Test
	public void canComputeFHTUsingJTransforms()
	{
		// Note: no need to test the correlation as the transformed data
		// is the same format as FHT so we just test that.

		final int size = 16;
		final int ex = 5, ey = 7;
		final int ox = 1, oy = 2;
		final FloatProcessor fp1 = createProcessor(size, ex, ey, 4, 4, null);
		final FloatProcessor fp2 = createProcessor(size, size / 2 + ox, size / 2 + oy, 4, 4, null);

		final float[] input1 = (float[]) fp1.getPixels();
		final float[] input2 = (float[]) fp2.getPixels();

		final FHT2 fht1 = new FHT2(input1.clone(), size, false);
		final FHT2 fht2 = new FHT2(input2.clone(), size, false);

		fht1.transform();
		fht2.transform();

		// Do the same with JTransforms
		final FloatDHT_2D dht = new FloatDHT_2D(size, size);

		dht.forward(input1);
		dht.forward(input2);

		Assertions.assertArrayEquals(fht1.getData(), input1, 1e-5f);
		Assertions.assertArrayEquals(fht2.getData(), input2, 1e-5f);
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

		@Override
		public int getSize()
		{
			return 1;
		}

		@Override
		public Object getData(int i)
		{
			return clone(data);
		}

		private float[][] clone(float[][] data)
		{
			final float[][] data2 = new float[data.length][];
			for (int i = 0; i < data.length; i++)
				data2[i] = data[i].clone();
			return data2;
		}

		@Override
		public Object run(Object data)
		{
			return run((float[][]) data);
		}

		abstract Object run(float[][] data);
	}

	private class NonDuplicatingFloatProcessor extends FloatProcessor
	{
		public NonDuplicatingFloatProcessor(int width, int height, float[] pixels)
		{
			super(width, height, pixels);
		}

		@Override
		public ImageProcessor duplicate()
		{
			return this;
		}

		@Override
		public ImageProcessor convertToFloat()
		{
			return this;
		}
	}

	private class IJFHTSpeedTask extends DHTSpeedTask
	{
		public IJFHTSpeedTask(int maxN, float[][] data)
		{
			super(FHT.class.getSimpleName(), maxN, data);
		}

		@Override
		Object run(float[][] data)
		{
			for (int i = 0; i < data.length; i += 2)
			{
				// Forward
				FHT fht = new FHT(new NonDuplicatingFloatProcessor(maxN, maxN, data[i]), false);
				fht.transform();
				// Reverse
				fht = new FHT(new NonDuplicatingFloatProcessor(maxN, maxN, data[i + 1]), true);
				fht.transform();
			}
			return null;
		}
	}

	private class IJFHT2SpeedTask extends DHTSpeedTask
	{
		FHT2 fht2;

		public IJFHT2SpeedTask(int maxN, float[][] data)
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

	@SpeedTag
	@SeededTest
	public void jTransforms2DDHTIsFasterThanFHT2(RandomSeed seed)
	{
		ExtraAssumptions.assumeSpeedTest();

		// Test the forward DHT of data. and reverse transform or the pre-computed correlation.

		final int size = 256;
		final int w = size / 4;
		final UniformRandomProvider r = TestSettings.getRandomGenerator(seed.getSeed());

		// Blob in the centre
		FloatProcessor fp = createProcessor(size, size / 2 - w / 2, size / 2 - w / 2, w, w, null);
		final FHT2 fht2 = new FHT2((float[]) fp.getPixels(), size, false);
		fht2.transform();
		fht2.initialiseFastMultiply();

		// Random blobs, original and correlated
		final int N = 40;
		final float[][] data = new float[N * 2][];
		final int lower = w;
		final int upper = size - w;
		final int range = upper - lower;
		for (int i = 0, j = 0; i < N; i++)
		{
			final int x = lower + r.nextInt(range);
			final int y = lower + r.nextInt(range);
			fp = createProcessor(size, x, y, w, w, r);
			final float[] pixels = (float[]) fp.getPixels();
			data[j++] = pixels.clone();
			final FHT2 fht1 = new FHT2(pixels, size, false);
			fht1.copyTables(fht2);
			fht2.transform();
			final float[] pixels2 = new float[pixels.length];
			fht2.conjugateMultiply(fht2, pixels2);
			data[j++] = pixels2;
		}

		//CommonUtils.setThreadsBeginN_1D_FFT_2Threads(Long.MAX_VALUE);
		//CommonUtils.setThreadsBeginN_1D_FFT_4Threads(Long.MAX_VALUE);
		CommonUtils.setThreadsBeginN_2D(Long.MAX_VALUE);

		final TimingService ts = new TimingService();
		ts.execute(new IJFHTSpeedTask(size, data));
		ts.execute(new IJFHT2SpeedTask(size, data));
		ts.execute(new JTransformsDHTSpeedTask(size, data));
		ts.repeat();
		if (logger.isLoggable(Level.INFO))
			ts.report();

		//Assertions.assertTrue(ts.get(-1).getMean() < ts.get(-2).getMean());

		final double t1 = ts.get(-1).getMean();
		final double t2 = ts.get(-2).getMean();
		TestLog.logTestResult(logger,t1 < t2, "%s %s => %s %s = %.2fx\n", ts.get(-2).getTask().getName(), t2,
				ts.get(-1).getTask().getName(), t1, t2 / t1);

	}
}
