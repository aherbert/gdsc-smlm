package gdsc.smlm.ij.utils;

import java.awt.Rectangle;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.ImageExtractor;
import gdsc.core.utils.SimpleArrayUtils;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

public class ImageConverterTest
{
	static byte[] bdata;
	static short[] sdata;
	static float[] fdata;
	final static int w = 200, h = 300;
	static
	{
		RandomGenerator r = new Well19937c(30051977);
		ByteProcessor bp = new ByteProcessor(w, h);
		bdata = (byte[]) bp.getPixels();
		sdata = new short[bdata.length];
		fdata = new float[bdata.length];
		for (int i = 0; i < bp.getPixelCount(); i++)
		{
			int value = r.nextInt(256);
			bp.set(i, value);
			fdata[i] = sdata[i] = (short) value;
		}
	}

	@Test
	public void canGetData()
	{
		Rectangle bounds = null;
		float[] fe = fdata;
		Assert.assertArrayEquals(fe, ImageConverter.getData(bdata, w, h, bounds, null), 0);
		Assert.assertArrayEquals(fe, ImageConverter.getData(sdata, w, h, bounds, null), 0);
		Assert.assertArrayEquals(fe, ImageConverter.getData(fdata, w, h, bounds, null), 0);
		// Check the double format
		double[] de = SimpleArrayUtils.toDouble(fe);
		Assert.assertArrayEquals(de, ImageConverter.getDoubleData(bdata, w, h, bounds, null), 0);
		Assert.assertArrayEquals(de, ImageConverter.getDoubleData(sdata, w, h, bounds, null), 0);
		Assert.assertArrayEquals(de, ImageConverter.getDoubleData(fdata, w, h, bounds, null), 0);
	}

	@Test
	public void canGetDataWithFullBounds()
	{
		Rectangle bounds = new Rectangle(0, 0, w, h);
		float[] fe = fdata;
		Assert.assertArrayEquals(fe, ImageConverter.getData(bdata, w, h, bounds, null), 0);
		Assert.assertArrayEquals(fe, ImageConverter.getData(sdata, w, h, bounds, null), 0);
		Assert.assertArrayEquals(fe, ImageConverter.getData(fdata, w, h, bounds, null), 0);
		// Check the double format
		double[] de = SimpleArrayUtils.toDouble(fe);
		Assert.assertArrayEquals(de, ImageConverter.getDoubleData(bdata, w, h, bounds, null), 0);
		Assert.assertArrayEquals(de, ImageConverter.getDoubleData(sdata, w, h, bounds, null), 0);
		Assert.assertArrayEquals(de, ImageConverter.getDoubleData(fdata, w, h, bounds, null), 0);
	}

	@Test
	public void canGetCropData()
	{
		RandomGenerator rand = new Well19937c(30051977);
		ImageExtractor ie = new ImageExtractor(fdata, w, h);
		for (int i = 0; i < 10; i++)
		{
			Rectangle bounds = ie.getBoxRegionBounds(10 + rand.nextInt(w - 20), 10 + rand.nextInt(h - 20),
					5 + rand.nextInt(5));
			FloatProcessor ip = new FloatProcessor(w, h, fdata.clone());
			ip.setRoi(bounds);
			float[] fe = (float[]) (ip.crop().getPixels());
			Assert.assertArrayEquals(fe, ImageConverter.getData(bdata, w, h, bounds, null), 0);
			Assert.assertArrayEquals(fe, ImageConverter.getData(sdata, w, h, bounds, null), 0);
			Assert.assertArrayEquals(fe, ImageConverter.getData(fdata, w, h, bounds, null), 0);
			// Check the double format
			double[] de = SimpleArrayUtils.toDouble(fe);
			Assert.assertArrayEquals(de, ImageConverter.getDoubleData(bdata, w, h, bounds, null), 0);
			Assert.assertArrayEquals(de, ImageConverter.getDoubleData(sdata, w, h, bounds, null), 0);
			Assert.assertArrayEquals(de, ImageConverter.getDoubleData(fdata, w, h, bounds, null), 0);
		}
	}
}
