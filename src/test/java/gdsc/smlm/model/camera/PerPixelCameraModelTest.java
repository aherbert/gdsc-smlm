package gdsc.smlm.model.camera;

import java.awt.Rectangle;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.ImageExtractor;
import ij.process.FloatProcessor;

public class PerPixelCameraModelTest
{
	static float[] bias, gain, variance, var_g2, image;
	static PerPixelCameraModel model;

	final static int w = 200, h = 300, size;
	static
	{
		RandomGenerator r = new Well19937c(30051977);
		size = w * h;
		bias = new float[size];
		gain = new float[size];
		variance = new float[size];
		var_g2 = new float[size];
		image = new float[size];
		for (int i = 0; i < size; i++)
		{
			bias[i] = r.nextFloat();
			gain[i] = 1f + r.nextFloat(); // Ensure positive
			variance[i] = r.nextFloat();
			image[i] = 100 * r.nextFloat();
			var_g2[i] = variance[i] / (gain[i] * gain[i]);
		}
		model = new PerPixelCameraModel(w, h, bias, gain, variance);
	}

	@Test
	public void canGetDataWithFullBounds()
	{
		Rectangle bounds = new Rectangle(0, 0, w, h);
		Assert.assertArrayEquals(bias, model.getBias(bounds), 0);
		Assert.assertArrayEquals(gain, model.getGain(bounds), 0);
		Assert.assertArrayEquals(var_g2, model.getNormalisedVariance(bounds), 0);
	}

	@Test
	public void canGetCropData()
	{
		RandomGenerator rand = new Well19937c(30051977);
		ImageExtractor ie = new ImageExtractor(bias, w, h);
		for (int i = 0; i < 10; i++)
		{
			Rectangle bounds = ie.getBoxRegionBounds(10 + rand.nextInt(w - 20), 10 + rand.nextInt(h - 20),
					5 + rand.nextInt(5));
			check(bias, bounds, model.getBias(bounds));
			check(gain, bounds, model.getGain(bounds));
			check(var_g2, bounds, model.getNormalisedVariance(bounds));
		}
	}

	private void check(float[] data, Rectangle bounds, float[] o)
	{
		FloatProcessor ip = new FloatProcessor(w, h, data.clone());
		ip.setRoi(bounds);
		float[] e = (float[]) (ip.crop().getPixels());
		Assert.assertArrayEquals(e, o, 0);
	}

	@Test
	public void canCropAndGetData()
	{
		RandomGenerator rand = new Well19937c(30051977);
		ImageExtractor ie = new ImageExtractor(bias, w, h);
		for (int i = 0; i < 10; i++)
		{
			Rectangle bounds = ie.getBoxRegionBounds(10 + rand.nextInt(w - 20), 10 + rand.nextInt(h - 20),
					5 + rand.nextInt(5));
			CameraModel model2 = model.crop(bounds);
			Assert.assertEquals(model2.getBounds(), bounds);
			check(bias, bounds, model2.getBias(bounds));
			check(gain, bounds, model2.getGain(bounds));
			check(var_g2, bounds, model2.getNormalisedVariance(bounds));
		}
	}
	
	@Test
	public void canConvertDataWithFullBounds()
	{
		checkConversion(new Rectangle(0, 0, w, h), model);
	}

	private void checkConversion(Rectangle bounds, CameraModel model)
	{
		FloatProcessor ip = new FloatProcessor(w, h, image.clone());
		ip.setRoi(bounds);
		float[] e = (float[]) (ip.crop().getPixels());
		float[] o = e.clone();
		float[] o2 = e.clone();

		ip.setPixels(bias);
		float[] bias = (float[]) (ip.crop().getPixels());

		for (int i = 0; i < e.length; i++)
			e[i] -= bias[i];
		model.removeBias(bounds, o);
		Assert.assertArrayEquals(e, o, 0);

		ip.setPixels(gain);
		float[] gain = (float[]) (ip.crop().getPixels());
		
		for (int i = 0; i < e.length; i++)
			e[i] /= gain[i];
		model.removeGain(bounds, o);
		Assert.assertArrayEquals(e, o, 0);

		o = o2;
		model.removeBiasAndGain(bounds, o);
		Assert.assertArrayEquals(e, o, 0);
	}

	@Test
	public void canConvertDataWithCropBounds()
	{
		RandomGenerator rand = new Well19937c(30051977);
		ImageExtractor ie = new ImageExtractor(bias, w, h);
		for (int j = 0; j < 10; j++)
		{
			Rectangle bounds = ie.getBoxRegionBounds(10 + rand.nextInt(w - 20), 10 + rand.nextInt(h - 20),
					5 + rand.nextInt(5));
			checkConversion(bounds, model);
		}
	}
	
	@Test
	public void canCropAndConvertDataWithCropBounds()
	{
		RandomGenerator rand = new Well19937c(30051977);
		ImageExtractor ie = new ImageExtractor(bias, w, h);
		for (int j = 0; j < 10; j++)
		{
			Rectangle bounds = ie.getBoxRegionBounds(10 + rand.nextInt(w - 20), 10 + rand.nextInt(h - 20),
					5 + rand.nextInt(5));
			checkConversion(bounds, model.crop(bounds));
		}
	}
}
