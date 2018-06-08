package gdsc.smlm.model.camera;

import java.awt.Rectangle;
import java.util.Arrays;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.ImageExtractor;
import gdsc.core.utils.Maths;
import gdsc.core.utils.SimpleArrayUtils;
import ij.process.FloatProcessor;

public class PerPixelCameraModelTest
{
	static float[] bias, gain, variance, var_g2, image;

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
	}

	@Test
	public void canGetDataWithFullBounds()
	{
		PerPixelCameraModel model = new PerPixelCameraModel(w, h, bias, gain, variance);
		Rectangle bounds = new Rectangle(0, 0, w, h);
		Assert.assertArrayEquals(bias, model.getBias(bounds), 0);
		Assert.assertArrayEquals(gain, model.getGain(bounds), 0);
		Assert.assertArrayEquals(variance, model.getVariance(bounds), 0);
		Assert.assertArrayEquals(var_g2, model.getNormalisedVariance(bounds), 0);
		Assert.assertArrayEquals(bias, model.getBias(), 0);
		Assert.assertArrayEquals(gain, model.getGain(), 0);
		Assert.assertArrayEquals(variance, model.getVariance(), 0);
		Assert.assertArrayEquals(var_g2, model.getNormalisedVariance(), 0);
	}

	@Test
	public void canGetCropData()
	{
		canGetCropData(true);
		canGetCropData(false);
	}

	private void canGetCropData(boolean initialise)
	{
		PerPixelCameraModel model = createModel(initialise);
		RandomGenerator rand = new Well19937c(30051977);
		ImageExtractor ie = new ImageExtractor(bias, w, h);
		for (int i = 0; i < 10; i++)
		{
			Rectangle bounds = ie.getBoxRegionBounds(10 + rand.nextInt(w - 20), 10 + rand.nextInt(h - 20),
					5 + rand.nextInt(5));
			check(bias, bounds, model.getBias(bounds));
			check(gain, bounds, model.getGain(bounds));
			check(variance, bounds, model.getVariance(bounds));
			check(var_g2, bounds, model.getNormalisedVariance(bounds));
		}
	}

	private PerPixelCameraModel createModel(boolean initialise)
	{
		PerPixelCameraModel model = new PerPixelCameraModel(w, h, bias, gain, variance);
		if (initialise)
			model.initialise();
		return model;
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
		canCropAndGetData(true);
		canCropAndGetData(false);
	}

	private void canCropAndGetData(boolean initialise)
	{
		PerPixelCameraModel model = createModel(initialise);
		RandomGenerator rand = new Well19937c(30051977);
		ImageExtractor ie = new ImageExtractor(bias, w, h);
		for (int i = 0; i < 10; i++)
		{
			Rectangle bounds = ie.getBoxRegionBounds(10 + rand.nextInt(w - 20), 10 + rand.nextInt(h - 20),
					5 + rand.nextInt(5));
			CameraModel model2 = model.crop(bounds, false);
			Assert.assertEquals(model2.getBounds(), bounds);
			check(bias, bounds, model2.getBias(bounds));
			check(gain, bounds, model2.getGain(bounds));
			check(variance, bounds, model2.getVariance(bounds));
			check(var_g2, bounds, model2.getNormalisedVariance(bounds));
		}
	}

	@Test
	public void canConvertDataWithFullBounds()
	{
		PerPixelCameraModel model = new PerPixelCameraModel(w, h, bias, gain, variance);
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
		PerPixelCameraModel model = new PerPixelCameraModel(w, h, bias, gain, variance);
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
		PerPixelCameraModel model = new PerPixelCameraModel(w, h, bias, gain, variance);
		RandomGenerator rand = new Well19937c(30051977);
		ImageExtractor ie = new ImageExtractor(bias, w, h);
		for (int j = 0; j < 10; j++)
		{
			Rectangle bounds = ie.getBoxRegionBounds(10 + rand.nextInt(w - 20), 10 + rand.nextInt(h - 20),
					5 + rand.nextInt(5));
			checkConversion(bounds, model.crop(bounds, false));
		}
	}

	@Test
	public void canGetWeightsWithPositiveVariance()
	{
		float[] var = variance.clone();
		for (int i = 0; i < var.length; i++)
			if (var[i] == 0)
				var[i] = 1;
		PerPixelCameraModel model = new PerPixelCameraModel(w, h, bias, gain, var);
		float[] w = model.getWeights(model.getBounds());
		float[] e = var;
		for (int i = 0; i < e.length; i++)
			e[i] = (float) (1.0 / e[i]);
		Assert.assertArrayEquals(e, w, 0);
	}

	@Test
	public void canGetWeightsWithAllZeroVariance()
	{
		float[] var = new float[variance.length];
		PerPixelCameraModel model = new PerPixelCameraModel(w, h, bias, gain, var);
		float[] w = model.getWeights(model.getBounds());
		float[] e = var;
		Arrays.fill(e, 1f);
		Assert.assertArrayEquals(e, w, 0);
	}

	@Test
	public void canGetWeightsWithZeroVariance()
	{
		float[] var = variance.clone();
		var[0] = 0;
		float min = SimpleArrayUtils.minAboveZero(var);
		PerPixelCameraModel model = new PerPixelCameraModel(w, h, bias, gain, var);
		float[] w = model.getWeights(model.getBounds());
		float[] e = var;
		for (int i = 0; i < e.length; i++)
			e[i] = (e[i] == 0) ? (float) (1.0 / min) : (float) (1.0 / e[i]);
		Assert.assertArrayEquals(e, w, 0);
	}

	@Test
	public void canGetMeanVariance()
	{
		canGetMeanVariance(true, false);
		canGetMeanVariance(false, false);
	}

	@Test
	public void canGetMeanNormalisedVariance()
	{
		canGetMeanVariance(true, true);
		canGetMeanVariance(false, true);
	}

	private void canGetMeanVariance(boolean initialise, boolean normalised)
	{
		PerPixelCameraModel model = createModel(initialise);
		RandomGenerator rand = new Well19937c(30051977);
		ImageExtractor ie = new ImageExtractor(bias, w, h);
		for (int i = 0; i < 10; i++)
		{
			Rectangle bounds = ie.getBoxRegionBounds(10 + rand.nextInt(w - 20), 10 + rand.nextInt(h - 20),
					5 + rand.nextInt(5));
			float[] v = (normalised) ? model.getNormalisedVariance(bounds) : model.getVariance(bounds);
			double e = Maths.sum(v) / v.length;
			double o = (normalised) ? model.getMeanNormalisedVariance(bounds) : model.getMeanVariance(bounds);
			Assert.assertEquals(e, o, 0);
		}
	}

	@Test
	public void canGetCoordinateData()
	{
		int ox = 2;
		int oy = 3;
		int w = 8;
		int h = 10;
		int size = w * h;
		float[] bias = Arrays.copyOf(PerPixelCameraModelTest.bias, size);
		float[] gain = Arrays.copyOf(PerPixelCameraModelTest.gain, size);
		float[] variance = Arrays.copyOf(PerPixelCameraModelTest.variance, size);
		PerPixelCameraModel model = new PerPixelCameraModel(ox, oy, w, h, bias, gain, variance);
		for (int y = 0, y1 = oy, i = 0; y < h; y++, y1++)
			for (int x = 0, x1 = ox; x < w; x++, x1++, i++)
			{
				Assert.assertEquals(bias[i], model.getBias(x1, y1), 0);
				Assert.assertEquals(gain[i], model.getGain(x1, y1), 0);
				Assert.assertEquals(variance[i], model.getVariance(x1, y1), 0);
				Assert.assertEquals(var_g2[i], model.getNormalisedVariance(x1, y1), 0);
			}
	}
}
