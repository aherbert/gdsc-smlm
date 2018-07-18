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
package uk.ac.sussex.gdsc.smlm.model.camera;

import java.awt.Rectangle;
import java.util.Arrays;

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import ij.process.FloatProcessor;
import uk.ac.sussex.gdsc.core.utils.ImageExtractor;
import uk.ac.sussex.gdsc.core.utils.Maths;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.test.TestSettings;

@SuppressWarnings({ "javadoc" })
public class PerPixelCameraModelTest
{
	static float[] bias, gain, variance, var_g2, image;

	final static int w = 113, h = 29, size;
	static
	{
		final RandomGenerator r = TestSettings.getRandomGenerator();
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
		final PerPixelCameraModel model = new PerPixelCameraModel(w, h, bias, gain, variance);
		final Rectangle bounds = new Rectangle(0, 0, w, h);
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

	private static void canGetCropData(boolean initialise)
	{
		final PerPixelCameraModel model = createModel(initialise);
		final RandomGenerator rand = TestSettings.getRandomGenerator();
		final ImageExtractor ie = new ImageExtractor(bias, w, h);
		for (int i = 0; i < 10; i++)
		{
			final Rectangle bounds = getBounds(rand, ie);
			check(bias, bounds, model.getBias(bounds));
			check(gain, bounds, model.getGain(bounds));
			check(variance, bounds, model.getVariance(bounds));
			check(var_g2, bounds, model.getNormalisedVariance(bounds));
		}
	}

	private static Rectangle getBounds(RandomGenerator rand, ImageExtractor ie)
	{
		final Rectangle bounds = ie.getBoxRegionBounds(5 + rand.nextInt(w - 10), 5 + rand.nextInt(h - 10),
				2 + rand.nextInt(3));
		return bounds;
	}

	private static PerPixelCameraModel createModel(boolean initialise)
	{
		final PerPixelCameraModel model = new PerPixelCameraModel(w, h, bias, gain, variance);
		if (initialise)
			model.initialise();
		return model;
	}

	private static void check(float[] data, Rectangle bounds, float[] o)
	{
		final FloatProcessor ip = new FloatProcessor(w, h, data.clone());
		ip.setRoi(bounds);
		final float[] e = (float[]) (ip.crop().getPixels());
		Assert.assertArrayEquals(e, o, 0);
	}

	@Test
	public void canCropAndGetData()
	{
		canCropAndGetData(true);
		canCropAndGetData(false);
	}

	private static void canCropAndGetData(boolean initialise)
	{
		final PerPixelCameraModel model = createModel(initialise);
		final RandomGenerator rand = TestSettings.getRandomGenerator();
		final ImageExtractor ie = new ImageExtractor(bias, w, h);
		for (int i = 0; i < 10; i++)
		{
			final Rectangle bounds = getBounds(rand, ie);
			final CameraModel model2 = model.crop(bounds, false);
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
		final PerPixelCameraModel model = new PerPixelCameraModel(w, h, bias, gain, variance);
		checkConversion(new Rectangle(0, 0, w, h), model);
	}

	private static void checkConversion(Rectangle bounds, CameraModel model)
	{
		final FloatProcessor ip = new FloatProcessor(w, h, image.clone());
		ip.setRoi(bounds);
		final float[] e = (float[]) (ip.crop().getPixels());
		float[] o = e.clone();
		final float[] o2 = e.clone();

		ip.setPixels(bias);
		final float[] bias = (float[]) (ip.crop().getPixels());

		for (int i = 0; i < e.length; i++)
			e[i] -= bias[i];
		model.removeBias(bounds, o);
		Assert.assertArrayEquals(e, o, 0);

		ip.setPixels(gain);
		final float[] gain = (float[]) (ip.crop().getPixels());

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
		final PerPixelCameraModel model = new PerPixelCameraModel(w, h, bias, gain, variance);
		final RandomGenerator rand = TestSettings.getRandomGenerator();
		final ImageExtractor ie = new ImageExtractor(bias, w, h);
		for (int j = 0; j < 10; j++)
		{
			final Rectangle bounds = getBounds(rand, ie);
			checkConversion(bounds, model);
		}
	}

	@Test
	public void canCropAndConvertDataWithCropBounds()
	{
		final PerPixelCameraModel model = new PerPixelCameraModel(w, h, bias, gain, variance);
		final RandomGenerator rand = TestSettings.getRandomGenerator();
		final ImageExtractor ie = new ImageExtractor(bias, w, h);
		for (int j = 0; j < 10; j++)
		{
			final Rectangle bounds = getBounds(rand, ie);
			checkConversion(bounds, model.crop(bounds, false));
		}
	}

	@Test
	public void canGetWeightsWithPositiveVariance()
	{
		final float[] var = variance.clone();
		for (int i = 0; i < var.length; i++)
			if (var[i] == 0)
				var[i] = 1;
		final PerPixelCameraModel model = new PerPixelCameraModel(w, h, bias, gain, var);
		final float[] w = model.getWeights(model.getBounds());
		final float[] e = var;
		for (int i = 0; i < e.length; i++)
			e[i] = (float) (1.0 / e[i]);
		Assert.assertArrayEquals(e, w, 0);
	}

	@Test
	public void canGetWeightsWithAllZeroVariance()
	{
		final float[] var = new float[variance.length];
		final PerPixelCameraModel model = new PerPixelCameraModel(w, h, bias, gain, var);
		final float[] w = model.getWeights(model.getBounds());
		final float[] e = var;
		Arrays.fill(e, 1f);
		Assert.assertArrayEquals(e, w, 0);
	}

	@Test
	public void canGetWeightsWithZeroVariance()
	{
		final float[] var = variance.clone();
		var[0] = 0;
		final float min = SimpleArrayUtils.minAboveZero(var);
		final PerPixelCameraModel model = new PerPixelCameraModel(w, h, bias, gain, var);
		final float[] w = model.getWeights(model.getBounds());
		final float[] e = var;
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

	private static void canGetMeanVariance(boolean initialise, boolean normalised)
	{
		final PerPixelCameraModel model = createModel(initialise);
		final RandomGenerator rand = TestSettings.getRandomGenerator();
		final ImageExtractor ie = new ImageExtractor(bias, w, h);
		for (int i = 0; i < 10; i++)
		{
			final Rectangle bounds = getBounds(rand, ie);
			final float[] v = (normalised) ? model.getNormalisedVariance(bounds) : model.getVariance(bounds);
			final double e = Maths.sum(v) / v.length;
			final double o = (normalised) ? model.getMeanNormalisedVariance(bounds) : model.getMeanVariance(bounds);
			Assert.assertEquals(e, o, 0);
		}
	}

	@Test
	public void canGetCoordinateData()
	{
		final int ox = 2;
		final int oy = 3;
		final int w = 8;
		final int h = 10;
		final int size = w * h;
		final float[] bias = Arrays.copyOf(PerPixelCameraModelTest.bias, size);
		final float[] gain = Arrays.copyOf(PerPixelCameraModelTest.gain, size);
		final float[] variance = Arrays.copyOf(PerPixelCameraModelTest.variance, size);
		final PerPixelCameraModel model = new PerPixelCameraModel(ox, oy, w, h, bias, gain, variance);
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
