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

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;

import ij.process.FloatProcessor;
import uk.ac.sussex.gdsc.core.utils.ImageExtractor;
import uk.ac.sussex.gdsc.core.utils.Maths;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.test.DataCache;
import uk.ac.sussex.gdsc.test.DataProvider;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;

@SuppressWarnings({ "javadoc" })
public class PerPixelCameraModelTest implements DataProvider<RandomSeed, Object>
{
	final static int w = 113, h = 29, size = w * h;

	private static class PerPixelCameraModelTestData
	{
		float[] bias, gain, variance, var_g2, image;
	}

	private static DataCache<RandomSeed, Object> dataCache = new DataCache<>();

	@Override
	public Object getData(RandomSeed seed)
	{
		final UniformRandomProvider r = TestSettings.getRandomGenerator(seed.getSeed());

		final PerPixelCameraModelTestData data = new PerPixelCameraModelTestData();
		data.bias = new float[size];
		data.gain = new float[size];
		data.variance = new float[size];
		data.var_g2 = new float[size];
		data.image = new float[size];
		for (int i = 0; i < size; i++)
		{
			data.bias[i] = r.nextFloat();
			data.gain[i] = 1f + r.nextFloat(); // Ensure positive
			data.variance[i] = r.nextFloat();
			data.image[i] = 100 * r.nextFloat();
			data.var_g2[i] = data.variance[i] / (data.gain[i] * data.gain[i]);
		}
		return data;
	}

	@SeededTest
	public void canGetDataWithFullBounds(RandomSeed seed)
	{
		final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache.getData(seed, this);
		final PerPixelCameraModel model = new PerPixelCameraModel(w, h, data.bias, data.gain, data.variance);
		final Rectangle bounds = new Rectangle(0, 0, w, h);
		Assertions.assertArrayEquals(data.bias, model.getBias(bounds));
		Assertions.assertArrayEquals(data.gain, model.getGain(bounds));
		Assertions.assertArrayEquals(data.variance, model.getVariance(bounds));
		Assertions.assertArrayEquals(data.var_g2, model.getNormalisedVariance(bounds));
		Assertions.assertArrayEquals(data.bias, model.getBias());
		Assertions.assertArrayEquals(data.gain, model.getGain());
		Assertions.assertArrayEquals(data.variance, model.getVariance());
		Assertions.assertArrayEquals(data.var_g2, model.getNormalisedVariance());
	}

	@SeededTest
	public void canGetCropData(RandomSeed seed)
	{
		canGetCropData(seed, true);
		canGetCropData(seed, false);
	}

	private void canGetCropData(RandomSeed seed, boolean initialise)
	{
		final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache.getData(seed, this);
		final PerPixelCameraModel model = createModel(data, initialise);
		final UniformRandomProvider rand = TestSettings.getRandomGenerator(seed.getSeed());
		final ImageExtractor ie = new ImageExtractor(data.bias, w, h);
		for (int i = 0; i < 10; i++)
		{
			final Rectangle bounds = getBounds(rand, ie);
			check(data.bias, bounds, model.getBias(bounds));
			check(data.gain, bounds, model.getGain(bounds));
			check(data.variance, bounds, model.getVariance(bounds));
			check(data.var_g2, bounds, model.getNormalisedVariance(bounds));
		}
	}

	private static Rectangle getBounds(UniformRandomProvider rand, ImageExtractor ie)
	{
		final Rectangle bounds = ie.getBoxRegionBounds(5 + rand.nextInt(w - 10), 5 + rand.nextInt(h - 10),
				2 + rand.nextInt(3));
		return bounds;
	}

	private static PerPixelCameraModel createModel(PerPixelCameraModelTestData data, boolean initialise)
	{
		final PerPixelCameraModel model = new PerPixelCameraModel(w, h, data.bias, data.gain, data.variance);
		if (initialise)
			model.initialise();
		return model;
	}

	private static void check(float[] data, Rectangle bounds, float[] o)
	{
		final FloatProcessor ip = new FloatProcessor(w, h, data.clone());
		ip.setRoi(bounds);
		final float[] e = (float[]) (ip.crop().getPixels());
		Assertions.assertArrayEquals(e, o);
	}

	@SeededTest
	public void canCropAndGetData(RandomSeed seed)
	{
		canCropAndGetData(seed, true);
		canCropAndGetData(seed, false);
	}

	private void canCropAndGetData(RandomSeed seed, boolean initialise)
	{
		final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache.getData(seed, this);
		final PerPixelCameraModel model = createModel(data, initialise);
		final UniformRandomProvider rand = TestSettings.getRandomGenerator(seed.getSeed());
		final ImageExtractor ie = new ImageExtractor(data.bias, w, h);
		for (int i = 0; i < 10; i++)
		{
			final Rectangle bounds = getBounds(rand, ie);
			final CameraModel model2 = model.crop(bounds, false);
			Assertions.assertEquals(model2.getBounds(), bounds);
			check(data.bias, bounds, model2.getBias(bounds));
			check(data.gain, bounds, model2.getGain(bounds));
			check(data.variance, bounds, model2.getVariance(bounds));
			check(data.var_g2, bounds, model2.getNormalisedVariance(bounds));
		}
	}

	@SeededTest
	public void canConvertDataWithFullBounds(RandomSeed seed)
	{
		final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache.getData(seed, this);
		final PerPixelCameraModel model = new PerPixelCameraModel(w, h, data.bias, data.gain, data.variance);
		checkConversion(data, new Rectangle(0, 0, w, h), model);
	}

	private static void checkConversion(PerPixelCameraModelTestData data, Rectangle bounds, CameraModel model)
	{
		final FloatProcessor ip = new FloatProcessor(w, h, data.image.clone());
		ip.setRoi(bounds);
		final float[] e = (float[]) (ip.crop().getPixels());
		float[] o = e.clone();
		final float[] o2 = e.clone();

		ip.setPixels(data.bias);
		final float[] bias = (float[]) (ip.crop().getPixels());

		for (int i = 0; i < e.length; i++)
			e[i] -= bias[i];
		model.removeBias(bounds, o);
		Assertions.assertArrayEquals(e, o);

		ip.setPixels(data.gain);
		final float[] gain = (float[]) (ip.crop().getPixels());

		for (int i = 0; i < e.length; i++)
			e[i] /= gain[i];
		model.removeGain(bounds, o);
		Assertions.assertArrayEquals(e, o);

		o = o2;
		model.removeBiasAndGain(bounds, o);
		Assertions.assertArrayEquals(e, o);
	}

	@SeededTest
	public void canConvertDataWithCropBounds(RandomSeed seed)
	{
		final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache.getData(seed, this);
		final PerPixelCameraModel model = new PerPixelCameraModel(w, h, data.bias, data.gain, data.variance);
		final UniformRandomProvider rand = TestSettings.getRandomGenerator(seed.getSeed());
		final ImageExtractor ie = new ImageExtractor(data.bias, w, h);
		for (int j = 0; j < 10; j++)
		{
			final Rectangle bounds = getBounds(rand, ie);
			checkConversion(data, bounds, model);
		}
	}

	@SeededTest
	public void canCropAndConvertDataWithCropBounds(RandomSeed seed)
	{
		final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache.getData(seed, this);
		final PerPixelCameraModel model = new PerPixelCameraModel(w, h, data.bias, data.gain, data.variance);
		final UniformRandomProvider rand = TestSettings.getRandomGenerator(seed.getSeed());
		final ImageExtractor ie = new ImageExtractor(data.bias, w, h);
		for (int j = 0; j < 10; j++)
		{
			final Rectangle bounds = getBounds(rand, ie);
			checkConversion(data, bounds, model.crop(bounds, false));
		}
	}

	@SeededTest
	public void canGetWeightsWithPositiveVariance(RandomSeed seed)
	{
		final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache.getData(seed, this);
		final float[] var = data.variance.clone();
		for (int i = 0; i < var.length; i++)
			if (var[i] == 0)
				var[i] = 1;
		final PerPixelCameraModel model = new PerPixelCameraModel(w, h, data.bias, data.gain, var);
		final float[] w = model.getWeights(model.getBounds());
		final float[] e = var;
		for (int i = 0; i < e.length; i++)
			e[i] = (float) (1.0 / e[i]);
		Assertions.assertArrayEquals(e, w);
	}

	@SeededTest
	public void canGetWeightsWithAllZeroVariance(RandomSeed seed)
	{
		final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache.getData(seed, this);
		final float[] var = new float[data.variance.length];
		final PerPixelCameraModel model = new PerPixelCameraModel(w, h, data.bias, data.gain, var);
		final float[] w = model.getWeights(model.getBounds());
		final float[] e = var;
		Arrays.fill(e, 1f);
		Assertions.assertArrayEquals(e, w);
	}

	@SeededTest
	public void canGetWeightsWithZeroVariance(RandomSeed seed)
	{
		final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache.getData(seed, this);
		final float[] var = data.variance.clone();
		var[0] = 0;
		final float min = SimpleArrayUtils.minAboveZero(var);
		final PerPixelCameraModel model = new PerPixelCameraModel(w, h, data.bias, data.gain, var);
		final float[] w = model.getWeights(model.getBounds());
		final float[] e = var;
		for (int i = 0; i < e.length; i++)
			e[i] = (e[i] == 0) ? (float) (1.0 / min) : (float) (1.0 / e[i]);
		Assertions.assertArrayEquals(e, w);
	}

	@SeededTest
	public void canGetMeanVariance(RandomSeed seed)
	{
		canGetMeanVariance(seed, true, false);
		canGetMeanVariance(seed, false, false);
	}

	@SeededTest
	public void canGetMeanNormalisedVariance(RandomSeed seed)
	{
		canGetMeanVariance(seed, true, true);
		canGetMeanVariance(seed, false, true);
	}

	private void canGetMeanVariance(RandomSeed seed, boolean initialise, boolean normalised)
	{
		final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache.getData(seed, this);
		final PerPixelCameraModel model = createModel(data, initialise);
		final UniformRandomProvider rand = TestSettings.getRandomGenerator(seed.getSeed());
		final ImageExtractor ie = new ImageExtractor(data.bias, w, h);
		for (int i = 0; i < 10; i++)
		{
			final Rectangle bounds = getBounds(rand, ie);
			final float[] v = (normalised) ? model.getNormalisedVariance(bounds) : model.getVariance(bounds);
			final double e = Maths.sum(v) / v.length;
			final double o = (normalised) ? model.getMeanNormalisedVariance(bounds) : model.getMeanVariance(bounds);
			Assertions.assertEquals(e, o);
		}
	}

	@SeededTest
	public void canGetCoordinateData(RandomSeed seed)
	{
		final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache.getData(seed, this);
		final int ox = 2;
		final int oy = 3;
		final int w = 8;
		final int h = 10;
		final int size = w * h;
		final float[] bias = Arrays.copyOf(data.bias, size);
		final float[] gain = Arrays.copyOf(data.gain, size);
		final float[] variance = Arrays.copyOf(data.variance, size);
		final PerPixelCameraModel model = new PerPixelCameraModel(ox, oy, w, h, bias, gain, variance);
		for (int y = 0, y1 = oy, i = 0; y < h; y++, y1++)
			for (int x = 0, x1 = ox; x < w; x++, x1++, i++)
			{
				Assertions.assertEquals(bias[i], model.getBias(x1, y1));
				Assertions.assertEquals(gain[i], model.getGain(x1, y1));
				Assertions.assertEquals(variance[i], model.getVariance(x1, y1));
				Assertions.assertEquals(data.var_g2[i], model.getNormalisedVariance(x1, y1));
			}
	}
}
