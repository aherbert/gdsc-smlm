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
package uk.ac.sussex.gdsc.smlm.ij.utils;

import java.awt.Rectangle;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;

import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import uk.ac.sussex.gdsc.core.utils.ImageExtractor;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.test.DataCache;
import uk.ac.sussex.gdsc.test.DataProvider;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;

@SuppressWarnings({ "javadoc" })
public class ImageConverterTest implements DataProvider<RandomSeed, Object>
{
	final static int w = 200, h = 300;

	private class ImageConverterTestData
	{
		byte[] bdata;
		short[] sdata;
		float[] fdata;
	}

	private static DataCache<RandomSeed, Object> dataCache = new DataCache<>();

	@Override
	public Object getData(RandomSeed seed)
	{
		final UniformRandomProvider r = TestSettings.getRandomGenerator(seed.getSeed());
		final ByteProcessor bp = new ByteProcessor(w, h);
		ImageConverterTestData data = new ImageConverterTestData();
		data.bdata = (byte[]) bp.getPixels();
		data.sdata = new short[data.bdata.length];
		data.fdata = new float[data.bdata.length];
		for (int i = 0; i < bp.getPixelCount(); i++)
		{
			final int value = r.nextInt(256);
			bp.set(i, value);
			data.fdata[i] = data.sdata[i] = (short) value;
		}
		return data;
	}

	@SeededTest
	public void canGetData(RandomSeed seed)
	{
		ImageConverterTestData data = (ImageConverterTestData) dataCache.getData(seed, this);
		byte[] bdata = data.bdata;
		short[] sdata = data.sdata;
		float[] fdata = data.fdata;
		final Rectangle bounds = null;
		final float[] fe = fdata;
		Assertions.assertArrayEquals(fe, IJImageConverter.getData(bdata, w, h, bounds, null));
		Assertions.assertArrayEquals(fe, IJImageConverter.getData(sdata, w, h, bounds, null));
		Assertions.assertArrayEquals(fe, IJImageConverter.getData(fdata, w, h, bounds, null));
		// Check the double format
		final double[] de = SimpleArrayUtils.toDouble(fe);
		Assertions.assertArrayEquals(de, IJImageConverter.getDoubleData(bdata, w, h, bounds, null));
		Assertions.assertArrayEquals(de, IJImageConverter.getDoubleData(sdata, w, h, bounds, null));
		Assertions.assertArrayEquals(de, IJImageConverter.getDoubleData(fdata, w, h, bounds, null));
	}

	@SeededTest
	public void canGetDataWithFullBounds(RandomSeed seed)
	{
		ImageConverterTestData data = (ImageConverterTestData) dataCache.getData(seed, this);
		byte[] bdata = data.bdata;
		short[] sdata = data.sdata;
		float[] fdata = data.fdata;
		final Rectangle bounds = new Rectangle(0, 0, w, h);
		final float[] fe = fdata;
		Assertions.assertArrayEquals(fe, IJImageConverter.getData(bdata, w, h, bounds, null));
		Assertions.assertArrayEquals(fe, IJImageConverter.getData(sdata, w, h, bounds, null));
		Assertions.assertArrayEquals(fe, IJImageConverter.getData(fdata, w, h, bounds, null));
		// Check the double format
		final double[] de = SimpleArrayUtils.toDouble(fe);
		Assertions.assertArrayEquals(de, IJImageConverter.getDoubleData(bdata, w, h, bounds, null));
		Assertions.assertArrayEquals(de, IJImageConverter.getDoubleData(sdata, w, h, bounds, null));
		Assertions.assertArrayEquals(de, IJImageConverter.getDoubleData(fdata, w, h, bounds, null));
	}

	@SeededTest
	public void canGetCropData(RandomSeed seed)
	{
		ImageConverterTestData data = (ImageConverterTestData) dataCache.getData(seed, this);
		byte[] bdata = data.bdata;
		short[] sdata = data.sdata;
		float[] fdata = data.fdata;
		final UniformRandomProvider rand = TestSettings.getRandomGenerator(seed.getSeed());
		final ImageExtractor ie = new ImageExtractor(fdata, w, h);
		for (int i = 0; i < 10; i++)
		{
			final Rectangle bounds = ie.getBoxRegionBounds(10 + rand.nextInt(w - 20), 10 + rand.nextInt(h - 20),
					5 + rand.nextInt(5));
			final FloatProcessor ip = new FloatProcessor(w, h, fdata.clone());
			ip.setRoi(bounds);
			final float[] fe = (float[]) (ip.crop().getPixels());
			Assertions.assertArrayEquals(fe, IJImageConverter.getData(bdata, w, h, bounds, null));
			Assertions.assertArrayEquals(fe, IJImageConverter.getData(sdata, w, h, bounds, null));
			Assertions.assertArrayEquals(fe, IJImageConverter.getData(fdata, w, h, bounds, null));
			// Check the double format
			final double[] de = SimpleArrayUtils.toDouble(fe);
			Assertions.assertArrayEquals(de, IJImageConverter.getDoubleData(bdata, w, h, bounds, null));
			Assertions.assertArrayEquals(de, IJImageConverter.getDoubleData(sdata, w, h, bounds, null));
			Assertions.assertArrayEquals(de, IJImageConverter.getDoubleData(fdata, w, h, bounds, null));
		}
	}
}
