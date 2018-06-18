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
package gdsc.smlm.ij.utils;

import java.awt.Rectangle;

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.ImageExtractor;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.test.TestSettings;
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
		RandomGenerator r = TestSettings.getRandomGenerator();
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
		Assert.assertArrayEquals(fe, IJImageConverter.getData(bdata, w, h, bounds, null), 0);
		Assert.assertArrayEquals(fe, IJImageConverter.getData(sdata, w, h, bounds, null), 0);
		Assert.assertArrayEquals(fe, IJImageConverter.getData(fdata, w, h, bounds, null), 0);
		// Check the double format
		double[] de = SimpleArrayUtils.toDouble(fe);
		Assert.assertArrayEquals(de, IJImageConverter.getDoubleData(bdata, w, h, bounds, null), 0);
		Assert.assertArrayEquals(de, IJImageConverter.getDoubleData(sdata, w, h, bounds, null), 0);
		Assert.assertArrayEquals(de, IJImageConverter.getDoubleData(fdata, w, h, bounds, null), 0);
	}

	@Test
	public void canGetDataWithFullBounds()
	{
		Rectangle bounds = new Rectangle(0, 0, w, h);
		float[] fe = fdata;
		Assert.assertArrayEquals(fe, IJImageConverter.getData(bdata, w, h, bounds, null), 0);
		Assert.assertArrayEquals(fe, IJImageConverter.getData(sdata, w, h, bounds, null), 0);
		Assert.assertArrayEquals(fe, IJImageConverter.getData(fdata, w, h, bounds, null), 0);
		// Check the double format
		double[] de = SimpleArrayUtils.toDouble(fe);
		Assert.assertArrayEquals(de, IJImageConverter.getDoubleData(bdata, w, h, bounds, null), 0);
		Assert.assertArrayEquals(de, IJImageConverter.getDoubleData(sdata, w, h, bounds, null), 0);
		Assert.assertArrayEquals(de, IJImageConverter.getDoubleData(fdata, w, h, bounds, null), 0);
	}

	@Test
	public void canGetCropData()
	{
		RandomGenerator rand = TestSettings.getRandomGenerator();
		ImageExtractor ie = new ImageExtractor(fdata, w, h);
		for (int i = 0; i < 10; i++)
		{
			Rectangle bounds = ie.getBoxRegionBounds(10 + rand.nextInt(w - 20), 10 + rand.nextInt(h - 20),
					5 + rand.nextInt(5));
			FloatProcessor ip = new FloatProcessor(w, h, fdata.clone());
			ip.setRoi(bounds);
			float[] fe = (float[]) (ip.crop().getPixels());
			Assert.assertArrayEquals(fe, IJImageConverter.getData(bdata, w, h, bounds, null), 0);
			Assert.assertArrayEquals(fe, IJImageConverter.getData(sdata, w, h, bounds, null), 0);
			Assert.assertArrayEquals(fe, IJImageConverter.getData(fdata, w, h, bounds, null), 0);
			// Check the double format
			double[] de = SimpleArrayUtils.toDouble(fe);
			Assert.assertArrayEquals(de, IJImageConverter.getDoubleData(bdata, w, h, bounds, null), 0);
			Assert.assertArrayEquals(de, IJImageConverter.getDoubleData(sdata, w, h, bounds, null), 0);
			Assert.assertArrayEquals(de, IJImageConverter.getDoubleData(fdata, w, h, bounds, null), 0);
		}
	}
}
