/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

import uk.ac.sussex.gdsc.core.utils.ImageExtractor;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;

import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;

import java.awt.Rectangle;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Function;

@SuppressWarnings({"javadoc"})
public class ImageConverterTest {
  private static ConcurrentHashMap<RandomSeed, Object> dataCache;

  @BeforeAll
  public static void beforeAll() {
    dataCache = new ConcurrentHashMap<>();
  }

  @AfterAll
  public static void afterAll() {
    dataCache.clear();
    dataCache = null;
  }

  static final int w = 200;
  static final int h = 300;

  private static class ImageConverterTestData {
    byte[] bdata;
    short[] sdata;
    float[] fdata;
  }

  private static Object createData(RandomSeed seed) {
    final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
    final ByteProcessor bp = new ByteProcessor(w, h);
    final ImageConverterTestData data = new ImageConverterTestData();
    data.bdata = (byte[]) bp.getPixels();
    data.sdata = new short[data.bdata.length];
    data.fdata = new float[data.bdata.length];
    for (int i = 0; i < bp.getPixelCount(); i++) {
      final int value = r.nextInt(256);
      bp.set(i, value);
      data.fdata[i] = data.sdata[i] = (short) value;
    }
    return data;
  }

  @SeededTest
  public void canGetData(RandomSeed seed) {
    final ImageConverterTestData data = (ImageConverterTestData) dataCache
        .computeIfAbsent(seed, ImageConverterTest::createData);
    final byte[] bdata = data.bdata;
    final short[] sdata = data.sdata;
    final float[] fdata = data.fdata;
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
  public void canGetDataWithFullBounds(RandomSeed seed) {
    final ImageConverterTestData data = (ImageConverterTestData) dataCache
        .computeIfAbsent(seed, ImageConverterTest::createData);
    final byte[] bdata = data.bdata;
    final short[] sdata = data.sdata;
    final float[] fdata = data.fdata;
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
  public void canGetCropData(RandomSeed seed) {
    final ImageConverterTestData data = (ImageConverterTestData) dataCache
        .computeIfAbsent(seed, ImageConverterTest::createData);
    final byte[] bdata = data.bdata;
    final short[] sdata = data.sdata;
    final float[] fdata = data.fdata;
    final UniformRandomProvider rand = RngUtils.create(seed.getSeedAsLong());
    final ImageExtractor ie = ImageExtractor.wrap(fdata, w, h);
    for (int i = 0; i < 10; i++) {
      final Rectangle bounds = ie.getBoxRegionBounds(10 + rand.nextInt(w - 20),
          10 + rand.nextInt(h - 20), 5 + rand.nextInt(5));
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
