/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

import ij.process.FloatProcessor;
import java.awt.Rectangle;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.ImageExtractor;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngFactory;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;

@SuppressWarnings({"javadoc"})
class PerPixelCameraModelTest {
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

  static final int w = 113;
  static final int h = 29;
  static final int size = w * h;

  private static class PerPixelCameraModelTestData {
    float[] bias;
    float[] gain;
    float[] variance;
    float[] varG2;
    float[] image;
  }

  private static Object createData(RandomSeed seed) {
    final UniformRandomProvider r = RngFactory.create(seed.get());

    final PerPixelCameraModelTestData data = new PerPixelCameraModelTestData();
    data.bias = new float[size];
    data.gain = new float[size];
    data.variance = new float[size];
    data.varG2 = new float[size];
    data.image = new float[size];
    for (int i = 0; i < size; i++) {
      data.bias[i] = r.nextFloat();
      data.gain[i] = 1f + r.nextFloat(); // Ensure positive
      data.variance[i] = r.nextFloat();
      data.image[i] = 100 * r.nextFloat();
      data.varG2[i] = data.variance[i] / (data.gain[i] * data.gain[i]);
    }
    return data;
  }

  @SeededTest
  void canGetDataWithFullBounds(RandomSeed seed) {
    final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache
        .computeIfAbsent(seed, PerPixelCameraModelTest::createData);
    final PerPixelCameraModel model =
        new PerPixelCameraModel(w, h, data.bias, data.gain, data.variance);
    final Rectangle bounds = new Rectangle(0, 0, w, h);
    Assertions.assertArrayEquals(data.bias, model.getBias(bounds));
    Assertions.assertArrayEquals(data.gain, model.getGain(bounds));
    Assertions.assertArrayEquals(data.variance, model.getVariance(bounds));
    Assertions.assertArrayEquals(data.varG2, model.getNormalisedVariance(bounds));
    Assertions.assertArrayEquals(data.bias, model.getBias());
    Assertions.assertArrayEquals(data.gain, model.getGain());
    Assertions.assertArrayEquals(data.variance, model.getVariance());
    Assertions.assertArrayEquals(data.varG2, model.getNormalisedVariance());
  }

  @SeededTest
  void canGetCropData(RandomSeed seed) {
    canGetCropData(seed, true);
    canGetCropData(seed, false);
  }

  private static void canGetCropData(RandomSeed seed, boolean initialise) {
    final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache
        .computeIfAbsent(seed, PerPixelCameraModelTest::createData);
    final PerPixelCameraModel model = createModel(data, initialise);
    final UniformRandomProvider rand = RngFactory.create(seed.get());
    final ImageExtractor ie = ImageExtractor.wrap(data.bias, w, h);
    for (int i = 0; i < 10; i++) {
      final Rectangle bounds = getBounds(rand, ie);
      check(data.bias, bounds, model.getBias(bounds));
      check(data.gain, bounds, model.getGain(bounds));
      check(data.variance, bounds, model.getVariance(bounds));
      check(data.varG2, bounds, model.getNormalisedVariance(bounds));
    }
  }

  private static Rectangle getBounds(UniformRandomProvider rand, ImageExtractor ie) {
    final Rectangle bounds = ie.getBoxRegionBounds(5 + rand.nextInt(w - 10),
        5 + rand.nextInt(h - 10), 2 + rand.nextInt(3));
    return bounds;
  }

  private static PerPixelCameraModel createModel(PerPixelCameraModelTestData data,
      boolean initialise) {
    final PerPixelCameraModel model =
        new PerPixelCameraModel(w, h, data.bias, data.gain, data.variance);
    if (initialise) {
      model.initialise();
    }
    return model;
  }

  private static void check(float[] data, Rectangle bounds, float[] observed) {
    final FloatProcessor ip = new FloatProcessor(w, h, data.clone());
    ip.setRoi(bounds);
    final float[] expected = (float[]) (ip.crop().getPixels());
    Assertions.assertArrayEquals(expected, observed);
  }

  @SeededTest
  void canCropAndGetData(RandomSeed seed) {
    canCropAndGetData(seed, true);
    canCropAndGetData(seed, false);
  }

  private static void canCropAndGetData(RandomSeed seed, boolean initialise) {
    final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache
        .computeIfAbsent(seed, PerPixelCameraModelTest::createData);
    final PerPixelCameraModel model = createModel(data, initialise);
    final UniformRandomProvider rand = RngFactory.create(seed.get());
    final ImageExtractor ie = ImageExtractor.wrap(data.bias, w, h);
    for (int i = 0; i < 10; i++) {
      final Rectangle bounds = getBounds(rand, ie);
      final CameraModel model2 = model.crop(bounds, false);
      Assertions.assertEquals(model2.getBounds(), bounds);
      check(data.bias, bounds, model2.getBias(bounds));
      check(data.gain, bounds, model2.getGain(bounds));
      check(data.variance, bounds, model2.getVariance(bounds));
      check(data.varG2, bounds, model2.getNormalisedVariance(bounds));
    }
  }

  @SeededTest
  void canConvertDataWithFullBounds(RandomSeed seed) {
    final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache
        .computeIfAbsent(seed, PerPixelCameraModelTest::createData);
    final PerPixelCameraModel model =
        new PerPixelCameraModel(w, h, data.bias, data.gain, data.variance);
    checkConversion(data, new Rectangle(0, 0, w, h), model);
  }

  private static void checkConversion(PerPixelCameraModelTestData data, Rectangle bounds,
      CameraModel model) {
    final FloatProcessor ip = new FloatProcessor(w, h, data.image.clone());
    ip.setRoi(bounds);
    final float[] e = (float[]) (ip.crop().getPixels());
    float[] o1 = e.clone();
    final float[] o2 = e.clone();

    ip.setPixels(data.bias);
    final float[] bias = (float[]) (ip.crop().getPixels());

    for (int i = 0; i < e.length; i++) {
      e[i] -= bias[i];
    }
    model.removeBias(bounds, o1);
    Assertions.assertArrayEquals(e, o1);

    ip.setPixels(data.gain);
    final float[] gain = (float[]) (ip.crop().getPixels());

    for (int i = 0; i < e.length; i++) {
      e[i] /= gain[i];
    }
    model.removeGain(bounds, o1);
    Assertions.assertArrayEquals(e, o1);

    o1 = o2;
    model.removeBiasAndGain(bounds, o1);
    Assertions.assertArrayEquals(e, o1);
  }

  @SeededTest
  void canConvertDataWithCropBounds(RandomSeed seed) {
    final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache
        .computeIfAbsent(seed, PerPixelCameraModelTest::createData);
    final PerPixelCameraModel model =
        new PerPixelCameraModel(w, h, data.bias, data.gain, data.variance);
    final UniformRandomProvider rand = RngFactory.create(seed.get());
    final ImageExtractor ie = ImageExtractor.wrap(data.bias, w, h);
    for (int j = 0; j < 10; j++) {
      final Rectangle bounds = getBounds(rand, ie);
      checkConversion(data, bounds, model);
    }
  }

  @SeededTest
  void canCropAndConvertDataWithCropBounds(RandomSeed seed) {
    final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache
        .computeIfAbsent(seed, PerPixelCameraModelTest::createData);
    final PerPixelCameraModel model =
        new PerPixelCameraModel(w, h, data.bias, data.gain, data.variance);
    final UniformRandomProvider rand = RngFactory.create(seed.get());
    final ImageExtractor ie = ImageExtractor.wrap(data.bias, w, h);
    for (int j = 0; j < 10; j++) {
      final Rectangle bounds = getBounds(rand, ie);
      checkConversion(data, bounds, model.crop(bounds, false));
    }
  }

  @SeededTest
  void canGetWeightsWithPositiveVariance(RandomSeed seed) {
    final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache
        .computeIfAbsent(seed, PerPixelCameraModelTest::createData);
    final float[] var = data.variance.clone();
    for (int i = 0; i < var.length; i++) {
      if (var[i] == 0) {
        var[i] = 1;
      }
    }
    final PerPixelCameraModel model = new PerPixelCameraModel(w, h, data.bias, data.gain, var);
    final float[] w = model.getWeights(model.getBounds());
    final float[] e = var;
    for (int i = 0; i < e.length; i++) {
      e[i] = (float) (1.0 / e[i]);
    }
    Assertions.assertArrayEquals(e, w);
  }

  @SeededTest
  void canGetWeightsWithAllZeroVariance(RandomSeed seed) {
    final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache
        .computeIfAbsent(seed, PerPixelCameraModelTest::createData);
    final float[] var = new float[data.variance.length];
    final PerPixelCameraModel model = new PerPixelCameraModel(w, h, data.bias, data.gain, var);
    final float[] w = model.getWeights(model.getBounds());
    final float[] e = var;
    Arrays.fill(e, 1f);
    Assertions.assertArrayEquals(e, w);
  }

  @SeededTest
  void canGetWeightsWithZeroVariance(RandomSeed seed) {
    final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache
        .computeIfAbsent(seed, PerPixelCameraModelTest::createData);
    final float[] var = data.variance.clone();
    var[0] = 0;
    final float min = SimpleArrayUtils.minAboveZero(var);
    final PerPixelCameraModel model = new PerPixelCameraModel(w, h, data.bias, data.gain, var);
    final float[] w = model.getWeights(model.getBounds());
    final float[] e = var;
    for (int i = 0; i < e.length; i++) {
      e[i] = (e[i] == 0) ? (float) (1.0 / min) : (float) (1.0 / e[i]);
    }
    Assertions.assertArrayEquals(e, w);
  }

  @SeededTest
  void canGetMeanVariance(RandomSeed seed) {
    canGetMeanVarianceData(seed, true, false);
    canGetMeanVarianceData(seed, false, false);
  }

  @SeededTest
  void canGetMeanNormalisedVariance(RandomSeed seed) {
    canGetMeanVarianceData(seed, true, true);
    canGetMeanVarianceData(seed, false, true);
  }

  private static void canGetMeanVarianceData(RandomSeed seed, boolean initialise,
      boolean normalised) {
    final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache
        .computeIfAbsent(seed, PerPixelCameraModelTest::createData);
    final PerPixelCameraModel model = createModel(data, initialise);
    final UniformRandomProvider rand = RngFactory.create(seed.get());
    final ImageExtractor ie = ImageExtractor.wrap(data.bias, w, h);
    for (int i = 0; i < 10; i++) {
      final Rectangle bounds = getBounds(rand, ie);
      final float[] v =
          (normalised) ? model.getNormalisedVariance(bounds) : model.getVariance(bounds);
      final double e = MathUtils.sum(v) / v.length;
      final double o =
          (normalised) ? model.getMeanNormalisedVariance(bounds) : model.getMeanVariance(bounds);
      Assertions.assertEquals(e, o);
    }
  }

  @SeededTest
  void canGetCoordinateData(RandomSeed seed) {
    final PerPixelCameraModelTestData data = (PerPixelCameraModelTestData) dataCache
        .computeIfAbsent(seed, PerPixelCameraModelTest::createData);
    final int ox = 2;
    final int oy = 3;
    final int w = 8;
    final int h = 10;
    final int size = w * h;
    final float[] bias = Arrays.copyOf(data.bias, size);
    final float[] gain = Arrays.copyOf(data.gain, size);
    final float[] variance = Arrays.copyOf(data.variance, size);
    final PerPixelCameraModel model = new PerPixelCameraModel(ox, oy, w, h, bias, gain, variance);
    for (int y = 0, y1 = oy, i = 0; y < h; y++, y1++) {
      for (int x = 0, x1 = ox; x < w; x++, x1++, i++) {
        Assertions.assertEquals(bias[i], model.getBias(x1, y1));
        Assertions.assertEquals(gain[i], model.getGain(x1, y1));
        Assertions.assertEquals(variance[i], model.getVariance(x1, y1));
        Assertions.assertEquals(data.varG2[i], model.getNormalisedVariance(x1, y1));
      }
    }
  }

  @Test
  void testLocate() {
    final int w = 10;
    final int h = 20;
    final float[] b = createFixedData(w, h, 100);
    final float[] g = createFixedData(w, h, 3);
    final float[] v = createFixedData(w, h, 7);
    final PerPixelCameraModel model =
        PerPixelCameraModel.create(new Rectangle(0, 0, w, h), b, g, v);

    Assertions.assertTrue(model.locate(new NullCameraModel()).isEmpty());

    final PerPixelCameraModel model2 = PerPixelCameraModel.create(new Rectangle(0, 0, w - 1, h - 1),
        createFixedData(w - 1, h - 1, 123), createFixedData(w - 1, h - 1, 4),
        createFixedData(w - 1, h - 1, 8));
    Assertions.assertTrue(model.locate(model2).isEmpty());

    float[] b2 = createFixedData(w + 1, h, b[0]);
    float[] g2 = createFixedData(w + 1, h, g[0]);
    float[] v2 = createFixedData(w + 1, h, v[0]);
    Assertions.assertTrue(model
        .locate(PerPixelCameraModel.create(new Rectangle(0, 0, w + 1, h), b2, g2, v2)).isEmpty());
    b2 = createFixedData(w, h + 1, b[0]);
    g2 = createFixedData(w, h + 1, g[0]);
    v2 = createFixedData(w, h + 1, v[0]);
    Assertions.assertTrue(model
        .locate(PerPixelCameraModel.create(new Rectangle(0, 0, w, h + 1), b2, g2, v2)).isEmpty());

    assertLocate(model, 1, 1, w - 1, h - 1, 0, 0, 1, 0, 0, 1, 1, 1);

    // Randomise
    final int n = b.length;
    final UniformRandomProvider rng = RngFactory.createWithFixedSeed();
    for (final float[] data : new float[][] {b, g, v}) {
      final float[] restore = data.clone();
      for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 100; j++) {
          data[rng.nextInt(n)] = rng.nextFloat();
        }
        assertLocate(model, 2, 3, 6, 13);
        assertLocate(model, 0, 0, w, h);
        System.arraycopy(restore, 0, data, 0, n);
      }
    }

    // Duplicate pixel
    final int x = 3;
    final int y = 7;
    final int ox = 1;
    final int oy = 10;
    final int i1 = y * w + x;
    final int i2 = (y + oy) * w + x + ox;
    b[i1] = ++b[i2];
    g[i1] = ++g[i2];
    v[i1] = ++v[i2];
    // Large patch matches correctly
    assertLocate(model, 2, 3, 7, 15);
    // Small second patch matches the first, including offsets
    assertLocate(model, 2 + ox, 3 + oy, 7, 5, 2, 3, 2 + ox, 3 + oy);
    assertLocate(model, 2 + ox - 1, 3 + oy, 7, 5, 2 - 1, 3, 2 + ox - 1, 3 + oy);
    assertLocate(model, 2 + ox, 3 + oy - 1, 7, 6, 2, 3 - 1, 2 + ox, 3 + oy - 1);
    assertLocate(model, 2 + ox - 1, 3 + oy - 1, 7, 6, 2 - 1, 3 - 1, 2 + ox - 1, 3 + oy - 1);
  }

  private static float[] createFixedData(int w, int h, float value) {
    final float[] d = new float[w * h];
    Arrays.fill(d, value);
    return d;
  }

  private static void assertLocate(PerPixelCameraModel model, int x, int y, int w, int h) {
    assertLocate(model, x, y, w, h, x, y);
  }

  private static void assertLocate(PerPixelCameraModel model, int x, int y, int w, int h,
      int... origin) {
    final Rectangle bounds = new Rectangle(x, y, w, h);
    final CameraModel crop = model.crop(bounds, false);
    final List<Rectangle> rectangles = model.locate(crop);
    Assertions.assertEquals(origin.length / 2, rectangles.size(), "Number of rectangles");
    for (int i = 0; i < origin.length; i += 2) {
      final Rectangle expectedBounds = new Rectangle(origin[i], origin[i + 1], w, h);
      Assertions.assertEquals(expectedBounds, rectangles.get(i / 2));
    }
  }
}
