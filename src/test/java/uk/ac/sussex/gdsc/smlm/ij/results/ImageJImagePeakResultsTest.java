/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.ij.results;

import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.awt.Rectangle;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.functions.IntArrayFormatSupplier;

/**
 * Test the IJImagePeakResults functionality.
 */
@SuppressWarnings({"javadoc"})
class ImageJImagePeakResultsTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(ImageJImagePeakResultsTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  static Calibration calibration;
  static PSF psf;

  static {
    psf = PsfHelper.create(PSFType.ONE_AXIS_GAUSSIAN_2D);
    final CalibrationWriter cw = new CalibrationWriter();
    cw.setDistanceUnit(DistanceUnit.PIXEL);
    cw.setIntensityUnit(IntensityUnit.COUNT);
    calibration = cw.getCalibration();
  }

  private static final String title = "Test";
  Rectangle bounds = new Rectangle(0, 0, 3, 5);

  @Test
  void canAddToSinglePixels() {
    final ImageJImagePeakResults r = new ImageJImagePeakResults(title, bounds, 1);
    final FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
    begin(r);
    add(fp, r, 1, 1, 1);
    add(fp, r, 1, 2, 4);
    add(fp, r, 0, 1, 2);
    r.end();
    final float[] expecteds = getImage(fp);
    final float[] actuals = getImage(r);
    Assertions.assertArrayEquals(expecteds, actuals);
  }

  @Test
  void canAddToSinglePixelsWithInvalidPositions() {
    final ImageJImagePeakResults r = new ImageJImagePeakResults(title, bounds, 1);
    final FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
    begin(r);
    add(fp, r, 1, 1, 1);
    add(fp, r, 1, 2, 4);
    add(fp, r, 0, 1, 2);
    for (final int x : new int[] {-1, 0, 1, bounds.width, bounds.width + 1}) {
      for (final int y : new int[] {-1, 0, 1, bounds.height, bounds.height + 1}) {
        add(fp, r, x, y, 1);
      }
    }
    r.end();
    final float[] expecteds = getImage(fp);
    final float[] actuals = getImage(r);
    Assertions.assertArrayEquals(expecteds, actuals);
  }

  @Test
  void canAddToZeroSizeImage() {
    final Rectangle zeroBounds = new Rectangle();
    canAddToZeroSizeImage(zeroBounds);
    zeroBounds.width = 1;
    canAddToZeroSizeImage(zeroBounds);
    zeroBounds.width = 0;
    zeroBounds.height = 1;
    canAddToZeroSizeImage(zeroBounds);
    zeroBounds.width = -1;
    zeroBounds.height = -1;
    canAddToZeroSizeImage(zeroBounds);
  }

  private static void canAddToZeroSizeImage(Rectangle bounds) {
    final ImageJImagePeakResults r = new ImageJImagePeakResults(title, bounds, 1);
    begin(r);
    for (final int x : new int[] {-1, 0, 1, bounds.width, bounds.width + 1}) {
      for (final int y : new int[] {-1, 0, 1, bounds.height, bounds.height + 1}) {
        r.add(x, y, 1);
      }
    }
    r.end();
    final float[] expecteds = new float[1];
    final float[] actuals = getImage(r);
    Assertions.assertArrayEquals(expecteds, actuals);
  }

  private static void add(FloatProcessor fp, ImageJImagePeakResults results, int x, int y,
      float value) {
    addValue(results, x, y, value);
    fp.putPixelValue(x, y, fp.getPixelValue(x, y) + value);
  }

  @Test
  void canInterpolateInMiddleOfPixel() {
    final ImageJImagePeakResults r = new ImageJImagePeakResults(title, bounds, 1);
    r.setDisplayFlags(ImageJImagePeakResults.DISPLAY_WEIGHTED);
    final FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
    begin(r);
    addValue(r, 1.5f, 1.5f, 1);
    fp.putPixelValue(1, 1, 1);
    r.end();
    final float[] expecteds = getImage(fp);
    final float[] actuals = getImage(r);
    Assertions.assertArrayEquals(expecteds, actuals);
  }

  @Test
  void canInterpolateDownInXAtPixelEdge() {
    final ImageJImagePeakResults r = new ImageJImagePeakResults(title, bounds, 1);
    r.setDisplayFlags(ImageJImagePeakResults.DISPLAY_WEIGHTED);
    final FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
    begin(r);
    addValue(r, 1f, 1.5f, 2);
    fp.putPixelValue(0, 1, 1);
    fp.putPixelValue(1, 1, 1);
    r.end();
    final float[] expecteds = getImage(fp);
    final float[] actuals = getImage(r);
    Assertions.assertArrayEquals(expecteds, actuals);
  }

  @Test
  void canInterpolateUpInXAtPixelEdge() {
    final ImageJImagePeakResults r = new ImageJImagePeakResults(title, bounds, 1);
    r.setDisplayFlags(ImageJImagePeakResults.DISPLAY_WEIGHTED);
    final FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
    begin(r);
    addValue(r, 2f, 1.5f, 2);
    fp.putPixelValue(1, 1, 1);
    fp.putPixelValue(2, 1, 1);
    r.end();
    final float[] expecteds = getImage(fp);
    final float[] actuals = getImage(r);
    Assertions.assertArrayEquals(expecteds, actuals);
  }

  @Test
  void canInterpolateDownInYAtPixelEdge() {
    final ImageJImagePeakResults r = new ImageJImagePeakResults(title, bounds, 1);
    r.setDisplayFlags(ImageJImagePeakResults.DISPLAY_WEIGHTED);
    final FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
    begin(r);
    addValue(r, 1.5f, 1f, 2);
    fp.putPixelValue(1, 0, 1);
    fp.putPixelValue(1, 1, 1);
    r.end();
    final float[] expecteds = getImage(fp);
    final float[] actuals = getImage(r);
    Assertions.assertArrayEquals(expecteds, actuals);
  }

  @Test
  void canInterpolateUpInYAtPixelEdge() {
    final ImageJImagePeakResults r = new ImageJImagePeakResults(title, bounds, 1);
    r.setDisplayFlags(ImageJImagePeakResults.DISPLAY_WEIGHTED);
    final FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
    begin(r);
    addValue(r, 1.5f, 2f, 2);
    fp.putPixelValue(1, 1, 1);
    fp.putPixelValue(1, 2, 1);
    r.end();
    final float[] expecteds = getImage(fp);
    final float[] actuals = getImage(r);
    Assertions.assertArrayEquals(expecteds, actuals);
  }

  @Test
  void canInterpolateDownInX() {
    final ImageJImagePeakResults r = new ImageJImagePeakResults(title, bounds, 1);
    r.setDisplayFlags(ImageJImagePeakResults.DISPLAY_WEIGHTED);
    final FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
    begin(r);
    addValue(r, 1.25f, 1.5f, 2);
    fp.putPixelValue(0, 1, 0.5f);
    fp.putPixelValue(1, 1, 1.5f);
    r.end();
    final float[] expecteds = getImage(fp);
    final float[] actuals = getImage(r);
    Assertions.assertArrayEquals(expecteds, actuals);
  }

  @Test
  void canInterpolateUpInX() {
    final ImageJImagePeakResults r = new ImageJImagePeakResults(title, bounds, 1);
    r.setDisplayFlags(ImageJImagePeakResults.DISPLAY_WEIGHTED);
    final FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
    begin(r);
    addValue(r, 1.75f, 1.5f, 2);
    fp.putPixelValue(1, 1, 1.5f);
    fp.putPixelValue(2, 1, 0.5f);
    r.end();
    final float[] expecteds = getImage(fp);
    final float[] actuals = getImage(r);
    Assertions.assertArrayEquals(expecteds, actuals);
  }

  @Test
  void canInterpolateDownInY() {
    final ImageJImagePeakResults r = new ImageJImagePeakResults(title, bounds, 1);
    r.setDisplayFlags(ImageJImagePeakResults.DISPLAY_WEIGHTED);
    final FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
    begin(r);
    addValue(r, 1.5f, 1.25f, 2);
    fp.putPixelValue(1, 0, 0.5f);
    fp.putPixelValue(1, 1, 1.5f);
    r.end();
    final float[] expecteds = getImage(fp);
    final float[] actuals = getImage(r);
    Assertions.assertArrayEquals(expecteds, actuals);
  }

  @Test
  void canInterpolateUpInY() {
    final ImageJImagePeakResults r = new ImageJImagePeakResults(title, bounds, 1);
    r.setDisplayFlags(ImageJImagePeakResults.DISPLAY_WEIGHTED);
    final FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
    begin(r);
    addValue(r, 1.5f, 1.75f, 2);
    fp.putPixelValue(1, 1, 1.5f);
    fp.putPixelValue(1, 2, 0.5f);
    r.end();
    final float[] expecteds = getImage(fp);
    final float[] actuals = getImage(r);
    Assertions.assertArrayEquals(expecteds, actuals);
  }

  @Test
  void canInterpolateDownInXyAtPixelEdge() {
    final ImageJImagePeakResults r = new ImageJImagePeakResults(title, bounds, 1);
    r.setDisplayFlags(ImageJImagePeakResults.DISPLAY_WEIGHTED);
    final FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
    begin(r);
    addValue(r, 1f, 1f, 4);
    fp.putPixelValue(0, 0, 1f);
    fp.putPixelValue(0, 1, 1f);
    fp.putPixelValue(1, 0, 1f);
    fp.putPixelValue(1, 1, 1f);
    r.end();
    final float[] expecteds = getImage(fp);
    final float[] actuals = getImage(r);
    Assertions.assertArrayEquals(expecteds, actuals);
  }

  @Test
  void canInterpolateUpInXyAtPixelEdge() {
    final ImageJImagePeakResults r = new ImageJImagePeakResults(title, bounds, 1);
    r.setDisplayFlags(ImageJImagePeakResults.DISPLAY_WEIGHTED);
    final FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
    begin(r);
    addValue(r, 2f, 2f, 4);
    fp.putPixelValue(1, 1, 1f);
    fp.putPixelValue(2, 1, 1f);
    fp.putPixelValue(1, 2, 1f);
    fp.putPixelValue(2, 2, 1f);
    r.end();
    final float[] expecteds = getImage(fp);
    final float[] actuals = getImage(r);
    Assertions.assertArrayEquals(expecteds, actuals);
  }

  @Test
  void canInterpolateDownInXy() {
    final ImageJImagePeakResults r = new ImageJImagePeakResults(title, bounds, 1);
    r.setDisplayFlags(ImageJImagePeakResults.DISPLAY_WEIGHTED);
    final FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
    begin(r);
    addValue(r, 1.25f, 1.25f, 1);
    fp.putPixelValue(0, 0, 0.25f * 0.25f);
    fp.putPixelValue(0, 1, 0.75f * 0.25f);
    fp.putPixelValue(1, 0, 0.75f * 0.25f);
    fp.putPixelValue(1, 1, 0.75f * 0.75f);
    r.end();
    final float[] expecteds = getImage(fp);
    final float[] actuals = getImage(r);
    Assertions.assertArrayEquals(expecteds, actuals);
  }

  @Test
  void canInterpolateUpInXy() {
    final ImageJImagePeakResults r = new ImageJImagePeakResults(title, bounds, 1);
    r.setDisplayFlags(ImageJImagePeakResults.DISPLAY_WEIGHTED);
    final FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
    begin(r);
    addValue(r, 1.75f, 1.75f, 1);
    fp.putPixelValue(1, 1, 0.75f * 0.75f);
    fp.putPixelValue(2, 1, 0.75f * 0.25f);
    fp.putPixelValue(1, 2, 0.75f * 0.25f);
    fp.putPixelValue(2, 2, 0.25f * 0.25f);
    r.end();
    final float[] expecteds = getImage(fp);
    final float[] actuals = getImage(r);
    Assertions.assertArrayEquals(expecteds, actuals);
  }

  @Test
  void noInterpolateDownInXAtImageEdge() {
    final ImageJImagePeakResults r = new ImageJImagePeakResults(title, bounds, 1);
    r.setDisplayFlags(ImageJImagePeakResults.DISPLAY_WEIGHTED);
    final FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
    begin(r);
    addValue(r, 0.5f, 1.5f, 2);
    fp.putPixelValue(0, 1, 2);
    r.end();
    final float[] expecteds = getImage(fp);
    final float[] actuals = getImage(r);
    Assertions.assertArrayEquals(expecteds, actuals);
  }

  @Test
  void noInterpolateUpInXAtImageEdge() {
    final ImageJImagePeakResults r = new ImageJImagePeakResults(title, bounds, 1);
    r.setDisplayFlags(ImageJImagePeakResults.DISPLAY_WEIGHTED);
    final FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
    begin(r);
    addValue(r, 2.5f, 1.5f, 2);
    fp.putPixelValue(2, 1, 2);
    r.end();
    final float[] expecteds = getImage(fp);
    final float[] actuals = getImage(r);
    Assertions.assertArrayEquals(expecteds, actuals);
  }

  @Test
  void noInterpolateDownInYAtImageEdge() {
    final ImageJImagePeakResults r = new ImageJImagePeakResults(title, bounds, 1);
    r.setDisplayFlags(ImageJImagePeakResults.DISPLAY_WEIGHTED);
    final FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
    begin(r);
    addValue(r, 1.5f, 0.5f, 2);
    fp.putPixelValue(1, 0, 2);
    r.end();
    final float[] expecteds = getImage(fp);
    final float[] actuals = getImage(r);
    Assertions.assertArrayEquals(expecteds, actuals);
  }

  @Test
  void noInterpolateUpInYAtImageEdge() {
    final ImageJImagePeakResults r = new ImageJImagePeakResults(title, bounds, 1);
    r.setDisplayFlags(ImageJImagePeakResults.DISPLAY_WEIGHTED);
    final FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
    begin(r);
    addValue(r, 1.5f, 2.5f, 2);
    fp.putPixelValue(1, 2, 2);
    r.end();
    final float[] expecteds = getImage(fp);
    final float[] actuals = getImage(r);
    Assertions.assertArrayEquals(expecteds, actuals);
  }

  @SeededTest
  void canAddUsingDifferentMethods(RandomSeed seed) {
    checkCanAddUsingDifferentMethods(seed, 0);
  }

  @SeededTest
  void canAddUsingDifferentMethodsEqualized(RandomSeed seed) {
    checkCanAddUsingDifferentMethods(seed, ImageJImagePeakResults.DISPLAY_EQUALIZED);
  }

  @SeededTest
  void canAddUsingDifferentMethodsEqualizedWeighted(RandomSeed seed) {
    checkCanAddUsingDifferentMethods(seed,
        ImageJImagePeakResults.DISPLAY_EQUALIZED | ImageJImagePeakResults.DISPLAY_WEIGHTED);
  }

  @SeededTest
  void canAddUsingDifferentMethodsMax(RandomSeed seed) {
    checkCanAddUsingDifferentMethods(seed, ImageJImagePeakResults.DISPLAY_MAX);
  }

  @SeededTest
  void canAddUsingDifferentMethodsReplace(RandomSeed seed) {
    checkCanAddUsingDifferentMethods(seed, ImageJImagePeakResults.DISPLAY_REPLACE);
  }

  @SeededTest
  void canAddUsingDifferentMethodsWeighted(RandomSeed seed) {
    checkCanAddUsingDifferentMethods(seed, ImageJImagePeakResults.DISPLAY_WEIGHTED);
  }

  private void checkCanAddUsingDifferentMethods(RandomSeed seed, int displayFlags) {
    final UniformRandomProvider rand = RngUtils.create(seed.getSeed());
    displayFlags |= ImageJImagePeakResults.DISPLAY_SIGNAL;

    final ImageJImagePeakResults[] r = new ImageJImagePeakResults[8];
    final PSF psf = PsfHelper.create(PSFType.ONE_AXIS_GAUSSIAN_2D);
    for (int i = 0; i < r.length; i++) {
      r[i] = new ImageJImagePeakResults(title + i, bounds, 1);
      r[i].setDisplayFlags(displayFlags);
      r[i].setPsf(psf);
      begin(r[i]);
    }

    final int size = 20;
    final int[] t = new int[size];
    final float[] x = new float[size];
    final float[] y = new float[size];
    final float[] v = new float[size];
    for (int i = 0; i < size; i++) {
      t[i] = i;
      x[i] = (rand.nextFloat() * bounds.width);
      y[i] = (rand.nextFloat() * bounds.height);
      v[i] = (rand.nextFloat());

      addPeakResult(r[0], x[i], y[i], v[i]);
      addPeakResult(r[1], t[i], x[i], y[i], v[i]);
      addValue(r[2], x[i], y[i], v[i]);
      addValue(r[3], t[i], x[i], y[i], v[i]);
    }

    addPeakResults(r[4], x, y, v);
    addPeakResults(r[5], t, x, y, v);
    addValues(r[6], x, y, v);
    addValues(r[7], t, x, y, v);

    final float[][] image = new float[r.length][];
    for (int i = 0; i < r.length; i++) {
      r[i].end();
      image[i] = getImage(r[i]);
      logger.log(TestLogUtils.getRecord(Level.FINE, "[%d] = %s", i, Arrays.toString(image[i])));
    }

    // Test single value adds
    float[] expecteds = image[0];
    IntArrayFormatSupplier msg = new IntArrayFormatSupplier("Single add image %d", 1);
    for (int i = 1; i < 4; i++) {
      final float[] actuals = image[i];
      Assertions.assertArrayEquals(expecteds, actuals, msg.set(0, i));
    }

    // Test multi value adds
    expecteds = image[4];
    msg = new IntArrayFormatSupplier("Multi add image %d", 1);
    for (int i = 5; i < image.length; i++) {
      final float[] actuals = image[i];
      Assertions.assertArrayEquals(expecteds, actuals, msg.set(0, i));
    }

    // Test they are roughly the same (differences occur due to floating point summation
    Assertions.assertArrayEquals(expecteds, image[0], 1e-5f, "Single != Multi");
  }

  private static void begin(ImageJImagePeakResults results) {
    results.setPsf(psf);
    results.setDisplayImage(false);
    results.setCalibration(calibration);
    results.begin();
  }

  private static void addPeakResult(ImageJImagePeakResults results, float x, float y, float value) {
    results.add(new PeakResult(x, y, value));
  }

  private static void addPeakResult(ImageJImagePeakResults results, int time, float x, float y,
      float value) {
    results.add(new PeakResult(time, 0, 0, 0, 0, 0, 0, createParams(x, y, value), null));
  }

  private static float[] createParams(float x, float y, float value) {
    return Gaussian2DPeakResultHelper.createOneAxisParams(0, value, x, y, 0, 1);
  }

  private static void addPeakResults(ImageJImagePeakResults results, float[] x, float[] y,
      float[] value) {
    final LocalList<PeakResult> list = new LocalList<>(x.length);
    for (int i = 0; i < x.length; i++) {
      list.add(new PeakResult(x[i], y[i], value[i]));
    }
    results.addAll(list);
  }

  private static void addPeakResults(ImageJImagePeakResults results, int[] time, float[] x,
      float[] y, float[] value) {
    final LocalList<PeakResult> list = new LocalList<>(x.length);
    for (int i = 0; i < x.length; i++) {
      list.add(new PeakResult(time[i], 0, 0, 0, 0, 0, 0, createParams(x[i], y[i], value[i]), null));
    }
    results.addAll(list);
  }

  private static void addValue(ImageJImagePeakResults results, float x, float y, float value) {
    results.add(x, y, value);
  }

  private static void addValue(ImageJImagePeakResults results, int time, float x, float y,
      float value) {
    results.add(time, x, y, value);
  }

  private static void addValues(ImageJImagePeakResults results, float[] x, float[] y,
      float[] value) {
    results.add(x, y, value);
  }

  private static void addValues(ImageJImagePeakResults results, int[] time, float[] x, float[] y,
      float[] value) {
    results.add(time, x, y, value);
  }

  private static float[] getImage(ImageJImagePeakResults results) {
    return getImage(results.getImagePlus().getProcessor());
  }

  private static float[] getImage(ImageProcessor ip) {
    return (float[]) ip.convertToFloat().getPixels();
  }
}
