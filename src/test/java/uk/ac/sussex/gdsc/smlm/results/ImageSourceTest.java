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

package uk.ac.sussex.gdsc.smlm.results;

import uk.ac.sussex.gdsc.core.utils.rng.RandomUtils;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

@SuppressWarnings({"javadoc", "unused"})
public class ImageSourceTest {
  @Test
  public void canConstructMemoryImageSource() {
    final int width = 5;
    final int height = 3;
    final int n = 15;
    final MemoryImageSource source =
        new MemoryImageSource(width, height, createData(width, height, n));
    Assertions.assertEquals(width, source.getWidth());
    Assertions.assertEquals(height, source.getHeight());
    Assertions.assertEquals(n, source.getFrames());
  }

  @Test
  public void nullDataThrowsConstructMemoryImageSource() {
    final int width = 5;
    final int height = 3;
    final float[][] data = null;
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      new MemoryImageSource(width, height, data);
    });
  }

  @Test
  public void nullArrayDataThrowsConstructMemoryImageSource() {
    final int width = 5;
    final int height = 3;
    final float[][] data = createData(width, height, 15);
    data[2] = null;
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      new MemoryImageSource(width, height, data);
    });
  }

  @Test
  public void invalidLengthArrayDataThrowsConstructMemoryImageSource() {
    final int width = 5;
    final int height = 3;
    final float[][] data = createData(width, height, 15);
    data[2] = new float[width * height + 1];
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      new MemoryImageSource(width, height, data);
    });
  }

  @Test
  public void invalidWidthThrowsConstructMemoryImageSource() {
    final int width = 5;
    final int height = 3;
    final float[][] data = createData(width, height, 15);
    int ok = 0;
    new MemoryImageSource(width, height, data);
    if (canConstruct(width + 1, height, data)) {
      ok++;
    }
    if (canConstruct(width - 1, height, data)) {
      ok++;
    }
    if (canConstruct(-1, height, data)) {
      ok++;
    }
    if (canConstruct(0, height, data)) {
      ok++;
    }
    Assertions.assertEquals(0, ok);
  }

  @Test
  public void invalidHeightThrowsConstructMemoryImageSource() {
    final int width = 5;
    final int height = 3;
    final float[][] data = createData(width, height, 15);
    int ok = 0;
    new MemoryImageSource(width, height, data);
    if (canConstruct(width, height + 1, data)) {
      ok++;
    }
    if (canConstruct(width, height - 1, data)) {
      ok++;
    }
    if (canConstruct(width, -1, data)) {
      ok++;
    }
    if (canConstruct(width, 0, data)) {
      ok++;
    }
    Assertions.assertEquals(0, ok);
  }

  private static boolean canConstruct(int width, int height, float[][] data) {
    try {
      new MemoryImageSource(width, height, data);
      return true;
    } catch (final RuntimeException ex) {
      return false;
    }
  }

  @Test
  public void memoryImageSourceCanReturnDataWithNext() {
    final int width = 5;
    final int height = 3;
    final float[][] data = createData(width, height, 15);
    final MemoryImageSource source = new MemoryImageSource(width, height, data);
    int index = 0;
    float[] next = null;
    Assertions.assertTrue(source.open());
    while ((next = source.next()) != null) {
      Assertions.assertArrayEquals(data[index++], next);
    }
    Assertions.assertEquals(index, data.length);
  }

  @SeededTest
  public void memoryImageSourceCanReturnDataWithGet(RandomSeed seed) {
    final int width = 5;
    final int height = 3;
    final float[][] data = createData(width, height, 15);
    final MemoryImageSource source = new MemoryImageSource(width, height, data);

    final int[] frames = new int[data.length];
    for (int i = 0; i < data.length; i++) {
      frames[i] = i + 1;
    }
    final UniformRandomProvider rg = RngUtils.create(seed.getSeed());
    RandomUtils.shuffle(frames, rg);

    Assertions.assertTrue(source.open());
    for (int i = 0; i < data.length; i++) {
      final int frame = frames[i];
      Assertions.assertTrue(source.isValid(frame));
      Assertions.assertArrayEquals(data[frame - 1], source.get(frame));
    }
    Assertions.assertFalse(source.isValid(0));
    Assertions.assertFalse(source.isValid(data.length + 1));
  }

  @Test
  public void memoryImageSourceCanReturnCroppedDataWithNext() {
    final int width = 5;
    final int height = 3;
    final Rectangle bounds = new Rectangle(2, 1, 3, 1);
    final float[][] data = createData(width, height, 15);
    final MemoryImageSource source = new MemoryImageSource(width, height, data);
    int index = 0;
    float[] next = null;
    Assertions.assertTrue(source.open());
    while ((next = source.next(bounds)) != null) {
      Assertions.assertEquals(bounds.width * bounds.height, next.length);
      Assertions.assertArrayEquals(crop(data[index], width, bounds), next);
      index++;
    }
    Assertions.assertEquals(index, data.length);
  }

  @SeededTest
  public void memoryImageSourceCanReturnCroppedDataWithGet(RandomSeed seed) {
    final int width = 5;
    final int height = 3;
    final Rectangle bounds = new Rectangle(2, 1, 3, 1);
    final float[][] data = createData(width, height, 15);
    final MemoryImageSource source = new MemoryImageSource(width, height, data);

    final int[] frames = new int[data.length];
    for (int i = 0; i < data.length; i++) {
      frames[i] = i + 1;
    }
    final UniformRandomProvider rg = RngUtils.create(seed.getSeed());
    RandomUtils.shuffle(frames, rg);

    Assertions.assertTrue(source.open());
    for (int i = 0; i < data.length; i++) {
      final int frame = frames[i];
      Assertions.assertTrue(source.isValid(frame));
      final float[] next = source.get(frame, bounds);
      Assertions.assertEquals(bounds.width * bounds.height, next.length);
      Assertions.assertArrayEquals(crop(data[frame - 1], width, bounds), next);
    }
    Assertions.assertFalse(source.isValid(0));
    Assertions.assertFalse(source.isValid(data.length + 1));
  }

  @Test
  public void memoryImageSourceThrowsWithInvalidBounds() {
    final int width = 5;
    final int height = 3;
    final float[][] data = createData(width, height, 15);
    final MemoryImageSource source = new MemoryImageSource(width, height, data);
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      source.next(new Rectangle(-1, 0, 1, 1));
    });
  }

  @Test
  public void canConstructInterlacedImageSource() {
    final int width = 5;
    final int height = 3;
    final int n = 15;
    final int start = 4;
    final int size = 2;
    final int skip = 1;
    final MemoryImageSource source =
        new MemoryImageSource(width, height, createData(width, height, n));
    final InterlacedImageSource iSource = new InterlacedImageSource(source, start, size, skip);
    Assertions.assertEquals(width, iSource.getWidth());
    Assertions.assertEquals(height, iSource.getHeight());
    Assertions.assertEquals(12 * 2 / 3, iSource.getFrames());
    Assertions.assertEquals(start, iSource.getStart());
    Assertions.assertEquals(size, iSource.getSize());
    Assertions.assertEquals(skip, iSource.getSkip());
  }

  @Test
  public void nullImageSourceThrowsConstructInterlacedImageSource() {
    final int start = 4;
    final int size = 2;
    final int skip = 1;
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      new InterlacedImageSource(null, start, size, skip);
    });
  }

  @Test
  public void invalidStartThrowsConstructInterlacedImageSource() {
    final int width = 5;
    final int height = 3;
    final int n = 15;
    final int start = 0;
    final int size = 2;
    final int skip = 1;
    final MemoryImageSource m = new MemoryImageSource(width, height, createData(width, height, n));
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      new InterlacedImageSource(m, start, size, skip);
    });
  }

  @Test
  public void invalidSizeThrowsConstructInterlacedImageSource() {
    final int width = 5;
    final int height = 3;
    final int n = 15;
    final int start = 4;
    final int size = 0;
    final int skip = 1;
    final MemoryImageSource m = new MemoryImageSource(width, height, createData(width, height, n));
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      new InterlacedImageSource(m, start, size, skip);
    });
  }

  @Test
  public void invalidSkipThrowsConstructInterlacedImageSource() {
    final int width = 5;
    final int height = 3;
    final int n = 15;
    final int start = 4;
    final int size = 2;
    final int skip = -1;
    final MemoryImageSource m = new MemoryImageSource(width, height, createData(width, height, n));
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      new InterlacedImageSource(m, start, size, skip);
    });
  }

  @Test
  public void interlacedImageSourceCanReturnDataWithNext() {
    final int width = 5;
    final int height = 3;
    final float[][] data = createData(width, height, 15);
    final int start = 4;
    final int size = 2;
    final int skip = 1;
    final ImageSource source =
        new InterlacedImageSource(new MemoryImageSource(width, height, data), start, size, skip);

    Assertions.assertTrue(source.open());

    final int[] expected = new int[] {4, 5, 7, 8, 10, 11, 13, 14};
    int index = 0;
    float[] next = null;
    while ((next = source.next()) != null) {
      final int startFrame = source.getStartFrameNumber();
      final int endFrame = source.getEndFrameNumber();
      Assertions.assertEquals(startFrame, endFrame, "Start and end frames do not match");
      Assertions.assertEquals(expected[index], startFrame);
      Assertions.assertArrayEquals(data[startFrame - 1], next);
      index++;
    }
    Assertions.assertEquals(index, source.getFrames());
  }

  @SeededTest
  public void interlacedImageSourceCanReturnDataWithGet(RandomSeed seed) {
    final int width = 5;
    final int height = 3;
    final float[][] data = createData(width, height, 15);
    final int start = 4;
    final int size = 2;
    final int skip = 1;
    final ImageSource source =
        new InterlacedImageSource(new MemoryImageSource(width, height, data), start, size, skip);

    final int[] frames = new int[data.length];
    for (int i = 0; i < data.length; i++) {
      frames[i] = i + 1;
    }
    final UniformRandomProvider rg = RngUtils.create(seed.getSeed());
    RandomUtils.shuffle(frames, rg);

    final int[] expected = new int[] {4, 5, 7, 8, 10, 11, 13, 14};
    Assertions.assertTrue(source.open());
    for (int i = 0; i < data.length; i++) {
      final int frame = frames[i];
      Assertions.assertTrue(source.isValid(frame));
      final float[] d = source.get(frame);
      if (isExpected(frame, expected)) {
        Assertions.assertArrayEquals(data[frame - 1], d);
      } else {
        Assertions.assertNull(d);
      }
    }
    Assertions.assertFalse(source.isValid(0));
    Assertions.assertFalse(source.isValid(data.length + 1));
  }

  @Test
  public void interlacedImageSourceCanReturnCroppedDataWithNext() {
    final int width = 5;
    final int height = 3;
    final float[][] data = createData(width, height, 15);
    final int start = 4;
    final int size = 2;
    final int skip = 1;
    final Rectangle bounds = new Rectangle(2, 1, 3, 1);
    final ImageSource source =
        new InterlacedImageSource(new MemoryImageSource(width, height, data), start, size, skip);

    Assertions.assertTrue(source.open());

    final int[] expected = new int[] {4, 5, 7, 8, 10, 11, 13, 14};
    int index = 0;
    float[] next = null;
    while ((next = source.next(bounds)) != null) {
      final int startFrame = source.getStartFrameNumber();
      final int endFrame = source.getEndFrameNumber();
      Assertions.assertEquals(startFrame, endFrame, "Start and end frames do not match");
      Assertions.assertEquals(expected[index], startFrame);
      Assertions.assertEquals(bounds.width * bounds.height, next.length);
      Assertions.assertArrayEquals(crop(data[startFrame - 1], width, bounds), next);
      index++;
    }
    Assertions.assertEquals(index, source.getFrames());
  }

  @SeededTest
  public void interlacedImageSourceCanReturnCroppedDataWithGet(RandomSeed seed) {
    final int width = 5;
    final int height = 3;
    final float[][] data = createData(width, height, 15);
    final int start = 4;
    final int size = 2;
    final int skip = 1;
    final Rectangle bounds = new Rectangle(2, 1, 3, 1);
    final ImageSource source =
        new InterlacedImageSource(new MemoryImageSource(width, height, data), start, size, skip);

    final int[] frames = new int[data.length];
    for (int i = 0; i < data.length; i++) {
      frames[i] = i + 1;
    }
    final UniformRandomProvider rg = RngUtils.create(seed.getSeed());
    RandomUtils.shuffle(frames, rg);

    final int[] expected = new int[] {4, 5, 7, 8, 10, 11, 13, 14};
    Assertions.assertTrue(source.open());
    for (int i = 0; i < data.length; i++) {
      final int frame = frames[i];
      Assertions.assertTrue(source.isValid(frame));
      final float[] next = source.get(frame, bounds);
      if (isExpected(frame, expected)) {
        Assertions.assertEquals(bounds.width * bounds.height, next.length);
        Assertions.assertArrayEquals(crop(data[frame - 1], width, bounds), next);
      } else {
        Assertions.assertNull(next);
      }
    }
    Assertions.assertFalse(source.isValid(0));
    Assertions.assertFalse(source.isValid(data.length + 1));
  }

  @Test
  public void canConstructAggregatedImageSource() {
    final int width = 5;
    final int height = 3;
    final int n = 15;
    final int aggregate = 3;
    final MemoryImageSource source =
        new MemoryImageSource(width, height, createData(width, height, n));
    final AggregatedImageSource iSource = new AggregatedImageSource(source, aggregate);
    Assertions.assertEquals(width, iSource.getWidth());
    Assertions.assertEquals(height, iSource.getHeight());
    Assertions.assertEquals(n / aggregate, iSource.getFrames());
    Assertions.assertEquals(aggregate, iSource.getAggregate());
  }

  @Test
  public void nullImageSourceThrowsConstructAggregatedImageSource() {
    final int aggregate = 3;
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      new AggregatedImageSource(null, aggregate);
    });
  }

  @Test
  public void invalidAggregateThrowsConstructAggregatedImageSource() {
    final int width = 5;
    final int height = 3;
    final int n = 15;
    final int aggregate = 1;
    final MemoryImageSource m = new MemoryImageSource(width, height, createData(width, height, n));
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      new AggregatedImageSource(m, aggregate);
    });
  }

  @Test
  public void aggregatedImageSourceCanReturnDataWithNext() {
    final int width = 5;
    final int height = 3;
    final int aggregate = 3;
    final float[][] data = createData(width, height, 15);
    final ImageSource source =
        new AggregatedImageSource(new MemoryImageSource(width, height, data), aggregate);

    Assertions.assertTrue(source.open());

    int index = 1;
    int ii = 0;
    float[] next = null;
    while ((next = source.next()) != null) {
      ii++;
      Assertions.assertEquals(index, source.getStartFrameNumber());
      Assertions.assertEquals(index + 2, source.getEndFrameNumber());
      final float[] all = combine(data[index - 1], data[index], data[index + 1]);
      Assertions.assertArrayEquals(all, next);
      index += 3;
    }
    Assertions.assertEquals(ii, source.getFrames());
  }

  @SeededTest
  public void aggregatedImageSourceCanReturnDataWithGet(RandomSeed seed) {
    final int width = 5;
    final int height = 3;
    final int aggregate = 3;
    final float[][] data = createData(width, height, 15);
    final ImageSource source =
        new AggregatedImageSource(new MemoryImageSource(width, height, data), aggregate);

    final int[] frames = new int[data.length / 3];
    for (int i = 0, frame = 1; i < frames.length; i++, frame += 3) {
      frames[i] = frame;
    }
    final UniformRandomProvider rg = RngUtils.create(seed.getSeed());
    RandomUtils.shuffle(frames, rg);

    Assertions.assertTrue(source.open());
    for (int i = 0; i < frames.length; i++) {
      final int frame = frames[i];
      Assertions.assertTrue(source.isValid(frame), () -> "Invalid frame " + frame);
      final float[] d = source.get(frame);
      Assertions.assertEquals(frame, source.getStartFrameNumber());
      Assertions.assertEquals(frame + 2, source.getEndFrameNumber());
      final float[] all = combine(data[frame - 1], data[frame], data[frame + 1]);
      Assertions.assertArrayEquals(all, d, () -> "Invalid frame data " + frame);
    }
    Assertions.assertFalse(source.isValid(0));
    Assertions.assertFalse(source.isValid(data.length + 1));
  }

  @Test
  public void aggregatedImageSourceCanReturnCroppedDataWithNext() {
    final int width = 5;
    final int height = 3;
    final int aggregate = 3;
    final float[][] data = createData(width, height, 15);
    final Rectangle bounds = new Rectangle(2, 1, 3, 1);
    final ImageSource source =
        new AggregatedImageSource(new MemoryImageSource(width, height, data), aggregate);

    Assertions.assertTrue(source.open());

    int index = 1;
    int ii = 0;
    float[] next = null;
    while ((next = source.next(bounds)) != null) {
      ii++;
      Assertions.assertEquals(index, source.getStartFrameNumber());
      Assertions.assertEquals(index + 2, source.getEndFrameNumber());
      final float[] all = combine(crop(data[index - 1], width, bounds),
          crop(data[index], width, bounds), crop(data[index + 1], width, bounds));
      Assertions.assertArrayEquals(all, next);
      index += 3;
    }
    Assertions.assertEquals(ii, source.getFrames());
  }

  @SeededTest
  public void aggregatedImageSourceCanReturnCroppedDataWithGet(RandomSeed seed) {
    final int width = 5;
    final int height = 3;
    final int aggregate = 3;
    final float[][] data = createData(width, height, 15);
    final Rectangle bounds = new Rectangle(2, 1, 3, 1);
    final ImageSource source =
        new AggregatedImageSource(new MemoryImageSource(width, height, data), aggregate);

    final int[] frames = new int[data.length / 3];
    for (int i = 0, frame = 1; i < frames.length; i++, frame += 3) {
      frames[i] = frame;
    }
    final UniformRandomProvider rg = RngUtils.create(seed.getSeed());
    RandomUtils.shuffle(frames, rg);

    Assertions.assertTrue(source.open());
    for (int i = 0; i < frames.length; i++) {
      final int frame = frames[i];
      Assertions.assertTrue(source.isValid(frame));
      final float[] d = source.get(frame, bounds);
      Assertions.assertEquals(frame, source.getStartFrameNumber());
      Assertions.assertEquals(frame + 2, source.getEndFrameNumber());
      final float[] all = combine(crop(data[frame - 1], width, bounds),
          crop(data[frame], width, bounds), crop(data[frame + 1], width, bounds));
      Assertions.assertArrayEquals(all, d, () -> "Invalid frame data " + frame);
    }
    Assertions.assertFalse(source.isValid(0));
    Assertions.assertFalse(source.isValid(data.length + 1));
  }

  @Test
  public void canConstructAggregatedInterlacedImageSource() {
    final int width = 5;
    final int height = 3;
    final int n = 15;
    final float[][] data = createData(width, height, n);
    final int aggregate = 3;
    final int start = 4;
    final int size = 2;
    final int skip = 1;
    final ImageSource source =
        new InterlacedImageSource(new MemoryImageSource(width, height, data), start, size, skip);
    final AggregatedImageSource iSource = new AggregatedImageSource(source, aggregate);
    Assertions.assertEquals(width, iSource.getWidth());
    Assertions.assertEquals(height, iSource.getHeight());
    Assertions.assertEquals(12 * 2 / 3, source.getFrames());
    Assertions.assertEquals((int) Math.ceil((double) source.getFrames() / aggregate),
        iSource.getFrames());
    Assertions.assertEquals(aggregate, iSource.getAggregate());
  }

  @Test
  public void constructInterlacedAggregatedImageSourceThrows() {
    final int width = 5;
    final int height = 3;
    final int n = 15;
    final float[][] data = createData(width, height, n);
    final int aggregate = 3;
    final int start = 4;
    final int size = 2;
    final int skip = 1;
    final ImageSource source =
        new AggregatedImageSource(new MemoryImageSource(width, height, data), aggregate);
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      new InterlacedImageSource(source, start, size, skip);
    });
  }

  @Test
  public void aggregatedInterlacedImageSourceCanReturnDataWithNext() {
    final int width = 5;
    final int height = 3;
    final int n = 15;
    final float[][] data = createData(width, height, n);
    final int aggregate = 3;
    final int start = 4;
    final int size = 2;
    final int skip = 1;
    final ImageSource source = new AggregatedImageSource(
        new InterlacedImageSource(new MemoryImageSource(width, height, data), start, size, skip),
        aggregate);

    Assertions.assertTrue(source.open());

    // Set the expected frames returned by the interlacing
    final int[] expected = new int[] {4, 5, 7, 8, 10, 11, 13, 14};

    int i1 = 0;
    int i2 = 0;
    float[] next = null;
    while ((next = source.next()) != null) {
      // Get the range for the data
      final int endE = FastMath.min(i1 + 2, expected.length - 1);
      final int startFrame = expected[i1];
      final int endFrame = expected[endE];

      // Check the correct range is returned
      Assertions.assertEquals(startFrame, source.getStartFrameNumber());
      Assertions.assertEquals(endFrame, source.getEndFrameNumber());

      // Check the data is collated correctly
      final float[] all = new float[data[0].length];
      for (int e = i1; e <= endE; e++) {
        final int frame = expected[e] - 1;
        for (int j = 0; j < all.length; j++) {
          all[j] += data[frame][j];
        }
      }
      Assertions.assertArrayEquals(all, next);
      i1 += 3;
      i2++;
    }
    Assertions.assertEquals(i2, source.getFrames());
  }

  @SeededTest
  public void aggregatedInterlacedImageSourceCanReturnDataWithGet(RandomSeed seed) {
    final int width = 5;
    final int height = 3;
    final int n = 15;
    final float[][] data = createData(width, height, n);
    final int aggregate = 3;
    final int start = 4;
    final int size = 2;
    final int skip = 1;
    final ImageSource source = new AggregatedImageSource(
        new InterlacedImageSource(new MemoryImageSource(width, height, data), start, size, skip),
        aggregate);

    // Set the expected frames returned by the interlacing
    final int[] expected = new int[] {4, 5, 7, 8, 10, 11, 13, 14};

    // Randomly pick points from the positions used by next()
    final int[] frames = new int[source.getFrames()];
    for (int i = 0, ii = 0; ii < expected.length; i++, ii += 3) {
      frames[i] = ii;
    }
    final UniformRandomProvider rg = RngUtils.create(seed.getSeed());
    RandomUtils.shuffle(frames, rg);

    Assertions.assertTrue(source.open());
    for (int i = 0; i < frames.length; i++) {
      // Get the range for the data
      int startE = frames[i];
      final int endE = FastMath.min(startE + 2, expected.length - 1);
      final int startFrame = expected[startE];
      final int endFrame = expected[endE];

      // Get the data
      final float[] d = source.get(startFrame);

      // Check the correct range is returned
      Assertions.assertEquals(startFrame, source.getStartFrameNumber());
      Assertions.assertEquals(endFrame, source.getEndFrameNumber());

      // Check the data is collated correctly
      final float[] all = new float[data[0].length];
      for (; startE <= endE; startE++) {
        final int frame = expected[startE] - 1;
        for (int j = 0; j < all.length; j++) {
          all[j] += data[frame][j];
        }
      }
      Assertions.assertArrayEquals(all, d);
    }

    // Check all the data is valid but skipped interlaced points return null
    for (int i = 0; i < data.length; i++) {
      final int frame = i + 1;
      Assertions.assertTrue(source.isValid(frame));
      if (!isExpected(frame, expected)) {
        Assertions.assertNull(source.get(frame));
      }
    }

    Assertions.assertFalse(source.isValid(0));
    Assertions.assertFalse(source.isValid(data.length + 1));
  }

  @Test
  public void canSerialiseMemoryImageSource() {
    final int width = 5;
    final int height = 3;
    final int n = 15;
    final String name = "canSerialiseMemoryImageSource";

    final MemoryImageSource source =
        new MemoryImageSource(width, height, createData(width, height, n));
    source.setName(name);
    source.setFreeMemoryOnClose(true);

    final String xml = source.toXml();
    // System.out.println(xml);

    final MemoryImageSource source2 = (MemoryImageSource) ImageSource.fromXml(xml);

    Assertions.assertEquals(width, source2.getWidth());
    Assertions.assertEquals(height, source2.getHeight());
    Assertions.assertEquals(n, source2.getFrames());
    Assertions.assertEquals(name, source2.getName());
    Assertions.assertEquals(true, source2.isFreeMemoryOnClose());

    float[] data;
    while ((data = source.next()) != null) {
      final float[] data2 = source2.next();
      Assertions.assertArrayEquals(data, data2, 1e-6f);
    }
  }

  /**
   * Create data using the specified dimensions and the number of frames. Each frame will have a
   * different base number and each index will be unique in the frame.
   *
   * @param width width
   * @param height height
   * @param n The number of frames
   * @return The data
   */
  private static float[][] createData(int width, int height, int n) {
    final List<float[]> data = new ArrayList<>(n);
    float base = 1;
    for (int i = 0; i < n; i++) {
      data.add(createData(width, height, base));
      base *= 2;
    }
    return data.toArray(new float[0][0]);
  }

  /**
   * Create data using the specified dimensions and the base number. Each index will be unique in
   * the frame.
   *
   * @param width width
   * @param height height
   * @param base the base number
   * @return The data
   */
  private static float[] createData(int width, int height, float base) {
    final float[] data = new float[width * height];
    Arrays.fill(data, base);
    for (int i = 0; i < data.length; i++) {
      data[i] += (double) i / data.length;
    }
    return data;
  }

  /**
   * Check if the frame is contained in the expected frames array.
   *
   * @param frame the frame
   * @param expected the expected frames
   * @return true if within the array
   */
  private static boolean isExpected(int frame, int[] expected) {
    for (final int e : expected) {
      if (e == frame) {
        return true;
      }
    }
    return false;
  }

  /**
   * Crop the data to the specified bounds.
   *
   * @param data the data
   * @param width the width
   * @param bounds the bounds
   * @return the cropped data
   */
  float[] crop(float[] data, int width, Rectangle bounds) {
    final float[] newData = new float[bounds.width * bounds.height];
    for (int y = bounds.y, ii = 0; y < bounds.y + bounds.height; y++) {
      for (int x = bounds.x; x < bounds.x + bounds.width; x++, ii++) {
        final int index = y * width + x;
        newData[ii] = data[index];
      }
    }
    return newData;
  }

  /**
   * Sum all the input arrays.
   *
   * @param arrays the arrays
   * @return The summed array
   */
  private static float[] combine(float[]... arrays) {
    final float[] all = Arrays.copyOf(arrays[0], arrays[0].length);
    for (int i = 1; i < arrays.length; i++) {
      for (int j = 0; j < all.length; j++) {
        all[j] += arrays[i][j];
      }
    }
    return all;
  }
}
