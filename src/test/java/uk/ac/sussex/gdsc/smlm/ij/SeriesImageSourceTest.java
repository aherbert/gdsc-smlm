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

package uk.ac.sussex.gdsc.smlm.ij;

import uk.ac.sussex.gdsc.core.utils.RandomUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.results.ImageSource.ReadHint;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;

import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileInfo;
import ij.io.TiffEncoder;
import ij.measure.Calibration;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

@SuppressWarnings({"javadoc"})
public class SeriesImageSourceTest {
  int w = 10, h = 5, d = 7;

  @Test
  public void canReadBigTIFFSequentially() throws IOException {
    canReadBigTIFFSequentially(false, true);
  }

  @Test
  public void canReadBigTIFFSequentiallyInMemory() throws IOException {
    canReadBigTIFFSequentially(true, true);
  }

  @SeededTest
  public void canReadBigTIFFSequentiallyBE() throws IOException {
    canReadBigTIFFSequentially(false, false);
  }

  @SeededTest
  public void canReadBigTIFFSequentiallyInMemoryBE() throws IOException {
    canReadBigTIFFSequentially(true, false);
  }

  private void canReadBigTIFFSequentially(boolean inMemory, boolean intelByteOrder)
      throws IOException {
    final int n = 2;
    final String[] filenames = createFilenames(n);
    final ImageStack[] stacks = createSeries(filenames, intelByteOrder);
    final SeriesImageSource source = new SeriesImageSource("Test", filenames);
    if (!inMemory) {
      source.setBufferLimit(0); // To force standard reading functionality
    }
    source.setReadHint(ReadHint.SEQUENTIAL);
    source.open();
    Assertions.assertEquals(w, source.getWidth());
    Assertions.assertEquals(h, source.getHeight());
    Assertions.assertEquals(d * n, source.getFrames());
    for (int i = 0; i < stacks.length; i++) {
      for (int j = 0; j < d; j++) {
        final float[] e = (float[]) stacks[i].getPixels(j + 1);
        final float[] o = source.next();
        Assertions.assertArrayEquals(e, o);
      }
    }
    Assertions.assertNull(source.next());
    source.close();
  }

  @SeededTest
  public void canReadBigTIFFNonSequentially(RandomSeed seed) throws IOException {
    canReadBigTIFFNonSequentially(seed, false, true);
  }

  @SeededTest
  public void canReadBigTIFFNonSequentiallyInMemory(RandomSeed seed) throws IOException {
    canReadBigTIFFNonSequentially(seed, true, true);
  }

  @SeededTest
  public void canReadBigTIFFNonSequentiallyBE(RandomSeed seed) throws IOException {
    canReadBigTIFFNonSequentially(seed, false, false);
  }

  @SeededTest
  public void canReadBigTIFFNonSequentiallyInMemoryBE(RandomSeed seed) throws IOException {
    canReadBigTIFFNonSequentially(seed, true, false);
  }

  private void canReadBigTIFFNonSequentially(RandomSeed seed, boolean inMemory,
      boolean intelByteOrder) throws IOException {
    final int n = 2;
    final String[] filenames = createFilenames(n);
    final ImageStack[] stacks = createSeries(filenames, intelByteOrder);
    final SeriesImageSource source = new SeriesImageSource("Test", filenames);
    if (!inMemory) {
      source.setBufferLimit(0); // To force standard reading functionality
    }
    source.setReadHint(ReadHint.NONSEQUENTIAL);
    source.open();
    Assertions.assertEquals(w, source.getWidth());
    Assertions.assertEquals(h, source.getHeight());
    Assertions.assertEquals(d * n, source.getFrames());
    final float[][] pixels = new float[n * d][];
    for (int i = 0, k = 0; i < stacks.length; i++) {
      for (int j = 0; j < d; j++) {
        pixels[k++] = (float[]) stacks[i].getPixels(j + 1);
      }
    }

    final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
    for (int i = 0; i < 3; i++) {
      final int[] random = RandomUtils.sample(pixels.length / 2, pixels.length, r);
      for (final int frame : random) {
        // logger.fine(FunctionUtils.getSupplier("[%d] frame = %d", i, frame);
        final float[] e = pixels[frame];
        final float[] o = source.get(frame + 1); // 1-base index on the frame
        Assertions.assertArrayEquals(e, o);
      }
    }
  }

  private String[] createFilenames(int n) throws IOException {
    final String[] filenames = new String[n];
    for (int i = 0; i < n; i++) {
      final File path = File.createTempFile(this.getClass().getSimpleName(), ".tif");
      path.deleteOnExit();
      filenames[i] = path.getCanonicalPath();
    }
    return filenames;
  }

  private ImageStack[] createSeries(String[] filenames, boolean intelByteOrder) throws IOException {
    final int n = filenames.length;
    final ImageStack[] stacks = new ImageStack[n];
    int index = 0;
    final int length = w * h;
    for (int i = 0; i < n; i++) {
      final ImageStack stack = new ImageStack(w, h);
      for (int j = 0; j < d; j++) {
        stack.addSlice(null, SimpleArrayUtils.newArray(length, index, 1f));
        index += length;
      }
      final ImagePlus imp = new ImagePlus(null, stack);
      // Add a calibration with origin
      final Calibration c = imp.getCalibration();
      c.xOrigin = 4;
      c.yOrigin = 5;
      saveAsTiff(imp, filenames[i], intelByteOrder);
      stacks[i] = stack;
    }
    return stacks;
  }

  private static void saveAsTiff(ImagePlus imp, String path, boolean intelByteOrder)
      throws IOException {
    // IJ.saveAsTiff(imp, path);

    final FileInfo fi = imp.getFileInfo();
    fi.nImages = imp.getStackSize();
    ij.Prefs.intelByteOrder = intelByteOrder;
    try (DataOutputStream out =
        new DataOutputStream(new BufferedOutputStream(new FileOutputStream(path)))) {
      final TiffEncoder file = new TiffEncoder(fi);
      file.write(out);
    }
  }
}
