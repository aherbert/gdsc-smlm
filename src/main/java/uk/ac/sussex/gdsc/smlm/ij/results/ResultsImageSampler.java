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

import gnu.trove.list.array.TLongArrayList;
import gnu.trove.map.hash.TLongObjectHashMap;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.process.ImageProcessor;
import java.awt.Rectangle;
import java.util.Arrays;
import java.util.Comparator;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.core.source64.SplitMix64;
import org.apache.commons.rng.simple.RandomSource;
import uk.ac.sussex.gdsc.core.ij.gui.OffsetPointRoi;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.rng.RandomUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Allows sampling sections from a source image for localisation results.
 */
public class ResultsImageSampler {
  private final MemoryPeakResults results;
  private final ImageStack stack;

  private final int lx;
  private final int ly;
  private final int xblocks;
  private final int yblocks;
  private final int xyblocks;

  /** The size for samples. */
  public final int size;

  /**
   * The max number of empty samples. Since they are empty it should not matter unless the noise
   * characteristics change over the image duration. Set to 0 to sample throughout the lifetime of
   * the localisation occurrences.
   */
  private int maxNumberOfEmptySamples = 500;

  private long[] no;
  private ResultsSample[] data;
  private int lower;
  private int upper;
  private final LocalList<ResultsSample> sampleList = new LocalList<>();
  private UniformRandomProvider rng = new SplitMix64(RandomSource.createLong());

  private static class PeakResultList {
    int size;
    PeakResult[] data;

    PeakResultList() {
      this(new PeakResult[1]);
    }

    PeakResultList(PeakResult[] data) {
      this.data = data;
    }

    void add(PeakResult peakResult) {
      data[size++] = peakResult;
      if (size == data.length) {
        final PeakResult[] data2 = new PeakResult[size * 2];
        System.arraycopy(data, 0, data2, 0, size);
        data = data2;
      }
    }

    PeakResult get(int index) {
      return data[index];
    }
  }

  private static class ResultsSample {
    long index;
    final PeakResultList list;

    static ResultsSample createEmpty(long index) {
      return new ResultsSample(index, new PeakResultList(null));
    }

    static ResultsSample create(long index) {
      return new ResultsSample(index, new PeakResultList());
    }

    ResultsSample(long index, PeakResultList data) {
      this.index = index;
      this.list = data;
    }

    public int size() {
      return list.size;
    }

    public void add(PeakResult peakResult) {
      list.add(peakResult);
    }
  }

  private static enum IndexComparator implements Comparator<ResultsSample> {
    /** An instance of the comparator. */
    INSTANCE;

    @Override
    public int compare(ResultsSample o1, ResultsSample o2) {
      return Long.compare(o1.index, o2.index);
    }
  }

  private static enum CountComparator implements Comparator<ResultsSample> {
    /** An instance of the comparator. */
    INSTANCE;

    @Override
    public int compare(ResultsSample o1, ResultsSample o2) {
      final int result = Integer.compare(o1.size(), o2.size());
      if (result == 0) {
        // Use index if the same count
        return Long.compare(o1.index, o2.index);
      }
      return result;
    }
  }

  private static enum ReverseCountComparator implements Comparator<ResultsSample> {
    /** An instance of the comparator. */
    INSTANCE;

    @Override
    public int compare(ResultsSample o1, ResultsSample o2) {
      final int result = Integer.compare(o2.size(), o1.size());
      if (result == 0) {
        // Use index if the same count
        return Long.compare(o1.index, o2.index);
      }
      return result;
    }
  }

  /**
   * Instantiates a new results image sampler.
   *
   * @param results the results (must be in pixel units)
   * @param stack the source stack for the results
   * @param size the size of the sample blocks
   */
  public ResultsImageSampler(MemoryPeakResults results, ImageStack stack, int size) {
    if (results.getDistanceUnit() != DistanceUnit.PIXEL) {
      throw new IllegalArgumentException("Results must be in pixel units");
    }

    this.results = results;
    this.stack = stack;
    this.size = size;

    final Rectangle bounds = results.getBounds(true);
    // Round the image dimensions to the nearest block interval
    lx = size * (bounds.x / size);
    ly = size * (bounds.y / size);
    final int ux = size * (int) Math.ceil((double) (bounds.x + bounds.width) / size);
    final int uy = size * (int) Math.ceil((double) (bounds.y + bounds.height) / size);
    xblocks = (ux - lx) / size;
    yblocks = (uy - ly) / size;
    xyblocks = xblocks * yblocks;
  }

  /**
   * Analyse the input data to allow sampling.
   *
   * @return true, if successful
   */
  public boolean analyse() {
    // Clear old samples
    data = null;
    no = null;

    if (results.isEmpty()) {
      return false;
    }

    createResultSamples();

    // Split the image into frames with zero localisations, low density, high density

    createEmptySampleIndices();

    // For the low and high sample we just split in half.
    // The data must be sorted by count.
    Arrays.sort(data, CountComparator.INSTANCE);
    lower = data.length / 2;
    upper = data.length - lower;

    return true;
  }

  private void createEmptySampleIndices() {
    // Create the empty blocks
    final long total = (long) stack.getSize() * xyblocks;
    final long empty = total - data.length;
    if (empty == 0) {
      // All indices are used
      no = new long[0];
    } else {
      // Sort by index
      Arrays.sort(data, IndexComparator.INSTANCE);

      // Randomly sample indices that are not used.
      if (maxNumberOfEmptySamples > 0) {
        // Just enumerate the first N. Since they are empty it should not matter
        // unless the noise characteristics change over the image duration.
        long emptyCandidate = 0;
        final long[] list = new long[(int) Math.min(empty, maxNumberOfEmptySamples)];
        int count = 0;
        OUTER: for (int i = 0; i < data.length; i++) {
          final long current = data[i].index;
          // If the current index is bigger than the candidate then it must be empty
          while (current > emptyCandidate) {
            // Add all those that are empty
            list[count++] = emptyCandidate++;
            if (count == list.length) {
              break OUTER;
            }
          }
          // Set the next candidate
          emptyCandidate = current + 1;
        }
        no = Arrays.copyOf(list, count);
      } else {
        // Sample throughout the localisation time course
        final TLongArrayList list = new TLongArrayList(data.length);
        if (empty < data.length) {
          // We can pick all the indices that are missing
          long emptyCandidate = 0;
          for (int i = 0; i < data.length; i++) {
            final long current = data[i].index;
            // If the current index is bigger than the candidate then it must be empty
            while (current > emptyCandidate) {
              // Add all those that are empty
              list.add(emptyCandidate++);
            }
            // Set the next candidate
            emptyCandidate = current + 1;
          }
        } else {
          // There are many empty blocks so just sample blocks
          // after those with localisations.
          long emptyCandidate = 1;
          for (int i = 0; i < data.length; i++) {
            final long current = data[i].index;
            // If the current index is bigger than the candidate then it must be empty
            if (current > emptyCandidate) {
              // Note: we only sample the next empty index after an index with data
              // This means the number of frames could be lower
              list.add(emptyCandidate);
            }
            // Set the next candidate
            emptyCandidate = current + 1;
          }
        }
        no = list.toArray();
      }
    }
  }

  /**
   * Creates the result samples. Do this by storing the coordinates at the region index.
   */
  private void createResultSamples() {
    final TLongObjectHashMap<ResultsSample> map = new TLongObjectHashMap<>(results.size());
    ResultsSample next = ResultsSample.create(-1);
    for (final PeakResult p : results.toArray()) {
      // Avoid invalid slices
      if (p.getFrame() < 1 || p.getFrame() > stack.getSize()) {
        continue;
      }

      final long index = getIndex(p.getXPosition(), p.getYPosition(), p.getFrame());

      ResultsSample current = map.putIfAbsent(index, next);
      if (current == null) {
        // If the return value is null then this is a new insertion.
        // Set the current value as the one we just added and create the next insertion object.
        current = next;
        current.index = index;
        next = ResultsSample.create(-1);
      }
      current.add(p);
    }

    // Create an array of all the sample entries.
    // This is used to sample regions by density.
    data = map.values(new ResultsSample[map.size()]);
  }

  /**
   * Gets the index of the block from the floating point coordinates.
   *
   * @param x the x
   * @param y the y
   * @param time the time frame
   * @return the index
   */
  private long getIndex(float x, float y, int time) {
    // Make frames start at zero for the index
    return (long) ((xyblocks) * (time - 1)) + xblocks * getY(y) + getX(x);
  }

  private int getX(float x) {
    return (int) ((x - lx) / size);
  }

  private int getY(float y) {
    return (int) ((y - ly) / size);
  }

  /**
   * Convert the single index into x,y,z coords, Input array must be length >= 3.
   *
   * @param index the index
   * @param xyz the xyz
   * @return The xyz array
   */
  private int[] getXyz(long index, int[] xyz) {
    xyz[2] = (int) (index / (xyblocks));
    final int mod = (int) (index % (xyblocks));
    xyz[1] = mod / xblocks;
    xyz[0] = mod % xblocks;

    // Convert back to real coords
    xyz[2]++; // Frames start at 1
    xyz[1] = (xyz[1] * size) + ly;
    xyz[0] = (xyz[0] * size) + lx;

    return xyz;
  }

  /**
   * Gets the results used to construct this instance.
   *
   * @return the results
   */
  public MemoryPeakResults getResults() {
    return results;
  }

  /**
   * Checks if is valid (i.e. samples can be obtained).
   *
   * @return true, if is valid
   */
  public boolean isValid() {
    return no != null;
  }

  /**
   * Gets the number of empty samples.
   *
   * @return the number of empty samples
   */
  public int getNumberOfEmptySamples() {
    return no.length;
  }

  /**
   * Gets the number of low density samples.
   *
   * @return the number of low density samples
   */
  public int getNumberOfLowDensitySamples() {
    return lower;
  }

  /**
   * Gets the number of high density samples.
   *
   * @return the number of high density samples
   */
  public int getNumberOfHighDensitySamples() {
    return upper;
  }

  /**
   * Gets the sample image. The image is a stack of the samples with an overlay of the localisation
   * positions. The info property is set with details of the localisations and the image is
   * calibrated.
   *
   * @param countNo the number of samples with no localisations
   * @param countLow the number of samples with low localisations
   * @param countHigh the number of samples with high localisations
   * @return the sample image (could be null if no samples were made)
   */
  public ImagePlus getSample(int countNo, int countLow, int countHigh) {
    final ImageStack out = new ImageStack(size, size);
    if (!isValid()) {
      return null;
    }

    sampleList.clear();

    // empty
    for (final int i : RandomUtils.sample(countNo, no.length, rng)) {
      sampleList.add(ResultsSample.createEmpty(no[i]));
    }
    // low
    for (final int i : RandomUtils.sample(countLow, lower, rng)) {
      sampleList.add(data[i]);
    }
    // high
    for (final int i : RandomUtils.sample(countHigh, upper, rng)) {
      sampleList.add(data[i + lower]);
    }

    if (sampleList.isEmpty()) {
      return null;
    }

    double nmPerPixel = 1;
    if (results.hasCalibration()) {
      final CalibrationReader calibration = results.getCalibrationReader();
      if (calibration.hasNmPerPixel()) {
        nmPerPixel = calibration.getNmPerPixel();
      }
    }

    // Sort descending by number in the frame
    final ResultsSample[] samples = sampleList.toArray(new ResultsSample[0]);
    Arrays.sort(samples, ReverseCountComparator.INSTANCE);

    final int[] xyz = new int[3];
    final Rectangle stackBounds = new Rectangle(stack.getWidth(), stack.getHeight());
    final Overlay overlay = new Overlay();
    final StringBuilder sb = new StringBuilder();
    if (nmPerPixel == 1) {
      sb.append("Sample X Y Z Signal\n");
    } else {
      sb.append("Sample X(nm) Y(nm) Z(nm) Signal\n");
    }

    for (final ResultsSample sample : samples) {
      getXyz(sample.index, xyz);

      // Construct the region to extract
      Rectangle target = new Rectangle(xyz[0], xyz[1], size, size);
      target = target.intersection(stackBounds);
      if (target.width == 0 || target.height == 0) {
        continue;
      }

      // Extract the frame
      final int slice = xyz[2];
      final ImageProcessor ip = stack.getProcessor(slice);

      // Cut out the desired pixels (some may be blank if the block overruns the source image)
      final ImageProcessor ip2 = ip.createProcessor(size, size);
      for (int y = 0; y < target.height; y++) {
        for (int x = 0, i = y * size, index = (y + target.y) * ip.getWidth() + target.x;
            x < target.width; x++, i++, index++) {
          ip2.setf(i, ip.getf(index));
        }
      }

      final int sampleSize = sample.size();
      if (sampleSize > 0) {
        final float[] ox = new float[sampleSize];
        final float[] oy = new float[sampleSize];
        final int position = out.getSize() + 1;
        // Create an ROI with the localisations
        for (int i = 0; i < sampleSize; i++) {
          final PeakResult p = sample.list.get(i);
          ox[i] = p.getXPosition() - xyz[0];
          oy[i] = p.getYPosition() - xyz[1];
          sb.append(position).append(' ');
          sb.append(MathUtils.rounded(ox[i] * nmPerPixel)).append(' ');
          sb.append(MathUtils.rounded(oy[i] * nmPerPixel)).append(' ');
          sb.append(MathUtils.rounded(p.getZPosition() * nmPerPixel)).append(' ');
          sb.append(MathUtils.rounded(p.getIntensity())).append('\n');
        }
        final PointRoi roi = new OffsetPointRoi(ox, oy, sampleSize);
        roi.setPosition(position);
        overlay.add(roi);
      }

      out.addSlice(String.format("Frame=%d @ %d,%d px (n=%d)", slice, xyz[0], xyz[1], sampleSize),
          ip2.getPixels());
    }

    if (out.getSize() == 0) {
      return null;
    }

    final ImagePlus imp = new ImagePlus("Sample", out);
    imp.setOverlay(overlay);
    // Note: Only the info property can be saved to a TIFF file
    imp.setProperty("Info", sb.toString());
    if (nmPerPixel != 1) {
      final ij.measure.Calibration cal = new ij.measure.Calibration();
      cal.setUnit("nm");
      cal.pixelHeight = cal.pixelWidth = nmPerPixel;
      imp.setCalibration(cal);
    }

    return imp;
  }

  /**
   * Set the random generator for use during sampling.
   *
   * @param rng the random generator to set (ignored if null)
   */
  public void setRandom(UniformRandomProvider rng) {
    if (rng != null) {
      this.rng = rng;
    }
  }

  /**
   * Gets the max number of empty samples.
   *
   * @return the max number of empty samples
   */
  public int getMaxNumberOfEmptySamples() {
    return maxNumberOfEmptySamples;
  }

  /**
   * Sets the max number of empty samples.
   *
   * <p>Since they are empty it should not matter unless the noise characteristics change over the
   * image duration. Set to 0 to sample throughout the lifetime of the localisation occurrences.
   *
   * @param maxNumberOfEmptySamples the new max number of empty samples
   */
  public void setMaxNumberOfEmptySamples(int maxNumberOfEmptySamples) {
    this.maxNumberOfEmptySamples = maxNumberOfEmptySamples;
  }
}
