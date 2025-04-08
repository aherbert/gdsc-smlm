/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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

import com.thoughtworks.xstream.annotations.XStreamOmitField;
import uk.ac.sussex.gdsc.core.utils.ValidationUtils;

/**
 * Wraps an image source and allows data to be interlaced within regular blocks.
 */
public class InterlacedImageSource extends ImageSource {
  private final int start;
  private final int size;
  private final int skip;
  private final ImageSource imageSource;

  // Record the number of frames returned from the current block
  @XStreamOmitField
  private int counter;

  /**
   * Create a new interlaced image source using the given image source.
   *
   * <p>Note: The input source cannot be aggregated as the data to interlace is assumed to be
   * contiguous from frame 1.
   *
   * @param imageSource The image source of interlaced data (must not be null or an
   *        AggregatedImageSource)
   * @param start The first frame that contains data
   * @param size The number of continuous frames containing data
   * @param skip The number of continuous frames to ignore before the next data
   * @throws IllegalArgumentException If the image is aggregated, or the interlace arguments
   *         are invalid
   */
  public InterlacedImageSource(ImageSource imageSource, int start, int size, int skip) {
    super("");
    ValidationUtils.checkNotNull(imageSource, "Image source must not be null");
    if (imageSource instanceof AggregatedImageSource) {
      throw new IllegalArgumentException("Image source must not be aggregated");
    }
    ValidationUtils.checkStrictlyPositive(start, "start frame");
    ValidationUtils.checkStrictlyPositive(size, "read size");
    ValidationUtils.checkPositive(skip, "skip length");
    setName(String.format("Interlaced (%d,%d,%d) %s", start, size, skip, imageSource.getName()));
    this.imageSource = imageSource;
    this.start = start;
    this.size = size;
    this.skip = skip;
  }

  @Override
  public int getXOrigin() {
    return imageSource.getXOrigin();
  }

  @Override
  public int getYOrigin() {
    return imageSource.getYOrigin();
  }

  @Override
  public int getWidth() {
    return imageSource.getWidth();
  }

  @Override
  public int getHeight() {
    return imageSource.getHeight();
  }

  @Override
  public int getFrames() {
    return (int) Math.ceil((imageSource.getFrames() - start + 1) * ((double) size / (size + skip)));
  }

  @Override
  public ImageSource getParent() {
    return imageSource;
  }

  @Override
  public ImageSource getOriginal() {
    return imageSource.getOriginal();
  }

  @Override
  protected boolean openSource() {
    return imageSource.openSource();
  }

  @Override
  protected void closeSource() {
    imageSource.closeSource();
  }

  @Override
  protected boolean initialiseSequentialRead() {
    if (imageSource.initialiseSequentialRead()) {
      imageSource.sequentialReadStatus = SequentialReadStatus.RUNNING;

      // Assume frame start at 1 and set the initial skip
      final int initialSkip = (start - 1);
      counter = -initialSkip;
      return true;
    }
    imageSource.sequentialReadStatus = SequentialReadStatus.CLOSED;
    return false;

  }

  @Override
  protected Object nextRawFrame() {
    // Skip frames until at the start of the next block
    while (counter < 0) {
      if (imageSource.nextRaw() == null) {
        return null;
      }
      counter++;
    }

    // Read the next frame in the current block
    final Object pixels = imageSource.nextRaw();

    // Check if this is the final frame in the current block
    if (++counter >= size) {
      counter = -skip;
    }
    // Set the frame to the last one read from the source
    setFrameNumber(imageSource.getStartFrameNumber(), imageSource.getEndFrameNumber());
    // System.out.printf("Interlaced %d-%d\n", getStartFrameNumber(), getEndFrameNumber());
    return pixels;
  }

  @Override
  protected Object getRawFrame(int frame) {
    if (frame < 1) {
      return null;
    }

    // Check if the frame is allowed:
    // Start
    // |
    // |----|Block|Skip|Block|Skip|Block|Skip
    if (frame < start) {
      return null;
    }
    final int frameInBlock = (frame - start) % (size + skip);
    if (frameInBlock >= size) {
      return null;
    }
    final Object pixels = imageSource.getRaw(frame);
    // Set the frame to the last one read from the source
    setFrameNumber(imageSource.getStartFrameNumber(), imageSource.getEndFrameNumber());
    return pixels;
  }

  /**
   * Checks if is valid.
   *
   * @param frame the frame
   * @return true, if is valid
   */
  @Override
  public boolean isValid(int frame) {
    return imageSource.isValid(frame);
  }

  @Override
  public String toString() {
    return String.format("%s (Interlaced %d,%d,%d)", imageSource.toString(), start, size, skip);
  }

  /**
   * Gets the first frame that contains data.
   *
   * @return The first frame that contains data.
   */
  public int getStart() {
    return start;
  }

  /**
   * Gets the number of continuous frames containing data.
   *
   * @return The number of continuous frames containing data.
   */
  public int getSize() {
    return size;
  }

  /**
   * Gets the number of continuous frames to ignore before the next data.
   *
   * @return The number of continuous frames to ignore before the next data.
   */
  public int getSkip() {
    return skip;
  }

  @Override
  public void setReadHint(ReadHint readHint) {
    imageSource.setReadHint(readHint);
  }

  @Override
  public ReadHint getReadHint() {
    return imageSource.getReadHint();
  }
}
