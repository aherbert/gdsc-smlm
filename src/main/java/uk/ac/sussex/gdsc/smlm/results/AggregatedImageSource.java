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

import com.thoughtworks.xstream.annotations.XStreamOmitField;

import java.awt.Rectangle;
import java.util.Arrays;

/**
 * Wraps an image source and allows aggregation of consecutive frames.
 */
public class AggregatedImageSource extends ImageSource {
  private final int aggregate;
  private final ImageSource imageSource;

  // Used for frame-based read
  @XStreamOmitField
  private int lastFrame;
  @XStreamOmitField
  private int lastStartFrame;
  @XStreamOmitField
  private int lastEndFrame;
  @XStreamOmitField
  private float[] lastImage;

  /**
   * Create a new aggregated image source using the given image source.
   *
   * @param imageSource The image source to aggregate (must not be null)
   * @param aggregate The number of frames to aggregate (must be above 1)
   * @throws IllegalArgumentException If the image is null or {@code aggregate} is not above 1.
   */
  public AggregatedImageSource(ImageSource imageSource, int aggregate) {
    super("");
    if (imageSource == null) {
      throw new IllegalArgumentException("Image source must not be null");
    }
    if (aggregate < 2) {
      throw new IllegalArgumentException(
          "The number of frames to aggregate (" + aggregate + ") must be above 1");
    }
    setName("Aggregated (" + aggregate + ") " + imageSource.getName());
    this.imageSource = imageSource;
    this.aggregate = aggregate;
  }

  /** {@inheritDoc} */
  @Override
  public int getXOrigin() {
    return imageSource.getXOrigin();
  }

  /** {@inheritDoc} */
  @Override
  public int getYOrigin() {
    return imageSource.getYOrigin();
  }

  /** {@inheritDoc} */
  @Override
  public int getWidth() {
    return imageSource.getWidth();
  }

  /** {@inheritDoc} */
  @Override
  public int getHeight() {
    return imageSource.getHeight();
  }

  /** {@inheritDoc} */
  @Override
  public int getFrames() {
    return (int) Math.ceil((double) imageSource.getFrames() / aggregate);
  }

  /** {@inheritDoc} */
  @Override
  public ImageSource getParent() {
    return imageSource;
  }

  /** {@inheritDoc} */
  @Override
  public ImageSource getOriginal() {
    return imageSource.getOriginal();
  }

  /** {@inheritDoc} */
  @Override
  protected boolean openSource() {
    return imageSource.openSource();
  }

  /** {@inheritDoc} */
  @Override
  protected void closeSource() {
    imageSource.closeSource();
  }

  /** {@inheritDoc} */
  @Override
  protected boolean initialiseSequentialRead() {
    if (imageSource.initialiseSequentialRead()) {
      imageSource.sequentialReadStatus = SequentialReadStatus.RUNNING;
      return true;
    }
    imageSource.sequentialReadStatus = SequentialReadStatus.CLOSED;
    return false;
  }

  /**
   * Overridden to allow the crop to be done before aggregation.
   */
  @Override
  public float[] next(Rectangle bounds) {
    if (!checkBounds(bounds)) {
      bounds = null;
    }

    if (sequentialReadStatus == SequentialReadStatus.READY) {
      if (initialiseSequentialRead()) {
        sequentialReadStatus = SequentialReadStatus.RUNNING;
      } else {
        sequentialReadStatus = SequentialReadStatus.CLOSED;
      }
    }
    if (sequentialReadStatus != SequentialReadStatus.RUNNING) {
      return null;
    }

    // Aggregate consecutive frames
    float[] image = imageSource.next(bounds);
    if (image != null) {
      // Ensure the original image is not updated by creating a copy
      image = Arrays.copyOf(image, image.length);

      final int start = imageSource.getStartFrameNumber();
      int end = imageSource.getEndFrameNumber();
      for (int n = 1; n < aggregate; n++) {
        final float[] image2 = imageSource.next(bounds);
        if (image2 == null) {
          break;
        }
        end = imageSource.getEndFrameNumber();
        for (int i = 0; i < image.length; i++) {
          image[i] += image2[i];
        }
      }
      // Ensure that the frame number is recorded
      setFrameNumber(start, end);
      // System.out.printf("Aggregated %d-%d\n", start, end);
    } else {
      setFrameNumber(0, 0);
      sequentialReadStatus = SequentialReadStatus.CLOSED;
    }
    return image;
  }

  /**
   * Overridden to allow the crop to be done before aggregation.
   */
  @Override
  public float[] get(int frame, Rectangle bounds) {
    if (frame < 1) {
      return null;
    }

    if (!checkBounds(bounds)) {
      bounds = null;
    }

    // Calculate if the cache is invalid
    if (frame != lastFrame || lastImage == null) {
      // Try and get the desired frame
      float[] image = imageSource.get(frame, bounds);
      if (image == null) {
        return null;
      }
      lastStartFrame = imageSource.getStartFrameNumber();
      lastEndFrame = imageSource.getEndFrameNumber();

      // Ensure the original image is not updated by creating a copy
      image = Arrays.copyOf(image, image.length);

      // Go forwards until the desired number of frames have been collated
      int collated = 1;
      int nextFrame = frame;
      while (collated < aggregate && imageSource.isValid(++nextFrame)) {
        final float[] image2 = imageSource.get(nextFrame, bounds);
        if (image2 != null) {
          lastEndFrame = imageSource.getEndFrameNumber();
          for (int i = 0; i < image.length; i++) {
            image[i] += image2[i];
          }
          collated++;
        }
      }
      // Cache the image
      lastImage = image;
      lastFrame = frame;
    }
    // Ensure that the frame number is recorded
    setFrameNumber(lastStartFrame, lastEndFrame);
    return lastImage;
  }

  /** {@inheritDoc} */
  @Override
  protected Object nextRawFrame() {
    // Aggregate frames consecutive frames
    float[] image = imageSource.next();
    if (image != null) {
      // Ensure the original image is not updated by creating a copy
      image = Arrays.copyOf(image, image.length);

      final int start = imageSource.getStartFrameNumber();
      int end = imageSource.getEndFrameNumber();
      for (int n = 1; n < aggregate; n++) {
        final float[] image2 = imageSource.next();
        if (image2 == null) {
          break;
        }
        end = imageSource.getEndFrameNumber();
        for (int i = 0; i < image.length; i++) {
          image[i] += image2[i];
        }
      }
      // Ensure that the frame number is recorded
      setFrameNumber(start, end);
      // System.out.printf("Aggregated %d-%d\n", start, end);
    }
    return image;
  }

  /** {@inheritDoc} */
  @Override
  protected Object getRawFrame(int frame) {
    if (frame < 1) {
      return null;
    }

    // Calculate if the cache is invalid
    if (frame != lastFrame || lastImage == null) {
      // Try and get the desired frame
      float[] image = imageSource.get(frame);
      if (image == null) {
        return null;
      }
      lastStartFrame = imageSource.getStartFrameNumber();
      lastEndFrame = imageSource.getEndFrameNumber();

      // Ensure the original image is not updated by creating a copy
      image = Arrays.copyOf(image, image.length);

      // Go forwards until the desired number of frames have been collated
      int collated = 1;
      int nextFrame = frame;
      while (collated < aggregate && imageSource.isValid(++nextFrame)) {
        final float[] image2 = imageSource.get(nextFrame);
        if (image2 != null) {
          lastEndFrame = imageSource.getEndFrameNumber();
          for (int i = 0; i < image.length; i++) {
            image[i] += image2[i];
          }
          collated++;
        }
      }
      // Cache the image
      lastImage = image;
      lastFrame = frame;
    }
    // Ensure that the frame number is recorded
    setFrameNumber(lastStartFrame, lastEndFrame);
    return lastImage;
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

  /**
   * Gets number of frames to aggregate.
   *
   * @return The number of frames to aggregate.
   */
  public int getAggregate() {
    return aggregate;
  }

  /** {@inheritDoc} */
  @Override
  public String toString() {
    return String.format("%s (Aggregate %d images)", imageSource.toString(), aggregate);
  }

  /** {@inheritDoc} */
  @Override
  public void setReadHint(ReadHint readHint) {
    imageSource.setReadHint(readHint);
  }

  /** {@inheritDoc} */
  @Override
  public ReadHint getReadHint() {
    return imageSource.getReadHint();
  }
}
