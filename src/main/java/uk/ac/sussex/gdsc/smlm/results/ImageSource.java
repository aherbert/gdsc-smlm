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

package uk.ac.sussex.gdsc.smlm.results;

import com.thoughtworks.xstream.XStream;
import com.thoughtworks.xstream.XStreamException;
import com.thoughtworks.xstream.annotations.XStreamOmitField;
import com.thoughtworks.xstream.io.xml.DomDriver;
import java.awt.Rectangle;
import java.util.logging.Level;
import java.util.logging.Logger;
import uk.ac.sussex.gdsc.core.annotation.Nullable;
import uk.ac.sussex.gdsc.smlm.utils.ImageConverter;
import uk.ac.sussex.gdsc.smlm.utils.XStreamUtils;

/**
 * Abstract base class for the image source for peak results.
 */
public abstract class ImageSource {
  @XStreamOmitField
  private static final ImageConverter IMAGE_CONVERTER = new ImageConverter();

  private String name;
  @XStreamOmitField
  private int startFrame;
  @XStreamOmitField
  private int endFrame;

  /** The x origin. */
  protected int xOrigin;

  /** The y origin. */
  protected int yOrigin;

  /** The width. */
  protected int width;

  /** The height. */
  protected int height;

  /** The frames. */
  protected int frames;
  /**
   * The instance that does the conversion. This can be manipulated for different RGB colour
   * conversions.
   */
  @XStreamOmitField
  private ImageConverter imageConverter = IMAGE_CONVERTER;

  /**
   * The status flag for sequential reading using the {@link ImageSource#next()} methods.
   */
  public enum SequentialReadStatus {
    /** Sequential reading is closed (no more data). */
    CLOSED,
    /** The image is open for sequential reading but no data has been read. */
    READY,
    /** Sequential reading is running (more data is available). */
    RUNNING;
  }

  /** The sequential read status. */
  @XStreamOmitField
  protected SequentialReadStatus sequentialReadStatus;

  /**
   * The hints for how the source will be used. This allows the source to optimise the opening
   * process to prepare the image.
   */
  public enum ReadHint {
    /** The source will be used only in sequential mode. */
    SEQUENTIAL,
    /** The source will be used only in non-sequential mode. */
    NONSEQUENTIAL,
    /** The source will be used in sequential and non-sequential mode. */
    BOTH;
  }

  @XStreamOmitField
  private ReadHint readHint = ReadHint.SEQUENTIAL;

  /**
   * Create the image source.
   *
   * @param name The name of the image source
   */
  public ImageSource(String name) {
    setName(name);
  }

  /**
   * Opens the source.
   *
   * @return true, if successful
   */
  public boolean open() {
    startFrame = endFrame = 0;
    // In the event of creation by reflection this may be null and will be reset to the default.
    setImageConverter(imageConverter);
    if (openSource()) {
      // Special flag for sequential reading
      sequentialReadStatus = SequentialReadStatus.READY;
      return true;
    }
    return false;
  }

  /**
   * Opens the source.
   *
   * @return true, if successful
   */
  protected abstract boolean openSource();

  /**
   * Closes the source.
   */
  public void close() {
    sequentialReadStatus = SequentialReadStatus.CLOSED;
    closeSource();
  }

  /**
   * Closes the source.
   */
  protected abstract void closeSource();

  /**
   * Gets the x origin of the image frame. This may be non-zero to specify a crop of an image frame.
   *
   * <p>Note that the origin is ignored by the method {@link #next(Rectangle)} and
   * {@link #get(int, Rectangle)} as these use a rectangle relative to the image source origin.
   *
   * @return the x origin
   */
  public int getXOrigin() {
    return xOrigin;
  }

  /**
   * Gets the y origin of the image frame. This may be non-zero to specify a crop of an image frame.
   *
   * <p>Note that the origin is ignored by the method {@link #next(Rectangle)} and
   * {@link #get(int, Rectangle)} as these use a rectangle relative to the image source origin.
   *
   * @return the y origin
   */
  public int getYOrigin() {
    return yOrigin;
  }

  /**
   * Get the width of the image frame. The frame returned by {@link #next()} will be equal to width
   * * height.
   *
   * @return the width
   */
  public int getWidth() {
    return width;
  }

  /**
   * Get the height of the image frame. The frame returned by {@link #next()} will be equal to width
   * * height.
   *
   * @return the height
   */
  public int getHeight() {
    return height;
  }

  /**
   * Gets the bounds. This just uses the X and Y origin and the width and height.
   *
   * @return the bounds
   */
  public Rectangle getBounds() {
    return new Rectangle(getXOrigin(), getYOrigin(), getWidth(), getHeight());
  }

  /**
   * Get the number of frames that can be extracted from the image source with calls to
   * {@link #next()}.
   *
   * @return The total number of frames
   */
  public int getFrames() {
    return frames;
  }

  /**
   * Get the start frame number of the source returned by the last call to {@link #get(int)} or
   * {@link #next()}.
   *
   * @return The start frame number of the latest block of data
   */
  public int getStartFrameNumber() {
    return startFrame;
  }

  /**
   * Get the end frame number of the source returned by the last call to {@link #get(int)} or
   * {@link #next()}.
   *
   * <p>This may be larger than the result returned by {@link #getFrames()} if the ImageSource is
   * selecting a subset of the possible frames.
   *
   * @return The end frame number of the latest block of data
   */
  public int getEndFrameNumber() {
    return endFrame;
  }

  /**
   * Set the current frame number(s) of the source returned by the last call to {@link #get(int)} or
   * {@link #next()}.
   *
   * <p>This should be called by subclasses that perform more complex frame manipulation than just
   * getting a single frame.
   *
   * @param startFrame the start frame of the current block of data
   * @param endFrame the end frame of the current block of data
   */
  protected void setFrameNumber(int startFrame, int endFrame) {
    this.startFrame = startFrame;
    this.endFrame = endFrame;
  }

  /**
   * Get the next frame. Return null if the frame is not available and set the current frame to
   * zero. The data is is packed in yx order: index = y * width + x;
   *
   * <p>Provides serial access to the data after a successful call to {@link #open()}.
   *
   * @return the next frame (or null if at the end)
   */
  public @Nullable float[] next() {
    final Object pixels = nextRaw();
    if (pixels != null) {
      return imageConverter.getData(pixels, getWidth(), getHeight(), null, null);
    }
    return null;
  }

  /**
   * Get the next frame. Return null if the frame is not available and set the current frame to
   * zero. The data is is packed in yx order: index = y * width + x;
   *
   * <p>Provides serial access to the data after a successful call to {@link #open()}
   *
   * <p>Note: The bounds are relative to the image source origin so that bounds.x + bounds.width
   * must be less or equal to than {@link #getWidth()}, similarly for height.
   *
   * @param bounds The bounding limits of the frame to extract
   * @return the next frame (or null if at the end)
   * @throws IllegalArgumentException if the bounds do not fit in the image
   */
  public @Nullable float[] next(Rectangle bounds) {
    if (!checkBounds(bounds)) {
      bounds = null;
    }
    final Object pixels = nextRaw();
    if (pixels != null) {
      return imageConverter.getData(pixels, getWidth(), getHeight(), bounds, null);
    }
    return null;
  }

  /**
   * Get the next frame of raw pixels. Return null if the frame is not available and set the current
   * frame to zero. The data is is packed in yx order: index = y * width + x;
   *
   * <p>Provides serial access to the data after a successful call to {@link #open()}
   *
   * @return the next frame (or null if at the end)
   */
  public @Nullable Object nextRaw() {
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

    startFrame = endFrame = (startFrame + 1);
    final Object data = nextRawFrame();
    if (data == null) {
      startFrame = endFrame = 0;
      sequentialReadStatus = SequentialReadStatus.CLOSED;
    }
    return data;
  }

  /**
   * Initialise for sequential read.
   *
   * @return true, if successful
   */
  protected abstract boolean initialiseSequentialRead();

  /**
   * Get the next frame of raw pixels. The data is is packed in yx order: index = y * width + x;
   *
   * <p>Must be implemented by sub-classes.
   *
   * @return the next frame (or null if at the end)
   */
  protected abstract @Nullable Object nextRawFrame();

  /**
   * Get a specific frame from the results. Return null if the frame is not available and set the
   * current frame to zero.
   *
   * <p>Provides random access to the data after a successful call to {@link #open()}. This
   * operation may be significantly slower than using {@link #next()} to read all the data.
   *
   * @param frame the frame
   * @return the frame (or null)
   */
  public @Nullable float[] get(int frame) {
    return get(frame, null);
  }

  /**
   * Get a specific frame from the results. Return null if the frame is not available and set the
   * current frame to zero.
   *
   * <p>Provides random access to the data after a successful call to {@link #open()}. This
   * operation may be significantly slower than using {@link #next()} to read all the data.
   *
   * <p>Note: The bounds are relative to the image source origin so that bounds.x + bounds.width
   * must be less or equal to than {@link #getWidth()}, similarly for height.
   *
   * @param frame the frame
   * @param bounds The bounding limits of the frame to extract
   * @return the frame (or null)
   * @throws IllegalArgumentException if the bounds do not fit in the image
   */
  public @Nullable float[] get(int frame, Rectangle bounds) {
    if (!checkBounds(bounds)) {
      bounds = null;
    }
    startFrame = endFrame = frame;
    final Object pixels = getRawFrame(frame);
    if (pixels != null) {
      return imageConverter.getData(pixels, getWidth(), getHeight(), bounds, null);
    }
    startFrame = endFrame = 0;
    return null;
  }

  /**
   * Get a specific frame of raw pixels from the results. Return null if the frame is not available
   * and set the current frame to zero.
   *
   * <p>Provides random access to the data after a successful call to {@link #open()}. This
   * operation may be significantly slower than using {@link #next()} to read all the data.
   *
   * @param frame the frame
   * @return the frame (or null)
   */
  public @Nullable Object getRaw(int frame) {
    startFrame = endFrame = frame;
    final Object data = getRawFrame(frame);
    if (data == null) {
      startFrame = endFrame = 0;
    }
    return data;
  }

  /**
   * Get a specific frame of raw pixels from the results. Return null if the frame is not available.
   *
   * <p>Must be implemented by sub-classes.
   *
   * @param frame the frame
   * @return The frame data
   */
  protected abstract @Nullable Object getRawFrame(int frame);

  /**
   * Get the name of the results source.
   *
   * @return the name
   */
  public String getName() {
    return name;
  }

  /**
   * Get the parent source upon which this source is based. The default is to return null.
   *
   * @return The parent source
   */
  public ImageSource getParent() {
    return null;
  }

  /**
   * Get the original source for the data provided. The default is to return this object.
   *
   * @return The original source
   */
  public ImageSource getOriginal() {
    return this;
  }

  /**
   * Set the name of the results source.
   *
   * @param name the new name
   */
  public void setName(String name) {
    if (name != null && name.length() > 0) {
      this.name = name;
    } else {
      this.name = "";
    }
  }

  /**
   * Return true if the frame is within the limits of the image source.
   *
   * <p>Note that the {@link #get(int)} method may still return null. This method can be used to
   * determine if the {@link #get(int)} method has skipped data, e.g. if interlaced, or if the data
   * has actually ended.
   *
   * @param frame the frame
   * @return true if valid
   */
  public abstract boolean isValid(int frame);

  @Override
  public String toString() {
    // Over-ride this to produce a nicer output description of the results source
    return String.format("%s [%d,%d:%dx%dx%d]", name, xOrigin, yOrigin, getWidth(), getHeight(),
        frames);
  }

  /**
   * Serialise to XML.
   *
   * @return An XML representation of this object
   * @see #fromXml(String)
   */
  public String toXml() {
    final XStream xs = new XStream(new DomDriver());
    try {
      XStream.setupDefaultSecurity(xs); // to be removed after 1.5
      xs.allowTypesByWildcard(new String[] {"uk.ac.sussex.gdsc.smlm.**"});
      xs.autodetectAnnotations(true);
      return xs.toXML(this);
    } catch (final XStreamException ex) {
      Logger.getLogger(getClass().getName()).log(Level.SEVERE, "Failed to serialise to XML", ex);
    }
    return "";
  }

  /**
   * Create the image source from serialised XML.
   *
   * @param xml the xml
   * @return the image source
   * @see #toXml()
   */
  public static ImageSource fromXml(String xml) {
    final XStream xs = new XStream(new DomDriver());
    try {
      XStream.setupDefaultSecurity(xs); // to be removed after 1.5
      xs.allowTypesByWildcard(new String[] {"uk.ac.sussex.gdsc.smlm.**"});
      xs.autodetectAnnotations(true);
      // Support package gdsc.smlm renamed to uk.ac.sussex.gdsc.smlm
      return (ImageSource) xs.fromXML(XStreamUtils.updateGdscPackageName(xml));
    } catch (final Exception ex) {
      Logger.getLogger(ImageSource.class.getName()).log(Level.SEVERE,
          "Failed to deserialise from XML", ex);
    }
    return null;
  }

  /**
   * Check if the bounds fit inside the image.
   *
   * @param bounds the bounds
   * @return True if the bounds are not null and are within the image, false if null or the bounds
   *         fit the image exactly
   * @throws IllegalArgumentException if the bounds do not fit in the image
   */
  public boolean checkBounds(Rectangle bounds) {
    return checkBounds(getWidth(), getHeight(), bounds);
  }

  /**
   * Check if the bounds fit inside the image.
   *
   * @param width the width
   * @param height the height
   * @param bounds the bounds
   * @return True if the bounds are not null and are within the image, false if null or the bounds
   *         fit the image exactly
   * @throws IllegalArgumentException if the bounds do not fit in the image
   */
  public static boolean checkBounds(int width, int height, Rectangle bounds) {
    if (bounds != null) {
      final int maxx = bounds.x + bounds.width;
      final int maxy = bounds.y + bounds.height;
      if (bounds.x < 0 || maxx > width || bounds.y < 0 || maxy > height) {
        throw new IllegalArgumentException("The bounds do not fit within the image");
      }
      return bounds.x != 0 || bounds.y != 0 || bounds.width != width || bounds.height != height;
    }
    return false;
  }

  /**
   * Gets the image converter.
   *
   * @return the image converter
   */
  public ImageConverter getImageConverter() {
    return imageConverter;
  }

  /**
   * Sets the image converter. If null then the default converter will be used.
   *
   * @param imageConverter the new image converter
   */
  public void setImageConverter(ImageConverter imageConverter) {
    this.imageConverter = (imageConverter == null) ? IMAGE_CONVERTER : imageConverter;
  }

  /**
   * Gets the sequential read status.
   *
   * @return the sequential read status
   */
  public SequentialReadStatus getSequentialReadStatus() {
    return sequentialReadStatus;
  }

  /**
   * Gets the read hint for how the source will be used.
   *
   * @return the read hint
   */
  public ReadHint getReadHint() {
    return readHint;
  }

  /**
   * Sets the read hint for how the source will be used.
   *
   * @param readHint the new read hint
   */
  public void setReadHint(ReadHint readHint) {
    this.readHint = readHint;
  }
}
