/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

import com.thoughtworks.xstream.annotations.XStreamOmitField;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.io.FileInfo;
import ij.measure.Calibration;
import ij.process.ImageProcessor;
import java.awt.Rectangle;
import uk.ac.sussex.gdsc.smlm.results.ImageSource;

/**
 * Represent an ImageJ image as a results source. Supports all greyscale images. Only processes
 * channel 0 of 32-bit colour images.
 */
// This is allowed to support backwards serialisation compatibility
// CHECKSTYLE.OFF: AbbreviationAsWordInName
public class IJImageSource extends ImageSource {
  // CHECKSTYLE.ON: AbbreviationAsWordInName

  @XStreamOmitField
  private int slice;
  private int singleFrame;
  private int extraFrames;
  @XStreamOmitField
  private Object[] imageArray;
  @XStreamOmitField
  private ImageStack imageStack;
  private String path;

  /**
   * Create a new image source using the name of the image file.
   *
   * @param name the name
   */
  public IJImageSource(String name) {
    super(name);
  }

  /**
   * Create a new image source using the path to the image file.
   *
   * @param name the name
   * @param path the path
   */
  public IJImageSource(String name, String path) {
    super(name);
    this.path = path;
  }

  /**
   * Create a new image source from an ImagePlus.
   *
   * @param imp the imp
   */
  public IJImageSource(ImagePlus imp) {
    super(imp.getTitle());
    initialise(imp);
  }

  private boolean initialise(ImagePlus imp) {
    imageArray = null;
    imageStack = null;
    if (imp == null) {
      return false;
    }
    final ImageStack s = imp.getImageStack();
    if (s.isVirtual()) {
      // We must use the image stack to get the image data for virtual images
      imageStack = s;
    } else {
      // We can access the image array directly
      imageArray = s.getImageArray();
    }
    width = imp.getWidth();
    height = imp.getHeight();
    // Store the number of valid frames
    if (singleFrame > 0) {
      frames = 1 + this.extraFrames;
    } else {
      frames = imp.getStackSize();
    }
    slice = 0;
    final FileInfo info = imp.getOriginalFileInfo();
    if (info != null) {
      path = info.directory + info.fileName;
    }
    final int[] origin = getOrigin(imp);
    setOrigin(origin[0], origin[1]);
    return true;
  }

  /**
   * Gets the origin from the Image calibration.
   *
   * @param imp the image
   * @return the origin
   * @throws IllegalArgumentException If the origin is not in integer pixel units
   */
  public static int[] getOrigin(ImagePlus imp) {
    final Calibration cal = imp.getLocalCalibration();
    if (cal != null && (cal.xOrigin != 0 || cal.yOrigin != 0)) {
      // Origin must be in integer pixels
      final double ox = Math.round(cal.xOrigin);
      final double oy = Math.round(cal.yOrigin);
      if (ox != cal.xOrigin || oy != cal.yOrigin) {
        throw new IllegalArgumentException("Origin must be in integer pixels");
      }
      // ImageJ has a negative origin to indicate 0,0 of the image
      // is greater than the actual origin. The ImageSource convention is
      // to have the origin represent the shift of the 0,0 pixel from the origin.
      return new int[] {(int) -ox, (int) -oy};
    }
    return new int[2];
  }

  /**
   * Gets the bounds from the Image calibration.
   *
   * @param imp the image
   * @return the bounds
   * @throws IllegalArgumentException If the origin is not in integer pixel units
   */
  public static Rectangle getBounds(ImagePlus imp) {
    final int[] origin = getOrigin(imp);
    return new Rectangle(origin[0], origin[1], imp.getWidth(), imp.getHeight());
  }

  /**
   * Create a new image source from an ImageProcessor.
   *
   * @param name the name
   * @param ip the ip
   */
  public IJImageSource(String name, ImageProcessor ip) {
    super(name);
    imageArray = new Object[] {ip.getPixels()};
    width = ip.getWidth();
    height = ip.getHeight();
    frames = 1;
  }

  /**
   * Create a new image source from an ImagePlus. Specify a single frame for processing.
   *
   * @param imp the imp
   * @param frame the frame
   */
  public IJImageSource(ImagePlus imp, int frame) {
    this(imp, frame, 0);
  }

  /**
   * Create a new image source from an ImagePlus. Specify a start frame and any additional frames
   * (after the start) for processing.
   *
   * @param imp the imp
   * @param startFrame The first frame to process
   * @param extraFrames The number of additional frames to process
   */
  public IJImageSource(ImagePlus imp, int startFrame, int extraFrames) {
    super(imp.getTitle());
    // Ensure only a single frame is processed
    singleFrame = startFrame;
    this.extraFrames = Math.max(0, extraFrames);
    initialise(imp);
  }

  /**
   * Sets the origin. This should be used if the source was a crop from a image camera sensor.
   *
   * @param x the x
   * @param y the y
   */
  public void setOrigin(int x, int y) {
    xOrigin = x;
    yOrigin = y;
  }

  /**
   * Gets the path to the original image file.
   *
   * @return the path (or null)
   */
  public String getPath() {
    return path;
  }

  @Override
  protected boolean openSource() {
    if (nullImageArray() && imageStack == null) {
      // Try and find the image using the name or path
      ImagePlus imp = null;
      if (getName() != null) {
        imp = WindowManager.getImage(getName());
      }

      if ((imp == null) && (path != null && path.length() > 0)) {
        // Try and open the original image from file
        imp = IJ.openImage(path);
        if (imp == null) {
          // Some readers return null and display the image, e.g. BioFormats.
          // Add code to handle this.
        } else if (getName() != null) {
          // Ensure the image has the correct name
          imp.setTitle(getName());
        }
      }
      return initialise(imp);
    }
    return true;
  }

  @Override
  public void closeSource() {
    imageArray = null;
    imageStack = null;
  }

  /**
   * Check if the image array or any part of it is null.
   *
   * @return True if the image array or any part of it is null.
   */
  private boolean nullImageArray() {
    if (imageArray == null) {
      return true;
    }
    // Check the image array. This is set to null when an image is closed by ImageJ
    // allowing us to detect when the image is still available
    for (final Object o : imageArray) {
      if (o == null) {
        return true;
      }
    }
    return false;
  }

  @Override
  protected boolean initialiseSequentialRead() {
    slice = 0;
    return true;
  }

  @Override
  protected Object nextRawFrame() {
    ++slice;
    if (singleFrame > 0) {
      // Return frames from the starting frame until the extra frames limit is reached
      return (slice - 1 <= extraFrames) ? get(singleFrame + slice - 1) : null;
    }
    return get(slice);
  }

  @Override
  protected Object getRawFrame(int frame) {
    if (imageArray != null) {
      if (frame > 0 && frame <= imageArray.length) {
        return imageArray[frame - 1];
      }
    } else if (imageStack != null
        // This is a virtual stack so access the image processor through the virtual stack object
        && frame > 0 && frame <= imageStack.getSize()) {
      return imageStack.getPixels(frame);
    }
    return null;
  }

  @Override
  public boolean isValid(int frame) {
    if (singleFrame > 0) {
      return (frame >= singleFrame && frame <= singleFrame + extraFrames);
    }
    return frame > 0 && frame <= frames;
  }

  @Override
  public String toString() {
    String string = super.toString();
    if (path != null) {
      string += String.format(" (%s)", path);
    }
    return string;
  }
}
