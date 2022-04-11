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
import ij.ImagePlus;
import ij.io.FileInfo;
import ij.io.Opener;
import java.awt.Rectangle;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import org.apache.commons.lang3.concurrent.ConcurrentRuntimeException;
import org.apache.commons.lang3.exception.ExceptionUtils;
import uk.ac.sussex.gdsc.core.annotation.Nullable;
import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.ij.ImageJPluginLoggerHelper;
import uk.ac.sussex.gdsc.core.ij.SeriesOpener;
import uk.ac.sussex.gdsc.core.ij.io.ByteArraySeekableStream;
import uk.ac.sussex.gdsc.core.ij.io.ExtendedFileInfo;
import uk.ac.sussex.gdsc.core.ij.io.FastImageReader;
import uk.ac.sussex.gdsc.core.ij.io.FastTiffDecoder;
import uk.ac.sussex.gdsc.core.ij.io.FastTiffDecoder.IndexMap;
import uk.ac.sussex.gdsc.core.ij.io.FastTiffDecoder.NumberOfImages;
import uk.ac.sussex.gdsc.core.ij.io.FileSeekableStream;
import uk.ac.sussex.gdsc.core.ij.io.SeekableStream;
import uk.ac.sussex.gdsc.core.logging.NullTrackProgress;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.logging.TrackProgress;
import uk.ac.sussex.gdsc.core.utils.FileUtils;
import uk.ac.sussex.gdsc.core.utils.concurrent.CloseableBlockingQueue;
import uk.ac.sussex.gdsc.smlm.results.ImageSource;

/**
 * Represent a series of TIFF image files as a results source. Supports all greyscale images. Only
 * processes channel 0 of 32-bit colour images.
 *
 * <p>Assumes that the width,height,depth dimensions of each file are the same. The depth for the
 * last image can be less (i.e. the last of the series) but the {@link #getFrames()} method will
 * return an incorrect value.
 */
public class SeriesImageSource extends ImageSource {
  private static final String UNEXPECTED_INTERRUPTION = "Unexpected interruption";

  /** The buffer limit for reading TIFF images into memory. Default = 30MB */
  private int bufferLimit = 31457280;

  /** The buffer limit for sequential reading of TIFF images. Default = 50MB */
  private long sequentialReadBufferLimit = 52428800L;

  /** The list of image filenames. */
  List<String> images;
  /**
   * Flag indicating if the image series contains only TIFF images. No other formats are currently
   * supported.
   */
  public final boolean isTiffSeries;

  /**
   * Contains the cumulative size of the TIFF image series. Used for a binary search to find the
   * image for the frame.
   */
  @XStreamOmitField
  int[] imageSize;

  /**
   * Used to cache data about the TIFF images.
   */
  @XStreamOmitField
  ImageData[] imageData;

  // Used for frame-based read
  @XStreamOmitField
  private Image lastImage;
  @XStreamOmitField
  private int lastImageId;

  private int numberOfThreads = 1;
  private int numberOfImages = 1;

  /** The track progress. */
  @XStreamOmitField
  TrackProgress trackProgress = NullTrackProgress.getInstance();

  // Used to process the files into images
  @XStreamOmitField
  private List<BaseWorker> workers;
  @XStreamOmitField
  private List<Thread> threads;

  /** The queue for in-memory buffered images awaiting TIFF decoding. */
  @XStreamOmitField
  CloseableBlockingQueue<NextSource> decodeQueue;
  /** The queue for in-memory buffered images awaiting TIFF reading. */
  @XStreamOmitField
  CloseableBlockingQueue<NextSource> readQueue;

  /** Used for sequential reading to queue the raw frames. */
  @XStreamOmitField
  CloseableBlockingQueue<Object> rawFramesQueue;

  /**
   * The first error that occurred during sequential reading. This is stored, the queue is shutdown
   * and then this can be thrown in the next() method.
   */
  @XStreamOmitField
  RuntimeException error;

  /**
   * Store details for an image.
   */
  private static abstract class Image {
    int width;
    int height;
    int size;

    Image(int width, int height, int size) {
      this.width = width;
      this.height = height;
      this.size = size;
    }

    /**
     * Gets the width.
     *
     * @return the width
     */
    int getWidth() {
      return width;
    }

    /**
     * Gets the height.
     *
     * @return the height
     */
    int getHeight() {
      return height;
    }

    /**
     * Gets the size.
     *
     * @return the size
     */
    int getSize() {
      return size;
    }

    /**
     * Gets the frame.
     *
     * @param index the frame number
     * @return the frame
     * @throws IOException Signals that an I/O exception has occurred.
     */
    abstract Object getFrame(int index) throws IOException;

    /**
     * Close the image.
     *
     * @param freeMemory Set to true to free any cached memory
     */
    abstract void close(boolean freeMemory);
  }

  /**
   * Support all images ImageJ can read using a fixed image array.
   */
  @SuppressWarnings("unused")
  private static class ArrayImage extends Image {
    final Object[] imageArray;

    ArrayImage() {
      super(0, 0, 0);
      imageArray = null;
    }

    ArrayImage(ImagePlus imp) {
      super(imp.getWidth(), imp.getHeight(), imp.getStackSize());
      this.imageArray = imp.getImageStack().getImageArray();
    }

    @Override
    Object getFrame(int index) {
      return imageArray[index];
    }

    @Override
    void close(boolean freeMemory) {
      // Do nothing
    }
  }

  /**
   * Special class adapted from ij.io.Opener.openTiffStack(...) and ij.io.FileOpener to handle large
   * contiguous TIFF images.
   *
   * <p>The methods that manipulate the internal state are synchronized to prevent multiple threads
   * causing read errors.
   */
  private static class TiffImage extends Image {
    /** The index map for OME-TIFF images, otherwise null. */
    final IndexMap indexMap;
    /** Cache of the file info for each image in the TIFF. */
    final ExtendedFileInfo[] info;
    /** The file info used by the FastImageReader. */
    final ExtendedFileInfo fi;
    final long bytesPerFrame;
    /** Flag indicating that the pixels are evenly spaced through the data from the first offset. */
    boolean contiguous;
    /** The image reader. */
    FastImageReader reader;
    /**
     * A reference to the primary seekable stream used to read the image pixels. This may be
     * buffered in memory as a ByteArraySeekableStream. If so it can be reused by another thread via
     * a copy.
     */
    SeekableStream ss;
    /**
     * A reference to a file seekable stream. This is used when reading metadata from a file based
     * image (i.e. the primary seekable stream is not held in memory).
     */
    FileSeekableStream fss;
    /**
     * The TIFF decoder. This holds a seekable stream but cannot be closed. The resource is
     * maintained separately.
     */
    FastTiffDecoder td;
    /** Flag indicating that the primary seekable stream is a ByteArraySeekableStream. */
    boolean inMemory;
    /** Flag indicating that no errors reading the image have occurred. */
    boolean canRead = true;
    /** The number of frames that have been read from the input stream. */
    int frameCount;
    /** Flag indicating that the Tiff info is complete. Relevant when opened using an index map */
    boolean complete;

    TiffImage(ExtendedFileInfo[] info, SeekableStream ss) {
      super(info[0].width, info[0].height, 0);
      indexMap = null;
      this.info = info;
      fi = info[0];

      // Only support certain types
      if (isSupported(fi.fileType)) {
        // We use the opened SeekableStream.
        // This may be in-memory data or may be from a random access file.
        if (ss != null) {
          try {
            ss.seek(0);
            this.ss = ss;

            // Store if the stream contains in-memory data
            this.inMemory = ss instanceof ByteArraySeekableStream;
          } catch (final IOException ex) {
            canRead = false;
            bytesPerFrame = 0;
            return;
          }
        }

        // Determine the number of images
        if (info.length > 1) {
          if (allSameSizeAndType(info)) {
            size = info.length;
          }
        } else {
          size = fi.nImages;
          contiguous = true;
        }
        if (size == 0) {
          canRead = false;
          bytesPerFrame = 0;
        } else {
          bytesPerFrame = getBytesPerFrame(fi.fileType);
          reader = new FastImageReader(fi);
        }
      } else {
        canRead = false;
        bytesPerFrame = 0;
      }
    }

    TiffImage() {
      super(0, 0, 0);
      indexMap = null;
      info = null;
      fi = null;
      bytesPerFrame = 0;
      canRead = false;
    }

    TiffImage(IndexMap indexMap, ExtendedFileInfo fi, SeekableStream ss) {
      super(fi.width, fi.height, 0);
      this.indexMap = indexMap;
      this.info = new ExtendedFileInfo[indexMap.getSize()];
      info[0] = fi;
      this.fi = fi;

      // Only support certain types
      if (isSupported(fi.fileType)) {
        // We use the opened SeekableStream.
        // This may be in-memory data or may be from a random access file.
        if (ss != null) {
          try {
            ss.seek(0);
            this.ss = ss;

            // Store if the stream contains in-memory data
            this.inMemory = ss instanceof ByteArraySeekableStream;
          } catch (final IOException ex) {
            canRead = false;
            bytesPerFrame = 0;
            return;
          }
        }

        // Determine number of images
        if (indexMap.getSize() > 1) {
          if (singlePlane()) {
            size = indexMap.getSize();
          }
        } else {
          size = indexMap.getSize();
          contiguous = true;
        }
        if (size == 0) {
          bytesPerFrame = getBytesPerFrame(fi.fileType);
          reader = new FastImageReader(fi);
        } else {
          canRead = false;
          bytesPerFrame = 0;
        }
      } else {
        canRead = false;
        bytesPerFrame = 0;
      }
    }

    private static boolean isSupported(int fileType) {
      switch (fileType) {
        // Greyscale images as we just want the raw pixels with no color model
        case FileInfo.GRAY8:
        case FileInfo.GRAY16_SIGNED:
        case FileInfo.GRAY16_UNSIGNED:
        case FileInfo.GRAY12_UNSIGNED:
        case FileInfo.GRAY32_INT:
        case FileInfo.GRAY32_UNSIGNED:
        case FileInfo.GRAY32_FLOAT:
        case FileInfo.GRAY24_UNSIGNED:
        case FileInfo.GRAY64_FLOAT:
          return true;
        default:
          return false;
      }
    }

    private int getBytesPerFrame(int fileType) {
      switch (fileType) {
        case FileInfo.GRAY8:
          return width * height;
        case FileInfo.GRAY16_SIGNED:
        case FileInfo.GRAY16_UNSIGNED:
          return 2 * width * height;
        case FileInfo.GRAY32_INT:
        case FileInfo.GRAY32_UNSIGNED:
        case FileInfo.GRAY32_FLOAT:
          return 4 * width * height;
        case FileInfo.GRAY64_FLOAT:
          return 8 * width * height;
        case FileInfo.GRAY24_UNSIGNED:
          return 3 * width * height;
        case FileInfo.GRAY12_UNSIGNED:
          return (int) (width * height * 1.5);
        default:
          return 0;
      }
    }

    /**
     * Check the file info are all the same size and type and have just 1 image.
     *
     * @param info the info
     * @return True if all the same size and type
     */
    final boolean allSameSizeAndType(ExtendedFileInfo[] info) {
      boolean ok = info[0].nImages == 1;
      contiguous = ok;
      final long startingOffset = info[0].getOffset();
      final long size = (long) info[0].width * info[0].height * info[0].getBytesPerPixel();
      for (int i = 1; i < info.length && ok; i++) {
        if (info[i].fileType != info[0].fileType || info[i].width != info[0].width
            || info[i].height != info[0].height || info[i].nImages != 1) {
          ok = false;
        }
        if (info[i].getOffset() != startingOffset + i * size) {
          contiguous = false;
        }
      }
      // We can read using only the first ExtendedFileInfo object if this is contiguous
      if (ok && contiguous) {
        info[0].gapBetweenImages = 0;
      }
      return ok;
    }

    /**
     * Check the file info is the same size and type as the first file info.
     *
     * @param info the info
     * @return True if the same size and type
     */
    boolean sameSizeAndType(ExtendedFileInfo info) {
      return info.fileType == fi.fileType && info.width == fi.width && info.height == fi.height
          && info.nImages == 1;
    }

    /**
     * Check if the index map has only a single plane. This is all we support at the moment.
     *
     * @return true, if successful
     */
    private boolean singlePlane() {
      if (!indexMap.isSingleChannel()) {
        return false;
      }

      // Assume this is stage positions
      if (!indexMap.isSinglePosition()) {
        return false;
      }

      // Only 1 of the following can be above 1.
      // If either is single then this is OK
      return indexMap.isSingleFrame() || indexMap.isSingleSlice();
    }

    /**
     * Gets the next frame. This is only supported if all IFDs have been read.
     *
     * @return the object
     * @throws IOException Signals that an I/O exception has occurred.
     */
    Object nextFrame() throws IOException {
      // Skip ahead
      long skip;

      if (contiguous) {
        // If sequential reading just skip the gap between frames
        if (frameCount == 0) {
          // The first frame we know the exact offset
          skip = fi.getOffset();
        } else {
          // If sequential reading just skip the gap between frames
          skip = fi.gapBetweenImages;
        }
      } else {
        // Adapted from ij.io.Opener.openTiffStack(...)

        // Each image offset is described by a separate ExtendedFileInfo object
        final ExtendedFileInfo fileInfo = getInfo(frameCount);
        skip = fileInfo.getOffset();
        this.fi.stripOffsets = fileInfo.stripOffsets;
        this.fi.stripLengths = fileInfo.stripLengths;

        if (frameCount != 0) {
          // We must subtract the current file location.
          // Assume the previous one is there as this is a sequential read
          skip -= (info[frameCount - 1].getOffset() + bytesPerFrame);
          if (skip < 0L) {
            canRead = false;
            throw new IllegalStateException("Bad TIFF offset " + skip);
          }
        }
      }

      // Store the number of frames that have been read
      frameCount++;

      FileUtils.skip(ss, skip);
      return readPixels();
    }

    private Object readPixels() throws IOException {
      if (inMemory) {
        return reader.readPixels((ByteArraySeekableStream) ss, 0);
      }
      return reader.readPixels(ss, 0);
    }

    @Override
    synchronized Object getFrame(int index) throws IOException {
      if (!openInputStream()) {
        throw new IllegalStateException("Cannot read the TIFF image");
      }

      // Skip ahead
      long offset;

      if (contiguous) {
        // Read using the first ExtendedFileInfo object

        // The first frame we know the exact offset
        offset = fi.getOffset();

        if (index != 0) {
          // Skip ahead
          offset += (bytesPerFrame + fi.gapBetweenImages) * index;
        }
      } else {
        // Adapted from ij.io.Opener.openTiffStack(...)

        // Each image offset is described by a separate ExtendedFileInfo object
        // We may have to read it first.
        offset = getInfo(index).getOffset();
        fi.stripOffsets = info[index].stripOffsets;
        fi.stripLengths = info[index].stripLengths;
      }

      // Store the number of frames that have been read
      frameCount = index + 1;

      ss.seek(offset);

      // Note:
      // The reader now throws exceptions.
      // If this occurs we could set the canRead flag to prevent further IO.
      // However do not do this so that an image can be scrolled up to the point at which
      // it cannot be read (e.g. in the case of TIFF file truncation during sequential writing).

      return readPixels();
    }

    /**
     * Gets the info of the specified index. Uses the cached FileInfo using the IndexMap.
     *
     * @param index the index
     * @return the info
     * @throws IOException Signals that an I/O exception has occurred.
     */
    private ExtendedFileInfo getInfo(int index) throws IOException {
      if (info[index] != null) {
        return info[index];
      }
      // We have to read the TIFF info
      if (td == null) {
        // We can clone the in-memory seekable stream
        if (inMemory) {
          td = FastTiffDecoder.create(((ByteArraySeekableStream) ss).copy(), fi.fileName);
        } else {
          // Open the file.
          // Potential resource leak so this second seekable stream is maintained for closing.
          fss = new FileSeekableStream(getFile());
          td = FastTiffDecoder.create(fss, fi.fileName);
        }
      }
      // If this throws a NullPointerException from the indexMap
      // then it will be handled in getFrame().
      info[index] = td.getTiffInfo(indexMap, index, true);
      if (!sameSizeAndType(info[index])) {
        throw new IOException("Not same size and type");
      }
      return info[index];
    }

    void readTiffInfo(ByteArraySeekableStream ss) throws IOException {
      if (isCompleteTiffInfo()) {
        return;
      }
      if (td == null) {
        td = FastTiffDecoder.create(ss.copy(), fi.fileName);
      }

      for (int i = 0; i < info.length; i++) {
        if (info[i] == null) {
          info[i] = td.getTiffInfo(indexMap, i, true);
          if (!sameSizeAndType(info[i])) {
            throw new IOException("Not same size and type");
          }
        }
      }
    }

    boolean isCompleteTiffInfo() {
      if (!complete) {
        for (final ExtendedFileInfo fi : info) {
          if (fi == null) {
            return false;
          }
        }
      }
      complete = true;
      return true;
    }

    boolean openInputStream() throws IOException {
      if (canRead && ss == null) {
        if (inMemory) {
          throw new IllegalStateException("No input file required for in-memory image");
        }

        try {
          final File f = getFile();
          if (validateFileInfo(f, fi)) {
            ss = new FileSeekableStream(f);
          }
        } finally {
          if (ss == null) {
            canRead = false;
          }
        }
      }
      return canRead;
    }

    File getFile() {
      if (fi.directory.length() > 0) {
        return new File(fi.directory, fi.fileName);
      }
      return new File(fi.fileName);
    }

    /**
     * Adapted from ij.io.FileOpener.validateFileInfo
     *
     * @param file the file
     * @param fi the file info
     * @return true, if successful
     */
    boolean validateFileInfo(File file, ExtendedFileInfo fi) {
      if (fi.width <= 0 || fi.height <= 0) {
        return false;
      }
      final long offset = fi.getOffset();
      if (offset >= 0 && offset < 1000L) {
        return true;
      }
      if (offset < 0L) {
        return false;
      }
      if (fi.fileType == FileInfo.BITMAP || fi.compression != FileInfo.COMPRESSION_NONE) {
        return true;
      }
      long size = (long) fi.width * fi.height * fi.getBytesPerPixel();
      size = fi.nImages > 1 ? size : size / 4;
      if (fi.height == 1) {
        size = 0; // allows plugins to read info of unknown length at end of file
      }
      return (offset + size <= file.length());
    }

    /**
     * Reset for sequential reading.
     *
     * @throws IOException Signals that an I/O exception has occurred.
     */
    void reset() throws IOException {
      frameCount = 0;
      if (ss == null) {
        openInputStream();
      } else {
        ss.seek(0);
      }
    }

    @Override
    public synchronized void close(boolean freeResources) {
      frameCount = 0;

      // This is done when sequentially reading so we clear the memory
      inMemory = false;

      // The dedicated tiff decoder is used for reading IFDs.
      // Ensure it is re-initialised.
      td = null;

      // Close a second seekable stream if present
      closeQuietly(fss);
      fss = null;

      // Close the primary seekable stream if present
      closeQuietly(ss);
      ss = null;
    }
  }

  /**
   * Close the stream quietly.
   *
   * @param ss the stream
   */
  static void closeQuietly(SeekableStream ss) {
    if (ss != null) {
      try {
        ss.close();
      } catch (final IOException ignored) {
        // Ignore
      }
    }
  }

  /**
   * Base class for workers.
   */
  private abstract static class BaseWorker implements Runnable {
    volatile boolean run = true;

    /**
     * Set the run flag to false.
     */
    void shutdown() {
      run = false;
    }
  }

  /**
   * Add source images to the queue to be read.
   */
  private class TiffWorker extends BaseWorker {
    @Override
    public void run() {
      FileSeekableStream ss = null;

      // All exceptions are caught so that the queue can be shutdown correctly
      try {
        for (int currentImage = 0; run && currentImage < images.size(); currentImage++) {
          final String path = images.get(currentImage);

          // Check the cache to avoid a re-read.
          TiffImage image = imageData[currentImage].tiffImage;
          if (image == null) {
            trackProgress.log("Reading TIFF info %s", path);

            // This is closed in a finally block if an error occurs
            // or by image.close(true) if the image was read.
            ss = createSeekableStream(path);

            // Reading the info can be done by the TiffImage dynamically
            // using a dedicated Tiff decoder.
            // Dynamically reading will not experience a big pause when
            // all IFDs are read on a large image (e.g. 5000 frames). This will
            // allow the source to continue queueing frames for processing
            // even on massive images.

            final FastTiffDecoder td = FastTiffDecoder.create(ss, path);
            td.setTrackProgress(trackProgress);
            final IndexMap indexMap = td.getIndexMap();

            // Check the image map is the correct size (only if we have sizes)
            if (indexMap != null
                && (imageSize == null || indexMap.getSize() == getImageSize(currentImage))) {
              // We need the first IFD to define the image pixel type and width/height
              final ExtendedFileInfo fi = td.getTiffInfo(indexMap, 0, true);
              // Do not pass the stream. A new one will be opened by the TiffImage
              // to read pixel data.
              image = new TiffImage(indexMap, fi, null);
              // Store the opened stream and the decoder.
              // These will be reused to read FileInfo for subsequent frames.
              image.fss = ss;
              image.td = td;
            } else {
              // We need to read all the IFDs

              // Reset after reading the index map
              td.reset();
              // This will throw if there is an error and the seekable stream
              // is closed in the finally block.
              final ExtendedFileInfo[] info = td.getTiffInfo(true);
              if (info == null) {
                setError(new DataException("No TIFF file info"));
                break;
              }

              // Use the same seekable stream to read pixel data.
              // The FileInfo should be complete and
              // the IFDs will not need to be read again.
              image = new TiffImage(info, ss);
            }

            storeTiffImage(currentImage, image);
          }

          // Check dimensions
          // Check it is the expected size
          if (image.getWidth() != getWidth() || image.getHeight() != getHeight()) {
            setError(new DataException("Dimension mismatch"));
            break;
          }
          if (imageSize != null && image.getSize() != getImageSize(currentImage)) {
            setError(new DataException("Unexpected image size"));
            break;
          }

          trackProgress.log("Reading TIFF %s", path);

          image.reset();

          try {
            // Read all the frames sequentially
            for (int i = 0; i < image.size; i++) {
              // This should throw if no pixels can be read
              final Object pixels = image.nextFrame();

              // This will block until the queue has capacity or is closed
              if (!rawFramesQueue.put(pixels)) {
                break;
              }
            }
          } finally {
            // Close the image and free resources
            image.close(true);
          }
        }
      } catch (final DataException ex) {
        setError(ex);
      } catch (final Exception ex) {
        setError(new DataException(ex));
      } finally {
        closeInputStream(ss);
      }

      // Close and clear if an error occurred
      rawFramesQueue.close((error != null));

      shutdown();
    }
  }

  /**
   * An in-memory source of image bytes.
   */
  private static class NextSource {
    final byte[] buffer;
    final int imageId;
    // ExtendedFileInfo[] info;
    Image image;

    NextSource(byte[] buffer, int imageId) {
      this.buffer = buffer;
      this.imageId = imageId;
    }

    // void setFileInfo(ExtendedFileInfo[] info)
    // {
    // this.info = info;
    // }

    void setImage(Image image) {
      this.image = image;
    }

    // NextSource()
    // {
    // this(null, -1);
    // }
  }

  /**
   * Read source image files into memory.
   */
  private class BufferWorker extends BaseWorker {
    @Override
    public void run() {
      FileSeekableStream fs = null;

      // All exceptions are caught so that the queue can be shutdown correctly
      try {
        for (int currentImage = 0; run && currentImage < images.size(); currentImage++) {
          final String path = images.get(currentImage);

          trackProgress.log("Reading TIFF into memory %s", path);

          // long start = System.nanoTime();
          final File file = new File(path);
          fs = new FileSeekableStream(file);

          // We already know they fit into memory (i.e. the size is not zero)
          final int size = (int) imageData[currentImage].fileSize;
          final byte[] buf = new byte[size];
          if (fs.readBytes(buf) != size) {
            setError(new DataException("Cannot buffer file into memory"));
            break;
          }
          try {
            fs.close();
          } finally {
            // Prevent another close attempt
            fs = null;
          }
          // System.out.printf("Read %d\n", (System.nanoTime() - start) / 1000000);

          // This may be closed upon error
          if (!decodeQueue.put(new NextSource(buf, currentImage))) {
            break;
          }
        }
      } catch (final DataException ex) {
        setError(ex);
      } catch (final Exception ex) {
        setError(new DataException(ex));
      } finally {
        closeInputStream(fs);

        if (error == null) {
          decodeQueue.close(false);
        } else {
          closeWorkflowQueues();
        }

        shutdown();
      }
    }
  }

  /**
   * Decode the TIFF IFDs. All IFDs are decoded on this thread so that the work of reading the TIFF
   * on the next thread is faster (i.e. avoiding dynamic IFD reading when reading the pixels).
   */
  private class DecodeWorker extends BaseWorker {
    @Override
    public void run() {
      try {
        while (run) {
          final NextSource nextSource = decodeQueue.take();
          if (nextSource == null || !run) {
            break;
          }
          final int currentImage = nextSource.imageId;
          if (currentImage == -1) {
            break;
          }

          @SuppressWarnings("resource")
          final ByteArraySeekableStream ss = ByteArraySeekableStream.wrap(nextSource.buffer);

          // Re-use the cache
          TiffImage image = imageData[currentImage].tiffImage;
          if (image == null) {
            final String path = images.get(currentImage);
            trackProgress.log("Reading TIFF info %s", path);

            final FastTiffDecoder td = FastTiffDecoder.create(ss, path);
            td.setTrackProgress(trackProgress);

            final ExtendedFileInfo[] info = td.getTiffInfo(true);

            if (info == null) {
              setError(new DataException("No TIFF file info"));
              break;
            }

            // Create as in-memory
            image = new TiffImage(info, ss);

            // Check dimensions
            // Check it is the expected size
            if (image.getWidth() != getWidth() || image.getHeight() != getHeight()) {
              setError(new DataException("Dimension mismatch"));
              break;
            }
            if (imageSize != null && image.getSize() != getImageSize(currentImage)) {
              setError(new DataException("Unexpected image size"));
              break;
            }

            storeTiffImage(currentImage, image);
          } else if (image.indexMap != null) {
            // Also re-read if the image was opened using an index map as
            // the tiff info may be incomplete.
            // This will always be the case
            // for the first image as that is opened in openSource().
            final String path = images.get(currentImage);
            trackProgress.log("Reading TIFF info %s", path);
            image.readTiffInfo(ss);
          } else {
            // Update to be in-memory. This will be used when reading the TIFF.
            // It will subsequently be closed to free-memory.
            image.close(true);
            image.inMemory = true;
            image.ss = ss;
            image.td = null;
          }

          // nextSource.setFileInfo(info);
          nextSource.setImage(image);

          // This may be closed upon error
          if (!readQueue.put(nextSource)) {
            break;
          }
        }
        // Suppress all exceptions. These are stored to be thrown by the nextRawFrame() method.
      } catch (final DataException ex) {
        setError(ex);
      } catch (final IOException ex) {
        setError(new UncheckedIOException(ex));
      } catch (final InterruptedException ex) {
        Thread.currentThread().interrupt();
        setError(new ConcurrentRuntimeException(UNEXPECTED_INTERRUPTION, ex));
      } catch (final Exception ex) {
        setError(new DataException(ex));
      } finally {
        if (error == null) {
          readQueue.close(false);
        } else {
          closeWorkflowQueues();
        }

        shutdown();
      }
    }
  }

  /**
   * Read source image files into memory.
   */
  private class ReadWorker extends BaseWorker {
    @Override
    public void run() {
      try {
        while (run) {
          final NextSource nextSource = readQueue.take();
          if (nextSource == null || !run) {
            break;
          }
          final int currentImage = nextSource.imageId;
          if (currentImage == -1) {
            break;
          }

          // It is assumed we are reading a TIFF series
          final TiffImage image = (TiffImage) nextSource.image;

          image.reset();

          try {
            // Read all the frames sequentially
            for (int i = 0; i < image.size; i++) {
              // This should throw if no pixels can be read
              final Object pixels = image.nextFrame();

              // This may be closed upon error
              if (!rawFramesQueue.put(pixels)) {
                break;
              }
            }
          } finally {
            // Close the image
            image.close(true);
          }
        }
      } catch (final DataException ex) {
        setError(ex);
      } catch (final IOException ex) {
        setError(new UncheckedIOException(ex));
      } catch (final InterruptedException ex) {
        Thread.currentThread().interrupt();
        setError(new ConcurrentRuntimeException(UNEXPECTED_INTERRUPTION, ex));
      } catch (final Exception ex) {
        setError(new DataException(ex));
      } finally {
        if (error == null) {
          rawFramesQueue.close(false);
        } else {
          closeWorkflowQueues();
        }

        shutdown();
      }
    }
  }

  /**
   * Hold image data.
   */
  private static class ImageData {
    long fileSize;
    TiffImage tiffImage;

    ImageData(long size) {
      fileSize = size;
    }
  }

  /**
   * Creates the seekable stream.
   *
   * @param path the path
   * @param size the size of the file
   * @return the seekable stream
   */
  private SeekableStream createSeekableStream(String path, long size) throws FileNotFoundException {
    final File file = new File(path);

    // Don't buffer massive images into memory
    if (belowBufferLimit(size)) {
      final SeekableStream ss = readByteArraySeekableStream(file, size);
      if (ss != null) {
        return ss;
      }
    }

    return createSeekableStream(path);
  }

  /**
   * Creates the seekable stream.
   *
   * @param path the path
   * @return the seekable stream
   * @throws FileNotFoundException the file not found exception
   */
  static FileSeekableStream createSeekableStream(String path) throws FileNotFoundException {
    return new FileSeekableStream(path);
  }

  /**
   * Close the input stream.
   *
   * @param is the input stream
   */
  static void closeInputStream(InputStream is) {
    if (is != null) {
      try {
        is.close();
      } catch (final IOException ignored) {
        // Ignore
      }
    }
  }

  private boolean belowBufferLimit(long size) {
    return (size > 0 && size <= bufferLimit);
  }

  private boolean belowBufferLimit() {
    for (final ImageData d : imageData) {
      if (!belowBufferLimit(d.fileSize)) {
        return false;
      }
    }
    return true;
  }

  private static long getSize(File file) {
    try {
      return file.length();
    } catch (final SecurityException ignored) {
      // No file size...
    }
    return 0;
  }

  private static SeekableStream readByteArraySeekableStream(File file, long size) {
    try (FileSeekableStream fs = new FileSeekableStream(file)) {
      final byte[] buf = new byte[(int) size];
      if (fs.readBytes(buf) == size) {
        return ByteArraySeekableStream.wrap(buf);
      }
    } catch (final IOException ex) {
      // This exception is not bubbled up.
      // At the moment if we ignore this then the ImageWorker will open the file
      // rather than process from the memory stream.
      ImageJPluginLoggerHelper.getLogger(SeriesImageSource.class.getName()).log(Level.WARNING,
          "Failed to read file into memory", ex);
    }
    return null;
  }

  /**
   * Create a new image source using the given image series.
   *
   * @param name the name
   * @param series the series
   */
  public SeriesImageSource(String name, SeriesOpener series) {
    super(name);
    if (series == null) {
      isTiffSeries = false;
    } else {
      final String[] names = series.getImageList();
      for (int i = 0; i < names.length; i++) {
        names[i] = new File(series.getPath(), names[i]).getPath();
      }
      isTiffSeries = isTiffSeries(names);
    }
  }

  /**
   * Create a new image source using the directory containing the images.
   *
   * <p>The directory is opened using {@link uk.ac.sussex.gdsc.core.ij.SeriesOpener }
   *
   * @param name the name
   * @param directory the directory
   */
  public SeriesImageSource(String name, String directory) {
    this(name, new SeriesOpener(directory));
  }

  /**
   * Create a new image source using the paths to the images.
   *
   * @param name the name
   * @param filenames the full path names of the image files
   */
  public SeriesImageSource(String name, String[] filenames) {
    super(name);
    isTiffSeries = isTiffSeries(filenames);
  }

  private boolean isTiffSeries(String[] filenames) {
    // Create this as it is needed for XStream serialisation
    images = new ArrayList<>();

    for (final String filename : filenames) {
      if (getFileType(filename) != Opener.TIFF) {
        // Only support TIFF images
        return false;
      }
    }

    // All images are TIFF so store the filenames
    images.addAll(Arrays.asList(filenames));

    return true;
  }

  /**
   * Attempts to determine the image file type by looking for 'magic numbers' and the file name
   * extension.
   *
   * <p>Copied from ij.io.Opener and removed all but the TIFF identification.
   */
  private static int getFileType(String filename) {
    final byte[] buf = new byte[132];
    int read;
    try (InputStream is = Files.newInputStream(Paths.get(filename))) {
      read = is.read(buf, 0, 132);
      if (read < 4) {
        return Opener.UNKNOWN;
      }
    } catch (final IOException ex) {
      return Opener.UNKNOWN;
    }

    final int b0 = buf[0] & 255;
    final int b1 = buf[1] & 255;
    final int b2 = buf[2] & 255;
    final int b3 = buf[3] & 255;

    // First check it is a possible TIFF:
    // Little-endian = II + 42 magic number
    // Big-endian = MM + 42 magic number
    final boolean littleEndian = (b0 == 73 && b1 == 73 && b2 == 42 && b3 == 0);
    final boolean bigEndian = (b0 == 77 && b1 == 77 && b2 == 0 && b3 == 42);
    if (!(littleEndian || bigEndian)) {
      return Opener.UNKNOWN;
    }

    // Don't check for OME-TIFF header.
    // This is actually slower than the few endsWith(...) checks below.

    // Rules out unsupported TIFF types

    // Combined TIFF and DICOM created by GE Senographe scanners
    if (read > 131 && buf[128] == 68 && buf[129] == 73 && buf[130] == 67 && buf[131] == 77) {
      return Opener.TIFF_AND_DICOM;
    }

    if (filename.endsWith(".lsm")) {
      return Opener.UNKNOWN; // The LSM Reader plugin opens these files
    }

    // This is only for little endian
    if (littleEndian && filename.endsWith(".flex")) {
      return Opener.UNKNOWN;
    }

    return Opener.TIFF;
  }

  /**
   * Initialise the TIFF image sizes and data structures.
   */
  private void initialise() {
    if (imageData == null) {
      trackProgress.status("Reading images sizes");
      final Ticker ticker = Ticker.createStarted(trackProgress, images.size(), false);

      // All images are TIFF. Get the size of each and count the total frames.
      imageData = new ImageData[images.size()];
      imageSize = new int[images.size()];
      final String[] names = new String[images.size()];
      frames = 0;
      int ok = 0;

      // We can guess for sequential read
      final boolean estimate = getReadHint() == ReadHint.SEQUENTIAL;
      boolean exact = true;

      for (int i = 0; i < names.length; i++) {
        final String path = images.get(i);

        final File file = new File(path);

        // Get the size of each file so we can determine if
        // they can fit into memory. We only use pre-loading for
        // sequential reading if all images fit into memory.
        final long size = getSize(file);

        try (SeekableStream ss = createSeekableStream(path)) {
          final FastTiffDecoder td = FastTiffDecoder.create(ss, path);

          NumberOfImages numberOfImages;
          if (estimate) {
            numberOfImages = td.getNumberOfImages(() -> size);
          } else {
            numberOfImages = td.getNumberOfImages();
          }
          if (numberOfImages.isExact()) {
            trackProgress.log("%s : images=%d (%d bytes)", path, numberOfImages.getImageCount(),
                size);
          } else {
            trackProgress.log("%s : images=%d (approx) (%d bytes)", path,
                numberOfImages.getImageCount(), size);
          }
          if (estimate) {
            // Track if this is exact
            exact = exact && numberOfImages.isExact();
          } else if (numberOfImages.getImageCount() <= 0) {
            // No TIFF images. This will break the non-sequential support
            // using the cumulative size array so remove the image.
            continue;
          }

          frames += numberOfImages.getImageCount();
          imageSize[ok] = frames;
          imageData[ok] = new ImageData(size);
          names[ok] = path;
          ok++;
        } catch (final IOException ex) {
          if (estimate) {
            // This is an untested method so log the error
            ImageJPluginLoggerHelper.getLogger(SeriesImageSource.class.getName()).log(Level.WARNING,
                "Failed to initialise size estimate", ex);
          }
        }
        ticker.tick();
      }

      trackProgress.status("");
      ticker.stop();

      if (ok < images.size()) {
        imageSize = Arrays.copyOf(imageSize, ok);
        imageData = Arrays.copyOf(imageData, ok);
        images.clear();
        images.addAll(Arrays.asList(names));
      }

      // No support for non-sequential access
      if (!exact) {
        imageSize = null;
      }
    }
  }

  /**
   * Gets the image size.
   *
   * @param index the index
   * @return the image size
   */
  int getImageSize(int index) {
    return (index == 0) ? imageSize[index] : imageSize[index] - imageSize[index - 1];
  }

  @Override
  protected boolean openSource() {
    // We now require a tiff series.
    // Object deserialisation of old data may have non tiff images so check.
    if (!isTiffSeries || images.isEmpty()) {
      return false;
    }

    initialise();
    if (frames == 0) {
      return false;
    }

    // Reset if already open
    // (Should we support only closing the sequential reading functionality)
    close();

    // Open the first TIFF image
    lastImage = openImage(0, images.get(0));
    lastImageId = 0;
    if (lastImage != null && lastImage.size > 0) {
      try {
        // Ignore the pixels returned. This will throw if they are null.
        lastImage.getFrame(0);
        setDimensions(lastImage.width, lastImage.height);

        // Attempt to get the origin if a MicroManager image
        final Rectangle roi = FastTiffDecoder.getOrigin(((TiffImage) lastImage).info[0]);
        if (roi != null && roi.width == getWidth() && roi.height == getHeight()) {
          setOrigin(roi.x, roi.y);
        }

        return true;
      } catch (final IOException ex) {
        // Q. Should the problem be reported. Currently if the exception is
        // not bubbled up then the source is not opened.
        if (trackProgress.isLog()) {
          trackProgress.log("Failed to open: %s", ExceptionUtils.getStackTrace(ex));
        }
      }
    }
    return false;
  }

  @Override
  protected synchronized void closeSource() {
    closeQueue();

    if (lastImage != null) {
      lastImage.close(true);
      lastImage = null;
    }
    lastImageId = 0;
  }

  /**
   * Creates a background thread to open the images sequentially.
   *
   * <p>If all the images are smaller than the buffer limit then a single thread is created to read
   * the raw bytes, one to decode the TIFF info and one to read the frames.
   *
   * <p>Otherwise a single thread is created to read each image in turn, decode the TIFF info and
   * read the frames.
   */
  private synchronized void createQueue() {
    // Q. What size is optimal?
    long bytesPerFrame = 1024L * 1024 * 2; // A typical unsigned short image
    // We should have successfully opened the first image and so find the pixel size
    if (imageData != null && imageData[0] != null && imageData[0].tiffImage != null
        && imageData[0].tiffImage.bytesPerFrame > 0) {
      bytesPerFrame = imageData[0].tiffImage.bytesPerFrame;
    }
    // Now create a queue to hold n images in memory
    final int n = Math.max(2, (int) Math.ceil(sequentialReadBufferLimit / (double) bytesPerFrame));
    rawFramesQueue = new CloseableBlockingQueue<>(n);

    if (belowBufferLimit() && images.size() > 1) {
      // For now just support a single thread for each of reading the
      // raw byte data, decoding the TIFF and reading the TIFF pixels.
      // We can control the number of images buffered into memory.
      final int nImages = numberOfImages;

      decodeQueue = new CloseableBlockingQueue<>(nImages);
      readQueue = new CloseableBlockingQueue<>(nImages);

      workers = new ArrayList<>(3);
      threads = new ArrayList<>(3);
      startWorker(new BufferWorker());
      startWorker(new DecodeWorker());
      startWorker(new ReadWorker());
    } else {
      // A single worker thread to read the series
      workers = new ArrayList<>(1);
      threads = new ArrayList<>(1);
      startWorker(new TiffWorker());
    }
  }

  private void startWorker(BaseWorker worker) {
    workers.add(worker);
    final Thread thread = new Thread(worker);
    threads.add(thread);
    thread.start();
  }

  /**
   * Close the background thread.
   */
  private synchronized void closeQueue() {
    if (threads != null) {
      // Signal the workers to stop
      for (final BaseWorker worker : workers) {
        worker.shutdown();
      }

      // Close the queues. This will wake any thread waiting for them to have capacity.
      if (decodeQueue != null) {
        decodeQueue.close(false);
        readQueue.close(false);
      }
      rawFramesQueue.close(false);

      // Join the threads and then set all to null
      for (final Thread thread : threads) {
        try {
          // This thread will be waiting on imageQueue.put() (which has been cleared)
          // or sourceQueue.take() which contains shutdown signals, it should be finished by now
          thread.join();
        } catch (final InterruptedException ex) {
          Thread.currentThread().interrupt();
          throw new ConcurrentRuntimeException(UNEXPECTED_INTERRUPTION, ex);
        }
      }

      threads = null;
      workers = null;
      decodeQueue = null;
      readQueue = null;

      // Do not set this to null as it is used in nextRawFrame()
      rawFramesQueue.close(false);
    }
  }

  /**
   * Close the queues for the in-memory workflow. This should be called on error as it shuts all the
   * queues.
   */
  void closeWorkflowQueues() {
    decodeQueue.close(true);
    readQueue.close(true);
    rawFramesQueue.close(true);
  }

  /**
   * Store the tiff image.
   *
   * @param imageId the image id
   * @param image the image
   */
  void storeTiffImage(int imageId, TiffImage image) {
    // This could be made optional to save memory
    imageData[imageId].tiffImage = image;
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

  private void setDimensions(int maxx, int maxy) {
    width = maxx;
    height = maxy;
  }

  @Override
  protected boolean initialiseSequentialRead() {
    // This should only be called once by the image source each time the series is opened.
    createQueue();
    // Reset the sequential read error
    error = null;
    return true;
  }

  /**
   * Sets the error.
   *
   * @param ex the new error
   */
  void setError(RuntimeException ex) {
    if (error != null) {
      ImageJPluginLoggerHelper.getDefaultLogger().log(Level.WARNING,
          "Encountered a second error during sequential read!");
    }
    error = ex;
  }

  /**
   * {@inheritDoc}
   *
   * @throws DataException If there was an error duing the sequential read
   */
  @Override
  protected Object nextRawFrame() {
    try {
      // We can take if closed
      final Object pixels = rawFramesQueue.take();
      if (pixels != null) {
        return pixels;
      }
      if (error != null) {
        throw error;
      }
    } catch (final InterruptedException ex) {
      Thread.currentThread().interrupt();
      throw new ConcurrentRuntimeException(UNEXPECTED_INTERRUPTION, ex);
    }

    // We are here because the pixels were null so shut down sequential reading
    rawFramesQueue.close(true);

    return null;
  }

  @Override
  protected Object getRawFrame(int frame) {
    if (imageSize == null || !isValid(frame)) {
      return null;
    }

    // Calculate the required image and slice
    int id = Arrays.binarySearch(imageSize, frame);
    if (id < 0) {
      id = -(id + 1);
    }
    // Note that frame is 1-based index and the slice is 0-based.
    final int slice = (id == 0) ? frame - 1 : frame - imageSize[id - 1] - 1;

    // This is not concurrent safe

    // Return from the cache if it exists
    if (id != lastImageId || lastImage == null) {
      if (lastImage != null) {
        lastImage.close(true);
      }

      // Used to cache the TIFF info
      lastImage = null;
      if (id < images.size()) {
        final String path = images.get(id);
        // Check the cache
        TiffImage tiffImage = imageData[id].tiffImage;
        if (tiffImage == null) {
          tiffImage = openImage(id, path);
        }

        lastImage = tiffImage;
      }
    }
    lastImageId = id;
    if (lastImage != null && slice < lastImage.size) {
      try {
        return lastImage.getFrame(slice);
      } catch (final IOException ex) {
        // Q. Should the problem be reported. Currently if the exception is
        // not bubbled up then the frame will be null.
        if (trackProgress.isLog()) {
          trackProgress.log("Failed to open frame %d: %s", frame, ExceptionUtils.getStackTrace(ex));
        }
      }
    }
    return null;
  }

  private TiffImage openImage(int id, String path) {
    TiffImage tiffImage = null;

    // We only need the meta data for the first image. We then assume check
    // all other images in the series match the pixel type and width of the first.
    final boolean pixelInfoOnly = id != 0;

    // Open using specialised TIFF reader for better non-sequential support
    SeekableStream ss = null;
    try {
      ss = createSeekableStream(path, imageData[id].fileSize);

      final FastTiffDecoder td = FastTiffDecoder.create(ss, path);
      td.setTrackProgress(trackProgress);

      // Try and use the index map
      final IndexMap indexMap = td.getIndexMap();

      // Check the image map is the correct size (only if we have sizes)
      if (indexMap != null && (imageSize == null || indexMap.getSize() == getImageSize(id))) {
        // We need the first IFD to define the image pixel type and width/height
        final ExtendedFileInfo fi = td.getTiffInfo(indexMap, 0, pixelInfoOnly);
        // A byte array seekable stream will ignore the close() method so we can re-use it
        tiffImage =
            new TiffImage(indexMap, fi, (ss instanceof ByteArraySeekableStream) ? ss : null);
      } else {
        // Reset after reading the index map
        td.reset();
        // Read all the IFDs
        final ExtendedFileInfo[] info = td.getTiffInfo(pixelInfoOnly);
        if (info != null) {
          // A byte array seekable stream will ignore the close() method so we can re-use it
          tiffImage = new TiffImage(info, (ss instanceof ByteArraySeekableStream) ? ss : null);
        }
      }
    } catch (final IOException ignored) {
      // Ignore
    } finally {
      closeInputStream(ss);
    }

    if (tiffImage == null) {
      // Prevent reading again. Skip the storeTiffImage(...) method
      // as that may be optional and we don't want to repeat the error.
      imageData[id].tiffImage = new TiffImage();
    } else {
      storeTiffImage(id, tiffImage);
    }

    return tiffImage;
  }

  @Override
  public boolean isValid(int frame) {
    return frame > 0 && frame <= frames;
  }

  @Override
  public String toString() {
    String string = super.toString();
    if (!images.isEmpty()) {
      string += String.format(" (%d images: %s ...)", images.size(), images.get(0));
    }
    return string;
  }

  /**
   * Gets the number of threads to use for opening images.
   *
   * @return The number of background threads to use for opening images
   * @deprecated Currently only 1 thread is used for opening images
   */
  @Deprecated
  public int getNumberOfThreads() {
    return numberOfThreads;
  }

  /**
   * Set the number of background threads to use for opening images.
   *
   * @param numberOfThreads The number of background threads to use for opening images
   * @deprecated Currently only 1 thread is used for opening images
   */
  @Deprecated
  public void setNumberOfThreads(int numberOfThreads) {
    this.numberOfThreads = Math.max(1, numberOfThreads);
  }

  /**
   * Gets the buffer limit for reading TIFF images into memory.
   *
   * @return the buffer limit
   */
  public int getBufferLimit() {
    return bufferLimit;
  }

  /**
   * Sets the buffer limit for reading TIFF images into memory.
   *
   * @param bufferLimit the new buffer limit
   */
  public void setBufferLimit(int bufferLimit) {
    this.bufferLimit = bufferLimit;
  }

  /**
   * Gets the buffer limit for sequential reading of a TIFF series. Reading will pause when the
   * output queue is full of images totalling this limit.
   *
   * @return the series buffer limit
   */
  public long getSequentialReadBufferLimit() {
    return sequentialReadBufferLimit;
  }

  /**
   * Sets the buffer limit for sequential reading of a TIFF series. Reading will pause when the
   * output queue is full of images totalling this limit.
   *
   * @param sequentialReadBufferLimit the new sequential read buffer limit
   */
  public void setSequentialReadBufferLimit(long sequentialReadBufferLimit) {
    this.sequentialReadBufferLimit = sequentialReadBufferLimit;
  }

  /**
   * Gets the number of images to buffer into memory.
   *
   * @return the number of images
   */
  public int getNumberOfImages() {
    return numberOfImages;
  }

  /**
   * Sets the number of images to buffer into memory.
   *
   * @param numberOfImages the new number of images
   */
  public void setNumberOfImages(int numberOfImages) {
    this.numberOfImages = Math.max(1, numberOfImages);
  }

  /**
   * Sets the track progress used for monitoring the progress of method execution.
   *
   * @param trackProgress the new track progress
   */
  public void setTrackProgress(TrackProgress trackProgress) {
    this.trackProgress = NullTrackProgress.createIfNull(trackProgress);
  }

  /**
   * Gets the a reference to the file info for image n in the series.
   *
   * @param index the image index
   * @return the file info
   */
  public @Nullable ExtendedFileInfo[] getFileInfo(int index) {
    if (imageData != null && index >= 0 && index < imageData.length) {
      final TiffImage image = imageData[index].tiffImage;
      if (image != null) {
        return image.info;
      }
    }
    return null;
  }
}
