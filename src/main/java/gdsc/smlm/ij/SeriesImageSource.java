package gdsc.smlm.ij;

import java.awt.Rectangle;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;

import com.thoughtworks.xstream.annotations.XStreamOmitField;

import gdsc.core.generics.CloseableBlockingQueue;
import gdsc.core.ij.SeriesOpener;
import gdsc.core.logging.NullTrackProgress;
import gdsc.core.logging.Ticker;
import gdsc.core.logging.TrackProgress;
import gdsc.smlm.results.ImageSource;
import ij.ImagePlus;
import ij.Prefs;
import ij.io.ByteArraySeekableStream;
import ij.io.ExtendedFileInfo;
import ij.io.FastTiffDecoder;
import ij.io.FastTiffDecoder.NumberOfImages;
import ij.io.FileSeekableStream;
import ij.io.ImageReader;
import ij.io.Opener;
import ij.io.SeekableStream;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Represent a series of TIFF image files as a results source. Supports all greyscale images. Only processes channel 0
 * of 32-bit colour images.
 * <p>
 * Assumes that the width,height,depth dimensions of each file are the same. The depth for the last image can be less
 * (i.e. the last of the series) but the {@link #getFrames()} method will return an incorrect value.
 */
public class SeriesImageSource extends ImageSource
{
	/** The buffer limit for reading TIFF images into memory. Default = 30MB */
	private int bufferLimit = 31457280;

	private abstract class Image
	{
		int width, height, size;

		Image(int width, int height, int size)
		{
			this.width = width;
			this.height = height;
			this.size = size;
		}

		int getWidth()
		{
			return width;
		}

		int getHeight()
		{
			return height;
		}

		int getSize()
		{
			return size;
		}

		abstract Object getFrame(int i);

		public void close(boolean freeMemory)
		{

		}
	}

	/**
	 * Support all images ImageJ can read using a fixed image array
	 */
	@SuppressWarnings("unused")
	private class ArrayImage extends Image
	{
		final Object[] imageArray;

		ArrayImage()
		{
			super(0, 0, 0);
			imageArray = null;
		}

		ArrayImage(ImagePlus imp)
		{
			super(imp.getWidth(), imp.getHeight(), imp.getStackSize());
			this.imageArray = imp.getImageStack().getImageArray();
		}

		@Override
		Object getFrame(int i)
		{
			return imageArray[i];
		}
	}

	/**
	 * Special class adapted from ij.io.Opener.openTiffStack(...) and ij.io.FileOpener to handle large contiguous TIFF
	 * images.
	 * <p>
	 * The methods that manipulate the internal state are synchronized to prevent multiple threads causing read errors.
	 */
	private class TiffImage extends Image
	{
		final ExtendedFileInfo[] info;
		ExtendedFileInfo fi; // Pointer to the file info used by the ImageReader
		boolean contiguous;
		final long bytesPerFrame;
		ImageReader reader = null;
		InputStream is = null;
		/**
		 * A reference to a seekable stream which may be buffered in memory.
		 */
		SeekableStream ss = null;
		boolean inMemory = false;
		/**
		 * Flag indicating that no errors reading the image have occurred
		 */
		boolean canRead = true;
		/**
		 * The number of frames that have been read from the input stream
		 */
		int frameCount = 0;

		TiffImage(ExtendedFileInfo[] info, SeekableStream ss)
		{
			super(info[0].width, info[0].height, 0);
			this.info = info;
			fi = info[0];

			// Only support certain types
			if (isSupported(fi.fileType))
			{
				// We use the opened SeekableStream.
				// This may be in-memory data or may be from a random access file.
				if (ss != null)
				{
					try
					{
						ss.seek(0);
						this.ss = ss;
						is = ss;

						// Store if the stream contains in-memory data						
						this.inMemory = ss instanceof ByteArraySeekableStream;
					}
					catch (IOException e)
					{
					}
				}

				// Determine the number of images
				if (info.length > 1)
				{
					if (allSameSizeAndType(info))
						size = info.length;
				}
				else
				{
					size = fi.nImages;
					contiguous = true;
				}
				bytesPerFrame = getBytesPerFrame(fi.fileType);
				reader = new ImageReader(fi);
			}
			else
			{
				canRead = false;
				bytesPerFrame = 0;
			}
		}

		TiffImage()
		{
			super(0, 0, 0);
			info = null;
			bytesPerFrame = 0;
			canRead = false;
		}

		private boolean isSupported(int fileType)
		{
			switch (fileType)
			{
				// Greyscale images as we just want the raw pixels with no color model
				case ExtendedFileInfo.GRAY8:
				case ExtendedFileInfo.GRAY16_SIGNED:
				case ExtendedFileInfo.GRAY16_UNSIGNED:
				case ExtendedFileInfo.GRAY12_UNSIGNED:
				case ExtendedFileInfo.GRAY32_INT:
				case ExtendedFileInfo.GRAY32_UNSIGNED:
				case ExtendedFileInfo.GRAY32_FLOAT:
				case ExtendedFileInfo.GRAY24_UNSIGNED:
				case ExtendedFileInfo.GRAY64_FLOAT:
					return true;
				default:
					return false;
			}
		}

		private int getBytesPerFrame(int fileType)
		{
			switch (fileType)
			{
				case ExtendedFileInfo.GRAY8:
					return width * height;
				case ExtendedFileInfo.GRAY16_SIGNED:
				case ExtendedFileInfo.GRAY16_UNSIGNED:
					return 2 * width * height;
				case ExtendedFileInfo.GRAY32_INT:
				case ExtendedFileInfo.GRAY32_UNSIGNED:
				case ExtendedFileInfo.GRAY32_FLOAT:
					return 4 * width * height;
				case ExtendedFileInfo.GRAY64_FLOAT:
					return 8 * width * height;
				case ExtendedFileInfo.GRAY24_UNSIGNED:
					return 3 * width * height;
				case ExtendedFileInfo.GRAY12_UNSIGNED:
					return (int) (width * height * 1.5);
				default:
					return 0;
			}
		}

		/**
		 * Check the file info are all the same size and type and have just 1 image
		 * 
		 * @param info
		 * @return True if all the same size and type
		 */
		boolean allSameSizeAndType(ExtendedFileInfo[] info)
		{
			boolean ok = info[0].nImages == 1;
			long startingOffset = info[0].getOffset();
			long size = info[0].width * info[0].height * info[0].getBytesPerPixel();
			for (int i = 1; i < info.length && ok; i++)
			{
				ok &= info[i].fileType == info[0].fileType && info[i].width == info[0].width &&
						info[i].height == info[0].height && info[i].nImages == 1;
				contiguous &= info[i].getOffset() == startingOffset + i * size;
			}
			if (ok)
			{
				// We can read using only the first ExtendedFileInfo object if this is contiguous
				if (contiguous)
					info[0].gapBetweenImages = 0;
			}
			return ok;
		}

		Object nextFrame()
		{
			try
			{
				// Skip ahead
				long skip;

				if (contiguous)
				{
					// If sequential reading just skip the gap between frames
					if (frameCount == 0)
					{
						// The first frame we know the exact offset
						skip = fi.getOffset();
					}
					else
					{
						// If sequential reading just skip the gap between frames
						skip = fi.gapBetweenImages;
					}
				}
				else
				{
					// Adapted from ij.io.Opener.openTiffStack(...)

					// Each image offset is described by a separate ExtendedFileInfo object
					skip = info[frameCount].getOffset();
					fi.stripOffsets = info[frameCount].stripOffsets;
					fi.stripLengths = info[frameCount].stripLengths;

					if (frameCount != 0)
					{
						// We must subtract the current file location.
						skip -= (info[frameCount - 1].getOffset() + bytesPerFrame);
						if (skip < 0L)
							throw new IllegalStateException("Bad TIFF offset " + skip);
					}
				}

				// Store the number of frames that have been read
				frameCount++;

				//long t = System.nanoTime();
				Object pixels = reader.readPixels(is, skip);
				//System.out.printf("IO Time = %f ms\n", (System.nanoTime()-t)/1e6);
				return pixels;
			}
			catch (Exception e)
			{
				System.out.println(e.toString());
				e.printStackTrace();
				canRead = false;
			}

			return null;
		}

		@Override
		synchronized Object getFrame(int i)
		{
			if (i < frameCount)
			{
				// Close will either seek to the start or free resources of a standard input stream.
				// Non-sequential access has poor support for a standard input stream as we just 
				// have to re-open the file from the beginning again.
				close(false);
			}

			if (!openInputStream())
				return null;

			try
			{

				// Skip ahead
				long skip;

				if (contiguous)
				{
					// Read using the first ExtendedFileInfo object

					if (i == 0)
					{
						// The first frame we know the exact offset
						skip = fi.getOffset();
					}
					else if (i == frameCount)
					{
						// If sequential reading just skip the gap between frames
						skip = fi.gapBetweenImages;
					}
					else
					{
						// Skipping ahead
						int nFrames = i - frameCount;
						skip = (bytesPerFrame + fi.gapBetweenImages) * nFrames;

						// Ensure we skip to the start position if nothing has been read
						if (frameCount == 0)
							skip += fi.getOffset();
					}
				}
				else
				{
					// Adapted from ij.io.Opener.openTiffStack(...)

					// Each image offset is described by a separate ExtendedFileInfo object
					skip = info[i].getOffset();
					fi.stripOffsets = info[i].stripOffsets;
					fi.stripLengths = info[i].stripLengths;

					if (frameCount != 0)
					{
						// We must subtract the current file location.
						skip -= (info[frameCount - 1].getOffset() + bytesPerFrame);
						if (skip < 0L)
							throw new IllegalStateException("Bad TIFF offset " + skip);
					}
				}

				// Store the number of frames that have been read
				frameCount = i + 1;

				//long t = System.nanoTime();
				Object pixels = reader.readPixels(is, skip);
				//System.out.printf("IO Time = %f ms\n", (System.nanoTime()-t)/1e6);
				return pixels;
			}
			catch (Exception e)
			{
				System.out.println(e.toString());
				e.printStackTrace();
				canRead = false;
			}

			return null;
		}

		public boolean openInputStream()
		{
			if (is == null && canRead)
			{
				try
				{
					is = createInputStream(fi);
				}
				catch (Exception e)
				{
					System.out.println(e.toString());
					e.printStackTrace();
				}
				finally
				{
					if (is == null)
					{
						canRead = false;
					}
				}
			}
			return canRead;
		}

		/** Returns an InputStream for the image described by this ExtendedFileInfo. */
		public FileInputStream createInputStream(ExtendedFileInfo fi) throws IOException
		{
			if (inMemory)
				throw new IllegalStateException("No input file for in-memory image");

			FileInputStream is = null;
			if (fi.directory.length() > 0 && !fi.directory.endsWith(Prefs.separator))
				fi.directory += Prefs.separator;
			File f = new File(fi.directory + fi.fileName);
			if (f == null || !f.exists() || f.isDirectory() || !validateFileInfo(f, fi))
				is = null;
			else
				is = new FileInputStream(f);
			return is;
		}

		/**
		 * Adapted from ij.io.FileOpener.validateFileInfo
		 */
		boolean validateFileInfo(File f, ExtendedFileInfo fi)
		{
			long offset = fi.getOffset();
			long length = 0;
			if (fi.width <= 0 || fi.height <= 0)
			{
				return false;
			}
			if (offset >= 0 && offset < 1000L)
				return true;
			if (offset < 0L)
			{
				return false;
			}
			if (fi.fileType == ExtendedFileInfo.BITMAP || fi.compression != ExtendedFileInfo.COMPRESSION_NONE)
				return true;
			length = f.length();
			long size = fi.width * fi.height * fi.getBytesPerPixel();
			size = fi.nImages > 1 ? size : size / 4;
			if (fi.height == 1)
				size = 0; // allows plugins to read info of unknown length at end of file
			if (offset + size > length)
			{
				return false;
			}
			return true;
		}

		public void reset()
		{
			if (frameCount != 0)
				close(false);
		}

		@Override
		synchronized public void close(boolean freeResources)
		{
			// Reset
			frameCount = 0;

			if (ss != null && !freeResources)
			{
				// We can seek to the start
				try
				{
					ss.seek(0);
					return;
				}
				catch (IOException e)
				{
					// Fall through to close the resources
				}
			}

			// This is done when sequentially reading so we clear the memory
			inMemory = false;
			ss = null;

			if (is != null)
			{
				try
				{
					is.close();
				}
				catch (IOException e)
				{
					// Ignore
				}
				// Reset
				is = null;
			}
		}
	}

	private abstract class BaseWorker implements Runnable
	{
		volatile boolean run = true;
	}

	/**
	 * Add source images to the queue to be read.
	 */
	private class TIFFWorker extends BaseWorker
	{
		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		public void run()
		{
			boolean error = false;
			try
			{
				for (int currentImage = 0; run && currentImage < images.size(); currentImage++)
				{
					final String path = images.get(currentImage);
					//System.out.println("Reading " + path);

					SeekableStream ss = createSeekableStream(path);
					if (ss == null)
					{
						error = true;
						break;
					}

					trackProgress.log("Reading TIFF info %s", path);

					// Check the cache to avoid a re-read. This will always be the case
					// for the first image as that is opened in openSource().
					TiffImage image = imageData[currentImage].tiffImage;
					if (image == null)
					{
						ExtendedFileInfo[] info = getTiffInfo(ss, path);
						if (info == null)
						{
							error = true;
							break;
						}
						image = new TiffImage(info, null);
						storeTiffImage(currentImage, image);
					}

					// Check dimensions
					// Check it is the expected size
					if (image.getWidth() != getWidth() || image.getHeight() != getHeight())
					{
						error = true;
						break;
					}
					if (imageSize != null && image.getSize() != getImageSize(currentImage))
					{
						error = true;
						break;
					}

					trackProgress.log("Reading TIFF %s", path);

					image.reset();
					if (!image.openInputStream())
					{
						error = true;
						break;
					}

					try
					{
						// Read all the frames sequentially
						for (int i = 0; i < image.size; i++)
						{
							Object pixels = image.nextFrame();
							if (pixels == null)
							{
								error = true;
								break;
							}

							// This will block until the queue has capacity or is closed
							if (!rawFramesQueue.putAndConfirm(pixels))
								break;
						}
					}
					finally
					{
						// Close the image
						image.close(true);
					}
				}
			}
			catch (InterruptedException e)
			{
				// This is from the queue put method, possibly an interrupt on the queue or thread? 
				System.out.println(e.toString());
				e.printStackTrace();
				error = true;
			}

			rawFramesQueue.close(error);

			run = false;
		}
	}

	private class NextSource
	{
		final byte[] buffer;
		final int imageId;
		//ExtendedFileInfo[] info;
		Image image;

		public NextSource(byte[] buffer, int imageId)
		{
			this.buffer = buffer;
			this.imageId = imageId;
		}

		//public void setFileInfo(ExtendedFileInfo[] info)
		//{
		//	this.info = info;
		//}

		public void setImage(Image image)
		{
			this.image = image;
		}

		//public NextSource()
		//{
		//	this(null, -1);
		//}
	}

	/**
	 * Read source image files into memory.
	 */
	private class BufferWorker extends BaseWorker
	{
		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		public void run()
		{
			boolean error = false;
			try
			{
				for (int currentImage = 0; run && currentImage < images.size(); currentImage++)
				{
					final String path = images.get(currentImage);
					//System.out.println("Reading " + path);

					File file = new File(path);
					FileInputStream fis = openStream(file);
					if (fis == null)
					{
						error = true;
						break;
					}

					try
					{
						// We already know they fit into memory (i.e. the size is not zero)
						int size = (int) imageData[currentImage].fileSize;
						byte[] buf = new byte[size];
						int read = fis.read(buf);
						if (read != size)
						{
							error = true;
							break;
						}

						// This may be closed upon error
						if (!decodeQueue.putAndConfirm(new NextSource(buf, currentImage)))
							break;
					}
					catch (IOException e)
					{
						System.out.println(e.toString());
						e.printStackTrace();
						error = true;
						break;
					}
					finally
					{
						try
						{
							fis.close();
						}
						catch (IOException e)
						{
							// Ignore
						}
					}
				}
			}
			catch (InterruptedException e)
			{
				// This is from the queue put method, possibly an interrupt on the queue or thread? 
				System.out.println(e.toString());
				e.printStackTrace();
				error = true;
			}

			if (error)
				closeWorkflowQueues();
			else
				decodeQueue.close(false);

			run = false;
		}
	}

	/**
	 * Read source image files into memory.
	 */
	private class DecodeWorker extends BaseWorker
	{
		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		public void run()
		{
			boolean error = false;
			try
			{
				while (run)
				{
					NextSource nextSource = decodeQueue.take();
					if (nextSource == null || !run)
						break;
					final int currentImage = nextSource.imageId;
					if (currentImage == -1)
						break;

					final String path = images.get(currentImage);
					//System.out.println("Reading " + path);

					SeekableStream ss = new ByteArraySeekableStream(nextSource.buffer);

					// Re-use the cache of the file info
					ExtendedFileInfo[] info = null;
					TiffImage image = imageData[currentImage].tiffImage;
					if (image == null)
					{
						info = getTiffInfo(ss, path);
						if (info == null)
						{
							error = true;
							break;
						}

						// Create as in-memory
						image = new TiffImage(info, ss);

						// Check dimensions
						// Check it is the expected size
						if (image.getWidth() != getWidth() || image.getHeight() != getHeight())
						{
							error = true;
							break;
						}
						if (imageSize != null && image.getSize() != getImageSize(currentImage))
						{
							error = true;
							break;
						}

						storeTiffImage(currentImage, image);
					}
					else
					{
						// Update to be in-memory. This will be used when reading the TIFF.
						// It will subsequently be closed to free-memory.
						image.inMemory = true;
						image.ss = ss;
						image.is = ss;
					}

					//nextSource.setFileInfo(info);
					nextSource.setImage(image);

					// This may be closed upon error
					if (!readQueue.putAndConfirm(nextSource))
						break;
				}
			}
			catch (InterruptedException e)
			{
				// This is from the queue put method, possibly an interrupt on the queue or thread? 
				System.out.println(e.toString());
				e.printStackTrace();
				error = true;
			}

			if (error)
				closeWorkflowQueues();
			else
				readQueue.close(false);

			run = false;
		}
	}

	/**
	 * Read source image files into memory.
	 */
	private class ReadWorker extends BaseWorker
	{
		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		public void run()
		{
			boolean error = false;
			try
			{
				while (run)
				{
					NextSource nextSource = readQueue.take();
					if (nextSource == null || !run)
						break;
					final int currentImage = nextSource.imageId;
					if (currentImage == -1)
						break;

					// It is assumed we are reading a TIFF series
					TiffImage image = (TiffImage) nextSource.image;

					image.reset();
					if (!image.openInputStream())
					{
						error = true;
						break;
					}

					try
					{
						// Read all the frames sequentially
						for (int i = 0; i < image.size; i++)
						{
							Object pixels = image.nextFrame();
							if (pixels == null)
							{
								error = true;
								break;
							}

							// This may be closed upon error
							if (!rawFramesQueue.putAndConfirm(pixels))
								break;
						}
					}
					finally
					{
						// Close the image
						image.close(true);
					}
				}
			}
			catch (InterruptedException e)
			{
				// This is from the queue put method, possibly an interrupt on the queue or thread? 
				System.out.println(e.toString());
				e.printStackTrace();
				error = true;
			}

			if (error)
				closeWorkflowQueues();
			else
				rawFramesQueue.close(false);

			run = false;
		}
	}

	//	/**
	//	 * Add source images to the queue to be read.
	//	 */
	//	private class SourceWorker extends BaseWorker
	//	{
	//		/*
	//		 * (non-Javadoc)
	//		 * 
	//		 * @see java.lang.Runnable#run()
	//		 */
	//		public void run()
	//		{
	//			try
	//			{
	//				// Read each image in sequence
	//				for (int currentImage = 0; run && currentImage < images.size(); currentImage++)
	//				{
	//					ExtendedFileInfo[] info = null;
	//					SeekableStream ss = null;
	//
	//					final String path = images.get(currentImage);
	//					//System.out.println("Reading " + images.get(currentImage));
	//
	//					// This may contain a custom SeekableStream that wraps an in-memory RandomAccessFile 
	//					ss = createSeekableStream(path);
	//
	//					if (ss != null)
	//					{
	//						TiffDecoder td = new CustomTiffDecoder(ss, path);
	//						if (logProgress)
	//						{
	//							long time = System.currentTimeMillis();
	//							if (time - lastTime > 500)
	//							{
	//								lastTime = time;
	//								IJ.log("Reading TIFF info " + path);
	//							}
	//						}
	//						try
	//						{
	//							//td.enableDebugging();
	//							// This will close the seekable stream.
	//							// If it is a custom in-memory object then the close call is ignored. 
	//							info = td.getTiffInfo();
	//							//if (info != null)
	//							//{
	//							//	System.out.println(info[0].debugInfo);
	//							//	
	//							//	// This contains OME TIFF metadata as serialised JSON when using 
	//							//	// MicroManager to save the TIFF
	//							//	System.out.println(info[0].info);
	//							//	
	//							//	ij.io.FileOpener fo = new ij.io.FileOpener(info[0]);
	//							//	java.util.Properties p = fo.decodeDescriptionString(info[0]);
	//							//	System.out.println(p);
	//							//	//double ox = 0, oy = 0;
	//							//	//if (p.containsKey("xorigin"))
	//							//	//	ox = Double.parseDouble(p.getProperty("xorigin"));
	//							//	//if (p.containsKey("yorigin"))
	//							//	//	oy = Double.parseDouble(p.getProperty("yorigin"));
	//							//	//// Should the origin be converted by the units?
	//							//	//setOrigin((int)ox, (int)oy);
	//							//}
	//						}
	//						catch (IOException e)
	//						{
	//							// TODO Auto-generated catch block
	//							e.printStackTrace();
	//						}
	//					}
	//
	//					sourceQueue.put(new NextSource(info, (inMemory[0]) ? ss : null, currentImage));
	//				}
	//			}
	//			catch (FileNotFoundException e)
	//			{
	//				e.printStackTrace();
	//			}
	//			catch (SecurityException e)
	//			{
	//
	//			}
	//			catch (InterruptedException e)
	//			{
	//				// This is from the queue put method, possibly an interrupt on the queue or thread? 
	//				System.out.println(e.toString()); e.printStackTrace();
	//			}
	//
	//			// Add jobs to shutdown all the workers
	//			try
	//			{
	//				for (int i = 0; run && i < workers.size(); i++)
	//					sourceQueue.put(new NextSource());
	//			}
	//			catch (InterruptedException e)
	//			{
	//				// This is from the queue put method, possibly an interrupt on the queue or thread? 
	//				// TODO - How should this be handled?
	//				System.out.println(e.toString()); e.printStackTrace();
	//			}
	//
	//			run = false;
	//		}
	//	}

	/**
	 * Creates the seekable stream.
	 *
	 * @param path
	 *            the path
	 * @param size
	 *            the size of the file
	 * @return the seekable stream
	 */
	private SeekableStream createSeekableStream(String path, long size)
	{
		File file = new File(path);

		// Don't buffer massive images into memory
		if (belowBufferLimit(size))
		{
			SeekableStream ss = readByteArraySeekableStream(file, size);
			if (ss != null)
				return ss;
		}

		try
		{
			return new FileSeekableStream(file);
		}
		catch (FileNotFoundException e)
		{
			e.printStackTrace();
		}
		catch (SecurityException e)
		{
		}

		return null;
	}

	/**
	 * Creates the seekable stream.
	 *
	 * @param path
	 *            the path
	 * @return the seekable stream
	 */
	private SeekableStream createSeekableStream(String path)
	{
		try
		{
			return new FileSeekableStream(path);
		}
		catch (FileNotFoundException e)
		{
			e.printStackTrace();
		}
		catch (SecurityException e)
		{
		}

		return null;
	}

	private ExtendedFileInfo[] getTiffInfo(SeekableStream ss, String path)
	{
		FastTiffDecoder td = new FastTiffDecoder(ss, path);
		td.setTrackProgress(trackProgress);
		try
		{
			//td.enableDebugging();
			// This will close the seekable stream.
			ExtendedFileInfo[] info = td.getTiffInfo();
			//if (info != null)
			//{
			//	System.out.println(info[0].debugInfo);
			//
			//	// This contains OME TIFF metadata as serialised JSON when using 
			//	// MicroManager to save the TIFF
			//	if (info[0].summaryMetaData != null)
			//		System.out.println(info[0].summaryMetaData.substring(0, Math.min(1000, info[0].summaryMetaData.length())));
			//	//	
			//	//	ij.io.FileOpener fo = new ij.io.FileOpener(info[0]);
			//	//	java.util.Properties p = fo.decodeDescriptionString(info[0]);
			//	//	System.out.println(p);
			//	//	//double ox = 0, oy = 0;
			//	//	//if (p.containsKey("xorigin"))
			//	//	//	ox = Double.parseDouble(p.getProperty("xorigin"));
			//	//	//if (p.containsKey("yorigin"))
			//	//	//	oy = Double.parseDouble(p.getProperty("yorigin"));
			//	//	//// Should the origin be converted by the units?
			//	//	//setOrigin((int)ox, (int)oy);
			//}
			return info;
		}
		catch (IOException e)
		{
			// This is from the TIFF decoder
			e.printStackTrace();

			try
			{
				// Close this if the TIFF decoder errored
				ss.close();
			}
			catch (IOException ex)
			{
			}
		}
		return null;
	}

	private boolean belowBufferLimit(long size)
	{
		return (size > 0 && size <= bufferLimit);
	}

	private static long getSize(File file)
	{
		try
		{
			return file.length();
		}
		catch (SecurityException e)
		{
		}
		return 0;
	}

	private SeekableStream readByteArraySeekableStream(File file, long size)
	{
		FileInputStream fis = openStream(file);
		if (fis != null)
		{
			try
			{
				byte[] buf = new byte[(int) size];
				int read = fis.read(buf);
				if (read == size)
					return new ByteArraySeekableStream(buf);
			}
			catch (IOException e)
			{
				// At the moment if we ignore this then the ImageWorker will open the file
				// rather than process from the memory stream.
				System.out.println(e.toString());
				e.printStackTrace();
			}
			finally
			{
				try
				{
					fis.close();
				}
				catch (IOException e)
				{
					// Ignore
				}
			}
		}
		return null;
	}

	private static FileInputStream openStream(File file)
	{
		try
		{
			return new FileInputStream(file);
		}
		catch (FileNotFoundException e)
		{
		}
		catch (SecurityException e)
		{
		}
		return null;
	}

	//	private class ImageWorker extends BaseWorker
	//	{
	//		final int id;
	//
	//		public ImageWorker(int id)
	//		{
	//			this.id = id;
	//		}
	//
	//		/*
	//		 * (non-Javadoc)
	//		 * 
	//		 * @see java.lang.Runnable#run()
	//		 */
	//		public void run()
	//		{
	//			try
	//			{
	//				while (run)
	//				{
	//					NextSource nextSource = sourceQueue.take();
	//					if (nextSource == null || !run)
	//						break;
	//					final int currentImage = nextSource.imageId;
	//					if (currentImage == -1)
	//						break;
	//
	//					// Open the image
	//					if (logProgress)
	//					{
	//						long time = System.currentTimeMillis();
	//						if (time - lastTime > 500)
	//						{
	//							lastTime = time;
	//							IJ.log("Opening " + images.get(currentImage));
	//						}
	//					}
	//
	//					// Only support TIFF images with pre-read ExtendedFileInfo
	//					Image image = null;
	//					if (nextSource.info != null)
	//					{
	//						image = new TiffImage(nextSource.info, nextSource.ss, nextSource.ss != null);
	//					}
	//					else
	//					{
	//						// An empty image
	//						image = new ArrayImage();
	//					}
	//
	//					//System.out.println(id + ": Opened " + images.get(currentImage));
	//
	//					if (image.getSize() != 0)
	//					{
	//						if (width == 0)
	//						{
	//							// Initialise dimensions on the first valid image
	//							setDimensions(image.getWidth(), image.getHeight());
	//						}
	//						else
	//						{
	//							// Check dimensions
	//							if (image.getWidth() != getWidth() || image.getHeight() != getHeight())
	//							{
	//								// Return no image data
	//								image = new ArrayImage();
	//							}
	//						}
	//					}
	//
	//					imageQueue.put(new NextImage(image, currentImage));
	//				}
	//			}
	//			catch (Exception e)
	//			{
	//				// TODO - handle appropriately 
	//				System.out.println(id + ": " + e.toString());
	//			}
	//			finally
	//			{
	//				run = false;
	//
	//				// Signal that no more images are available
	//				imageQueue.offer(new NextImage());
	//			}
	//		}
	//	}

	private ArrayList<String> images;
	/**
	 * Flag indicating if the image series contains only TIFF images. No other formats are currently supported.
	 */
	public final boolean isTiffSeries;

	/**
	 * Hold image data
	 */
	private class ImageData
	{
		long fileSize;
		TiffImage tiffImage;

		public ImageData(long size)
		{
			fileSize = size;
		}
	}

	/**
	 * Contains the cumulative size of the TIFF image series. Used for a binary search to find the image for the
	 * frame.
	 */
	@XStreamOmitField
	private int[] imageSize;

	/**
	 * Used to cache data about the TIFF images
	 */
	@XStreamOmitField
	private ImageData[] imageData;

	//	// Used for sequential read
	//	@XStreamOmitField
	//	private Image image = null;
	//	@XStreamOmitField
	//	private int nextImageId;
	//	@XStreamOmitField
	//	private NextSource[] nextImages;
	//	@XStreamOmitField
	//	private int currentSlice;

	// Used for frame-based read
	@XStreamOmitField
	private Image lastImage = null;
	@XStreamOmitField
	private int lastImageId;

	@XStreamOmitField
	private long lastTime = 0;
	private int numberOfThreads = 1;
	private int numberOfImages = 1;

	@XStreamOmitField
	private TrackProgress trackProgress = NullTrackProgress.INSTANCE;

	// Used to process the files into images
	@XStreamOmitField
	private ArrayList<BaseWorker> workers = null;
	@XStreamOmitField
	private ArrayList<Thread> threads = null;

	/** The queue for in-memory buffered images awaiting TIFF decoding. */
	@XStreamOmitField
	private CloseableBlockingQueue<NextSource> decodeQueue;
	/** The queue for in-memory buffered images awaiting TIFF reading. */
	@XStreamOmitField
	private CloseableBlockingQueue<NextSource> readQueue;

	/** Used for sequential reading to queue the raw frames */
	@XStreamOmitField
	private CloseableBlockingQueue<Object> rawFramesQueue = null;

	/**
	 * Create a new image source using the given image series
	 * 
	 * @param name
	 * @param path
	 */
	public SeriesImageSource(String name, SeriesOpener series)
	{
		super(name);
		if (series != null)
		{
			String[] names = series.getImageList();
			for (int i = 0; i < names.length; i++)
			{
				names[i] = new File(series.getPath(), names[i]).getPath();
			}
			isTiffSeries = isTiffSeries(names);
		}
		else
		{
			isTiffSeries = false;
		}
	}

	/**
	 * Create a new image source using the directory containing the images.
	 * <p>
	 * The directory is opened using {@link gdsc.core.ij.SeriesOpener }
	 * 
	 * @param name
	 * @param path
	 */
	public SeriesImageSource(String name, String directory)
	{
		this(name, new SeriesOpener(directory, false, 1));
	}

	/**
	 * Create a new image source using the paths to the images
	 * 
	 * @param name
	 * @param filenames
	 *            the full path names of the image files
	 */
	public SeriesImageSource(String name, String[] filenames)
	{
		super(name);
		isTiffSeries = isTiffSeries(filenames);
	}

	private boolean isTiffSeries(String[] filenames)
	{
		// Create this as it is needed for XStream serialisation
		images = new ArrayList<String>();

		for (int i = 0; i < filenames.length; i++)
		{
			int fileType = getFileType(filenames[i]);
			if (fileType != Opener.TIFF)
				// Only support TIFF images
				return false;
		}

		// All images are TIFF so store the filenames
		for (int i = 0; i < filenames.length; i++)
		{
			images.add(filenames[i]);
		}

		return true;
	}

	/**
	 * Attempts to determine the image file type by looking for
	 * 'magic numbers' and the file name extension.
	 * <p>
	 * Copied from ij.io.Opener and removed all but the TIFF identification.
	 */
	private int getFileType(String filename)
	{
		File file = new File(filename);
		InputStream is;
		byte[] buf = new byte[132];
		try
		{
			is = new FileInputStream(file);
			is.read(buf, 0, 132);
			is.close();
		}
		catch (IOException e)
		{
			return Opener.UNKNOWN;
		}

		int b0 = buf[0] & 255, b1 = buf[1] & 255, b2 = buf[2] & 255, b3 = buf[3] & 255;
		//IJ.log("getFileType: "+ name+" "+b0+" "+b1+" "+b2+" "+b3);

		// Combined TIFF and DICOM created by GE Senographe scanners
		if (buf[128] == 68 && buf[129] == 73 && buf[130] == 67 && buf[131] == 77 &&
				((b0 == 73 && b1 == 73) || (b0 == 77 && b1 == 77)))
			return Opener.TIFF_AND_DICOM;

		// Big-endian TIFF ("MM")
		String name = file.getName();
		if (name.endsWith(".lsm"))
			return Opener.UNKNOWN; // The LSM Reader plugin opens these files
		if (b0 == 73 && b1 == 73 && b2 == 42 && b3 == 0 && !(name.endsWith(".flex")))
			return Opener.TIFF;

		// Little-endian TIFF ("II")
		if (b0 == 77 && b1 == 77 && b2 == 0 && b3 == 42)
			return Opener.TIFF;

		return Opener.UNKNOWN;
	}

	/**
	 * Initialise the TIFF image sizes and data structures.
	 */
	private void initialise()
	{
		if (imageData == null)
		{
			trackProgress.status("Reading images sizes");
			Ticker ticker = Ticker.createStarted(trackProgress, images.size(), false);

			// All images are TIFF. Get the size of each and count the total frames.
			imageData = new ImageData[images.size()];
			imageSize = new int[images.size()];
			String[] names = new String[images.size()];
			frames = 0;
			int ok = 0;

			// We can guess for sequential read
			boolean estimate = getReadHint() == ReadHint.SEQUENTIAL;
			boolean exact = true;

			for (int i = 0; i < names.length; i++)
			{
				String path = images.get(i);

				SeekableStream ss = null;
				try
				{
					File file = new File(path);

					// Get the size of each file so we can determine if 
					// they can fit into memory. We only use pre-loading for
					// sequential reading if all images fit into memory.
					long size = getSize(file);

					//System.out.printf("%s = %d bytes\n", path, size);

					ss = createSeekableStream(path);
					FastTiffDecoder td = new FastTiffDecoder(ss, path);

					NumberOfImages nImages = td.getNumberOfImages(estimate);
					if (nImages.exact)
						trackProgress.log("%s : images=%d (%d bytes)", path, nImages.numberOfImages, size);
					else
						trackProgress.log("%s : images=%d (approx) (%d bytes)", path, nImages.numberOfImages, size);
					if (estimate)
					{
						// Track if this is exact
						exact = exact && nImages.exact;
					}
					else if (nImages.numberOfImages <= 0)
					{
						// No TIFF images. This will break the non-sequential support
						// using the cumulative size array so remove the image. 
						continue;
					}

					frames += nImages.numberOfImages;
					imageSize[ok] = frames;
					imageData[ok] = new ImageData(size);
					names[ok] = path;
					ok++;
				}
				catch (Throwable e)
				{
					if (estimate)
						// This is an utested method so log the error
						e.printStackTrace();
				}
				finally
				{
					if (ss != null)
					{
						try
						{
							ss.close();
						}
						catch (IOException e1)
						{
						}
					}
				}
				ticker.tick();
			}

			trackProgress.status("");
			ticker.stop();

			if (ok < images.size())
			{
				imageSize = Arrays.copyOf(imageSize, ok);
				imageData = Arrays.copyOf(imageData, ok);
				images.clear();
				images.addAll(Arrays.asList(names));
			}

			// No support for non-sequential access
			if (!exact)
				imageSize = null;
		}
	}

	private int getImageSize(int i)
	{
		return (i == 0) ? imageSize[i] : imageSize[i] - imageSize[i - 1];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#openSource()
	 */
	@Override
	protected boolean openSource()
	{
		// We now require a tiff series. 
		// Object deserialisation of old data may have non tiff images so check. 
		if (!isTiffSeries || images.isEmpty())
			return false;

		initialise();
		if (frames == 0)
			return false;

		// Reset if already open 
		// (Should we support only closing the sequential reading functionality) 
		close();

		// Open the first TIFF image
		lastImage = openImage(0, images.get(0));
		lastImageId = 0;
		if (lastImage != null && lastImage.size > 0)
		{
			Object pixels = lastImage.getFrame(0);
			if (pixels != null)
			{
				setDimensions(imageData[0].tiffImage.width, imageData[0].tiffImage.height);

				// Attempt to get the origin if a MicroManager image
				Rectangle roi = FastTiffDecoder.getOrigin(imageData[0].tiffImage.info[0]);
				if (roi != null && roi.width == getWidth() && roi.height == getHeight())
					setOrigin(roi.x, roi.y);

				return true;
			}
		}
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#closeSource()
	 */
	protected void closeSource()
	{
		if (threads != null)
			closeQueue();
		//setDimensions(0, 0, 0);
		//if (image != null)
		//	image.close(true);

		if (lastImage != null)
			lastImage.close(true);
		//image = 
		lastImage = null;
		//nextImageId = currentSlice = 
		lastImageId = 0;
	}

	//	private Image getNextImage()
	//	{
	//		image = null;
	//		currentSlice = 0;
	//		if (nextImageId < nextImages.length)
	//		{
	//			try
	//			{
	//				for (;;)
	//				{
	//					// Images may be out of order due to multiple thread processing.
	//					// Check if we have processed the next image we need.
	//					NextImage next = nextImages[nextImageId];
	//					if (next != null)
	//					{
	//						// Clear memory
	//						nextImages[nextImageId] = null;
	//
	//						nextImageId++;
	//
	//						// Check if there is image data. It may be null if the image was invalid
	//						if (next.image != null)
	//						{
	//							image = next.image;
	//
	//							//System.out.println("Found image: " + images.get(next.image));
	//
	//							// Fill cache
	//							lastImage = image;
	//							lastImageId = next.imageId;
	//							// The cache is used for non-sequential reading. To prevent 
	//							// memory usage during sequential reading only cache the first 
	//							// one as this is generated when the source is opened and it is 
	//							// unclear if the source is to be used non-sequentially. 
	//							if (lastImageId == 0 && image instanceof TiffImage)
	//								storeTiffImage(lastImageId, (TiffImage) image);
	//							return image;
	//						}
	//						else
	//						{
	//							// Process the next stored image
	//							continue;
	//						}
	//					}
	//
	//					// We are still awaiting the next image.
	//					// Get the images processed by the worker threads.
	//					next = imageQueue.poll();
	//
	//					// If there is nothing then check if any workers are alive
	//					if (next == null && workersRunning())
	//					{
	//						// This will block until something produces an image
	//						next = imageQueue.take();
	//					}
	//
	//					if (next != null)
	//					{
	//						// -1 is used when the worker has finished
	//						if (next.imageId != -1)
	//						{
	//							// Valid image so store it
	//							nextImages[next.imageId] = next;
	//						}
	//
	//						continue;
	//					}
	//
	//					// Nothing is alive producing images so break
	//					break;
	//				}
	//			}
	//			catch (InterruptedException e)
	//			{
	//
	//			}
	//		}
	//		return image;
	//
	//    	private boolean workersRunning()
	//    	{
	//    		for (BaseWorker worker : workers)
	//    			if (worker.run)
	//    				return true;
	//    		return false;
	//    	}
	//	}

	/**
	 * Creates a background thread to open the images sequentially.
	 * <p>
	 * If all the images are smaller than the buffer limit then a single thread is created to read the raw bytes,
	 * one to
	 * decode the TIFF info and one to read the frames.
	 * <p>
	 * Otherwise a single thread is created to read each image in turn, decode the TIFF info and read the frames.
	 */
	private synchronized void createQueue()
	{
		// Q. What size is optimal?
		rawFramesQueue = new CloseableBlockingQueue<Object>(Runtime.getRuntime().availableProcessors() * 2);

		if (belowBufferLimit() && images.size() > 1)
		{
			// A list of images that may have been read
			//nextImages = new NextSource[images.size()];

			// For now just support a single thread for reading the 
			// raw byte data, decoding and read the TIFF. We can control 
			// the number of images buffered into memory.
			final int nImages = numberOfImages;

			decodeQueue = new CloseableBlockingQueue<NextSource>(nImages);
			readQueue = new CloseableBlockingQueue<NextSource>(nImages);

			workers = new ArrayList<BaseWorker>(3);
			threads = new ArrayList<Thread>(3);
			startWorker(new BufferWorker());
			startWorker(new DecodeWorker());
			startWorker(new ReadWorker());
		}
		else
		{
			// A single worker thread to read the series
			workers = new ArrayList<BaseWorker>(1);
			threads = new ArrayList<Thread>(1);
			startWorker(new TIFFWorker());
		}
	}

	private void startWorker(BaseWorker worker)
	{
		workers.add(worker);
		Thread thread = new Thread(worker);
		threads.add(thread);
		thread.start();
	}

	private boolean belowBufferLimit()
	{
		for (int i = 0; i < imageData.length; i++)
			if (!belowBufferLimit(imageData[i].fileSize))
				return false;
		return true;
	}

	/**
	 * Close the background thread
	 */
	private synchronized void closeQueue()
	{
		if (threads != null)
		{
			// Signal the workers to stop
			for (BaseWorker worker : workers)
				worker.run = false;

			// Close the queues. This will wake any thread waiting for them to have capacity.
			if (decodeQueue != null)
			{
				decodeQueue.close(false);
				readQueue.close(false);
			}
			rawFramesQueue.close(false);

			// Join the threads and then set all to null
			for (Thread thread : threads)
			{
				try
				{
					// This thread will be waiting on imageQueue.put() (which has been cleared)
					// or sourceQueue.take() which contains shutdown signals, it should be finished by now
					thread.join();
				}
				catch (InterruptedException e)
				{
					// Ignore thread errors
					//e.printStackTrace();
				}
			}
			threads.clear();
			threads = null;
			workers.clear();
			workers = null;

			if (decodeQueue != null)
			{
				decodeQueue.close(false);
				readQueue.close(false);
				decodeQueue = null;
				readQueue = null;
			}

			// Do not set this to null as it is used in nextRawFrame()
			rawFramesQueue.close(false);
		}
	}

	/**
	 * Close the queues for the in-memory workflow. This should be called on error as it shuts all the queues.
	 */
	private void closeWorkflowQueues()
	{
		decodeQueue.close(true);
		readQueue.close(true);
		rawFramesQueue.close(true);
	}

	private void storeTiffImage(int imageId, TiffImage image)
	{
		// This could be made optional to save memory
		imageData[imageId].tiffImage = image;
	}

	/**
	 * Sets the origin. This should be used if the source was a crop from a image camera sensor.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 */
	public void setOrigin(int x, int y)
	{
		xOrigin = x;
		yOrigin = y;
	}

	private void setDimensions(int maxx, int maxy)
	{
		width = maxx;
		height = maxy;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#initialiseSequentialReading()
	 */
	@Override
	protected boolean initialiseSequentialRead()
	{
		// This should only be called once by the image source each time the series is opened.
		createQueue();
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#nextRawFrame()
	 */
	@Override
	protected Object nextRawFrame()
	{
		try
		{
			// We can take if closed
			Object pixels = rawFramesQueue.take();
			if (pixels != null)
			{
				return pixels;
			}
		}
		catch (IllegalStateException e)
		{
			// We do not expect these as we use the queue without exceptions when closed
			e.printStackTrace();
		}
		catch (InterruptedException e)
		{
			// Maybe if we forced the threads to stop
			e.printStackTrace();
		}

		// We are here because the pixels were null so shut down sequential reading 
		rawFramesQueue.close(true);

		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#getRawFrame(int)
	 */
	@Override
	protected Object getRawFrame(int frame)
	{
		if (imageSize == null || !isValid(frame))
			return null;

		// Calculate the required image and slice
		int id = Arrays.binarySearch(imageSize, frame);
		if (id < 0)
			id = -(id + 1);
		// Note that frame is 1-based index and the slice is 0-based.
		int slice = (id == 0) ? frame - 1 : frame - imageSize[id - 1] - 1;

		// Return from the cache if it exists
		if (id != lastImageId || lastImage == null)
		{
			if (lastImage != null)
				lastImage.close(true);

			// Used to cache the TIFF info
			lastImage = null;
			if (id < images.size())
			{
				String path = images.get(id);
				if (isTiffSeries)
				{
					// Check the cache
					TiffImage tiffImage = imageData[id].tiffImage;
					if (tiffImage == null)
					{
						tiffImage = openImage(id, path);
					}

					lastImage = tiffImage;

					//try
					//{
					//	if (lastImage == null || lastImage.getSize() == 0)
					//	{
					//		// Not supported - Fall back to IJ objects
					//		ImagePlus imp = IJ.openImage(path);
					//		if (imp != null)
					//		{
					//			lastImage = new ArrayImage(imp);
					//		}
					//	}
					//}
					//catch (Throwable e)
					//{
					//	System.out.println(e.toString()); e.printStackTrace();
					//}
				}
			}
		}
		lastImageId = id;
		if (lastImage != null)
		{
			if (slice < lastImage.size)
			{
				return lastImage.getFrame(slice);
			}
		}
		return null;
	}

	private TiffImage openImage(int id, String path)
	{
		TiffImage tiffImage = null;

		// Open using specialised TIFF reader for better non-sequential support
		SeekableStream ss = createSeekableStream(path, imageData[id].fileSize);
		if (ss != null)
		{
			ExtendedFileInfo[] info = getTiffInfo(ss, path);
			try
			{
				ss.close();
			}
			catch (IOException ioe)
			{
			}

			if (info != null)
			{

				// Store the opened stream as we will use it
				tiffImage = new TiffImage(info, ss);

				storeTiffImage(id, tiffImage);
			}
		}

		if (tiffImage == null)
		{
			// Prevent reading again. Skip the storeTiffImage(...) method 
			// as that may be optional and we don't want to repeat the error.
			imageData[id].tiffImage = new TiffImage();
		}

		return tiffImage;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#isValid(int)
	 */
	@Override
	public boolean isValid(int frame)
	{
		return frame > 0 && frame <= frames;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#toString()
	 */
	@Override
	public String toString()
	{
		String s = super.toString();
		if (!images.isEmpty())
			s += String.format(" (%d images: %s ...)", images.size(), images.get(0));
		return s;
	}

	/**
	 * @return The number of background threads to use for opening images
	 * @deprecated Currently only 1 thread is used for opening images
	 */
	@Deprecated
	public int getNumberOfThreads()
	{
		return numberOfThreads;
	}

	/**
	 * Set the number of background threads to use for opening images.
	 *
	 * @param numberOfThreads
	 *            The number of background threads to use for opening images
	 * @deprecated Currently only 1 thread is used for opening images
	 */
	@Deprecated
	public void setNumberOfThreads(int numberOfThreads)
	{
		this.numberOfThreads = Math.max(1, numberOfThreads);
	}

	/**
	 * Gets the buffer limit for reading TIFF images into memory.
	 *
	 * @return the buffer limit
	 */
	public int getBufferLimit()
	{
		return bufferLimit;
	}

	/**
	 * Sets the buffer limit for reading TIFF images into memory.
	 *
	 * @param bufferLimit
	 *            the new buffer limit
	 */
	public void setBufferLimit(int bufferLimit)
	{
		this.bufferLimit = bufferLimit;
	}

	/**
	 * Gets the number of images to buffer into memory.
	 *
	 * @return the number of images
	 */
	public int getNumberOfImages()
	{
		return numberOfImages;
	}

	/**
	 * Sets the number of images to buffer into memory.
	 *
	 * @param numberOfImages
	 *            the new number of images
	 */
	public void setNumberOfImages(int numberOfImages)
	{
		this.numberOfImages = Math.max(1, numberOfImages);
	}

	/**
	 * Sets the track progress used for monitoring the progress of method execution.
	 *
	 * @param p
	 *            the new track progress
	 */
	public void setTrackProgress(TrackProgress p)
	{
		this.trackProgress = NullTrackProgress.createIfNull(p);
	}

	/**
	 * Gets the a reference to the file info for image n in the series.
	 *
	 * @param index
	 *            the image index
	 * @return the file info
	 */
	public ExtendedFileInfo[] getFileInfo(int index)
	{
		if (imageData != null && index >= 0 && index < imageData.length)
		{
			TiffImage image = imageData[index].tiffImage;
			if (image != null)
			{
				return image.info;
			}
		}
		return null;
	}
}
