package gdsc.smlm.ij;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import com.thoughtworks.xstream.annotations.XStreamOmitField;

import gdsc.core.generics.CloseableBlockingQueue;
import gdsc.core.ij.SeriesOpener;
import gdsc.smlm.results.ImageSource;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.io.CustomTiffDecoder;
import ij.io.FileInfo;
import ij.io.ImageReader;
import ij.io.Opener;
import ij.io.RandomAccessStream;
import ij.io.TiffDecoder;

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
		final FileInfo[] info;
		FileInfo fi; // Pointer to the file info used by the ImageReader
		boolean contiguous;
		final long bytesPerFrame;
		ImageReader reader = null;
		InputStream is = null;
		/**
		 * A reference to a stream with an underlying ByteArrayRandomAccessFile. This is effectively the file contents
		 * held in memory.
		 */
		RandomAccessStream ras = null;
		boolean inMemory = false;
		/**
		 * Flag indicating that no errors reading the image have occurred
		 */
		boolean canRead = true;
		/**
		 * The number of frames that have been read from the input stream
		 */
		int frameCount = 0;

		TiffImage(FileInfo[] info, RandomAccessStream ras, boolean inMemory)
		{
			super(info[0].width, info[0].height, 0);
			this.info = info;
			fi = info[0];

			// Only support certain types
			if (isSupported(fi.fileType))
			{
				// We use the opened RandomAccessStream.
				// This may be in-memory data or may be from a random access file.
				if (ras != null)
				{
					try
					{
						ras.seek(0);
						this.ras = ras;
						is = ras;

						// Store if the stream contains in-memory data						
						this.inMemory = inMemory;
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

		private int getBytesPerFrame(int fileType)
		{
			switch (fileType)
			{
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
		 * Check the file info are all the same size and type and have just 1 image
		 * 
		 * @param info
		 * @return True if all the same size and type
		 */
		boolean allSameSizeAndType(FileInfo[] info)
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
				// We can read using only the first FileInfo object if this is contiguous
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

					// Each image offset is described by a separate FileInfo object
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
					// Read using the first FileInfo object

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

					// Each image offset is described by a separate FileInfo object
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

		/** Returns an InputStream for the image described by this FileInfo. */
		public FileInputStream createInputStream(FileInfo fi) throws IOException
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
		boolean validateFileInfo(File f, FileInfo fi)
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
			if (fi.fileType == FileInfo.BITMAP || fi.compression != FileInfo.COMPRESSION_NONE)
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

			if (ras != null && !freeResources)
			{
				// We can seek to the start
				try
				{
					ras.seek(0);
					return;
				}
				catch (IOException e)
				{
					// Fall through to close the resources
				}
			}

			// This is done when sequentially reading so we clear the memory
			inMemory = false;
			ras = null;

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

	/**
	 * Extend RandomAccessFile to use a byte buffer held in memory. The source file is not manipulated but used to allow
	 * construction.
	 */
	public class ByteArrayRandomAccessFile extends RandomAccessFile
	{
		int p = 0;
		byte[] bytes;
		final int length;

		public ByteArrayRandomAccessFile(byte[] bytes, File file) throws IOException
		{
			super(file, "r");
			this.bytes = bytes;
			length = bytes.length;

			// Close the file again!
			super.close();
		}

		// Override the methods that are not final

		@Override
		public int read() throws IOException
		{
			if (p < length)
				return bytes[p++] & 0xff;
			return -1;
		}

		@Override
		public int read(byte[] b, int off, int len) throws IOException
		{
			if (p < length)
			{
				if (len > 0)
				{
					int size = (p + len <= length) ? len : length - p;
					System.arraycopy(bytes, p, b, off, size);
					p += size;
					return size;
				}
				throw new IOException("No length specified");
			}
			return -1;
		}

		@Override
		public void write(byte[] b) throws IOException
		{
			throw new IOException("Not supported");
		}

		@Override
		public void write(byte[] b, int off, int len) throws IOException
		{
			throw new IOException("Not supported");
		}

		@Override
		public long getFilePointer() throws IOException
		{
			return p;
		}

		@Override
		public void seek(long pos) throws IOException
		{
			if (pos < 0)
				throw new IOException("Negative position");
			// Allow seek to the end
			p = (pos > length) ? length : (int) pos;
		}

		@Override
		public long length() throws IOException
		{
			return length;
		}

		@Override
		public void setLength(long newLength) throws IOException
		{
			throw new IOException("Not supported");
		}

		@Override
		public void close() throws IOException
		{
			// Do nothing
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
					//System.out.println("Reading " + images.get(currentImage));

					RandomAccessStream ras = createRandomAccessStream(path);
					if (ras == null)
					{
						error = true;
						break;
					}

					if (logProgress)
					{
						long time = System.currentTimeMillis();
						if (time - lastTime > 500)
						{
							lastTime = time;
							IJ.log("Reading TIFF info " + path);
						}
					}

					// Check the cache to avoid a re-read. This will always be the case
					// for the first image as that is opened in openSource().
					TiffImage image = imageData[currentImage].tiffImage;
					if (image == null)
					{
						FileInfo[] info = getTiffInfo(ras, path);
						if (info == null)
						{
							error = true;
							break;
						}
						image = new TiffImage(info, null, false);
						storeTiffImage(currentImage, image);
					}

					// Check it is the expected size
					if (image.size != getImageSize(currentImage))
					{
						error = true;
						break;
					}

					if (logProgress)
					{
						long time = System.currentTimeMillis();
						if (time - lastTime > 500)
						{
							lastTime = time;
							IJ.log("Reading TIFF " + path);
						}
					}

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
							rawFrames.put(pixels);
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
				error = true;
			}

			rawFrames.close(error);

			run = false;
		}

	}

	private class NextSource
	{
		final byte[] buffer;
		final int imageId;
		FileInfo[] info;
		Image image;

		public NextSource(byte[] buffer, int imageId)
		{
			this.buffer = buffer;
			this.imageId = imageId;
		}

		public void setFileInfo(FileInfo[] info)
		{
			this.info = info;
		}

		public void setImage(Image image)
		{

		}

		public NextSource()
		{
			this(null, -1);
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
	//					FileInfo[] info = null;
	//					RandomAccessStream ras = null;
	//
	//					final String path = images.get(currentImage);
	//					//System.out.println("Reading " + images.get(currentImage));
	//
	//					// This may contain a custom RandomAccessStream that wraps an in-memory RandomAccessFile 
	//					ras = createRandomAccessStream(path);
	//
	//					if (ras != null)
	//					{
	//						TiffDecoder td = new CustomTiffDecoder(ras, path);
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
	//							// This will close the random access stream.
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
	//					sourceQueue.put(new NextSource(info, (inMemory[0]) ? ras : null, currentImage));
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
	//				System.out.println(e.toString());
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
	//				System.out.println(e.toString());
	//			}
	//
	//			run = false;
	//		}
	//	}

	/**
	 * Creates the random access stream.
	 *
	 * @param path
	 *            the path
	 * @param inMemory
	 *            the in memory
	 * @param size
	 *            the size of the file
	 * @return the random access stream
	 */
	private RandomAccessStream createRandomAccessStream(String path, boolean[] inMemory, long size)
	{
		File file = new File(path);

		// Don't buffer massive images into memory
		if (belowBufferLimit(size))
		{
			RandomAccessStream ras = readBytes(file, size);
			if (ras != null)
			{
				inMemory[0] = true;
				return ras;
			}
		}

		try
		{
			return new RandomAccessStream(new RandomAccessFile(file, "r"));
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
	 * Creates the random access stream.
	 *
	 * @param path
	 *            the path
	 * @return the random access stream
	 */
	private RandomAccessStream createRandomAccessStream(String path)
	{
		File file = new File(path);
		try
		{
			return new RandomAccessStream(new RandomAccessFile(file, "r"));
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

	private FileInfo[] getTiffInfo(RandomAccessStream ras, String path)
	{
		TiffDecoder td = new CustomTiffDecoder(ras, path);
		try
		{
			//td.enableDebugging();
			// This will close the random access stream.
			FileInfo[] info = td.getTiffInfo();
			//if (info != null)
			//{
			//	System.out.println(info[0].debugInfo);
			//	
			//	// This contains OME TIFF metadata as serialised JSON when using 
			//	// MicroManager to save the TIFF
			//	System.out.println(info[0].info);
			//	
			//	ij.io.FileOpener fo = new ij.io.FileOpener(info[0]);
			//	java.util.Properties p = fo.decodeDescriptionString(info[0]);
			//	System.out.println(p);
			//	//double ox = 0, oy = 0;
			//	//if (p.containsKey("xorigin"))
			//	//	ox = Double.parseDouble(p.getProperty("xorigin"));
			//	//if (p.containsKey("yorigin"))
			//	//	oy = Double.parseDouble(p.getProperty("yorigin"));
			//	//// Should the origin be converted by the units?
			//	//setOrigin((int)ox, (int)oy);
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
				ras.close();
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

	private RandomAccessStream readBytes(File file, long size)
	{
		FileInputStream fis = openStream(file);
		if (fis != null)
		{
			try
			{
				byte[] buf = new byte[(int) size];
				int read = fis.read(buf);
				if (read == size)
					return new RandomAccessStream(new ByteArrayRandomAccessFile(buf, file));
			}
			catch (IOException e)
			{
				// At the moment if we ignore this then the ImageWorker will open the file
				// rather than process from the memory stream.
				System.out.println(e.toString());
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
	//					// Only support TIFF images with pre-read FileInfo
	//					Image image = null;
	//					if (nextSource.info != null)
	//					{
	//						image = new TiffImage(nextSource.info, nextSource.ras, nextSource.ras != null);
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

	// Used for sequential read
	@XStreamOmitField
	private Image image = null;
	@XStreamOmitField
	private int nextImageId;
	@XStreamOmitField
	private NextSource[] nextImages;
	@XStreamOmitField
	private int currentSlice;
	@XStreamOmitField
	private byte sequentialReadStatus;
	
	private static final byte CLOSED = 0;
	private static final byte OPEN = 1;
	private static final byte RUNNING = 2;

	// Used for frame-based read
	@XStreamOmitField
	private Image lastImage = null;
	@XStreamOmitField
	private int lastImageId;

	@XStreamOmitField
	private boolean logProgress = false;
	@XStreamOmitField
	private long lastTime = 0;
	private int numberOfThreads = 1;

	//	@XStreamOmitField
	//	private ArrayBlockingQueue<NextSource> sourceQueue;
	//	@XStreamOmitField
	//	private ArrayBlockingQueue<NextImage> imageQueue;

	// Used to process the files into images
	@XStreamOmitField
	private ArrayList<BaseWorker> workers = null;
	@XStreamOmitField
	private ArrayList<Thread> threads = null;

	/**
	 * Used for sequential reading to queue the raw frames
	 */
	@XStreamOmitField
	private CloseableBlockingQueue<Object> rawFrames = null;

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
		if (imageSize == null)
		{
			// All images are TIFF. Get the size of each and count the total frames.
			imageSize = new int[images.size()];
			imageData = new ImageData[images.size()];
			String[] names = new String[images.size()];
			frames = 0;
			int ok = 0;

			for (int i = 0; i < names.length; i++)
			{
				String path = images.get(i);

				RandomAccessStream ras = null;
				try
				{
					File file = new File(path);

					// Get the size of each file so we can determine if 
					// they can fit into memory. We only use pre-loading for
					// sequential reading if all images fit into memory.
					long size = getSize(file);

					ras = new RandomAccessStream(new RandomAccessFile(file, "r"));
					CustomTiffDecoder td = new CustomTiffDecoder(ras, path);

					int n = td.getNumberOfImages();
					System.out.printf("%s = %d\n", path, n);
					if (n <= 0)
					{
						// No TIFF images. This will break the non-sequential support
						// using the cumulative size array so remove the image. 
						continue;
					}

					frames += n;
					imageSize[ok] = frames;
					imageData[ok] = new ImageData(size);
					names[ok] = path;
					ok++;
				}
				catch (Throwable e)
				{
				}
				finally
				{
					if (ras != null)
					{
						try
						{
							ras.close();
						}
						catch (IOException e1)
						{
						}
					}
				}
			}

			if (ok < images.size())
			{
				imageSize = Arrays.copyOf(imageSize, ok);
				imageData = Arrays.copyOf(imageData, ok);
				images.clear();
				images.addAll(Arrays.asList(names));
			}
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
	public boolean openSource()
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
		Object pixels = getRawFrame(1);
		if (pixels != null)
		{
			setDimensions(imageData[0].tiffImage.width, imageData[0].tiffImage.height);
			// Special flag for sequential reading
			sequentialReadStatus = OPEN;
			return true;
		}
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#close()
	 */
	public void close()
	{
		if (threads != null)
			closeQueue();
		//setDimensions(0, 0, 0);
		if (image != null)
			image.close(true);
		if (lastImage != null)
			lastImage.close(true);
		image = lastImage = null;
		nextImageId = currentSlice = lastImageId = 0;
		sequentialReadStatus = CLOSED;
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
		if (sequentialReadStatus != OPEN)
			return;
		sequentialReadStatus = RUNNING;
			
		// Q. What size is optimal?
		rawFrames = new CloseableBlockingQueue<Object>(Runtime.getRuntime().availableProcessors() * 2);

		if (false && belowBufferLimit() && images.size() > 1)
		{
			// A list of images that may have been read
			nextImages = new NextSource[images.size()];
			final int nThreads = numberOfThreads;

			//			// A blocking queue is used so that the threads do not read too many images in advance
			//			sourceQueue = new ArrayBlockingQueue<NextSource>(nThreads);
			//
			//			// Start a thread to queue up the images
			//			sourceWorker = new SourceWorker();
			//			sourceThread = new Thread(sourceWorker);
			//			sourceThread.start();
			//
			//			// A blocking queue is used so that the threads do not read too many images in advance
			//			imageQueue = new ArrayBlockingQueue<NextImage>(nThreads + 2);
			//			workers = new ArrayList<BaseWorker>(nThreads);
			//			threads = new ArrayList<Thread>(nThreads);
			//			for (int i = 0; i < nThreads; i++)
			//			{
			//				ImageWorker worker = new ImageWorker(i + 1);
			//				workers.add(worker);
			//				Thread thread = new Thread(worker);
			//				threads.add(thread);
			//				thread.start();
			//			}
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

			//			// Prevent processing more source images
			//			sourceQueue.clear();
			//
			//			// Ensure any images already waiting on a blocked queue can be added 
			//			imageQueue.clear();
			//
			//			// Send shutdown signals to anything for another source to process
			//			for (int i = 0; i < workers.size(); i++)
			//				sourceQueue.offer(new NextSource());

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

			// Do not set this to null as it is used in nextRawFrame()
			rawFrames.close(false);

			//sourceThread = null;
			//imageQueue = null;
			//sourceQueue = null;
		}
	}

	private void storeTiffImage(int imageId, TiffImage image)
	{
		// This could be made optional to save memory
		imageData[imageId].tiffImage = image;
	}

	private boolean workersRunning()
	{
		for (BaseWorker worker : workers)
			if (worker.run)
				return true;
		return false;
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
	 * @see gdsc.smlm.results.ImageSource#nextRawFrame()
	 */
	@Override
	protected Object nextRawFrame()
	{
		if (sequentialReadStatus == OPEN)
			createQueue();

		if (sequentialReadStatus != RUNNING)
			return null;

		try
		{
			// We can take if closed
			Object pixels = rawFrames.take();
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
		sequentialReadStatus = CLOSED;
		rawFrames.close(true);
		rawFrames = null;

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
						// Open using specialised TIFF reader for better non-sequential support
						RandomAccessStream ras = null;
						try
						{
							boolean[] inMemory = new boolean[1];
							ras = createRandomAccessStream(path, inMemory, imageData[id].fileSize);
							CustomTiffDecoder td = new CustomTiffDecoder(ras, path);
							FileInfo[] info = td.getTiffInfo();

							// Store the opened stream as we will use it
							tiffImage = new TiffImage(info, ras, inMemory[0]);

							storeTiffImage(id, tiffImage);
						}
						catch (Throwable e)
						{
							System.out.println(e.toString());

							// Close resources
							if (ras != null)
							{
								try
								{
									ras.close();
								}
								catch (IOException ioe)
								{
								}
							}

							// Prevent reading again. Skip the storeTiffImage(...) method 
							// as that may be optional and we don't want to repeat the error.
							imageData[id].tiffImage = new TiffImage();
						}
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
					//	System.out.println(e.toString());
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
	 * @return If true send progress to the ImageJ log
	 */
	public boolean isLogProgress()
	{
		return logProgress;
	}

	/**
	 * @param logProgress
	 *            Set to true to send progress to the ImageJ log
	 */
	public void setLogProgress(boolean logProgress)
	{
		this.logProgress = logProgress;
	}

	/**
	 * @return The number of background threads to use for opening images
	 */
	public int getNumberOfThreads()
	{
		return numberOfThreads;
	}

	/**
	 * Set the number of background threads to use for opening images
	 * 
	 * @param numberOfThreads
	 *            The number of background threads to use for opening images
	 */
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
}
