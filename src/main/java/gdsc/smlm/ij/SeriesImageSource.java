package gdsc.smlm.ij;

import java.awt.Rectangle;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ArrayBlockingQueue;

import com.thoughtworks.xstream.annotations.XStreamOmitField;

import gdsc.core.ij.SeriesOpener;
import gdsc.core.ij.Utils;
import gdsc.smlm.results.ImageSource;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.io.FileInfo;
import ij.io.ImageReader;
import ij.io.Opener;
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
 * Represent a series of image files as a results source. Supports all greyscale images. Only processes channel 0 of
 * 32-bit colour images.
 * <p>
 * Assumes that the width,height,depth dimensions of each file are the same. The depth for the last image can be less
 * (i.e. the last of the series) but the {@link #getFrames()} method will return an incorrect value.
 * <p>
 * Support for the {@link #get(int)} and {@link #get(int, Rectangle)} methods is only provided for TIFF images if the
 * images are stacks.
 */
public class SeriesImageSource extends ImageSource
{
	/** The buffer limit for reading TIFF images into memory. Default = 30MB */
	private int bufferLimit = 31457280;

	private class NextSource
	{
		final InputStream is;
		final FileInfo[] info;
		final int image;

		public NextSource(InputStream is, FileInfo[] info, int image)
		{
			this.is = is;
			this.info = info;
			this.image = image;
		}

		public NextSource()
		{
			this(null, null, -1);
		}
	}

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

		public void close()
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
		 * Flag indicating that no errors reading the image have occurred
		 */
		boolean canRead = true;
		/**
		 * The number of frames that have been read from the input stream
		 */
		int frameCount = 0;

		TiffImage(FileInfo[] info)
		{
			super(info[0].width, info[0].height, 0);
			this.info = info;
			fi = info[0];

			// Only support certain types
			if (isSupported(fi.fileType))
			{
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

		@Override
		Object getFrame(int i)
		{
			if (i < frameCount)
			{
				// Non-sequential access has poor support as we close the file, re-open and 
				// skip to the given frame
				close();
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

					// Each image is described by a separate FileInfo object
					if (i == 0)
					{
						// The first frame we know the exact offset
						skip = fi.getOffset();
					}
					else
					{
						// Skipping ahead from the current file location
						// (including the bytes that were read last time) 
						skip = info[i].getOffset() - (info[frameCount - 1].getOffset() + bytesPerFrame);
						if (skip < 0L)
							throw new IllegalStateException("Bad TIFF offset " + skip);
						fi.stripOffsets = info[i].stripOffsets;
						fi.stripLengths = info[i].stripLengths;
					}
				}

				// Store the number of frames that have been read
				frameCount = i + 1;

				return reader.readPixels(is, skip);
			}
			catch (Exception e)
			{
				System.out.println(e.toString());
				canRead = false;
			}

			return null;
		}

		private boolean openInputStream()
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
		public InputStream createInputStream(FileInfo fi) throws IOException
		{
			// Close the old input
			close();

			InputStream is = null;
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

		@Override
		public void close()
		{
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
				frameCount = 0;
			}
		}
	}

	private class NextImage
	{
		final Image image;
		final int imageId;

		public NextImage(Image image, int imageId)
		{
			this.image = image;
			this.imageId = imageId;
		}

		public NextImage()
		{
			this(null, -1);
		}
	}

	/**
	 * Add source images to the queue to be read.
	 * <p>
	 * If the series is all TIFF images then we can open the file and read it into memory.
	 */
	private class SourceWorker implements Runnable
	{
		volatile boolean run = true;

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		public void run()
		{
			try
			{
				for (int currentImage = 0; run && currentImage < images.size(); currentImage++)
				{
					InputStream is = null;
					FileInfo[] info = null;
					// For a TIFF series ImageJ can open as input stream. We support this by using this
					// thread to read the file from disk sequentially and cache into memory. The images
					// can then be opened by multiple threads without IO contention.
					if (isTiffSeries)
					{
						final String path = images.get(currentImage);
						FileInputStream fis = null;
						try
						{
							//System.out.println("Reading " + images.get(currentImage));

							fis = new FileInputStream(path);
							// Don't buffer massive images into memory
							if (fis.getChannel().size() <= bufferLimit)
							{
								byte[] buf = new byte[fis.available()];
								int read = fis.read(buf);
								is = new ByteArrayInputStream(buf, 0, read);
							}
							else
							{
								// This is a big image. We support reading the TIFF sequentially.
								TiffDecoder td = new TiffDecoder(fis, path);
								info = td.getTiffInfo();
							}
						}
						catch (Throwable e)
						{
							// TODO - handle appropriately
							// At the moment if we ignore this then the ImageWorker will open the file
							// rather than process from the memory stream.
							System.out.println(e.toString());
						}
						finally
						{
							if (fis != null)
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

					sourceQueue.put(new NextSource(is, info, currentImage));
				}
			}
			catch (InterruptedException e)
			{
				// This is from the queue put method, possibly an interrupt on the queue or thread? 
				System.out.println(e.toString());
			}

			// Add jobs to shutdown all the workers
			try
			{
				for (int i = 0; run && i < workers.size(); i++)
					sourceQueue.put(new NextSource());
			}
			catch (InterruptedException e)
			{
				// This is from the queue put method, possibly an interrupt on the queue or thread? 
				// TODO - How should this be handled?
				System.out.println(e.toString());
			}

			run = false;
		}
	}

	private class ImageWorker implements Runnable
	{
		final int id;
		volatile boolean run = true;

		public ImageWorker(int id)
		{
			this.id = id;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		public void run()
		{
			try
			{
				while (run)
				{
					NextSource nextSource = sourceQueue.take();
					if (nextSource == null || !run)
						break;
					final int currentImage = nextSource.image;
					if (currentImage == -1)
						break;

					// Open the image
					Opener opener = new Opener();
					opener.setSilentMode(true);
					Utils.setShowProgress(false);
					if (logProgress)
					{
						long time = System.currentTimeMillis();
						if (time - lastTime > 500)
						{
							lastTime = time;
							IJ.log("Opening " + images.get(currentImage));
						}
					}
					Image image = null;
					ImagePlus imp = null;
					// The TIFF info is used when a very large TIFF image
					if (nextSource.info != null)
					{
						image = new TiffImage(nextSource.info);
					}
					else if (nextSource.is != null)
					{
						//System.out.println(id + ": Processing " + images.get(currentImage));
						imp = opener.openTiff(nextSource.is, images.get(currentImage));
					}
					else
					{
						//System.out.println(id + ": Opening " + images.get(currentImage));
						imp = opener.openImage(images.get(currentImage));
					}
					Utils.setShowProgress(true);

					//System.out.println(id + ": Opened " + images.get(currentImage));

					if (image == null)
					{
						if (imp != null)
						{
							image = new ArrayImage(imp);
						}
						else
						{
							image = new ArrayImage();
						}
					}

					if (image.getSize() != 0)
					{
						if (width == 0)
						{
							// Initialise dimensions on the first valid image
							setDimensions(image.getWidth(), image.getHeight(), image.getSize());
						}
						else
						{
							// Check dimensions
							if (image.getWidth() != getWidth() || image.getHeight() != getHeight())
							{
								// Return no image data
								image = new ArrayImage();
							}
						}
					}

					imageQueue.put(new NextImage(image, currentImage));
				}
			}
			catch (Exception e)
			{
				// TODO - handle appropriately 
				System.out.println(id + ": " + e.toString());
			}
			finally
			{
				run = false;

				// Signal that no more images are available
				imageQueue.offer(new NextImage());
			}
		}
	}

	private ArrayList<String> images;
	public final boolean isTiffSeries;

	/**
	 * Used to cache the TIFF info for non-sequential reading
	 */
	@XStreamOmitField
	private TiffImage[] tiffImages;

	@XStreamOmitField
	private int maxz;

	// Used for sequential read
	@XStreamOmitField
	private Image image = null;
	@XStreamOmitField
	private int nextImageId;
	@XStreamOmitField
	private NextImage[] nextImages;
	@XStreamOmitField
	private int currentSlice;

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

	@XStreamOmitField
	private ArrayBlockingQueue<NextSource> sourceQueue;
	@XStreamOmitField
	private ArrayBlockingQueue<NextImage> imageQueue;

	// Used to queue the files from disk. This is sequential so only one thread is required.
	@XStreamOmitField
	private SourceWorker sourceWorker = null;
	@XStreamOmitField
	private Thread sourceThread = null;

	// Used to process the files into images
	@XStreamOmitField
	private ArrayList<ImageWorker> workers = null;
	@XStreamOmitField
	private ArrayList<Thread> threads = null;

	/**
	 * Create a new image source using the given image series
	 * 
	 * @param name
	 * @param path
	 */
	public SeriesImageSource(String name, SeriesOpener series)
	{
		super(name);
		images = new ArrayList<String>();
		if (series != null)
		{
			for (String imageName : series.getImageList())
			{
				images.add(new File(series.getPath(), imageName).getPath());
			}
		}
		isTiffSeries = isTiffSeries();
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
		images = new ArrayList<String>();
		images.addAll(Arrays.asList(filenames));
		isTiffSeries = isTiffSeries();
	}

	private boolean isTiffSeries()
	{
		Opener opener = new Opener();
		for (String path : images)
			if (opener.getFileType(path) != Opener.TIFF)
				return false;
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#openSource()
	 */
	@Override
	public boolean openSource()
	{
		// reset
		close();

		if (images.isEmpty())
			return false;

		// Create the queue for loading the images
		createQueue();

		if (isTiffSeries && (tiffImages == null || tiffImages.length != images.size()))
			tiffImages = new TiffImage[images.size()];

		return getNextImage() != null;
	}

	/**
	 * Creates a background thread to open the images sequentially
	 */
	private void createQueue()
	{
		nextImages = new NextImage[images.size()];
		final int nThreads = numberOfThreads;

		// A blocking queue is used so that the threads do not read too many images in advance
		sourceQueue = new ArrayBlockingQueue<NextSource>(nThreads);

		// Start a thread to queue up the images
		sourceWorker = new SourceWorker();
		sourceThread = new Thread(sourceWorker);
		sourceThread.start();

		// A blocking queue is used so that the threads do not read too many images in advance
		imageQueue = new ArrayBlockingQueue<NextImage>(nThreads + 2);
		workers = new ArrayList<ImageWorker>(nThreads);
		threads = new ArrayList<Thread>(nThreads);
		for (int i = 0; i < nThreads; i++)
		{
			ImageWorker worker = new ImageWorker(i + 1);
			workers.add(worker);
			Thread thread = new Thread(worker);
			threads.add(thread);
			thread.start();
		}
	}

	/**
	 * Close the background thread
	 */
	private void closeQueue()
	{
		if (threads != null)
		{
			// Signal the workers to stop
			sourceWorker.run = false;
			for (ImageWorker worker : workers)
				worker.run = false;

			// Prevent processing more source images
			sourceQueue.clear();

			// Ensure any images already waiting on a blocked queue can be added 
			imageQueue.clear();

			// Send shutdown signals to anything for another source to process
			for (int i = 0; i < workers.size(); i++)
				sourceQueue.offer(new NextSource());

			// Join the threads and then set all to null
			try
			{
				// This thread will be waiting on sourceQueue.put() so it should be finished by now
				sourceThread.join();
			}
			catch (InterruptedException e)
			{
				// Ignore thread errors
				//e.printStackTrace();
			}

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
			sourceThread = null;
			threads.clear();
			threads = null;
			workers.clear();
			workers = null;
			imageQueue = null;
			sourceQueue = null;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#close()
	 */
	public void close()
	{
		setDimensions(0, 0, 0);
		image = lastImage = null;
		nextImageId = currentSlice = lastImageId = 0;
		closeQueue();
	}

	private Image getNextImage()
	{
		image = null;
		currentSlice = 0;
		if (nextImageId < nextImages.length)
		{
			try
			{
				for (;;)
				{
					// Images may be out of order due to multiple thread processing.
					// Check if we have processed the next image we need.
					NextImage next = nextImages[nextImageId];
					if (next != null)
					{
						// Clear memory
						nextImages[nextImageId] = null;

						nextImageId++;

						// Check if there is image data. It may be null if the image was invalid
						if (next.image != null)
						{
							image = next.image;

							//System.out.println("Found image: " + images.get(next.image));

							// Fill cache
							lastImage = image;
							lastImageId = next.imageId;
							// The cache is used for non-sequential reading. To prevent 
							// memory usage during sequential reading only cache the first 
							// one as this is generated when the source is opened and it is 
							// unclear if the source is to be used non-sequentially. 
							if (isTiffSeries && lastImageId == 0 && image instanceof TiffImage)
								storeTiffImage(lastImageId, (TiffImage) image);
							return image;
						}
						else
						{
							// Process the next stored image
							continue;
						}
					}

					// We are still awaiting the next image.
					// Get the images processed by the worker threads.
					next = imageQueue.poll();

					// If there is nothing then check if any workers are alive
					if (next == null && workersRunning())
					{
						// This will block until something produces an image
						next = imageQueue.take();
					}

					if (next != null)
					{
						// -1 is used when the worker has finished
						if (next.imageId != -1)
						{
							// Valid image so store it
							nextImages[next.imageId] = next;
						}

						continue;
					}

					// Nothing is alive producing images so break
					break;
				}
			}
			catch (InterruptedException e)
			{

			}
		}
		return image;
	}

	private void storeTiffImage(int imageId, TiffImage image)
	{
		// This could be made optional to save memory
		tiffImages[imageId] = (TiffImage) image;
	}

	private boolean workersRunning()
	{
		for (ImageWorker worker : workers)
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

	private void setDimensions(int maxx, int maxy, int maxz)
	{
		width = maxx;
		height = maxy;
		this.maxz = maxz;
		// This will be wrong if the stacks are different sizes or images are missing
		frames = maxz * images.size();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#nextRawFrame()
	 */
	@Override
	protected Object nextRawFrame()
	{
		// Rolling access
		if (image != null)
		{
			// Check if all frames have been accessed in the current image
			if (currentSlice >= image.size)
			{
				image.close();
				// If no more images then return null
				if (getNextImage() == null)
					return null;
			}
			return image.getFrame(currentSlice++);
		}
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
		if (maxz == 0 || frame < 1)
			return null;

		// Calculate the required image and slice
		int id = (frame - 1) / maxz;
		int slice = (frame - 1) % maxz;

		// Return from the cache if it exists
		if (id != lastImageId || lastImage == null)
		{
			if (lastImage != null)
				lastImage.close();

			// Used to cache the TIFF info
			lastImage = null;
			if (id < images.size())
			{
				String path = images.get(id);
				if (isTiffSeries)
				{
					// Check the cache
					TiffImage tiffImage = tiffImages[id];
					if (tiffImage == null)
					{
						// Open using specialised TIFF reader for better non-sequential support
						FileInputStream fis = null;
						try
						{
							fis = new FileInputStream(path);
							TiffDecoder td = new TiffDecoder(fis, path);
							FileInfo[] info = td.getTiffInfo();
							tiffImage = new TiffImage(info);

							storeTiffImage(id, tiffImage);
						}
						catch (Throwable e)
						{
							System.out.println(e.toString());
							// Prevent reading again. Skip the storeTiffImage(...) method 
							// as that may be optional and we don't want to repeat the error.
							tiffImages[id] = new TiffImage();
						}
						finally
						{
							if (fis != null)
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

					lastImage = tiffImage;
					try
					{
						if (lastImage == null || lastImage.getSize() == 0)
						{
							// Not supported - Fall back to IJ objects
							ImagePlus imp = IJ.openImage(path);
							if (imp != null)
							{
								lastImage = new ArrayImage(imp);
							}
						}
					}
					catch (Throwable e)
					{
						System.out.println(e.toString());
					}
				}
				else
				{
					ImagePlus imp = IJ.openImage(path);
					if (imp != null)
					{
						lastImage = new ArrayImage(imp);
					}
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
