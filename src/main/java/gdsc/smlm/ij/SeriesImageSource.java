package gdsc.smlm.ij;

import java.awt.Rectangle;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ArrayBlockingQueue;

import com.thoughtworks.xstream.annotations.XStreamOmitField;

import gdsc.core.ij.SeriesOpener;
import gdsc.core.ij.Utils;
import gdsc.smlm.ij.utils.ImageConverter;
import gdsc.smlm.results.ImageSource;
import ij.IJ;
import ij.ImagePlus;
import ij.io.Opener;

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
	private class NextSource
	{
		final InputStream is;
		final int image;

		public NextSource(InputStream is, int image)
		{
			this.is = is;
			this.image = image;
		}

		public NextSource()
		{
			this(null, -1);
		}
	}

	private class NextImage
	{
		final Object[] imageArray;
		final int imageSize;
		final int image;

		public NextImage(Object[] imageArray, int imageSize, int image)
		{
			this.imageArray = imageArray;
			this.imageSize = imageSize;
			this.image = image;
		}

		public NextImage()
		{
			this(null, 0, -1);
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
					// For a TIFF series ImageJ can open as input stream. We support this by using this
					// thread to read the file from disk sequentially and cache into memory. The images
					// can then be opened by multiple threads without IO contention.
					if (isTiffSeries)
					{
						final String path = images.get(currentImage);
						try
						{
							//System.out.println("Reading " + images.get(currentImage));

							FileInputStream fis = new FileInputStream(path);
							byte[] buf = new byte[fis.available()];
							int read = fis.read(buf);
							fis.close();
							is = new ByteArrayInputStream(buf, 0, read);
						}
						catch (Exception e)
						{
							// TODO - handle appropriately
							// At the moment if we ignore this then the ImageWorker will open the file
							// rather than process from the memory stream.
							System.out.println(e.toString());
						}
					}

					sourceQueue.put(new NextSource(is, currentImage));
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
					ImagePlus imp;
					if (nextSource.is != null)
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

					Object[] imageArray = null;
					int currentImageSize = 0;
					if (imp != null)
					{
						imageArray = imp.getImageStack().getImageArray();
						currentImageSize = imp.getStackSize();

						// Initialise dimensions on the first valid image
						if (width == 0)
							setDimensions(imp.getWidth(), imp.getHeight(), currentImageSize);
					}
					imageQueue.put(new NextImage(imageArray, currentImageSize, currentImage));
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

	@XStreamOmitField
	private int maxz;

	// Used for sequential read
	@XStreamOmitField
	private Object[] imageArray = null;
	@XStreamOmitField
	private int nextImage;
	@XStreamOmitField
	private NextImage[] nextImages;
	@XStreamOmitField
	private int currentSlice;
	@XStreamOmitField
	private int currentImageSize;

	// Used for frame-based read
	@XStreamOmitField
	private Object[] lastImageArray = null;
	@XStreamOmitField
	private int lastImage;
	@XStreamOmitField
	private int lastImageSize;

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
		imageArray = lastImageArray = null;
		nextImage = currentSlice = currentImageSize = lastImage = lastImageSize = 0;
		closeQueue();
	}

	private Object[] getNextImage()
	{
		imageArray = null;
		currentSlice = currentImageSize = 0;
		if (nextImage < nextImages.length)
		{
			try
			{
				for (;;)
				{
					// Images may be out of order due to multiple thread processing.
					// Check if we have processed the next image we need.
					NextImage next = nextImages[nextImage];
					if (next != null)
					{
						// Clear memory
						nextImages[nextImage] = null;

						nextImage++;

						// Check if there is image data. It may be null if the image was invalid
						if (next.imageArray != null)
						{
							imageArray = next.imageArray;
							currentImageSize = next.imageSize;

							//System.out.println("Found image: " + images.get(next.image));

							// Fill cache
							lastImageArray = imageArray;
							lastImageSize = currentImageSize;
							lastImage = next.image;
							return imageArray;
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
						if (next.image != -1)
						{
							// Valid image so store it
							nextImages[next.image] = next;
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
		return imageArray;
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
	 * @see gdsc.smlm.results.ImageSource#nextFrame(java.awt.Rectangle)
	 */
	@Override
	protected float[] nextFrame(Rectangle bounds)
	{
		// Rolling access
		if (imageArray != null)
		{
			// Check if all frames have been accessed in the current image
			if (currentSlice >= currentImageSize)
			{
				// If no more images then return null
				if (getNextImage() == null)
					return null;
			}
			return ImageConverter.getData(imageArray[currentSlice++], width, height, bounds, null);
		}
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#getFrame(int, java.awt.Rectangle)
	 */
	@Override
	protected float[] getFrame(int frame, Rectangle bounds)
	{
		if (maxz == 0 || frame < 1)
			return null;

		// Calculate the required image and slice
		int image = (frame - 1) / maxz;
		int slice = (frame - 1) % maxz;

		// Return from the cache if it exists
		if (image != lastImage || lastImageArray == null)
		{
			lastImageArray = null;
			lastImageSize = 0;
			if (image < images.size())
			{
				ImagePlus imp = IJ.openImage(images.get(image++));
				if (imp != null)
				{
					lastImageArray = imp.getImageStack().getImageArray();
					lastImageSize = imp.getStackSize();
				}
			}
		}
		lastImage = image;
		if (lastImageArray != null)
		{
			if (slice < lastImageSize)
			{
				return ImageConverter.getData(lastImageArray[slice], width, height, bounds, null);
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
}
