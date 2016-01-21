package gdsc.smlm.ij;

import gdsc.smlm.ij.utils.ImageConverter;
import gdsc.smlm.ij.utils.SeriesOpener;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.results.ImageSource;
import ij.IJ;
import ij.ImagePlus;
import ij.io.Opener;

import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.thoughtworks.xstream.XStream;

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
 * 32-bit
 * colour images.
 * <p>
 * Assumes that the width,height,depth dimensions of each file are the same. The depth for the last image can be less
 * (i.e. the last of the series) but the {@link #getFrames()} method will return an incorrect value.
 * <p>
 * Support for the {@link #get(int)} and {@link #get(int, Rectangle)} methods is only provided for TIFF images if the
 * images are stacks.
 */
public class SeriesImageSource extends ImageSource
{
	private class NextImage
	{
		final Object[] imageArray;
		final int imageSize;
		final int image;

		public NextImage(Object[] imageArray, int currentImageSize, int currentImage)
		{
			this.imageArray = imageArray;
			this.imageSize = currentImageSize;
			this.image = currentImage;
		}
	}

	private class ImageWorker implements Runnable
	{
		volatile boolean run = true;

		@Override
		public void run()
		{
			try
			{
				while (run && !nextQueue.isEmpty())
				{
					Integer i = nextQueue.poll();
					if (i == null)
						break;
					final int currentImage = i.intValue();

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
					//System.out.println("Opening " + images.get(currentImage));
					ImagePlus imp = opener.openImage(images.get(currentImage));
					Utils.setShowProgress(true);

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
				System.out.println(e.toString());
			}
			finally
			{
				run = false;

				// Signal that no more images are available
				try
				{
					imageQueue.put(new NextImage(null, 0, -1));
				}
				catch (InterruptedException e)
				{
					// Ignore
				}
			}
		}
	}

	private ArrayList<String> images;

	private int maxz;

	// Used for sequential read
	private Object[] imageArray = null;
	private int nextImage;
	private NextImage[] nextImages;
	private int currentSlice;
	private int currentImageSize;

	// Used for frame-based read
	private Object[] lastImageArray = null;
	private int lastImage;
	private int lastImageSize;

	private boolean logProgress = false;
	private long lastTime = 0;
	private int numberOfThreads = 1;

	private ConcurrentLinkedQueue<Integer> nextQueue;
	private ArrayBlockingQueue<NextImage> imageQueue;
	private ArrayList<ImageWorker> workers = null;
	private ArrayList<Thread> threads = null;

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ResultsSource#init(com.thoughtworks.xstream.XStream)
	 */
	@Override
	public void init(XStream xs)
	{
		super.init(xs);
		xs.omitField(SeriesImageSource.class, "maxz");
		xs.omitField(SeriesImageSource.class, "imageArray");
		xs.omitField(SeriesImageSource.class, "nextImage");
		xs.omitField(SeriesImageSource.class, "nextImages");
		xs.omitField(SeriesImageSource.class, "currentImageSize");
		xs.omitField(SeriesImageSource.class, "lastImageArray");
		xs.omitField(SeriesImageSource.class, "lastImage");
		xs.omitField(SeriesImageSource.class, "lastImageSize");
		xs.omitField(SeriesImageSource.class, "logProgress");
		xs.omitField(SeriesImageSource.class, "lastTime");
		xs.omitField(SeriesImageSource.class, "nextQueue");
		xs.omitField(SeriesImageSource.class, "imageQueue");
		xs.omitField(SeriesImageSource.class, "workers");
		xs.omitField(SeriesImageSource.class, "threads");
	}

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
				images.add(series.getPath() + imageName);
			}
		}
	}

	/**
	 * Create a new image source using the directory containing the images.
	 * <p>
	 * The directory is opened using {@link gdsc.smlm.ij.utils.SeriesOpener }
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
		nextQueue = new ConcurrentLinkedQueue<Integer>();
		for (int i = 0; i < images.size(); i++)
			nextQueue.add(i);
		final int nThreads = numberOfThreads;
		// A blocking queue is used so that the threads do not read too many images in advance
		imageQueue = new ArrayBlockingQueue<NextImage>(nThreads + 2);
		workers = new ArrayList<ImageWorker>(nThreads);
		threads = new ArrayList<Thread>(nThreads);
		for (int i = 0; i < nThreads; i++)
		{
			ImageWorker worker = new ImageWorker();
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
			// Signal the worker to stop
			for (ImageWorker worker : workers)
				worker.run = false;

			// Prevent processing more images
			nextQueue.clear();

			// Ensure any images already waiting on a blocked queue can be added 
			imageQueue.clear();

			// Join the thread and then set all to null
			for (Thread thread : threads)
			{
				try
				{
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
			imageQueue = null;
			nextQueue = null;
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

							//System.out.printf("Found image %d: %s\n", nextImage - 1, images.get(nextImage - 1));

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
