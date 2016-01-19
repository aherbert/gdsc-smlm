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
		Object[] imageArray = null;
		int currentImageSize;
		int currentImage;

		public NextImage(Object[] imageArray, int currentImageSize, int currentImage)
		{
			this.imageArray = imageArray;
			this.currentImageSize = currentImageSize;
			this.currentImage = currentImage;
		}
	}

	private class ImageWorker implements Runnable
	{
		boolean run = true;

		@Override
		public void run()
		{
			try
			{
				for (int currentImage=0; run && currentImage < images.size(); currentImage++)
				{
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
						
						imageQueue.put(new NextImage(imageArray, currentImageSize, currentImage));
					}
				}
				
				// Signal that no more images are available
				imageQueue.put(new NextImage(null, 0, 0));
			}
			catch (Exception e)
			{
				System.out.println(e.toString());
				throw new RuntimeException(e);
			}
			finally
			{
			}
		}
	}

	private ArrayList<String> images;

	private int maxz;

	// Used for sequential read
	private Object[] imageArray = null;
	private int currentSlice;
	private int currentImageSize;

	// Used for frame-based read
	private Object[] lastImageArray = null;
	private int lastImage;
	private int lastImageSize;

	private boolean logProgress = false;
	private long lastTime = 0;

	private ArrayBlockingQueue<NextImage> imageQueue;
	private ImageWorker worker = null;
	private Thread thread = null;

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
		xs.omitField(SeriesImageSource.class, "currentSlice");
		xs.omitField(SeriesImageSource.class, "currentImageSize");
		xs.omitField(SeriesImageSource.class, "lastImageArray");
		xs.omitField(SeriesImageSource.class, "lastImage");
		xs.omitField(SeriesImageSource.class, "lastImageSize");
		xs.omitField(SeriesImageSource.class, "logProgress");
		xs.omitField(SeriesImageSource.class, "lastTime");
		xs.omitField(SeriesImageSource.class, "worker");
		xs.omitField(SeriesImageSource.class, "thread");
		xs.omitField(SeriesImageSource.class, "imageQueue");
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
		this(name, new SeriesOpener(directory, false));
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
		imageQueue = new ArrayBlockingQueue<NextImage>(2);
		worker = new ImageWorker();
		thread = new Thread(worker);
		thread.start();
	}

	/**
	 * Close the background thread
	 */
	private void closeQueue()
	{
		if (thread != null)
		{
			// Signal the worker to stop
			worker.run = false;
			
			// Ensure any images already waiting on a blocked queue can be added 
			imageQueue.clear();

			// Join the thread and then set all to null
			try
			{
				thread.join();
			}
			catch (InterruptedException e)
			{
				// Ignore thread errors
				//e.printStackTrace();
			}
			thread = null;
			worker = null;
			imageQueue = null;
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
		currentSlice = currentImageSize = lastImage = lastImageSize = 0;
		closeQueue();
	}

	private Object[] getNextImage()
	{
		imageArray = null;
		currentSlice = currentImageSize = 0;
		try
		{
			final NextImage next = imageQueue.take();

			if (next != null)
			{
				imageArray = next.imageArray;
				currentImageSize = next.currentImageSize;

				// Fill cache
				lastImageArray = imageArray;
				lastImageSize = currentImageSize;
				lastImage = next.currentImage;
			}
		}
		catch (InterruptedException e)
		{

		}
		return imageArray;
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
}
