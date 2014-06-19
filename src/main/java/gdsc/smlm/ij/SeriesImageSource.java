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
	private ArrayList<String> images;

	private int maxz;

	// Used for sequential read
	private int currentSlice;
	private int currentImage;
	private Object[] imageArray = null;
	private int currentImageSize;

	// Used for frame-based read
	private int lastImage;
	private Object[] lastImageArray = null;
	private int lastImageSize;
	
	private boolean logProgress = false;

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
		xs.omitField(SeriesImageSource.class, "currentSlice");
		xs.omitField(SeriesImageSource.class, "currentImage");
		xs.omitField(SeriesImageSource.class, "currentImageSize");
		xs.omitField(SeriesImageSource.class, "imageArray");
		xs.omitField(SeriesImageSource.class, "lastImage");
		xs.omitField(SeriesImageSource.class, "lastImageArray");
		xs.omitField(SeriesImageSource.class, "logProgress");
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
	 * @see gdsc.smlm.results.ResultsSource#open()
	 */
	@Override
	public boolean open()
	{
		// reset
		setDimensions(0, 0, 0);
		imageArray = lastImageArray = null;
		currentImageSize = lastImageSize = 0;

		if (images.isEmpty())
			return false;

		// Open the first image to initialise dimensions
		currentImage = 0;
		return getNextImage() != null;
	}

	private Object[] getNextImage()
	{
		currentImageSize = currentSlice = 0;
		if (currentImage >= images.size())
			return null;
		
		// Disable the progress bar when opening files
		Opener opener = new Opener();
		opener.setSilentMode(true);
		Utils.setShowProgress(false);
		if (logProgress)
			IJ.log("Opening " + images.get(currentImage));
		ImagePlus imp = opener.openImage(images.get(currentImage++));
		Utils.setShowProgress(true);
		
		//ImagePlus imp = IJ.openImage(images.get(currentImage++));
		
		if (imp == null)
			return null;
		imageArray = imp.getImageStack().getImageArray();
		currentImageSize = imp.getStackSize();

		// Initialise on the first image
		if (currentImage == 1)
			setDimensions(imp.getWidth(), imp.getHeight(), currentImageSize);

		// Fill cache if empty
		if (lastImageArray == null)
		{
			lastImageArray = imageArray;
			lastImageSize = currentImageSize;
			lastImage = currentImage - 1;
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
	 * @see gdsc.smlm.results.ResultsSource#next(java.awt.Rectangle)
	 */
	@Override
	public float[] next(Rectangle bounds)
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
	 * @see gdsc.smlm.results.ResultsSource#get(int, java.awt.Rectangle)
	 */
	@Override
	public float[] get(int frame, Rectangle bounds)
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
	 * @param logProgress Set to true to send progress to the ImageJ log
	 */
	public void setLogProgress(boolean logProgress)
	{
		this.logProgress = logProgress;
	}
}
