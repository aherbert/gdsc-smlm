package gdsc.smlm.ij;

import gdsc.smlm.ij.utils.ImageConverter;
import gdsc.smlm.results.ImageSource;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.io.FileInfo;
import ij.process.ImageProcessor;

import java.awt.Rectangle;

import org.apache.commons.math3.util.FastMath;

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
 * Represent an ImageJ image as a results source. Supports all greyscale images. Only processes channel 0 of 32-bit
 * colour images.
 */
public class IJImageSource extends ImageSource
{
	private int slice;
	private int singleFrame = 0;
	private int extraFrames = 0;
	private Object[] imageArray;
	private ImageStack imageStack;
	private String path;

	/**
	 * Create a new image source using the path to the image file
	 * 
	 * @param name
	 * @param path
	 */
	public IJImageSource(String name, String path)
	{
		super(name);
		this.path = path;
	}

	/**
	 * Create a new image source from an ImagePlus
	 * 
	 * @param imp
	 */
	public IJImageSource(ImagePlus imp)
	{
		super(imp.getTitle());
		initialise(imp);
	}

	private boolean initialise(ImagePlus imp)
	{
		imageArray = null;
		imageStack = null;
		if (imp == null)
			return false;
		ImageStack s = imp.getImageStack();
		if (s.isVirtual())
		{
			// We must use the image stack to get the image data for virtual images
			imageStack = s;
		}
		else
		{
			// We can access the image array directly
			imageArray = s.getImageArray();
		}
		width = imp.getWidth();
		height = imp.getHeight();
		// Store the number of valid frames
		if (singleFrame > 0)
			frames = 1 + this.extraFrames;
		else
			frames = imp.getStackSize();			
		slice = 0;
		FileInfo info = imp.getOriginalFileInfo();
		if (info != null)
			path = info.directory + info.fileName;
		return true;
	}

	/**
	 * Create a new image source from an ImageProcessor
	 * 
	 * @param name
	 * @param ip
	 */
	public IJImageSource(String name, ImageProcessor ip)
	{
		super(name);
		imageArray = new Object[] { ip.getPixels() };
		width = ip.getWidth();
		height = ip.getHeight();
		frames = 1;
	}

	/**
	 * Create a new image source from an ImagePlus. Specify a single frame for processing.
	 * 
	 * @param imp
	 * @param frame
	 */
	public IJImageSource(ImagePlus imp, int frame)
	{
		this(imp, frame, 0);
	}

	/**
	 * Create a new image source from an ImagePlus. Specify a start frame and any additional frames (after the start)
	 * for processing.
	 * 
	 * @param imp
	 * @param startFrame
	 *            The first frame to process
	 * @param extraFrames
	 *            The number of additional frames to process
	 */
	public IJImageSource(ImagePlus imp, int startFrame, int extraFrames)
	{
		super(imp.getTitle());
		// Ensure only a single frame is processed
		singleFrame = startFrame;
		this.extraFrames = FastMath.max(0, extraFrames);
		initialise(imp);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#openSource()
	 */
	@Override
	protected boolean openSource()
	{
		if (nullImageArray() && imageStack == null)
		{
			// Try and find the image using the name or path
			ImagePlus imp = null;
			if (getName() != null)
				imp = WindowManager.getImage(getName());

			if (imp == null)
			{
				// Try and open the original image from file
				if (path != null && path.length() > 0)
				{
					imp = IJ.openImage(path);
					if (imp == null)
					{
						// Some readers return null and display the image, e.g. BioFormats.
						// Add code to handle this.
					}
					else
					{
						// Ensure the image has the correct name
						if (getName() != null)
							imp.setTitle(getName());
					}
				}
			}
			return initialise(imp);
		}
		slice = 0;
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#close()
	 */
	@Override
	public void close()
	{
		imageArray = null;
		imageStack = null;
	}

	/**
	 * @return True is the image array or any part of it is null
	 */
	private boolean nullImageArray()
	{
		if (imageArray == null)
			return true;
		// Check the image array. This is set to null when an image is closed by ImageJ
		// allowing us to detect when the image is still available
		for (Object o : imageArray)
			if (o == null)
				return true;
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#nextFrame(java.awt.Rectangle)
	 */
	@Override
	protected float[] nextFrame(Rectangle bounds)
	{
		++slice;
		if (singleFrame > 0)
		{
			// Return frames from the starting frame until the extra frames limit is reached
			return (slice - 1 <= extraFrames) ? get(singleFrame + slice - 1, bounds) : null;
		}
		return get(slice, bounds);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#getFrame(int, java.awt.Rectangle)
	 */
	@Override
	protected float[] getFrame(int frame, Rectangle bounds)
	{
		if (imageArray != null)
		{
			if (frame > 0 && frame <= imageArray.length)
			{
				return ImageConverter.getData(imageArray[frame - 1], width, height, bounds, null);
			}
		}
		else if (imageStack != null)
		{
			// This is a virtual stack so access the image processor through the virtual stack object
			if (frame > 0 && frame <= imageStack.getSize())
			{
				return ImageConverter.getData(imageStack.getPixels(frame), width, height, bounds, null);
			}
		}
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ResultsSource#init(com.thoughtworks.xstream.XStream)
	 */
	@Override
	public void init(XStream xs)
	{
		super.init(xs);
		xs.omitField(IJImageSource.class, "slice");
		// Q. Should this be omitted? It may be useful when reloading results to restore the 
		// single frame read property of the next() function
		//xs.omitField(IJImageSource.class, "singleFrame");
		xs.omitField(IJImageSource.class, "imageArray");
		xs.omitField(IJImageSource.class, "imageStack");
	}

	@Override
	public boolean isValid(int frame)
	{
		if (singleFrame > 0)
		{
			return (frame >= singleFrame && frame <= singleFrame + extraFrames);
		}
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
		if (path != null)
			s += String.format(" (%s)", path);
		return s;
	}
}
