package gdsc.smlm.ij;

import gdsc.smlm.results.ImageSource;

import java.awt.Rectangle;

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
 * Wraps an image source and allows aggregation of consecutive frames.
 */
public class AggregatedImageSource extends ImageSource
{
	private int aggregate;
	private ImageSource imageSource;

	// Used for frame-based read
	private int lastFrame;
	private float[] lastImage = null;

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ResultsSource#init(com.thoughtworks.xstream.XStream)
	 */
	@Override
	public void init(XStream xs)
	{
		super.init(xs);
		if (imageSource != null)
			imageSource.init(xs);
		//---
		// TODO - This does not work as the fields are defined in the super class ImageSource.
		// using ImageSource.class to omit the fields would mean they are also omitted from the 
		// imageSource variable. So for now we get an output that has the width,height,frames fields
		// twice in the XML, once for each class (AggregatedImageSource and the ImageSource field).
		//xs.omitField(getClass(), "width");
		//xs.omitField(getClass(), "height");
		//xs.omitField(getClass(), "frames");
		//---
		xs.omitField(AggregatedImageSource.class, "lastFrame");
		xs.omitField(AggregatedImageSource.class, "lastImage");
	}

	/**
	 * Create a new aggregated image source using the given image source
	 * 
	 * @param imageSource
	 *            The image source to aggregate (must not be null)
	 * @param aggregate
	 *            The number of frames to aggregate (must be above 1)
	 */
	public AggregatedImageSource(ImageSource imageSource, int aggregate)
	{
		super("");
		if (imageSource == null)
			throw new IllegalArgumentException("Image source must not be null");
		if (aggregate < 2)
			throw new IllegalArgumentException("The number of frames to aggregate ({0}) must be above 1");
		setName("Aggregated (" + aggregate + ") " + imageSource.getName());
		this.imageSource = imageSource;
		this.aggregate = aggregate;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#getWidth()
	 */
	@Override
	public int getWidth()
	{
		return imageSource.getWidth();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#getHeight()
	 */
	@Override
	public int getHeight()
	{
		return imageSource.getHeight();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#getFrames()
	 */
	@Override
	public int getFrames()
	{
		return (int) Math.ceil((double) imageSource.getFrames() / aggregate);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ResultsSource#open()
	 */
	@Override
	public boolean open()
	{
		return imageSource.open();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ResultsSource#next(java.awt.Rectangle)
	 */
	@Override
	public float[] next(Rectangle bounds)
	{
		// Aggregate frames consecutive frames
		float[] image = imageSource.next(bounds);
		if (image != null)
		{
			for (int n = 1; n < aggregate; n++)
			{
				float[] image2 = imageSource.next(bounds);
				if (image2 == null)
					break;
				for (int i = 0; i < image.length; i++)
					image[i] += image2[i];
			}
		}
		return image;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ResultsSource#get(int, java.awt.Rectangle)
	 */
	@Override
	public float[] get(int frame, Rectangle bounds)
	{
		if (frame < 1)
			return null;

		// Calculate if the cache is invalid
		if (frame != lastFrame || lastImage == null)
		{
			// Calculate the required images
			int firstRealFrame = (frame - 1) * aggregate + 1;

			float[] image = imageSource.get(firstRealFrame, bounds);
			if (image != null)
			{
				for (int n = 1; n < aggregate; n++)
				{
					float[] image2 = imageSource.get(firstRealFrame + n, bounds);
					if (image2 == null)
						break;
					for (int i = 0; i < image.length; i++)
						image[i] += image2[i];
				}
			}
			// Cache the image
			lastImage = image;
			lastFrame = frame;
		}
		return lastImage;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#toString()
	 */
	@Override
	public String toString()
	{
		return String.format("%s (Aggregate %d images)", imageSource.toString(), aggregate);
	}
}
