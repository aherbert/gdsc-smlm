package gdsc.smlm.results;

import java.awt.Rectangle;

import com.thoughtworks.xstream.annotations.XStreamOmitField;

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
 * Wraps an image source and allows data to be interlaced within regular blocks.
 */
public class InterlacedImageSource extends ImageSource
{
	private final int start, size, skip;
	private final ImageSource imageSource;

	// Record the number of frames returned from the current block
	@XStreamOmitField
	private int counter;

	/**
	 * Create a new interlaced image source using the given image source
	 * <p>
	 * Note: The input source cannot be aggregated as the data to interlace is assumed to be contiguous from frame 1
	 * 
	 * @param imageSource
	 *            The image source of interlaced data (must not be null or an AggregatedImageSource)
	 * @param start
	 *            The first frame that contains data
	 * @param size
	 *            The number of continuous frames containing data
	 * @param skip
	 *            The number of continuous frames to ignore before the next data
	 */
	public InterlacedImageSource(ImageSource imageSource, int start, int size, int skip)
	{
		super("");
		if (imageSource == null)
			throw new IllegalArgumentException("Image source must not be null");
		if (imageSource instanceof AggregatedImageSource)
			throw new IllegalArgumentException("Image source must not be aggregated");
		if (start < 1)
			throw new IllegalArgumentException("The start frame must be 1 or above");
		if (size < 1)
			throw new IllegalArgumentException("The read size must be 1 or above");
		if (skip < 0)
			throw new IllegalArgumentException("The skip length must be 0 or above");
		setName(String.format("Interlaced (%d,%d,%d) %s", start, size, skip, imageSource.getName()));
		this.imageSource = imageSource;
		this.start = start;
		this.size = size;
		this.skip = skip;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#getXOrigin()
	 */
	@Override
	public int getXOrigin()
	{
		return imageSource.getXOrigin();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#getYOrigin()
	 */
	@Override
	public int getYOrigin()
	{
		return imageSource.getYOrigin();
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
		return (int) Math.ceil((imageSource.getFrames() - start + 1) * ((double) size / (size + skip)));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#getParent()
	 */
	@Override
	public ImageSource getParent()
	{
		return imageSource;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#getOriginal()
	 */
	@Override
	public ImageSource getOriginal()
	{
		return imageSource.getOriginal();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ResultsSource#openSource()
	 */
	@Override
	protected boolean openSource()
	{
		// Assume frame start at 1 and set the intial skip		
		final int initialSkip = (start - 1);
		counter = -initialSkip;
		return imageSource.openSource();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ResultsSource#close()
	 */
	@Override
	public void close()
	{
		imageSource.close();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#nextFrame(java.awt.Rectangle)
	 */
	@Override
	protected float[] nextFrame(Rectangle bounds)
	{
		// Skip frames until at the start of the next block
		while (counter < 0)
		{
			if (imageSource.next(bounds) == null)
				return null;
			counter++;
		}

		// Read the next frame in the current block
		final float[] image = imageSource.next(bounds);

		// Check if this is the final frame in the current block
		if (++counter >= size)
		{
			counter = -skip;
		}
		// Set the frame to the last one read from the source
		setFrameNumber(imageSource.getStartFrameNumber(), imageSource.getEndFrameNumber());
		//System.out.printf("Interlaced %d-%d\n", getStartFrameNumber(), getEndFrameNumber());
		return image;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#getFrame(int, java.awt.Rectangle)
	 */
	@Override
	protected float[] getFrame(int frame, Rectangle bounds)
	{
		if (frame < 1)
			return null;

		// Check if the frame is allowed:
		//    Start
		//      |
		// |----|Block|Skip|Block|Skip|Block|Skip
		if (frame < start)
		{
			return null;
		}
		int frameInBlock = (frame - start) % (size + skip);
		if (frameInBlock >= size)
		{
			return null;
		}
		return imageSource.get(frame, bounds);
	}

	/**
	 * @param frame
	 * @return
	 */
	@Override
	public boolean isValid(int frame)
	{
		return imageSource.isValid(frame);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#toString()
	 */
	@Override
	public String toString()
	{
		return String.format("%s (Interlaced %d,%d,%d)", imageSource.toString(), start, size, skip);
	}

	/**
	 * @return The first frame that contains data
	 */
	public int getStart()
	{
		return start;
	}

	/**
	 * @return The number of continuous frames containing data
	 */
	public int getSize()
	{
		return size;
	}

	/**
	 * @return The number of continuous frames to ignore before the next data
	 */
	public int getSkip()
	{
		return skip;
	}
}
