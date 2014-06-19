package gdsc.smlm.results;

import gdsc.smlm.utils.XmlUtils;

import java.awt.Rectangle;

import com.thoughtworks.xstream.XStream;
import com.thoughtworks.xstream.XStreamException;
import com.thoughtworks.xstream.io.xml.DomDriver;

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
 * Abstract base class for the image source for peak results.
 */
public abstract class ImageSource
{
	private String name;
	protected int width, height, frames;

	/**
	 * Create the image source
	 * @param name The name of the image source
	 */
	public ImageSource(String name)
	{
		setName(name);
	}

	/**
	 * Opens the source
	 */
	public abstract boolean open();

	/**
	 * Get the width of the results. The frame returned by {@link #next()} will be equal to width * height.
	 * 
	 * @return
	 */
	public int getWidth()
	{
		return width;
	}

	/**
	 * Get the height of the results. The frame returned by {@link #next()} will be equal to width * height.
	 * 
	 * @return
	 */
	public int getHeight()
	{
		return height;
	}

	/**
	 * Get the number of frames
	 * 
	 * @return
	 */
	public int getFrames()
	{
		return frames;
	}

	/**
	 * Get the next frame. The data is is packed in yx order: index = y * width + x;
	 * <p>
	 * Provides serial access to the data after a successful call to {@link #open()}
	 * 
	 * @return the next frame (or null if at the end)
	 */
	public float[] next()
	{
		return next(null);
	}

	/**
	 * Get the next frame. The data is is packed in yx order: index = y * width + x;
	 * <p>
	 * Provides serial access to the data after a successful call to {@link #open()}
	 * <p>
	 * Note: bounds.x + bounds.width must be less or equal to than {@link #getWidth()}, similarly for height.
	 * 
	 * @param bounds
	 *            The bounding limits of the frame to extract
	 * @return the next frame (or null if at the end)
	 */
	public abstract float[] next(Rectangle bounds);

	/**
	 * Get a specific frame from the results. Return null if the frame is not available.
	 * <p>
	 * Provides random access to the data after a successful call to {@link #open()}. This operation may be
	 * significantly slower than using {@link #next()} to read all the data.
	 * 
	 * @param frame
	 * @return the frame (or null)
	 */
	public float[] get(int frame)
	{
		return get(frame, null);
	}

	/**
	 * Get a specific frame from the results. Return null if the frame is not available.
	 * <p>
	 * Provides random access to the data after a successful call to {@link #open()}. This operation may be
	 * significantly slower than using {@link #next()} to read all the data.
	 * <p>
	 * Note: bounds.x + bounds.width must be less or equal to than {@link #getWidth()}, similarly for height.
	 * 
	 * @param frame
	 * @param bounds
	 *            The bounding limits of the frame to extract
	 * @return the frame (or null)
	 */
	public abstract float[] get(int frame, Rectangle bounds);

	/**
	 * Get the name of the results source
	 * 
	 * @return
	 */
	public String getName()
	{
		return name;
	}

	/**
	 * Set the name of the results source
	 * 
	 * @param name
	 */
	public void setName(String name)
	{
		if (name != null && name.length() > 0)
			this.name = name;
		else
			this.name = "";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		// Over-ride this to produce a nicer output description of the results source
		return String.format("%s [%dx%dx%d]", name, width, height, frames);
	}

	/**
	 * @return An XML representation of this object
	 */
	public String toXML()
	{
		XStream xs = new XStream(new DomDriver());
		try
		{
			init(xs);
			return xs.toXML(this);
		}
		catch (XStreamException ex)
		{
			//ex.printStackTrace();
		}
		return "";
	}

	public static ImageSource fromXML(String xml)
	{
		return (ImageSource) XmlUtils.fromXML(xml);
	}

	/**
	 * Initialise the XStream class to be used for XML conversion. Provides the opportunity for a 
	 * class to omit fields from the XStream serialisation.
	 * 
	 * @param xs
	 */
	public void init(XStream xs)
	{
	}
}
