package gdsc.smlm.results;

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

import java.awt.Rectangle;
import java.util.Collection;

/**
 * Abstract base class for peak results.
 */
public abstract class AbstractPeakResults implements PeakResults
{
	protected ImageSource source = null;
	protected Rectangle bounds = null;
	protected Calibration calibration = null;
	protected String configuration = "";
	protected String name = "";

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#begin()
	 */
	public abstract void begin();

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#add(int, int, int, float, double, float, float[], float[])
	 */
	public abstract void add(int peak, int origX, int origY, float origValue, double chiSquared, float noise,
			float[] params, float[] paramsStdDev);

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#addAll(java.util.Collection)
	 */
	public abstract void addAll(Collection<PeakResult> results);

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#size()
	 */
	public abstract int size();

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#end()
	 */
	public abstract void end();

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#isActive()
	 */
	abstract public boolean isActive();

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#setSource(java.lang.String)
	 */
	public void setSource(ImageSource source)
	{
		this.source = source;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#getSource()
	 */
	public ImageSource getSource()
	{
		return source;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#setBounds(java.lang.String)
	 */
	public void setBounds(Rectangle bounds)
	{
		this.bounds = bounds;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#getBounds()
	 */
	public Rectangle getBounds()
	{
		return bounds;
	}

	public String getBoundsString()
	{
		if (bounds != null)
		{
			return String.format("x%d y%d w%d h%d", bounds.x, bounds.y, bounds.width, bounds.height);
		}
		return "";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#setCalibration(java.lang.String)
	 */
	public void setCalibration(Calibration calibration)
	{
		this.calibration = calibration;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#getCalibration()
	 */
	public Calibration getCalibration()
	{
		return calibration;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#setConfiguration(java.lang.String)
	 */
	public void setConfiguration(String configuration)
	{
		this.configuration = configuration;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#getConfiguration()
	 */
	public String getConfiguration()
	{
		return configuration;
	}


	/**
	 * @return The name of the results set (or the source if empty)
	 */
	public String getName()
	{
		if (name.length() > 0)
			return name;
		return (getSource() != null) ? getSource().getName() : "";
	}

	/**
	 * @param name
	 *            The name of the results set
	 */
	public void setName(String name)
	{
		if (name == null)
			this.name = "";
		else
			this.name = name;
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#copySettings(gdsc.utils.fitting.results.PeakResults)
	 */
	public void copySettings(PeakResults peakResults)
	{
		this.setSource(peakResults.getSource());
		this.setBounds(peakResults.getBounds());
		this.setCalibration(peakResults.getCalibration());
		this.setConfiguration(peakResults.getConfiguration());
		this.setName(peakResults.getName());
	}
}
