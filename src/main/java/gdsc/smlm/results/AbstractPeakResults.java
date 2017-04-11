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
	public static final double DEFAULT_NM_PER_PIXEL = 100;
	public static final double DEFAULT_GAIN = 1;
	public static final boolean DEFAULT_EMCCD = true;

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

	/**
	 * Adds the result.
	 * <p>
	 * Convenience method that just calls {@link #add(int, int, int, float, double, float, float[], float[])} passing
	 * the properties from the PeakResult. Any other properties in derived classes of PeakResult will be lost.
	 *
	 * @param result
	 *            the result
	 */
	public void add(PeakResult result)
	{
		add(result.getFrame(), result.origX, result.origY, result.origValue, result.error, result.noise, result.params,
				result.paramsStdDev);
	}

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

	/**
	 * Get the nm-per-pixel from the calibration, or if not available, return the {@link #DEFAULT_NM_PER_PIXEL}
	 * 
	 * @return the nmPerPixel
	 */
	public double getNmPerPixel()
	{
		return (calibration != null && calibration.hasNmPerPixel()) ? calibration.getNmPerPixel()
				: DEFAULT_NM_PER_PIXEL;
	}

	/**
	 * Get the gain from the calibration, or if not available, return the {@link #DEFAULT_GAIN}
	 * 
	 * @return the gain
	 */
	public double getGain()
	{
		return (calibration != null && calibration.hasGain()) ? calibration.getGain() : DEFAULT_GAIN;
	}

	/**
	 * Get the EMCCD flag from the calibration, or if not available, return the {@link #DEFAULT_EMCCD}
	 * 
	 * @return the gain
	 */
	public boolean isEMCCD()
	{
		return (calibration != null && calibration.hasEMCCD()) ? calibration.isEmCCD() : DEFAULT_EMCCD;
	}

	/**
	 * Checks if the results have a valid calibration to compute the localisation precision. This requires the pixel
	 * size and camera details (gain and CCD type).
	 *
	 * @return true, if is calibrated for precision
	 */
	public boolean isCalibratedForPrecision()
	{
		if (calibration != null)
		{
			return calibration.hasNmPerPixel() && calibration.hasGain() && calibration.hasEMCCD();
		}
		return false;
	}

	/**
	 * Checks if the results have a basic calibration. This requires the pixel
	 * size and camera gain.
	 *
	 * @return true, if is calibrated
	 */
	public boolean isCalibrated()
	{
		if (calibration != null)
		{
			return calibration.hasNmPerPixel() && calibration.hasGain();
		}
		return false;
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
