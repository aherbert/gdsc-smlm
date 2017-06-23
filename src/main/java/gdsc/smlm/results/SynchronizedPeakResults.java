package gdsc.smlm.results;

import java.awt.Rectangle;
import java.util.Collection;

import gdsc.smlm.data.config.PSFConfig.PSF;
import gdsc.smlm.data.config.CalibrationConfig.Calibration;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Wraps a peak results with synchronized methods.
 */
public class SynchronizedPeakResults implements ThreadSafePeakResults
{
	private final PeakResults r;
	private final Object lock = new Object();

	/**
	 * Instantiates a new synchronized peak results.
	 *
	 * @param peakResults
	 *            the peak results
	 * @throws IllegalArgumentException
	 *             if the results are null
	 */
	public SynchronizedPeakResults(PeakResults peakResults)
	{
		if (peakResults == null)
			throw new IllegalArgumentException("PeakResults must not be null");
		this.r = peakResults;
	}

	/**
	 * Creates a PeakResults object that is synchronized if the thread count is above 1, otherwise the input results are
	 * returned.
	 *
	 * @param peakResults
	 *            the peak results
	 * @param threadCount
	 *            the thread count
	 * @return the peak results
	 * @throws IllegalArgumentException
	 *             if the results are null
	 */
	public static PeakResults create(PeakResults peakResults, int threadCount)
	{
		return (threadCount <= 1) ? peakResults : new SynchronizedPeakResults(peakResults);
	}

	//@formatter:off
	
	public void begin()
	{
		synchronized (lock)	{ r.begin(); }
	}

	public void add(int peak, int origX, int origY, float origValue, double error, float noise, float[] params,
			float[] paramsStdDev)
	{
		synchronized (lock)	{ r.add(peak, origX, origY, origValue, error, noise, params, paramsStdDev); }
	}

	public void add(PeakResult result)
	{
		synchronized (lock)	{ r.add(result); }
	}
	
	public void addAll(Collection<PeakResult> results)
	{
		synchronized (lock)	{ r.addAll(results); }
	}
	
	public void addAll(PeakResult[] results)
	{
		synchronized (lock)	{ r.addAll(results); }
	}

	public int size()
	{
		synchronized (lock)	{ return r.size(); }
	}

	public void end()
	{
		synchronized (lock)	{ r.end(); }
	}

	public boolean isActive()
	{
		synchronized (lock)	{ return r.isActive(); }
	}

	public void setSource(ImageSource source)
	{
		synchronized (lock)	{ r.setSource(source); }
	}

	public ImageSource getSource()
	{
		synchronized (lock)	{ return r.getSource(); }
	}

	public void setBounds(Rectangle bounds)
	{
		synchronized (lock)	{ r.setBounds(bounds); }
	}

	public Rectangle getBounds()
	{
		synchronized (lock)	{ return r.getBounds(); }
	}

	public void setCalibration(Calibration calibration)
	{
		synchronized (lock)	{ r.setCalibration(calibration); }
	}

	public Calibration getCalibration()
	{
		synchronized (lock)	{ return r.getCalibration(); }
	}

	public void setPSF(PSF psf)
	{
		synchronized (lock)	{ r.setPSF(psf); }
	}

	public PSF getPSF()
	{
		synchronized (lock)	{ return r.getPSF(); }
	}

	public void setConfiguration(String configuration)
	{
		synchronized (lock)	{ r.setConfiguration(configuration); }
	}

	public String getConfiguration()
	{
		synchronized (lock)	{ return r.getConfiguration(); }
	}

	public String getName()
	{
		synchronized (lock)	{ return r.getName(); }
	}

	public void setName(String name)
	{
		synchronized (lock)	{ r.setName(name); }
	}

	public void copySettings(PeakResults peakResults)
	{
		synchronized (lock)	{ r.copySettings(peakResults); }
	}
	
	//@formatter:on
}
