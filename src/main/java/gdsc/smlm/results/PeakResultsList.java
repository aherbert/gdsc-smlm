package gdsc.smlm.results;

import java.awt.Rectangle;
import java.util.Collection;

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

import java.util.LinkedList;
import java.util.List;

import gdsc.smlm.data.config.CalibrationConfig.Calibration;
import gdsc.smlm.data.config.PSFConfig.PSF;

/**
 * Wrapper class to output to multiple results destinations
 */
public class PeakResultsList extends AbstractPeakResults implements PeakResults
{
	private List<PeakResults> results = new LinkedList<PeakResults>();

	/**
	 * Add a result format to the output. If a PeakResultsList is passed then it will be
	 * separated into the child PeakResults instances. This will break the size() function
	 * of any input PeakResultsList since only the children will remain within this list.
	 * <p>
	 * Sets the settings (source and configuration) of the child to the same as this list
	 * 
	 * @param peakResults
	 */
	public void addOutput(PeakResults peakResults)
	{
		if (peakResults instanceof PeakResultsList)
		{
			for (PeakResults r : ((PeakResultsList) peakResults).results)
				addOutput(r);
		}
		else
		{
			peakResults.copySettings(this);
			results.add(peakResults);
		}
	}

	/**
	 * @return The number of outputs contained in the list
	 */
	public int numberOfOutputs()
	{
		return results.size();
	}

	/**
	 * @return The outputs
	 */
	public PeakResults[] toArray()
	{
		return results.toArray(new PeakResults[results.size()]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#begin()
	 */
	public void begin()
	{
		for (PeakResults peakResults : results)
			peakResults.begin();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#add(int, int, int, float, double, float, float[], float[])
	 */
	public void add(int peak, int origX, int origY, float origValue, double error, float noise, float[] params,
			float[] paramsStdDev)
	{
		for (PeakResults peakResults : results)
			peakResults.add(peak, origX, origY, origValue, error, noise, params, paramsStdDev);
	}

	public void add(PeakResult result)
	{
		for (PeakResults peakResults : results)
			peakResults.add(result);
	}

	public void addAll(PeakResult[] results)
	{
		for (PeakResults peakResults : this.results)
			peakResults.addAll(results);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#size()
	 */
	public int size()
	{
		return (results.isEmpty()) ? 0 : results.get(0).size();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#end()
	 */
	public void end()
	{
		for (PeakResults peakResults : results)
			peakResults.end();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#isActive()
	 */
	public boolean isActive()
	{
		for (PeakResults peakResults : this.results)
			if (peakResults.isActive())
				return true;
		return false;
	}

	/**
	 * Checks all the results in the list. If any are not thread safe then they are wrapped with a
	 * SynchronizedPeakResults container.
	 *
	 * @return the thread safe list
	 */
	public PeakResultsList getThreadSafeList()
	{
		PeakResultsList newList = new PeakResultsList();
		for (PeakResults peakResults : this.results)
		{
			if (!(peakResults instanceof ThreadSafePeakResults))
			{
				peakResults = new SynchronizedPeakResults(peakResults);
			}
			newList.addOutput(peakResults);
		}
		return newList;
	}

	// Pass through all the modifications to the list objects

	public void addAll(Collection<PeakResult> results)
	{
		super.addAll(results);
		for (PeakResults peakResults : this.results)
			peakResults.addAll(results);
	}

	public void setSource(ImageSource source)
	{
		super.setSource(source);
		for (PeakResults peakResults : results)
			peakResults.setSource(source);
	}

	public void setBounds(Rectangle bounds)
	{
		super.setBounds(bounds);
		for (PeakResults peakResults : results)
			peakResults.setBounds(bounds);
	}

	public void setCalibration(Calibration calibration)
	{
		super.setCalibration(calibration);
		for (PeakResults peakResults : results)
			peakResults.setCalibration(calibration);
	}

	public void setPSF(PSF psf)
	{
		super.setPSF(psf);
		for (PeakResults peakResults : results)
			peakResults.setPSF(psf);
	}

	public void setConfiguration(String configuration)
	{
		super.setConfiguration(configuration);
		for (PeakResults peakResults : results)
			peakResults.setConfiguration(configuration);
	}

	public void setName(String name)
	{
		super.setName(name);
		for (PeakResults peakResults : results)
			peakResults.setName(name);
	}

	public void copySettings(PeakResults results)
	{
		super.copySettings(results);
		for (PeakResults peakResults : this.results)
			peakResults.copySettings(results);
	}
}
