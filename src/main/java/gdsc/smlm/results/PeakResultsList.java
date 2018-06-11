/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 * 
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.results;

import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import gdsc.smlm.data.config.CalibrationProtos.Calibration;
import gdsc.smlm.data.config.PSFProtos.PSF;

/**
 * Wrapper class to output to multiple results destinations
 */
public class PeakResultsList extends AbstractPeakResults implements PeakResults
{
	private List<PeakResults> results = new ArrayList<PeakResults>();

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
	 * Gets the output.
	 *
	 * @param index
	 *            the index
	 * @return the output
	 */
	public PeakResults getOutput(int index)
	{
		return results.get(index);
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
	 * @see gdsc.smlm.results.PeakResults#add(int, int, int, float, double, float, float, float[], float[])
	 */
	public void add(int peak, int origX, int origY, float origValue, double error, float noise, float meanIntensity,
			float[] params, float[] paramsStdDev)
	{
		for (PeakResults peakResults : results)
			peakResults.add(peak, origX, origY, origValue, error, noise, meanIntensity, params, paramsStdDev);
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

	public void addAll(Collection<PeakResult> results)
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
		newList.copySettings(this);
		for (PeakResults peakResults : this.results)
		{
			// This will copy the settings
			//newList.addOutput(peakResults);

			// This assumes the settings are OK, i.e. the result was added
			// using addOutput(...). 
			newList.results.add(SynchronizedPeakResults.create(peakResults));
		}
		return newList;
	}

	// Pass through all the modifications to the list objects as well as this 
	// so that any new list objects can copy the settings.

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
