package gdsc.smlm.ij.settings;

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

import gdsc.smlm.results.Calibration;

import java.util.ArrayList;

/**
 * Contain the configuration file settings for the batch fitting plugin
 */
public class BatchSettings
{
	public ArrayList<String> images = new ArrayList<String>();
	public ArrayList<ParameterSettings> parameters = new ArrayList<ParameterSettings>();
	public String resultsDirectory = null;
	public boolean runPeakFit = true; 
	private Calibration calibration = new Calibration();
	
	/**
	 * @return the calibration
	 */
	public Calibration getCalibration()
	{
		if (calibration == null)
			calibration = new Calibration();
		return calibration;
	}
	/**
	 * @param calibration the calibration to set
	 */
	public void setCalibration(Calibration calibration)
	{
		this.calibration = calibration;
	}
}
