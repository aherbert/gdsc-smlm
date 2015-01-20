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

/**
 * Contain the settings for the PSF Estimator
 */
public class PSFEstimatorSettings
{
	public int numberOfPeaks = 1000;
	public double pValue = 0.01;
	public boolean updatePreferences = true;
	public boolean debugPSFEstimator = false;
	public boolean iterate = true;
	public boolean showHistograms = false;
	public int histogramBins = 100;
}
