package gdsc.smlm.ij.results;

import gdsc.smlm.results.AbstractPeakResults;
import gdsc.smlm.results.Calibration;
import gdsc.smlm.results.ThreadSafePeakResults;

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
 * Wrap the fit results to provide convenience methods
 */
public abstract class IJAbstractPeakResults extends AbstractPeakResults implements ThreadSafePeakResults
{
	public void setCalibration(double nmPerPixel, double gain)
	{
		Calibration calibration = getCalibration();
		if (calibration == null)
			calibration = new Calibration();
		calibration.setNmPerPixel(nmPerPixel);
		calibration.setGain(gain);
		setCalibration(calibration);
	}
}
