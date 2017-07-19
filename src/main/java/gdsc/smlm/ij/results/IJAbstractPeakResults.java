package gdsc.smlm.ij.results;

import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.results.AbstractPeakResults;
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
		CalibrationWriter cw = getCalibrationWriterSafe();
		cw.setNmPerPixel(nmPerPixel);
		cw.setCountPerPhoton(gain);
		setCalibration(cw.getCalibration());
	}
}
