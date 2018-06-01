package gdsc.smlm.results.data;

import gdsc.smlm.results.Gaussian2DPeakResultHelper;
import gdsc.smlm.results.PeakResult;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Gets the SNR from the PeakResult assuming a Gaussian 2D PSF. This is ratio of the average signal value to the
 * standard deviation of the background. The result must have the standard deviation for
 * each dimension in the first two additional parameters of the PeakResult parameter array.
 * <p>
 * Assumes that the mean signal is the total signal within 1 standard deviation of the centre divided by the elliptical
 * area of the Gaussian.
 */
public class Gaussian2DPeakResultSNRData extends PeakResultDataFloat
{
	final int i = PeakResult.STANDARD_PARAMETERS;
	final int j = i + 1;

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultData#getValue(gdsc.smlm.results.PeakResult)
	 */
	public Float getValue(PeakResult result)
	{
		return new Float(Gaussian2DPeakResultHelper.getSNR1(result.getSignal(), result.getParameter(i),
				result.getParameter(j), result.getNoise()));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultData#getValueName()
	 */
	public String getValueName()
	{
		return "Gaussian2D SNR";
	}
}