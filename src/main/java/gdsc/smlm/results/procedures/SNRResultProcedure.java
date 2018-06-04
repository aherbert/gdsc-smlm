package gdsc.smlm.results.procedures;

import gdsc.core.data.DataException;
import gdsc.smlm.results.MemoryPeakResults;
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
 * Contains functionality to obtain the Signal-to-Noise Ratio (SNR) for results.
 */
//@formatter:off
public class SNRResultProcedure extends AbstractResultProcedure implements 
	PeakResultProcedure
//@formatter:on
{
	/** The Signal-to-Noise Ratio (SNR). */
	public float[] snr;

	/**
	 * Instantiates a new precision result procedure.
	 *
	 * @param results
	 *            the results
	 * @throws DataException
	 *             if the results have no noise
	 */
	public SNRResultProcedure(MemoryPeakResults results) throws DataException
	{
		super(results);
		if (!results.hasNoise())
			throw new DataException("Results do not have noise");
		if (!results.hasMeanIntensity())
			throw new DataException("Results do not have mean intensity");
	}

	/**
	 * Gets the SNR for the results.
	 * <p>
	 * The SNR is computed using the mean signal divided by the noise.
	 *
	 * @return the snr
	 */
	public float[] getSNR()
	{
		i = 0;
		snr = allocate(snr);
		results.forEach((PeakResultProcedure) this);
		return snr;
	}

	@Override
	public void execute(PeakResult peakResult)
	{
		this.snr[i++] = peakResult.getSNR();
	}
}