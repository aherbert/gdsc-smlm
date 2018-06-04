package gdsc.smlm.results.procedures;

import gdsc.core.data.DataException;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResultData;
import gdsc.smlm.results.data.PeakResultDataFloat;

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

	private PeakResultData<Float> mean;

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
	}

	/**
	 * Gets the SNR for the results.
	 * <p>
	 * An attempt is made to compute the SNR using the average signal divided by the noise. If the average signal cannot
	 * be computed then the total signal is used, effectively making this represent the entire result intensity in a 1x1
	 * area.
	 * 
	 * @return the snr
	 */
	public float[] getSNR()
	{
		PeakResultData<Float> mean;
		if (PSFHelper.isGaussian2D(results.getPSF()))
		{
			int[] indices = PSFHelper.getGaussian2DWxWyIndices(results.getPSF());
			final int isx = indices[0];
			final int isy = indices[1];
			mean = new PeakResultDataFloat()
			{
				@Override
				public Float getValue(PeakResult result)
				{
					return new Float(Gaussian2DPeakResultHelper.getMeanSignalUsingP05(result.getSignal(),
							result.getParameter(isx), result.getParameter(isy)));
				}
			};
		}
		else
		{
			// This should be extended to support other PSFs.
			//
			// One mechanism would be to have the PSF contain a parameter that describes 
			// the area of the localisation. The PSF may not be entirely captured by
			// that area. In this case it could be called the normalisation area and 
			// it represents the area of the PSF that contains a set fraction of the signal.
			// E.g. in the case of a Gaussian 2D function then the area of pi*sx*sy 
			// contains 0.6827 of the signal.
			//
			// If this is to be stored then it may be better to store the mean signal
			// as a parameter. This could be done during fitting by evaluating the PSF, 
			// sorting the values, and then counting the number of values to achieve
			// a fraction of 1. This is the area and can be used to compute the 
			// mean signal = [cumulative signal] / area.
			// This can be done using a partial sort with a heap.

			mean = new PeakResultDataFloat()
			{
				@Override
				public Float getValue(PeakResult result)
				{
					return new Float(result.getSignal());
				}
			};
		}
		return getSNR(mean);
	}

	/**
	 * Gets the SNR for the results using the given method to compute the mean signal.
	 * <p>
	 * The SNR is computed using the mean signal divided by the noise.
	 *
	 * @param mean
	 *            the method to compute the mean signal
	 * @return the snr
	 */
	public float[] getSNR(PeakResultData<Float> mean)
	{
		this.mean = mean;
		i = 0;
		snr = allocate(snr);
		results.forEach((PeakResultProcedure) this);
		return snr;
	}

	@Override
	public void execute(PeakResult peakResult)
	{
		this.snr[i++] = mean.getValue(peakResult) / peakResult.getNoise();
	}
}