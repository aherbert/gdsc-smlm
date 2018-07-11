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

import java.nio.ByteBuffer;
import java.security.MessageDigest;

import gdsc.core.utils.Digest;

/**
 * Provide digest functionality for ImageJ images to digest the pixels array
 */
public class PeakResultDigest
{
	/** The expected data bytes without the parameters */
	private static final int EXPECTED_DATA_BYTES = 48;

	private MessageDigest digest;
	// Allocate assuming 8 parameters and deviations
	private ByteBuffer buffer = ByteBuffer.allocate(EXPECTED_DATA_BYTES + 4 * 2 * 8);

	/**
	 * Instantiates a new IJ digest.
	 */
	public PeakResultDigest()
	{
		this(Digest.MD5);
	}

	/**
	 * Instantiates a new IJ digest.
	 *
	 * @param algorithm
	 *            the algorithm
	 */
	public PeakResultDigest(String algorithm)
	{
		digest = Digest.getDigest(algorithm);
	}

	/**
	 * Reset the digest.
	 */
	public void reset()
	{
		digest.reset();
	}

	/**
	 * Update the digest with the peak result.
	 *
	 * @param peakResult
	 *            the peak result
	 * @return the string
	 */
	public void update(PeakResult peakResult)
	{
		// Check buffer size
		final int n = peakResult.getNumberOfParameters();
		int required = n * 4;
		if (peakResult.hasParameterDeviations())
			required *= 2;

		if (buffer.capacity() < EXPECTED_DATA_BYTES + required)
			buffer = ByteBuffer.allocate(EXPECTED_DATA_BYTES + required);
		else
			buffer.clear();

		// Add standard data
		buffer.putInt(peakResult.getFrame()); // 4
		buffer.putInt(peakResult.getOrigX()); // 8
		buffer.putInt(peakResult.getOrigY()); // 12
		buffer.putFloat(peakResult.getOrigValue()); // 16
		buffer.putDouble(peakResult.getError()); // 24
		buffer.putFloat(peakResult.getNoise()); // 28
		buffer.putFloat(peakResult.getMeanIntensity()); // 32

		// Optional data
		if (peakResult.hasId())
			buffer.putInt(peakResult.getId()); // 36
		if (peakResult.hasEndFrame())
			buffer.putInt(peakResult.getEndFrame()); // 40
		if (peakResult.hasPrecision())
			buffer.putDouble(peakResult.getPrecision()); // 48

		for (int i = 0; i < n; i++)
		{
			buffer.putFloat(peakResult.getParameter(i));
		}
		if (peakResult.hasParameterDeviations())
		{
			for (int i = 0; i < n; i++)
			{
				buffer.putFloat(peakResult.getParameterDeviation(i));
			}
		}

		buffer.flip();

		digest.update(buffer);
	}

	/**
	 * Get the digest and reset.
	 *
	 * @return the hex string
	 */
	public String digest()
	{
		return Digest.toHex(digest.digest());
	}

	/**
	 * Get a snapshot of the current digest. This is done by cloning the digest. If this is not support then return
	 * null.
	 *
	 * @return the hex string (or null)
	 */
	public String snapshot()
	{
		try
		{
			MessageDigest d = (MessageDigest) digest.clone();
			return Digest.toHex(d.digest());
		}
		catch (CloneNotSupportedException e)
		{
			return null;
		}
	}
}
