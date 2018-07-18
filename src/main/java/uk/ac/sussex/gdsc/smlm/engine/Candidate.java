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
package uk.ac.sussex.gdsc.smlm.engine;

import uk.ac.sussex.gdsc.smlm.filters.Spot;

/**
 * Extends the Spot class with fields used during fitting.
 */
class Candidate extends Spot
{

	/** The index. */
	final public int index;

	/** Flag to indicate if the candidate has been fit. */
	public boolean fit = false;

	// Results of fitting

	/** The params. */
	public float[] params;

	/** The param deviations. */
	public float[] paramDevs;

	/** The error. */
	public double error;

	/** The noise. */
	public float noise;

	/** The mean intensity. */
	public float meanIntensity;

	/** The precision. */
	public double precision;

	/**
	 * Instantiates a new candidate.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param intensity
	 *            the intensity
	 * @param score
	 *            the score
	 * @param index
	 *            the index
	 */
	public Candidate(int x, int y, float intensity, float score, int index)
	{
		super(x, y, intensity, score);
		this.index = index;
	}

	/**
	 * Instantiates a new candidate.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param index
	 *            the index
	 * @param params
	 *            the params
	 * @param paramDevs
	 *            the param deviations
	 * @param error
	 *            the error
	 * @param noise
	 *            the noise
	 * @param meanIntensity
	 *            the mean intensity
	 * @param valid
	 *            the valid
	 */
	public Candidate(int x, int y, int index, float[] params, float[] paramDevs, double error, float noise,
			float meanIntensity, boolean valid)
	{
		super(x, y, 0, 0);
		this.index = index;
		this.params = params;
		this.paramDevs = paramDevs;
		this.error = error;
		this.noise = noise;
		this.meanIntensity = meanIntensity;
		this.fit = valid;
	}

	/**
	 * Instantiates a new candidate.
	 *
	 * @param spot
	 *            the spot
	 * @param index
	 *            the index
	 */
	public Candidate(Spot spot, int index)
	{
		super(spot.x, spot.y, spot.intensity, spot.getScore());
		this.index = index;
	}

	/**
	 * Creates a fitted candidate with fitted parameters.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param index
	 *            the index
	 * @param params
	 *            the params
	 * @param paramDevs
	 *            the param deviations
	 * @param error
	 *            the error
	 * @param noise
	 *            the noise
	 * @param meanIntensity
	 *            the mean intensity
	 * @param valid
	 *            the valid
	 * @return the candidate
	 */
	public Candidate createFitted(int x, int y, int index, float[] params, float[] paramDevs, double error, float noise,
			float meanIntensity, boolean valid)
	{
		final Candidate c = new Candidate(x, y, intensity, getScore(), index);
		c.params = params;
		c.paramDevs = paramDevs;
		c.error = error;
		c.noise = noise;
		c.meanIntensity = meanIntensity;
		c.fit = valid;
		return c;
	}
}
