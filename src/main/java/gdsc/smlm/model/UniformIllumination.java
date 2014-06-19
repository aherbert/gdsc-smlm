package gdsc.smlm.model;

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
 * Specifies the same illumination for any position
 */
public class UniformIllumination implements SpatialIllumination
{
	private double photons;
	private double pulsePhotons;
	private int pulseInterval;

	/**
	 * @param photons
	 *            The number of photons in a time frame
	 */
	public UniformIllumination(double photons)
	{
		this(photons, 0, 0);
	}

	/**
	 * @param photons
	 *            The number of photons in a time frame
	 * @param pulsePhotons
	 *            The number of photons in a pulse
	 * @param pulseInterval
	 *            The interval between pulses (t=1 is the first pulse). Must be above 1.
	 */
	public UniformIllumination(double photons, double pulsePhotons, int pulseInterval)
	{
		this.photons = photons;
		this.pulsePhotons = pulsePhotons;
		this.pulseInterval = pulseInterval;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialIllumination#getPhotons(double[])
	 */
	public double getPhotons(double[] xyz)
	{
		return photons;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialIllumination#getPulsedPhotons(double[], int)
	 */
	public double[] getPulsedPhotons(double[] xyz, int t)
	{

		if (pulseInterval > 1)
		{
			return new double[] { (t % pulseInterval == 1) ? pulsePhotons : 0, photons };
		}
		return new double[] { 0, photons };
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.SpatialIllumination#getAveragePhotons()
	 */
	public double getAveragePhotons()
	{
		if (pulseInterval > 1)
		{
			return photons + pulsePhotons / pulseInterval;
		}
		return photons;
	}
}
