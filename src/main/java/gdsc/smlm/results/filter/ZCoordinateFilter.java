package gdsc.smlm.results.filter;

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

import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;

/**
 * Filter results using a z-coordinate range
 */
public class ZCoordinateFilter extends DirectFilter
{
	// Assuming units are in pixels with 100nm/px set the default range as +/- 1000nm
	public static final double DEFAULT_INCREMENT = 0.1;
	public static final double DEFAULT_RANGE = 10;
	public static final double UPPER_LIMIT = 50; // This may need to be changed

	@XStreamAsAttribute
	final float minZ;
	@XStreamAsAttribute
	final float maxZ;

	public ZCoordinateFilter(float minZ, float maxZ)
	{
		if (maxZ < minZ)
		{
			float f = maxZ;
			maxZ = minZ;
			minZ = f;
		}
		this.minZ = minZ;
		this.maxZ = maxZ;
	}

	@Override
	protected String generateName()
	{
		return "Z " + minZ + "-" + maxZ;
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
	}

	@Override
	public boolean accept(PeakResult peak)
	{
		return peak.getZPosition() >= minZ && peak.getZPosition() <= maxZ;
	}

	public int getValidationFlags()
	{
		return V_Z;
	}

	@Override
	public int validate(final PreprocessedPeakResult peak)
	{
		if (peak.getZ() < minZ || peak.getZ() > maxZ)
			return V_Z;
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getDescription()
	 */
	@Override
	public String getDescription()
	{
		return "Filter results using a z-coordinate range.";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getNumberOfParameters()
	 */
	@Override
	public int getNumberOfParameters()
	{
		return 2;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getParameterValueInternal(int)
	 */
	@Override
	protected double getParameterValueInternal(int index)
	{
		switch (index)
		{
			case 0:
				return minZ;
			default:
				return maxZ;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getParameterIncrement(int)
	 */
	@Override
	public double getParameterIncrement(int index)
	{
		checkIndex(index);
		return DEFAULT_INCREMENT;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getDisabledParameterValue(int)
	 */
	@Override
	public double getDisabledParameterValue(int index)
	{
		checkIndex(index);
		switch (index)
		{
			case 0:
				return Double.NEGATIVE_INFINITY;
			default:
				return Double.POSITIVE_INFINITY;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getParameterType(int)
	 */
	@Override
	public ParameterType getParameterType(int index)
	{
		checkIndex(index);
		switch (index)
		{
			case 0:
				return ParameterType.MIN_Z;
			default:
				return ParameterType.MAX_Z;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#adjustParameter(int, double)
	 */
	@Override
	public Filter adjustParameter(int index, double delta)
	{
		checkIndex(index);
		switch (index)
		{
			case 0:
				return new ZCoordinateFilter(updateParameter(minZ, delta, DEFAULT_RANGE), maxZ);
			default:
				return new ZCoordinateFilter(minZ, updateParameter(maxZ, delta, DEFAULT_RANGE));
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#create(double[])
	 */
	@Override
	public Filter create(double... parameters)
	{
		return new ZCoordinateFilter((float) parameters[0], (float) parameters[1]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#weakestParameters(double[])
	 */
	@Override
	public void weakestParameters(double[] parameters)
	{
		setMin(parameters, 0, minZ);
		setMax(parameters, 1, maxZ);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.DirectFilter#lowerBoundOrientation(int)
	 */
	@Override
	public int lowerBoundOrientation(int index)
	{
		return (index == 1) ? 1 : -1;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#length()
	 */
	public int length()
	{
		return 2;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#sequence()
	 */
	public double[] sequence()
	{
		// Ignore the mode parameters
		return new double[] { minZ, maxZ };
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#mutationStepRange()
	 */
	public double[] mutationStepRange()
	{
		return new double[] { DEFAULT_RANGE, DEFAULT_RANGE };
	}
}