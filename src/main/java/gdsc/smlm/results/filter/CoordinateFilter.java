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
 * Filter results using a coordinate range
 */
public class CoordinateFilter extends Filter
{
	static double DEFAULT_RANGE = 1;

	@XStreamAsAttribute
	final float minX;
	@XStreamAsAttribute
	final float maxX;
	@XStreamAsAttribute
	final float minY;
	@XStreamAsAttribute
	final float maxY;

	public CoordinateFilter(float minX, float maxX, float minY, float maxY)
	{
		if (maxX < minX)
		{
			float f = maxX;
			maxX = minX;
			minX = f;
		}
		this.minX = minX;
		this.maxX = maxX;
		if (maxY < minY)
		{
			float f = maxY;
			maxY = minY;
			minY = f;
		}
		this.minY = minY;
		this.maxY = maxY;
	}

	@Override
	protected String generateName()
	{
		return "X " + minX + "-" + maxX + ", Y " + minY + "-" + maxY;
	}

	@Override
	protected String generateType()
	{
		return "Coordinate";
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
	}

	@Override
	public boolean accept(PeakResult peak)
	{
		return peak.getXPosition() >= minX && peak.getXPosition() <= maxX && peak.getYPosition() >= minY &&
				peak.getYPosition() <= maxY;
	}

	@Override
	public double getNumericalValue()
	{
		return minX;
	}

	@Override
	public String getNumericalValueName()
	{
		return "X";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getDescription()
	 */
	@Override
	public String getDescription()
	{
		return "Filter results using a coordinate range.";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getNumberOfParameters()
	 */
	@Override
	public int getNumberOfParameters()
	{
		return 4;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getParameterValue(int)
	 */
	@Override
	public double getParameterValue(int index)
	{
		checkIndex(index);
		switch (index)
		{
			case 0:
				return minX;
			case 1:
				return maxX;
			case 2:
				return minY;
			default:
				return maxY;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getParameterName(int)
	 */
	@Override
	public String getParameterName(int index)
	{
		checkIndex(index);
		switch (index)
		{
			case 0:
				return "Min X";
			case 1:
				return "Max X";
			case 2:
				return "Min Y";
			default:
				return "Max Y";
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
				return new CoordinateFilter(updateParameter(minX, delta, DEFAULT_RANGE), maxX, minY, maxY);
			case 1:
				return new CoordinateFilter(minX, updateParameter(maxX, delta, DEFAULT_RANGE), minY, maxY);
			case 2:
				return new CoordinateFilter(minX, maxX, updateParameter(minY, delta, DEFAULT_RANGE), maxY);
			default:
				return new CoordinateFilter(minX, maxX, minY, updateParameter(maxY, delta, DEFAULT_RANGE));
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
		return new CoordinateFilter((float) parameters[0], (float) parameters[1], (float) parameters[2],
				(float) parameters[3]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#weakestParameters(double[])
	 */
	@Override
	public void weakestParameters(double[] parameters)
	{
		setMin(parameters, 0, minX);
		setMax(parameters, 1, maxX);
		setMin(parameters, 2, minY);
		setMax(parameters, 3, maxY);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#length()
	 */
	@Override
	public int length()
	{
		return 4;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#sequence()
	 */
	@Override
	public double[] sequence()
	{
		// Ignore the mode parameters
		return new double[] { minX, maxX, minY, maxY };
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#mutationStepRange()
	 */
	@Override
	public double[] mutationStepRange()
	{
		return new double[] { DEFAULT_RANGE, DEFAULT_RANGE, DEFAULT_RANGE, DEFAULT_RANGE };
	}
}