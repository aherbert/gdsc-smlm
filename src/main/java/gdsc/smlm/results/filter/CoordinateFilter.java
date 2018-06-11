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
package gdsc.smlm.results.filter;

import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;

/**
 * Filter results using a coordinate range
 */
public class CoordinateFilter extends DirectFilter
{
	public static final double DEFAULT_INCREMENT = 0.01;
	public static final double DEFAULT_RANGE = 1;

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
	public int getValidationFlags()
	{
		return V_X | V_Y;
	}

	@Override
	public int validate(final PreprocessedPeakResult peak)
	{
		if (peak.getX() < minX || peak.getX() > maxX)
			return V_X;
		if (peak.getY() < minY || peak.getY() > maxY)
			return V_Y;
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
	 * @see gdsc.smlm.results.filter.Filter#getParameterValueInternal(int)
	 */
	@Override
	protected double getParameterValueInternal(int index)
	{
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
			case 1:
				return Double.POSITIVE_INFINITY;
			case 2:
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
				return ParameterType.MIN_X;
			case 1:
				return ParameterType.MAX_X;
			case 2:
				return ParameterType.MIN_Y;
			default:
				return ParameterType.MAX_Y;
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
	 * @see gdsc.smlm.results.filter.DirectFilter#lowerBoundOrientation(int)
	 */
	@Override
	public int lowerBoundOrientation(int index)
	{
		return (index == 1 || index == 3) ? 1 : -1;
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
