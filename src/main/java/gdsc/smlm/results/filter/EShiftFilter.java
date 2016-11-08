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

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

/**
 * Filter results using a X/Y coordinate shift as a Euclidian distance. This filter requires that the result X and Y
 * coordinates are reported relative to their initial positions.
 */
public class EShiftFilter extends DirectFilter implements IMultiFilter
{
	public static final double DEFAULT_INCREMENT = 0.05;
	public static final double DEFAULT_RANGE = 10;
	public static final double UPPER_LIMIT = 5;

	@XStreamAsAttribute
	final double eshift;
	@XStreamOmitField
	float eoffset;
	@XStreamOmitField
	float eshift2;
	@XStreamOmitField
	boolean shiftEnabled;

	public EShiftFilter(double eshift)
	{
		this.eshift = Math.max(0, eshift);
	}

	@Override
	protected String generateName()
	{
		return "EShift " + eshift;
	}

	@Override
	protected String generateType()
	{
		return "EShift";
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		// Set the shift limit
		eoffset = Float.POSITIVE_INFINITY;
		Pattern pattern = Pattern.compile("initialSD0>([\\d\\.]+)");
		Matcher match = pattern.matcher(peakResults.getConfiguration());
		if (match.find())
		{
			// Convert to squared distance
			eoffset = getUpperSquaredLimit(Double.parseDouble(match.group(1)) * eshift);
		}
	}

	@Override
	public void setup()
	{
		setup(true);
	}

	@Override
	public void setup(int flags)
	{
		setup(!areSet(flags, DirectFilter.NO_SHIFT));
	}

	private void setup(final boolean shiftEnabled)
	{
		this.shiftEnabled = shiftEnabled;
		if (shiftEnabled)
		{
			eshift2 = getUpperSquaredLimit(eshift);
		}
	}

	@Override
	public boolean accept(PeakResult peak)
	{
		final float dx = peak.getXPosition();
		final float dy = peak.getYPosition();
		return dx * dx + dy * dy <= eoffset;
	}

	@Override
	public int validate(final PreprocessedPeakResult peak)
	{
		if (shiftEnabled)
		{
			if (peak.getXRelativeShift2() + peak.getYRelativeShift2() > eshift2)
				return V_X_RELATIVE_SHIFT | V_Y_RELATIVE_SHIFT;
		}
		return 0;
	}

	@Override
	public double getNumericalValue()
	{
		return eshift;
	}

	@Override
	public String getNumericalValueName()
	{
		return "EShift";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getDescription()
	 */
	@Override
	public String getDescription()
	{
		return "Filter results using a Euclidian shift factor. (Euclidian shift is relative to initial peak width.)";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getNumberOfParameters()
	 */
	@Override
	public int getNumberOfParameters()
	{
		return 1;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getParameterValueInternal(int)
	 */
	@Override
	protected double getParameterValueInternal(int index)
	{
		return eshift;
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
	 * @see gdsc.smlm.results.filter.Filter#getParameterName(int)
	 */
	@Override
	public String getParameterName(int index)
	{
		checkIndex(index);
		return "Shift";
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
		return new EShiftFilter(updateParameter(eshift, delta, DEFAULT_RANGE));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#create(double[])
	 */
	@Override
	public Filter create(double... parameters)
	{
		return new EShiftFilter(parameters[0]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#weakestParameters(double[])
	 */
	@Override
	public void weakestParameters(double[] parameters)
	{
		setMax(parameters, 0, eshift);
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#length()
	 */
	public int length()
	{
		return 1;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#upperLimit()
	 */
	@Override
	public double[] upperLimit()
	{
		return new double[] { UPPER_LIMIT };
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#sequence()
	 */
	public double[] sequence()
	{
		return new double[] { eshift };
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#mutationStepRange()
	 */
	public double[] mutationStepRange()
	{
		return new double[] { DEFAULT_RANGE };
	}

	public double getSignal()
	{
		return 0;
	}

	public double getSNR()
	{
		return 0;
	}

	public double getMinWidth()
	{
		return 0;
	}

	public double getMaxWidth()
	{
		return 0;
	}

	public double getShift()
	{
		return 0;
	}

	public double getEShift()
	{
		return eshift;
	}

	public double getPrecision()
	{
		return 0;
	}

	public boolean isPrecisionUsesLocalBackground()
	{
		return false;
	}
}