package gdsc.smlm.results.filter;

import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.match.MatchResult;

import java.util.List;

import com.thoughtworks.xstream.annotations.XStreamOmitField;

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
 * Filter a set of peak results into accepted/rejected.
 */
public abstract class Filter implements Comparable<Filter>
{
	@XStreamOmitField
	private String name;
	@XStreamOmitField
	private String type;

	/**
	 * Generate the name of the filter using the filter settings
	 * 
	 * @return The name of the filter
	 */
	protected abstract String generateName();

	/**
	 * Generate the type of the filter using the filter settings
	 * 
	 * @return The type of the filter
	 */
	protected abstract String generateType();

	/**
	 * Filter the results
	 * 
	 * @param results
	 * @return the filtered results
	 */
	public MemoryPeakResults filter(MemoryPeakResults results)
	{
		MemoryPeakResults newResults = new MemoryPeakResults();
		newResults.copySettings(results);
		setup(results);
		for (PeakResult peak : results.getResults())
		{
			if (accept(peak))
				newResults.add(peak);
		}
		return newResults;
	}

	/**
	 * Filter the results and return the performance score. Allows benchmarking the filter by marking the results as
	 * true or false.
	 * <p>
	 * Any input PeakResult with an original value that is not zero will be treated as a true result, all other results
	 * are false. The filter is run and the results are marked as true positive, false negative and false positive.
	 * 
	 * @param resultsList
	 *            a list of results to analyse
	 * @return the score
	 */
	public MatchResult score(List<MemoryPeakResults> resultsList)
	{
		int tp = 0, fp = 0, fn = 0;
		for (MemoryPeakResults peakResults : resultsList)
		{
			setup(peakResults);

			for (PeakResult peak : peakResults.getResults())
			{
				boolean isTrue = peak.origValue != 0;
				boolean isPositive = accept(peak);
				if (isTrue)
				{
					if (isPositive)
						tp++; // true positive
					else
						fn++; // false negative
				}
				else if (isPositive)
					fp++; // false positive
				// Default: true negative
			}
		}
		return new MatchResult(tp, fp, fn, 0);
	}

	/**
	 * Called before the accept method is called for each peak in the results. Allows pre-processing of the results.
	 * 
	 * @param peakResults
	 */
	public abstract void setup(MemoryPeakResults peakResults);

	/**
	 * Called for each peak in the results that are filtered.
	 * 
	 * @param peak
	 * @return true if the peak should be accepted, otherwise false to reject.
	 */
	public abstract boolean accept(PeakResult peak);

	/**
	 * @return The numerical value of the filter. Used for plotting value against performance score.
	 */
	public abstract double getNumericalValue();

	/**
	 * @return The name of the numerical value of the filter. Used for plotting value against performance score.
	 */
	public abstract String getNumericalValueName();

	/**
	 * @return the name
	 */
	public String getName()
	{
		if (name == null)
			name = generateName();
		return name;
	}

	/**
	 * @return the type (excluding any parameter values)
	 */
	public String getType()
	{
		if (type == null)
			type = generateType();
		return type;
	}

	/**
	 * @return Describes the functionality of the filter
	 */
	public abstract String getDescription();

	/**
	 * @return An XML representation of this object
	 */
	public String toXML()
	{
		return XStreamWrapper.toXML(this);
	}

	/**
	 * Create the filter from the XML representation
	 * 
	 * @param xml
	 * @return the filter
	 */
	public static Filter fromXML(String xml)
	{
		try
		{
			return (Filter) XStreamWrapper.fromXML(xml);
		}
		catch (ClassCastException ex)
		{
			//ex.printStackTrace();
		}
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(Filter o)
	{
		// Null to end of list
		if (o == null)
			return -1;
		final double v1 = getNumericalValue();
		final double v2 = o.getNumericalValue();
		if (v1 < v2)
			return -1;
		if (v1 > v2)
			return 1;
		return 0;
	}

	/**
	 * @return The number of parameters for the filter
	 */
	public abstract int getNumberOfParameters();

	protected void checkIndex(final int index)
	{
		if (index < 0 || index >= getNumberOfParameters())
			throw new IndexOutOfBoundsException("Index must be >= 0 and < " + getNumberOfParameters());
	}

	/**
	 * @param index
	 * @return The value of the specified parameter
	 */
	public abstract double getParameterValue(int index);

	/**
	 * @param index
	 * @return The name of the specified parameter
	 */
	public abstract String getParameterName(int index);

	/**
	 * Create a new filter by adjusting the specified parameter.
	 * <p>
	 * A positive delta will adjust the parameter to be larger. A negative delta will adjust the parameter to be
	 * smaller. The adjustment is relative to the parameter value, e.g. 0.1 is 10%.
	 * 
	 * @param index
	 *            The parameter index
	 * @param delta
	 *            The amount to adjust the parameter
	 * @return The new filter
	 */
	public abstract Filter adjustParameter(int index, double delta);
	
	/**
	 * Adjust the specified parameter value.
	 * <p>
	 * A positive delta will adjust the parameter to be larger. A negative delta will adjust the parameter to be
	 * smaller. The adjustment is relative to the parameter value, e.g. 0.1 is 10%.
	 * 
	 * @param value
	 * @param delta
	 * @return
	 */
	protected double updateParameter(double value, double delta)
	{
		return value + value * delta;
	}

	/**
	 * Adjust the specified parameter value.
	 * <p>
	 * A positive delta will adjust the parameter to be larger. A negative delta will adjust the parameter to be
	 * smaller. The adjustment is relative to the parameter value, e.g. 0.1 is 10%.
	 * 
	 * @param value
	 * @param delta
	 * @return
	 */
	protected float updateParameter(float value, double delta)
	{
		return (float) (value + value * delta);
	}
	
	/**
	 * Adjust the specified parameter value.
	 * <p>
	 * A positive delta will adjust the parameter to be larger. A negative delta will adjust the parameter to be
	 * smaller. The adjustment is relative to the parameter value, e.g. 0.1 is 10%. The adjustment is rounded
	 * up to the next valid integer to ensure a new parameter value is created.
	 * 
	 * @param value
	 * @param delta
	 * @return
	 */
	protected int updateParameter(int value, double delta)
	{
		final int update = (int)Math.ceil(value * Math.abs(delta));
		if (delta < 0)
			return value - update;
		return value + update;
	}
}