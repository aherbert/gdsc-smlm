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

/**
 * Filter results using the combination of two filters
 */
public abstract class CombinedFilter extends Filter
{
	protected Filter filter1, filter2;

	public CombinedFilter(Filter filter1, Filter filter2)
	{
		this.filter1 = filter1;
		this.filter2 = filter2;
	}

	@Override
	protected String generateName()
	{
		StringBuilder sb = new StringBuilder();
		addText(sb, filter1, filter1.getName());
		sb.append(" ").append(getOperator()).append(" ");
		addText(sb, filter2, filter2.getName());
		return sb.toString();
	}

	private void addText(StringBuilder sb, Filter f, String text)
	{
		if (f instanceof CombinedFilter)
			sb.append("(");
		sb.append(text);
		if (f instanceof CombinedFilter)
			sb.append(")");
	}

	@Override
	protected String generateType()
	{
		StringBuilder sb = new StringBuilder();
		addText(sb, filter1, filter1.getType());
		sb.append(" ").append(getOperator()).append(" ");
		addText(sb, filter2, filter2.getType());
		return sb.toString();
	}

	/**
	 * Get the string representation of the operator used to combine the two filters. This is used in the filter name.
	 * 
	 * @return The operator
	 */
	protected abstract String getOperator();

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#setup(gdsc.smlm.results.MemoryPeakResults)
	 */
	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		filter1.setup(peakResults);
		filter2.setup(peakResults);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getNumericalValue()
	 */
	@Override
	public double getNumericalValue()
	{
		return filter1.getNumericalValue();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getNumericalValueName()
	 */
	@Override
	public String getNumericalValueName()
	{
		return filter1.getNumericalValueName();
	}

	@Override
	public int getNumberOfParameters()
	{
		return filter1.getNumberOfParameters() + filter2.getNumberOfParameters();
	}

	@Override
	public double getParameterValue(int index)
	{
		checkIndex(index);
		if (index < filter1.getNumberOfParameters())
			return filter1.getParameterValue(index);
		return filter2.getParameterValue(index - filter1.getNumberOfParameters());
	}

	@Override
	public String getParameterName(int index)
	{
		checkIndex(index);
		if (index < filter1.getNumberOfParameters())
			return filter1.getParameterName(index);
		return filter2.getParameterName(index - filter1.getNumberOfParameters());
	}
	
	@Override
	public Filter adjustParameter(int index, double delta)
	{
		checkIndex(index);
		Filter f1 = filter1;
		Filter f2 = filter2;
		if (index < filter1.getNumberOfParameters())
			f1 = filter1.adjustParameter(index, delta);
		else
			f2 = filter2.adjustParameter(index - filter1.getNumberOfParameters(), delta);
		return createFilter(f1, f2);
	}

	/**
	 * Create a new combined filter from the two input filters 
	 * @param f1
	 * @param f2
	 * @return
	 */
	protected abstract Filter createFilter(Filter f1, Filter f2);
}