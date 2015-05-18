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
import gdsc.smlm.results.Trace;
import gdsc.smlm.results.TraceManager;

import java.util.HashSet;
import java.util.Set;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

/**
 * Filter results that can be traced over time frames.
 */
public class TraceFilter extends Filter
{
	static double DEFAULT_DISTANCE_RANGE = 2;
	static int DEFAULT_TIME_RANGE = 10;
	
	@XStreamAsAttribute
	final double d;
	@XStreamAsAttribute
	final int t;
	@XStreamOmitField
	Set<PeakResult> ok;

	public TraceFilter(double d, int t)
	{
		this.d = Math.max(0, d);
		this.t = Math.max(0, t);
	}

	@Override
	protected String generateName()
	{
		return String.format("Trace d=%.2f, t=%d", d, t);
	}

	@Override
	protected String generateType()
	{
		return "Trace";
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		ok = new HashSet<PeakResult>();

		// Trace molecules. Anything that is part of a trace is OK
		TraceManager tm = new TraceManager(peakResults);
		tm.traceMolecules(d, t);
		Trace[] traces = tm.getTraces();
		for (Trace trace : traces)
		{
			if (trace.size() > 1)
			{
				for (PeakResult result : trace.getPoints())
				{
					ok.add(result);
				}
			}
		}
	}

	/**
	 * @throws NullPointerException
	 *             if not first initialised with a call to {@link #setup(MemoryPeakResults)}
	 * @see gdsc.smlm.results.filter.Filter#accept(gdsc.smlm.results.PeakResult)
	 */
	@Override
	public boolean accept(PeakResult peak) throws NullPointerException
	{
		return ok.contains(peak);
	}

	@Override
	public double getNumericalValue()
	{
		return t;
	}

	@Override
	public String getNumericalValueName()
	{
		return "Time";
	}

	@Override
	public String getDescription()
	{
		return "Filter results that can be traced over time frames.";
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
	 * @see gdsc.smlm.results.filter.Filter#getParameterValue(int)
	 */
	@Override
	public double getParameterValue(int index)
	{
		checkIndex(index);
		switch (index)
		{
			case 0:
				return d;
			default:
				return t;
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
				return "d-threshold";
			default:
				return "t-threshold";
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
				return new TraceFilter(updateParameter(d, delta, DEFAULT_DISTANCE_RANGE), t);
			default:
				return new TraceFilter(d, updateParameter(t, delta, DEFAULT_TIME_RANGE));
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
		return new TraceFilter(parameters[0], (int) parameters[1]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#weakestParameters(double[])
	 */
	@Override
	public void weakestParameters(double[] parameters)
	{
		setMax(parameters, 0, d);
		setMax(parameters, 1, t);
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
		return new double[] { d, t };
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#mutationStepRange()
	 */
	public double[] mutationStepRange()
	{
		return new double[] { DEFAULT_DISTANCE_RANGE, DEFAULT_TIME_RANGE };
	}
}