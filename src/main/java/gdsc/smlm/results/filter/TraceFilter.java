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

import java.util.HashSet;
import java.util.Set;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.Trace;
import gdsc.smlm.results.TraceManager;

/**
 * Filter results that can be traced over time frames.
 */
public class TraceFilter extends Filter
{
	private static final double DEFAULT_DISTANCE_INCREMENT = 0.05;
	private static final int DEFAULT_TIME_INCREMENT = 1;
	private static final double DEFAULT_DISTANCE_RANGE = 2;
	private static final int DEFAULT_TIME_RANGE = 10;

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
				for (int i = 0; i < trace.size(); i++)
				{
					ok.add(trace.get(i));
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
		return ParameterType.TIME_THRESHOLD.toString();
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
	 * @see gdsc.smlm.results.filter.Filter#getParameterValueInternal(int)
	 */
	@Override
	protected double getParameterValueInternal(int index)
	{
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
	 * @see gdsc.smlm.results.filter.Filter#getParameterIncrement(int)
	 */
	@Override
	public double getParameterIncrement(int index)
	{
		checkIndex(index);
		switch (index)
		{
			case 0:
				return DEFAULT_DISTANCE_INCREMENT;
			default:
				return DEFAULT_TIME_INCREMENT;
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
				return ParameterType.DISTANCE_THRESHOLD;
			default:
				return ParameterType.TIME_THRESHOLD;
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
	 * @see gdsc.smlm.ga.Chromosome#mutationStepRange()
	 */
	@Override
	public double[] mutationStepRange()
	{
		return new double[] { DEFAULT_DISTANCE_RANGE, DEFAULT_TIME_RANGE };
	}
}
