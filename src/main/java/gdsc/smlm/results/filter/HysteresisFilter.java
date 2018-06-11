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
import java.util.LinkedList;
import java.util.Set;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

import gdsc.core.utils.NotImplementedException;
import gdsc.smlm.data.config.CalibrationReader;
import gdsc.smlm.ga.Chromosome;
import gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.Trace;
import gdsc.smlm.results.TraceManager;
import gdsc.smlm.results.TraceManager.TraceMode;
import gdsc.smlm.results.procedures.PeakResultProcedure;

/**
 * Filter results using a precision threshold. Any results below the lower
 * precision limit are included. Any results above the upper precision limit are
 * excluded. Any results between the limits are included only if they can be
 * traced through time, optionally via other candidates, to a valid result.
 */
public abstract class HysteresisFilter extends Filter
{
	public static final double DEFAULT_ABSOLUTE_DISTANCE_INCREMENT = 5;
	public static final double DEFAULT_RELATIVE_DISTANCE_INCREMENT = 0.05;
	public static final double DEFAULT_SECONDS_TIME_INCREMENT = 0.05;
	public static final double DEFAULT_FRAMES_TIME_INCREMENT = 1;
	public static final double DEFAULT_ABSOLUTE_DISTANCE_RANGE = 200;
	public static final double DEFAULT_RELATIVE_DISTANCE_RANGE = 1;
	public static final double DEFAULT_SECONDS_TIME_RANGE = 5;
	public static final double DEFAULT_FRAMES_TIME_RANGE = 10;

	@XStreamAsAttribute
	final double searchDistance;
	@XStreamAsAttribute
	final int searchDistanceMode;
	@XStreamAsAttribute
	final double timeThreshold;
	@XStreamAsAttribute
	final int timeThresholdMode;
	@XStreamOmitField
	Set<PeakResult> ok;

	protected enum PeakStatus
	{
		OK, CANDIDATE, REJECT
	}

	/**
	 * @param searchDistance
	 * @param searchDistanceMode
	 *            0 = relative to the precision of the candidates; 1 = Absolute (in nm)
	 * @param timeThreshold
	 * @param timeThresholdMode
	 *            0 = frames; 1 = seconds
	 */
	public HysteresisFilter(double searchDistance, int searchDistanceMode, double timeThreshold, int timeThresholdMode)
	{
		this.searchDistance = Math.max(0, searchDistance);
		this.searchDistanceMode = searchDistanceMode;
		this.timeThreshold = Math.max(0, timeThreshold);
		this.timeThresholdMode = timeThresholdMode;
	}

	/**
	 * @return The name of the configured search distance mode
	 */
	public String getSearchName()
	{
		switch (searchDistanceMode)
		{
			case 1:
				return "Absolute";

			case 0:
			default:
				return "Candidate precision";
		}
	}

	protected double getDefaultSearchRange()
	{
		switch (searchDistanceMode)
		{
			case 1:
				return DEFAULT_ABSOLUTE_DISTANCE_RANGE;

			case 0:
			default:
				return DEFAULT_RELATIVE_DISTANCE_RANGE;
		}
	}

	/**
	 * @return The name of the configured time threshold mode
	 */
	public String getTimeName()
	{
		switch (timeThresholdMode)
		{
			case 1:
				return "Seconds";

			case 0:
			default:
				return "Frames";
		}
	}

	protected double getDefaultTimeRange()
	{
		switch (timeThresholdMode)
		{
			case 1:
				return DEFAULT_SECONDS_TIME_RANGE;

			case 0:
			default:
				return DEFAULT_FRAMES_TIME_RANGE;
		}
	}

	protected String getTraceParameters()
	{
		return String.format("@%.2f %s, %.2f %s", searchDistance, getSearchName(), timeThreshold, getTimeName());
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
				return searchDistance;
			case 1:
				return searchDistanceMode;
			case 2:
				return timeThreshold;
			default:
				return timeThresholdMode;
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
				return (searchDistanceMode == 1) ? DEFAULT_RELATIVE_DISTANCE_INCREMENT
						: DEFAULT_ABSOLUTE_DISTANCE_INCREMENT;
			case 1:
				return 1;
			case 2:
				return (timeThresholdMode == 1) ? DEFAULT_SECONDS_TIME_INCREMENT : DEFAULT_FRAMES_TIME_INCREMENT;
			default:
				return 1;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getDisabledParameterValue(int)
	 */
	@Override
	public double getDisabledParameterValue(int index)
	{
		throw new NotImplementedException("Parameters in hysteresis filters cannot be disabled");
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
			case 1:
				return ParameterType.DISTANCE_THRESHOLD_MODE;
			case 2:
				return ParameterType.TIME_THRESHOLD;
			default:
				return ParameterType.TIME_THRESHOLD_MODE;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#weakestParameters(double[])
	 */
	@Override
	public void weakestParameters(double[] parameters)
	{
		setMax(parameters, 0, searchDistance);
		setMax(parameters, 2, timeThreshold);
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		ok = new HashSet<PeakResult>();

		// Create a set of candidates and valid peaks
		final MemoryPeakResults traceResults = new MemoryPeakResults();

		// Initialise peaks to check
		final LinkedList<PeakResult> candidates = new LinkedList<PeakResult>();
		peakResults.forEach(new PeakResultProcedure()
		{
			@Override
			public void execute(PeakResult result)
			{
				switch (getStatus(result))
				{
					case OK:
						ok.add(result);
						traceResults.add(result);
						break;
					case CANDIDATE:
						candidates.add(result);
						traceResults.add(result);
						break;
					default:
						break;
				}
			}
		});

		if (candidates.isEmpty())
		{
			// No candidates for tracing so just return
			return;
		}

		double distanceThreshold;
		switch (searchDistanceMode)
		{
			case 1:
				distanceThreshold = searchDistance / peakResults.getNmPerPixel();
				break;

			case 0:
			default:
				distanceThreshold = getSearchDistanceUsingCandidates(peakResults, candidates);
		}

		if (distanceThreshold <= 0)
			return;

		int myTimeThreshold;
		switch (timeThresholdMode)
		{
			case 1:
				myTimeThreshold = 1;
				if (peakResults.hasCalibration())
				{
					CalibrationReader cr = peakResults.getCalibrationReader();
					double et = cr.getExposureTime();
					if (et > 0)
						myTimeThreshold = (int) Math.round((this.timeThreshold / et));
				}
				else

					break;

			case 0:
			default:
				myTimeThreshold = (int) this.timeThreshold;
		}

		if (myTimeThreshold <= 0)
			return;

		// Trace through candidates
		TraceManager tm = new TraceManager(traceResults);
		tm.setTraceMode(TraceMode.LATEST_FORERUNNER);
		tm.traceMolecules(distanceThreshold, myTimeThreshold);
		Trace[] traces = tm.getTraces();

		for (Trace trace : traces)
		{
			if (trace.size() > 1)
			{
				// Check if the trace touches a valid point
				boolean isOk = false;
				for (int i = 0; i < trace.size(); i++)
				{
					if (ok.contains(trace.get(i)))
					{
						isOk = true;
						break;
					}
				}
				// Add the entire trace to the OK points
				if (isOk)
				{
					for (int i = 0; i < trace.size(); i++)
					{
						ok.add(trace.get(i));
					}
				}
			}
		}
	}

	/**
	 * Find average precision of the candidates and use it for the search
	 * distance
	 * 
	 * @param peakResults
	 * @param candidates
	 * @return
	 */
	private double getSearchDistanceUsingCandidates(MemoryPeakResults peakResults, LinkedList<PeakResult> candidates)
	{
		Gaussian2DPeakResultCalculator calculator = Gaussian2DPeakResultHelper.create(peakResults.getPSF(),
				peakResults.getCalibration(), Gaussian2DPeakResultHelper.LSE_PRECISION);
		double sum = 0;
		for (PeakResult peakResult : candidates)
		{
			sum += calculator.getLSEPrecision(peakResult.getParameters(), peakResult.getNoise());
		}
		final double nmPerPixel = peakResults.getNmPerPixel();
		double distanceThreshold = (sum / candidates.size()) * searchDistance / nmPerPixel;
		return distanceThreshold;
	}

	protected abstract PeakStatus getStatus(PeakResult result);

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

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#end()
	 */
	@Override
	public void end()
	{
		ok.clear();
		ok = null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getDescription()
	 */
	@Override
	public String getDescription()
	{
		return "Any results between the limits (candidates) are included only if they can be traced " +
				"through time, potentially via other candidates, to a valid result. The distance used for " +
				"tracing is the search distance multiplied by the average precision of the candidates.";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#subsetWithFailCount()
	 */
	@Override
	public boolean subsetWithFailCount()
	{
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#newChromosome(double[])
	 */
	@Override
	public Chromosome<FilterScore> newChromosome(double[] sequence)
	{
		// Hysteresis filters remove their search and time mode parameters in their Chromosome sequence
		// so add it back
		double[] parameters = new double[sequence.length];
		parameters[0] = sequence[0];
		parameters[1] = searchDistanceMode;
		parameters[2] = sequence[1];
		parameters[3] = timeThresholdMode;
		System.arraycopy(sequence, 2, parameters, 4, sequence.length - 2);
		return create(parameters);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getChromosomeParameters()
	 */
	@Override
	public int[] getChromosomeParameters()
	{
		// Hysteresis filters remove their search and time mode parameters in their Chromosome sequence
		// Skip the search mode [param 1]
		// Skip the time mode [param 3]
		int[] indices = new int[length()];
		indices[1] = 2;
		for (int i = 2; i < indices.length; i++)
			indices[i] = i + 2;
		return indices;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#length()
	 */
	@Override
	public int length()
	{
		// Hysteresis filters remove their search and time mode parameters in their Chromosome sequence
		return getNumberOfParameters() - 2;
	}

	@Override
	public double[] sequence()
	{
		// Remind derived classes to implement this.
		//throw new NotImplementedException();
		// Implement a default version using the results of getChromosomeParameters() and getParameters().
		double[] sequence = new double[length()];
		double[] params = getParameters();
		int[] indices = getChromosomeParameters();
		for (int i = 0; i < indices.length; i++)
			sequence[i] = params[indices[i]];
		return sequence;
	}
}
