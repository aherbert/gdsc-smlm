package gdsc.smlm.results.filter;

import gdsc.smlm.ga.Chromosome;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.match.ClassificationResult;
import gdsc.smlm.results.match.FractionClassificationResult;

import java.util.List;

import org.apache.commons.math3.util.FastMath;

import com.thoughtworks.xstream.annotations.XStreamOmitField;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2015 Alex Herbert
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
public abstract class Filter implements Comparable<Filter>, Chromosome
{
	@XStreamOmitField
	private String name;
	@XStreamOmitField
	private String type;
	@XStreamOmitField
	private double fitness;

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
		end();
		return newResults;
	}

	/**
	 * Filter the results
	 * <p>
	 * The number of consecutive rejections are counted per frame. When the configured number of failures is reached all
	 * remaining results for the frame are rejected. This assumes the results are ordered by the frame.
	 * 
	 * @param results
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @return the filtered results
	 */
	public MemoryPeakResults filter(MemoryPeakResults results, int failures)
	{
		MemoryPeakResults newResults = new MemoryPeakResults();
		newResults.copySettings(results);
		setup(results);
		int frame = -1;
		int failCount = 0;
		for (PeakResult peak : results.getResults())
		{
			if (frame != peak.peak)
			{
				frame = peak.peak;
				failCount = 0;
			}

			// Reject all peaks if we have exceeded the fail count
			final boolean isPositive;
			if (failCount > failures)
			{
				isPositive = false;
			}
			else
			{
				// Otherwise assess the peak
				isPositive = accept(peak);
			}

			if (isPositive)
			{
				failCount = 0;
				newResults.add(peak);
			}
			else
			{
				failCount++;
			}
		}
		end();
		return newResults;
	}

	/**
	 * Filter the results
	 * <p>
	 * Input PeakResults must be allocated a score for true positive, false positive, true negative and false negative
	 * (accessed via the object property get methods). The filter is run and results that pass accumulate scores for
	 * true positive and false positive, otherwise the scores are accumulated for true negative and false negative. The
	 * simplest scoring scheme is to mark valid results as tp=fn=1 and fp=tn=0 and invalid results the opposite.
	 * <p>
	 * The number of consecutive rejections are counted per frame. When the configured number of failures is reached all
	 * remaining results for the frame are rejected. This assumes the results are ordered by the frame.
	 * <p>
	 * The number of failures before each peak is stored in the origX property of the PeakResult.
	 * 
	 * @param results
	 * @param score
	 *            If not null will be populated with the fraction score [ tp, fp, tn, fn ]
	 * @return the filtered results
	 */
	public MemoryPeakResults filterSubset(MemoryPeakResults results, double[] score)
	{
		MemoryPeakResults newResults = new MemoryPeakResults();
		newResults.copySettings(results);
		setup(results);
		int frame = -1;
		int failCount = 0;
		double fp = 0, fn = 0;
		double tp = 0, tn = 0;
		for (PeakResult peak : results.getResults())
		{
			if (frame != peak.peak)
			{
				frame = peak.peak;
				failCount = 0;
			}

			// Reject all peaks if we have exceeded the fail count
			final boolean isPositive = accept(peak);

			if (isPositive)
			{
				peak.origX = failCount;
				failCount = 0;
				newResults.add(peak);
			}
			else
			{
				failCount++;
			}

			if (isPositive)
			{
				tp += peak.getTruePositiveScore();
				fp += peak.getFalsePositiveScore();
			}
			else
			{
				fn += peak.getFalseNegativeScore();
				tn += peak.getTrueNegativeScore();
			}
		}
		end();

		if (score != null && score.length > 3)
		{
			score[0] = tp;
			score[1] = fp;
			score[2] = tn;
			score[3] = fn;
		}

		return newResults;
	}

	/**
	 * Filter the results
	 * <p>
	 * Input PeakResults must be allocated a score for true positive, false positive, true negative and false negative
	 * (accessed via the object property get methods). The filter is run and results that pass accumulate scores for
	 * true positive and false positive, otherwise the scores are accumulated for true negative and false negative. The
	 * simplest scoring scheme is to mark valid results as tp=fn=1 and fp=tn=0 and invalid results the opposite.
	 * <p>
	 * The number of consecutive rejections are counted per frame. When the configured number of failures is reached all
	 * remaining results for the frame are rejected. This assumes the results are ordered by the frame.
	 * <p>
	 * The number of failures before each peak is stored in the origX property of the PeakResult.
	 * 
	 * @param results
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @param score
	 *            If not null will be populated with the fraction score [ tp, fp, tn, fn ]
	 * @return the filtered results
	 */
	public MemoryPeakResults filterSubset(MemoryPeakResults results, int failures, double[] score)
	{
		MemoryPeakResults newResults = new MemoryPeakResults();
		newResults.copySettings(results);
		setup(results);
		int frame = -1;
		int failCount = 0;
		double fp = 0, fn = 0;
		double tp = 0, tn = 0;
		for (PeakResult peak : results.getResults())
		{
			if (frame != peak.peak)
			{
				frame = peak.peak;
				failCount = 0;
			}

			// Reject all peaks if we have exceeded the fail count
			final boolean isPositive;
			if (failCount > failures)
			{
				isPositive = false;
			}
			else
			{
				// Otherwise assess the peak
				isPositive = accept(peak);
			}

			if (isPositive)
			{
				peak.origX = failCount;
				failCount = 0;
				newResults.add(peak);
			}
			else
			{
				failCount++;
			}

			if (isPositive)
			{
				tp += peak.getTruePositiveScore();
				fp += peak.getFalsePositiveScore();
			}
			else
			{
				fn += peak.getFalseNegativeScore();
				tn += peak.getTrueNegativeScore();
			}
		}
		end();

		if (score != null && score.length > 3)
		{
			score[0] = tp;
			score[1] = fp;
			score[2] = tn;
			score[3] = fn;
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
	public ClassificationResult score(List<MemoryPeakResults> resultsList)
	{
		int tp = 0, fp = 0, tn = 0, fn = 0;
		for (MemoryPeakResults peakResults : resultsList)
		{
			setup(peakResults);
			for (PeakResult peak : peakResults.getResults())
			{
				final boolean isTrue = peak.origValue != 0;
				boolean isPositive = accept(peak);
				if (isTrue)
				{
					if (isPositive)
						tp++; // true positive
					else
						fn++; // false negative
				}
				else
				{
					if (isPositive)
						fp++; // false positive
					else
						tn++; // true negative
				}
			}
			end();
		}
		return new ClassificationResult(tp, fp, tn, fn);
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
	 * @param tn
	 *            The initial true negatives (used when the results have been pre-filtered)
	 * @param fn
	 *            The initial false negatives (used when the results have been pre-filtered)
	 * @return
	 */
	public ClassificationResult score(List<MemoryPeakResults> resultsList, int tn, int fn)
	{
		int tp = 0, fp = 0;
		for (MemoryPeakResults peakResults : resultsList)
		{
			setup(peakResults);
			for (PeakResult peak : peakResults.getResults())
			{
				final boolean isTrue = peak.origValue != 0;
				boolean isPositive = accept(peak);
				if (isTrue)
				{
					if (isPositive)
						tp++; // true positive
					else
						fn++; // false negative
				}
				else
				{
					if (isPositive)
						fp++; // false positive
					else
						tn++; // true negative
				}
			}
			end();
		}
		return new ClassificationResult(tp, fp, tn, fn);
	}

	/**
	 * Filter the results and return the performance score. Allows benchmarking the filter by marking the results as
	 * true or false.
	 * <p>
	 * Any input PeakResult with an original value that is not zero will be treated as a true result, all other results
	 * are false. The filter is run and the results are marked as true positive, false negative and false positive.
	 * <p>
	 * The number of consecutive rejections are counted per frame. When the configured number of failures is reached all
	 * remaining results for the frame are rejected. This assumes the results are ordered by the frame.
	 * 
	 * @param resultsList
	 *            a list of results to analyse
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @return the score
	 */
	public ClassificationResult score(List<MemoryPeakResults> resultsList, int failures)
	{
		int tp = 0, fp = 0, tn = 0, fn = 0;
		for (MemoryPeakResults peakResults : resultsList)
		{
			setup(peakResults);

			int frame = -1;
			int failCount = 0;
			for (PeakResult peak : peakResults.getResults())
			{
				final boolean isTrue = peak.origValue != 0;

				// Reset fail count for new frames
				if (frame != peak.peak)
				{
					frame = peak.peak;
					failCount = 0;
				}

				// Reject all peaks if we have exceeded the fail count
				final boolean isPositive;
				if (failCount > failures)
				{
					isPositive = false;
				}
				else
				{
					// Otherwise assess the peak
					isPositive = accept(peak);
				}

				if (isPositive)
				{
					failCount = 0;
				}
				else
				{
					failCount++;
				}

				if (isTrue)
				{
					if (isPositive)
						tp++; // true positive
					else
						fn++; // false negative
				}
				else
				{
					if (isPositive)
						fp++; // false positive
					else
						tn++; // true negative
				}
			}
			end();
		}
		return new ClassificationResult(tp, fp, tn, fn);
	}

	/**
	 * Filter the results and return the performance score. Allows benchmarking the filter by marking the results as
	 * true or false.
	 * <p>
	 * Any input PeakResult with an original value that is not zero will be treated as a true result, all other results
	 * are false. The filter is run and the results are marked as true positive, false negative and false positive.
	 * <p>
	 * The number of consecutive rejections are counted per frame. When the configured number of failures is reached all
	 * remaining results for the frame are rejected. This assumes the results are ordered by the frame.
	 * <p>
	 * Note that this method is to be used to score a subset that was generated using
	 * {@link #filterSubset(MemoryPeakResults, int)} since the number of consecutive failures before each peak are
	 * expected to be stored in the origX property.
	 * 
	 * @param resultsList
	 *            a list of results to analyse
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @param tn
	 *            The initial true negatives (used when the results have been pre-filtered)
	 * @param fn
	 *            The initial false negatives (used when the results have been pre-filtered)
	 * @return the score
	 */
	public ClassificationResult scoreSubset(List<MemoryPeakResults> resultsList, int failures, int tn, int fn)
	{
		int tp = 0, fp = 0;
		for (MemoryPeakResults peakResults : resultsList)
		{
			setup(peakResults);

			int frame = -1;
			int failCount = 0;
			for (PeakResult peak : peakResults.getResults())
			{
				final boolean isTrue = peak.origValue != 0;

				// Reset fail count for new frames
				if (frame != peak.peak)
				{
					frame = peak.peak;
					failCount = 0;
				}

				failCount += peak.origX;

				// Reject all peaks if we have exceeded the fail count
				final boolean isPositive;
				if (failCount > failures)
				{
					isPositive = false;
				}
				else
				{
					// Otherwise assess the peak
					isPositive = accept(peak);
				}

				if (isPositive)
				{
					failCount = 0;
				}
				else
				{
					failCount++;
				}

				if (isTrue)
				{
					if (isPositive)
						tp++; // true positive
					else
						fn++; // false negative
				}
				else
				{
					if (isPositive)
						fp++; // false positive
					else
						tn++; // true negative
				}
			}
			end();
		}
		return new ClassificationResult(tp, fp, tn, fn);
	}

	/**
	 * Filter the results and return the performance score. Allows benchmarking the filter by marking the results as
	 * true or false.
	 * <p>
	 * Input PeakResults must be allocated a score for true positive, false positive, true negative and false negative
	 * (accessed via the object property get methods). The filter is run and results that pass accumulate scores for
	 * true positive and false positive, otherwise the scores are accumulated for true negative and false negative. The
	 * simplest scoring scheme is to mark valid results as tp=fn=1 and fp=tn=0 and invalid results the opposite.
	 * <p>
	 * The number of consecutive rejections are counted per frame. When the configured number of failures is reached all
	 * remaining results for the frame are rejected. This assumes the results are ordered by the frame.
	 * 
	 * @param resultsList
	 *            a list of results to analyse
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @return the score
	 */
	public FractionClassificationResult fractionScore(List<MemoryPeakResults> resultsList, int failures)
	{
		double fp = 0, fn = 0;
		double tp = 0, tn = 0;
		for (MemoryPeakResults peakResults : resultsList)
		{
			setup(peakResults);

			int frame = -1;
			int failCount = 0;
			for (PeakResult peak : peakResults.getResults())
			{
				// Reset fail count for new frames
				if (frame != peak.peak)
				{
					frame = peak.peak;
					failCount = 0;
				}

				// Reject all peaks if we have exceeded the fail count
				final boolean isPositive;
				if (failCount > failures)
				{
					isPositive = false;
				}
				else
				{
					// Otherwise assess the peak
					isPositive = accept(peak);
				}

				if (isPositive)
				{
					failCount = 0;
				}
				else
				{
					failCount++;
				}

				if (isPositive)
				{
					tp += peak.getTruePositiveScore();
					fp += peak.getFalsePositiveScore();
				}
				else
				{
					fn += peak.getFalseNegativeScore();
					tn += peak.getTrueNegativeScore();
				}
			}
			end();
		}
		return new FractionClassificationResult(tp, fp, tn, fn);
	}

	/**
	 * Filter the results and return the performance score. Allows benchmarking the filter by marking the results as
	 * true or false.
	 * <p>
	 * Input PeakResults must be allocated a score for true positive, false positive, true negative and false negative
	 * (accessed via the object property get methods). The filter is run and results that pass accumulate scores for
	 * true positive and false positive, otherwise the scores are accumulated for true negative and false negative. The
	 * simplest scoring scheme is to mark valid results as tp=fn=1 and fp=tn=0 and invalid results the opposite.
	 * <p>
	 * The number of consecutive rejections are counted per frame. When the configured number of failures is reached all
	 * remaining results for the frame are rejected. This assumes the results are ordered by the frame.
	 * <p>
	 * Note that this method is to be used to score a subset that was generated using
	 * {@link #filterSubset(MemoryPeakResults, int)} since the number of consecutive failures before each peak are
	 * expected to be stored in the origX property.
	 * 
	 * @param resultsList
	 *            a list of results to analyse
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @param tn
	 *            The initial true negatives (used when the results have been pre-filtered)
	 * @param fn
	 *            The initial false negatives (used when the results have been pre-filtered)
	 * @return the score
	 */
	public FractionClassificationResult fractionScoreSubset(List<MemoryPeakResults> resultsList, int failures,
			double tn, double fn)
	{
		double fp = 0;
		double tp = 0;

		for (MemoryPeakResults peakResults : resultsList)
		{
			setup(peakResults);

			int frame = -1;
			int failCount = 0;
			for (PeakResult peak : peakResults.getResults())
			{
				// Reset fail count for new frames
				if (frame != peak.peak)
				{
					frame = peak.peak;
					failCount = 0;
				}

				failCount += peak.origX;

				// Reject all peaks if we have exceeded the fail count
				final boolean isPositive;
				if (failCount > failures)
				{
					isPositive = false;
				}
				else
				{
					// Otherwise assess the peak
					isPositive = accept(peak);
				}

				if (isPositive)
				{
					failCount = 0;
				}
				else
				{
					failCount++;
				}

				if (isPositive)
				{
					tp += peak.getTruePositiveScore();
					fp += peak.getFalsePositiveScore();
				}
				else
				{
					fn += peak.getFalseNegativeScore();
					tn += peak.getTrueNegativeScore();
				}
			}
			end();
		}
		return new FractionClassificationResult(tp, fp, tn, fn);
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
	 * Called after the accept method has been called for each peak in the results. Allows memory clean-up of the
	 * results.
	 */
	public void end()
	{
	}

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

		// Use all the parameters
		if (getNumberOfParameters() == o.getNumberOfParameters())
		{
			//for (int i=getNumberOfParameters(); i-- > 0; )
			for (int i = 0; i < getNumberOfParameters(); i++)
			{
				final double d1 = getParameterValue(i);
				final double d2 = o.getParameterValue(i);
				if (d1 < d2)
					return -1;
				if (d1 > d2)
					return 1;
			}
		}

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
	 * <p>
	 * Filters can adjust the parameter by a different amount, e.g. by the delta multiplied by a range expected to
	 * change the filter performance. This may be relevant in the case where the value is presently zero since no
	 * relative change is possible.
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
	 * @param defaultRange
	 *            The default range to apply the delta to in the case where the value is zero and no relative adjustment
	 *            is possible.
	 * @return
	 */
	protected double updateParameter(double value, double delta, double defaultRange)
	{
		if (value != 0)
			return (float) (value + value * delta);
		return (float) (value + defaultRange * delta);
	}

	/**
	 * Adjust the specified parameter value.
	 * <p>
	 * A positive delta will adjust the parameter to be larger. A negative delta will adjust the parameter to be
	 * smaller. The adjustment is relative to the parameter value, e.g. 0.1 is 10%.
	 * 
	 * @param value
	 * @param delta
	 * @param defaultRange
	 *            The default range to apply the delta to in the case where the value is zero and no relative adjustment
	 *            is possible.
	 * @return
	 */
	protected float updateParameter(float value, double delta, double defaultRange)
	{
		if (value != 0)
			return (float) (value + value * delta);
		return (float) (value + defaultRange * delta);
	}

	/**
	 * Adjust the specified parameter value.
	 * <p>
	 * A positive delta will adjust the parameter to be larger. A negative delta will adjust the parameter to be
	 * smaller. The adjustment is relative to the parameter value, e.g. 0.1 is 10%. The adjustment is rounded up to the
	 * next valid integer to ensure a new parameter value is created.
	 * 
	 * @param value
	 * @param delta
	 * @param defaultRange
	 *            The default range to apply the delta to in the case where the value is zero and no relative adjustment
	 *            is possible.
	 * @return
	 */
	protected int updateParameter(int value, double delta, int defaultRange)
	{
		final int update;
		if (value != 0)
			update = (int) Math.ceil(value * Math.abs(delta));
		else
			update = (int) Math.ceil(defaultRange * Math.abs(delta));
		if (delta < 0)
			return value - update;
		return value + update;
	}

	protected void setMin(double[] parameters, int index, double value)
	{
		if (parameters[index] > value)
			parameters[index] = value;
	}

	protected void setMax(double[] parameters, int index, double value)
	{
		if (parameters[index] < value)
			parameters[index] = value;
	}

	/**
	 * Create a new filter with the specified parameters
	 * 
	 * @param parameters
	 * @return A new filter
	 */
	public abstract Filter create(double... parameters);

	/**
	 * Update the input array if the Filter's parameters are weaker. This method can be used to find the weakest
	 * parameters across a set of filters of the same type. The weakest filter can then be used to create a subset of
	 * pre-filtered results to use for testing the filter set.
	 * 
	 * @param parameters
	 *            The parameters
	 */
	public abstract void weakestParameters(double[] parameters);

	/**
	 * Some filters requires all the data in a subset for scoring analysis. Others can create a subset using the fail
	 * count parameter for a smaller subset that will evaluate faster. This method returns true if the subset can be
	 * created using the fail count parameter that will be used to score the subset.
	 * 
	 * @return True if the {@link #filterSubset(MemoryPeakResults, int, double[])} is valid
	 */
	public boolean subsetWithFailCount()
	{
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#newChromosome(double[])
	 */
	public Chromosome newChromosome(double[] sequence)
	{
		return create(sequence);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#lowerLimit()
	 */
	public double[] lowerLimit()
	{
		// Set zero as the lower limit
		return new double[length()];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#upperLimit()
	 */
	public double[] upperLimit()
	{
		// No need for upper limits on filters
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#setFitness(double)
	 */
	public void setFitness(double fitness)
	{
		this.fitness = fitness;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#getFitness()
	 */
	public double getFitness()
	{
		return fitness;
	}

	/**
	 * Return the Manhattan (city-block) distance between two chromosomes. This measure is intended to return if the
	 * sequences are the same (zero distance) or not). It is not intended for use in distance analysis.
	 * 
	 * @see gdsc.smlm.ga.Chromosome#distance(gdsc.smlm.ga.Chromosome)
	 */
	public double distance(Chromosome other)
	{
		// NOTE: If the distance is required for a certain type of analysis then this could be done
		// using injection of an interface for calculating the distance.

		final int n = FastMath.min(length(), other.length());
		double[] s1 = sequence();
		double[] s2 = other.sequence();
		double d = 0;
		for (int i = 0; i < n; i++)
			d += Math.abs(s1[i] - s2[i]);
		return d;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#equals(gdsc.smlm.ga.Chromosome)
	 */
	public boolean equals(Chromosome other)
	{
		if (length() != other.length())
			return false;
		final int n = length();
		double[] s1 = sequence();
		double[] s2 = other.sequence();
		for (int i = 0; i < n; i++)
			if (s1[i] != s2[i])
				return false;
		return true;
	}

	/**
	 * Get the indices of the parameters that are included in the Chromosome interface. This can be used to look up the
	 * name of the parameter using {@link #getParameterName(int)}.
	 * 
	 * @return The indices of the parameters that are included in the Chromosome interface
	 */
	public int[] getChromosomeParameters()
	{
		// Assume all the parameters are included in the Chromosome
		return Utils.newArray(getNumberOfParameters(), 0, 1);
	}
}