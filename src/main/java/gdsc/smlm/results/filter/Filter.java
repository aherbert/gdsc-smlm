package gdsc.smlm.results.filter;

import gdsc.smlm.ga.Chromosome;
import gdsc.core.ij.Utils;
import gdsc.core.match.ClassificationResult;
import gdsc.core.match.FractionClassificationResult;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.MemoryPeakResults;

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
public abstract class Filter implements Comparable<Filter>, Chromosome<FilterScore>, Cloneable
{
	@XStreamOmitField
	private String name;
	@XStreamOmitField
	private String type;
	@XStreamOmitField
	private FilterScore fitness;

	/**
	 * Generate the name of the filter using the filter settings (defaults to the first parameter)
	 * 
	 * @return The name of the filter
	 */
	protected String generateName()
	{
		return getParameterName(0) + " " + getParameterValue(0);
	}

	/**
	 * Generate the type of the filter using the filter settings (default to the class name with 'Filter' removed)
	 * 
	 * @return The type of the filter
	 */
	protected String generateType()
	{
		return this.getClass().getSimpleName().replaceAll("Filter", "");
	}

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
	 * The number of consecutive rejections are counted per frame. When the configured number of failures is reached all
	 * remaining results for the frame are rejected. This assumes the results are ordered by the frame.
	 * <p>
	 * Note that this method is to be used to score a set of results that may have been extracted from a larger set
	 * since the number of consecutive failures before each peak are expected to be stored in the origY property. Set
	 * this to zero and the results should be identical to {@link #filter(MemoryPeakResults, int)}
	 * 
	 * @param results
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @return the filtered results
	 */
	public MemoryPeakResults filter2(MemoryPeakResults results, int failures)
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

			failCount += peak.origY;

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
	 *            If not null will be populated with the fraction score [ tp, fp, tn, fn, p, n ]
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
		int p = 0;
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
				p++;
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

		if (score != null && score.length > 5)
		{
			score[0] = tp;
			score[1] = fp;
			score[2] = tn;
			score[3] = fn;
			score[4] = p;
			score[5] = results.size() - p;
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
	 * Note that this method is to be used to score a set of results that may have been extracted from a larger set
	 * since the number of consecutive failures before each peak are expected to be stored in the origY property. Set
	 * this to zero and the results should be identical to {@link #filterSubset(MemoryPeakResults, double[])}.
	 * <p>
	 * The number of failures before each peak is stored in the origX property of the PeakResult.
	 * 
	 * @param results
	 * @param score
	 *            If not null will be populated with the fraction score [ tp, fp, tn, fn, p, n ]
	 * @return the filtered results
	 */
	public MemoryPeakResults filterSubset2(MemoryPeakResults results, double[] score)
	{
		MemoryPeakResults newResults = new MemoryPeakResults();
		newResults.copySettings(results);
		setup(results);
		int frame = -1;
		int failCount = 0;
		double fp = 0, fn = 0;
		double tp = 0, tn = 0;
		int p = 0;
		for (PeakResult peak : results.getResults())
		{
			if (frame != peak.peak)
			{
				frame = peak.peak;
				failCount = 0;
			}

			failCount += peak.origY;

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
				p++;
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

		if (score != null && score.length > 5)
		{
			score[0] = tp;
			score[1] = fp;
			score[2] = tn;
			score[3] = fn;
			score[4] = p;
			score[5] = results.size() - p;
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
	 *            If not null will be populated with the fraction score [ tp, fp, tn, fn, p, n ]
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

		if (score != null && score.length > 5)
		{
			score[0] = tp;
			score[1] = fp;
			score[2] = tn;
			score[3] = fn;
			score[4] = newResults.size();
			score[5] = results.size() - newResults.size();
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
	 * Note that this method is to be used to score a set of results that may have been extracted from a larger set
	 * since the number of consecutive failures before each peak are expected to be stored in the origY property. Set
	 * this to zero and the results should be identical to {@link #filterSubset(MemoryPeakResults, int, double[])}.
	 * <p>
	 * The number of failures before each peak is stored in the origX property of the PeakResult.
	 * 
	 * @param results
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @param score
	 *            If not null will be populated with the fraction score [ tp, fp, tn, fn, p, n ]
	 * @return the filtered results
	 */
	public MemoryPeakResults filterSubset2(MemoryPeakResults results, int failures, double[] score)
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

			failCount += peak.origY;

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

		if (score != null && score.length > 5)
		{
			score[0] = tp;
			score[1] = fp;
			score[2] = tn;
			score[3] = fn;
			score[4] = newResults.size();
			score[5] = results.size() - newResults.size();
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
		int p = 0, n = 0;
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
					p++;
					tp += peak.getTruePositiveScore();
					fp += peak.getFalsePositiveScore();
				}
				else
				{
					fn += peak.getFalseNegativeScore();
					tn += peak.getTrueNegativeScore();
				}
			}
			n += peakResults.size();
			end();
		}
		n -= p;
		return new FractionClassificationResult(tp, fp, tn, fn, p, n);
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
	 * Note that this method is to be used to score a set of results that may have been extracted from a larger set
	 * since the number of consecutive failures before each peak are expected to be stored in the origY property. Set
	 * this to zero and the results should be identical to {@link #fractionScore(List, int)}.
	 * 
	 * @param resultsList
	 *            a list of results to analyse
	 * @param failures
	 *            the number of failures to allow per frame before all peaks are rejected
	 * @return the score
	 */
	public FractionClassificationResult fractionScore2(List<MemoryPeakResults> resultsList, int failures)
	{
		int p = 0, n = 0;
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

				failCount += peak.origY;

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
					p++;
					tp += peak.getTruePositiveScore();
					fp += peak.getFalsePositiveScore();
				}
				else
				{
					fn += peak.getFalseNegativeScore();
					tn += peak.getTrueNegativeScore();
				}
			}
			n += peakResults.size();
			end();
		}
		n -= p;
		return new FractionClassificationResult(tp, fp, tn, fn, p, n);
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
	 * @param n
	 *            The initial negatives (used when the results have been pre-filtered)
	 * @return the score
	 */
	public FractionClassificationResult fractionScoreSubset(List<MemoryPeakResults> resultsList, int failures,
			double tn, double fn, int n)
	{
		int p = 0;
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
					p++;
					tp += peak.getTruePositiveScore();
					fp += peak.getFalsePositiveScore();
				}
				else
				{
					fn += peak.getFalseNegativeScore();
					tn += peak.getTrueNegativeScore();
				}
			}
			n += peakResults.size();
			end();
		}
		n -= p;
		return new FractionClassificationResult(tp, fp, tn, fn, p, n);
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
	 * The numerical value of the filter (defaults to the first parameter)
	 * 
	 * @return The numerical value of the filter. Used for plotting value against performance score.
	 */
	public double getNumericalValue()
	{
		return getParameterValue(0);
	}

	/**
	 * The name of the numerical value of the filter (defaults to the first parameter)
	 * 
	 * @return The name of the numerical value of the filter. Used for plotting value against performance score.
	 */
	public String getNumericalValueName()
	{
		return getParameterName(0);
	}

	/**
	 * @return the name (including any parameter values)
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
			Filter f = (Filter) XStreamWrapper.fromXML(xml);
			f.initialiseState();
			return f;
		}
		catch (ClassCastException ex)
		{
			//ex.printStackTrace();
		}
		return null;
	}

	/**
	 * Run after the filter is deserialised using XStream or cloned. Overrride this method if the filter has state that
	 * requires resetting.
	 */
	protected void initialiseState()
	{

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
		final int size = getNumberOfParameters();
		if (size == o.getNumberOfParameters())
		{
			for (int i = 0; i < size; i++)
			{
				final double d1 = getParameterValueInternal(i);
				final double d2 = o.getParameterValueInternal(i);
				if (d1 < d2)
					return -1;
				if (d1 > d2)
					return 1;
			}
		}

		return 0;
	}

	/**
	 * Compare to the other filter, count the number of weakest parameters. If negative then this filter has more weak
	 * parameters. If positive then this filter has less weak parameters. If the same or the number of parameters do not
	 * match then return 0. If the other filter is null return -1.
	 * 
	 * @param o
	 *            The other filter
	 * @return the count difference
	 */
	public int weakest(Filter o)
	{
		// Null to end of list
		if (o == null)
			return -1;

		// Use all the parameters
		int i = getNumberOfParameters();
		if (i == o.getNumberOfParameters())
		{
			// Extract the parameters
			final double[] p1 = getParameters();
			final double[] p2 = o.getParameters();

			// Find the weakest

			final double[] weakest = p1.clone();
			o.weakestParameters(weakest);
			// Count the number of weakest
			int c = 0;
			while (i-- > 0)
			{
				if (p1[i] != p2[i])
				{
					if (p1[i] == weakest[i])
						--c;
					else
						++c;
				}
			}
			return c;
		}

		return 0;
	}

	/**
	 * Compare to the other filter, count the number of weakest parameters. If negative then this filter has more weak
	 * parameters. If positive then this filter has less weak parameters.
	 * <p>
	 * This method does not check for null or if the other filter has a different number of parameters.
	 * 
	 * @param o
	 *            The other filter
	 * @return the count difference
	 */
	public int weakestUnsafe(Filter o)
	{
		// Use all the parameters
		int i = getNumberOfParameters();

		// Extract the parameters
		final double[] p1 = getParameters();
		final double[] p2 = o.getParameters();

		// Find the weakest

		final double[] weakest = p1.clone();
		o.weakestParameters(weakest);
		// Count the number of weakest
		int c = 0;
		while (i-- > 0)
		{
			if (p1[i] != p2[i])
			{
				if (p1[i] == weakest[i])
					--c;
				else
					++c;
			}
		}
		return c;
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
	 * Get the parameter value.
	 * 
	 * @param index
	 * @return The value of the specified parameter
	 */
	public double getParameterValue(int index)
	{
		checkIndex(index);
		return getParameterValueInternal(index);
	}

	/**
	 * Get the parameter value. The index should always be between 0 and {@link #getNumberOfParameters()}
	 * 
	 * @param index
	 * @return The value of the specified parameter
	 */
	protected abstract double getParameterValueInternal(int index);

	/**
	 * Gets the parameters as an array.
	 *
	 * @return the parameters
	 */
	public double[] getParameters()
	{
		final int n = getNumberOfParameters();
		final double[] p = new double[n];
		for (int i = 0; i < n; i++)
			p[i] = getParameterValueInternal(i);
		return p;
	}

	/**
	 * Get the recommended minimum amount by which to increment the parameter
	 * 
	 * @param index
	 * @return The increment value of the specified parameter
	 */
	public abstract double getParameterIncrement(int index);

	/**
	 * Return a value to use to disable the parameter
	 * <p>
	 * Override this method if zero does not disable the parameter
	 * 
	 * @param index
	 * @return The disabled value of the specified parameter
	 */
	public double getDisabledParameterValue(int index)
	{
		checkIndex(index);
		return 0;
	}

	/**
	 * @param index
	 * @return The name of the specified parameter
	 */
	public String getParameterName(int index)
	{
		return getParameterType(index).toString();
	}

	/**
	 * @param index
	 * @return The type of the specified parameter
	 */
	public abstract ParameterType getParameterType(int index);

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
			return (value + value * delta);
		return (value + defaultRange * delta);
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
	 * Creates a new filter with only the specified parameters enabled.
	 *
	 * @param enable
	 *            the enabled flags
	 * @return the filter
	 */
	public Filter create(boolean[] enable)
	{
		if (enable == null || enable.length != getNumberOfParameters())
			throw new IllegalArgumentException(
					"Enable array must match the number of parameters: " + getNumberOfParameters());
		final double[] p = new double[enable.length];
		for (int i = 0; i < p.length; i++)
			p[i] = (enable[i]) ? getParameterValueInternal(i) : getDisabledParameterValue(i);
		return create(p);
	}

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
	 * Compare the two values and return a sort result for the minimum of the two
	 *
	 * @param value1
	 *            the value 1
	 * @param value2
	 *            the value 2
	 * @return the result (-1 is value1 is lower, 0 is equal, 1 is value2 is lower)
	 */
	public static int compareMin(double value1, double value2)
	{
		if (value1 < value2)
			return -1;
		if (value1 > value2)
			return 1;
		return 0;
	}

	/**
	 * Compare the two values and return a sort result for the maximum of the two
	 *
	 * @param value1
	 *            the value 1
	 * @param value2
	 *            the value 2
	 * @return the result (-1 is value1 is higher, 0 is equal, 1 is value2 is higher)
	 */
	public static int compareMax(double value1, double value2)
	{
		if (value1 < value2)
			return 1;
		if (value1 > value2)
			return -1;
		return 0;
	}

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
	 * @see gdsc.smlm.ga.Chromosome#length()
	 */
	public int length()
	{
		// Assume all the parameters are included in the Chromosome
		return getNumberOfParameters();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#sequence()
	 */
	public double[] sequence()
	{
		// Assume all the parameters are included in the Chromosome
		return getParameters();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#newChromosome(double[])
	 */
	public Chromosome<FilterScore> newChromosome(double[] sequence)
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
	public void setFitness(FilterScore fitness)
	{
		this.fitness = fitness;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#getFitness()
	 */
	public FilterScore getFitness()
	{
		return fitness;
	}

	/**
	 * Return the Manhattan (city-block) distance between two chromosomes. This measure is intended to return if the
	 * sequences are the same (zero distance) or not). It is not intended for use in distance analysis.
	 * 
	 * @see gdsc.smlm.ga.Chromosome#distance(gdsc.smlm.ga.Chromosome)
	 */
	public double distance(Chromosome<FilterScore> other)
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
	public boolean equals(Chromosome<FilterScore> other)
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

	/**
	 * Return the value or Float.POSITIVE_INFINITY if value is not positive
	 * 
	 * @param value
	 * @return The limit
	 */
	public static float getUpperLimit(double value)
	{
		if (value > 0)
			return (float) value;
		else
			return Float.POSITIVE_INFINITY;
	}

	/**
	 * Return the value squared or Float.POSITIVE_INFINITY if value is not positive
	 * 
	 * @param value
	 * @return The squared limit
	 */
	public static float getUpperSquaredLimit(double value)
	{
		if (value > 0)
			return (float) (value * value);
		else
			return Float.POSITIVE_INFINITY;
	}

	/**
	 * Return the value or Double.POSITIVE_INFINITY if value is not positive
	 * 
	 * @param value
	 * @return The limit
	 */
	public static double getDUpperLimit(double value)
	{
		if (value > 0)
			return value;
		else
			return Double.POSITIVE_INFINITY;
	}

	/**
	 * Return the value squared or Double.POSITIVE_INFINITY if value is not positive
	 * 
	 * @param value
	 * @return The squared limit
	 */
	public static double getDUpperSquaredLimit(double value)
	{
		if (value > 0)
			return value * value;
		else
			return Double.POSITIVE_INFINITY;
	}

	/**
	 * @return The filter type
	 */
	public FilterType getFilterType()
	{
		return FilterType.STANDARD;
	}

	@Override
	public Filter clone()
	{
		try
		{
			Filter f = (Filter) super.clone();
			f.initialiseState();
			return f;
		}
		catch (CloneNotSupportedException e)
		{
			return null;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj)
	{
		if (obj == null)
		{
			return false;
		}
		if (obj == this)
		{
			return true;
		}
		if (!(obj instanceof Filter))
		{
			return false;
		}
		final Filter other = (Filter) obj;
		final int size = getNumberOfParameters();
		if (size != other.getNumberOfParameters())
		{
			return false;
		}
		// Check the types are the same before a parameter comparison
		if (!this.getType().equals(other.getType()))
		{
			return false;
		}
		for (int i = 0; i < size; i++)
		{
			final double d1 = getParameterValueInternal(i);
			final double d2 = other.getParameterValueInternal(i);
			if (d1 != d2)
				return false;
		}
		return true;
	}
}