package gdsc.smlm.results.match;

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
 * Class to store the result of a binary scoring analysis.
 * 
 * Can calculate the F-score statistic with a given beta weighting between the precision and recall.
 * 
 * @see http://en.wikipedia.org/wiki/Precision_and_recall#F-measure
 */
public class MatchResult
{
	private final int tp;
	private final int fp;
	private final int fn;
	private final double precision;
	private final double recall;
	private final double jaccard;
	private final double rmsd;

	/**
	 * @param tp
	 *            The number of true positives
	 * @param fp
	 *            The number of false positives
	 * @param fn
	 *            The number of false negatives
	 * @param rmsd
	 *            The root mean squared distance between true positives
	 */
	public MatchResult(int tp, int fp, int fn, double rmsd)
	{
		this.tp = tp;
		this.fp = fp;
		this.fn = fn;
		this.rmsd = rmsd;

		precision = divide(tp, tp + fp);
		recall = divide(tp, tp + fn);
		jaccard = divide(tp, tp + fp + fn);
	}

	private static double divide(final double numerator, final int denominator)
	{
		if (denominator == 0)
			return 0;
		return numerator / denominator;
	}

	/**
	 * Return the F-Score statistic, a weighted combination of the precision and recall
	 * 
	 * @param precision
	 * @param recall
	 * @param beta
	 *            The weight
	 * @return The F-Score
	 */
	public static double calculateFScore(double precision, double recall, double beta)
	{
		final double b2 = beta * beta;
		final double f = ((1.0 + b2) * precision * recall) / (b2 * precision + recall);
		return (Double.isNaN(f) ? 0 : f);
	}

	/**
	 * Return the F-Score statistic, a weighted combination of the precision and recall
	 * 
	 * @param beta
	 *            The weight
	 * @return The F-Score
	 */
	public double getFScore(double beta)
	{
		return calculateFScore(precision, recall, beta);
	}

	/**
	 * Return the F1-Score statistic, a equal weighted combination of the precision and recall
	 * 
	 * @return The F1-Score
	 */
	public double getF1Score()
	{
		final double f = (2 * precision * recall) / (precision + recall);
		return (Double.isNaN(f) ? 0 : f);
	}

	/**
	 * @return the n
	 */
	public int getNumberPredicted()
	{
		return tp + fp;
	}

	/**
	 * @return the number of actual
	 */
	public int getNumberActual()
	{
		return tp + fn;
	}

	/**
	 * @return the tp
	 */
	public int getTruePositives()
	{
		return tp;
	}

	/**
	 * @return the fp
	 */
	public int getFalsePositives()
	{
		return fp;
	}

	/**
	 * @return the fn
	 */
	public int getFalseNegatives()
	{
		return fn;
	}

	/**
	 * @return the precision
	 */
	public double getPrecision()
	{
		return precision;
	}

	/**
	 * @return the recall
	 */
	public double getRecall()
	{
		return recall;
	}

	/**
	 * @return the Jaccard index (defined as the size of the intersection divided by the size of the union of the sample
	 *         sets)
	 */
	public double getJaccard()
	{
		return jaccard;
	}

	/**
	 * @return the root mean squared distance between true positives
	 */
	public double getRMSD()
	{
		return rmsd;
	}
}
