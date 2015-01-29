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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

/**
 * Calculates the match between a set of predicted points and the actual points.
 */
public class MatchCalculator
{
	/**
	 * Calculate the match results for the given actual and predicted points.
	 * Points that are within the distance threshold are identified as a match.
	 * The number of true positives, false positives and false negatives are calculated.
	 * 
	 * @param actualPoints
	 * @param predictedPoints
	 * @param dThreshold
	 *            The distance threshold
	 * @return The match results
	 */
	public static MatchResult analyseResults2D(Coordinate[] actualPoints, Coordinate[] predictedPoints,
			double dThreshold)
	{
		return analyseResults2D(actualPoints, predictedPoints, dThreshold, null, null, null, null);
	}

	/**
	 * Calculate the match results for the given actual and predicted points.
	 * Points that are within the distance threshold are identified as a match.
	 * The number of true positives, false positives and false negatives are calculated.
	 * 
	 * @param actualPoints
	 * @param predictedPoints
	 * @param dThreshold
	 *            The distance threshold
	 * @param TP
	 *            True Positives
	 * @param FP
	 *            False Positives
	 * @param FN
	 *            False Negatives
	 * @return The match results
	 */
	public static MatchResult analyseResults2D(Coordinate[] actualPoints, Coordinate[] predictedPoints,
			double dThreshold, List<Coordinate> TP, List<Coordinate> FP, List<Coordinate> FN)
	{
		return analyseResults2D(actualPoints, predictedPoints, dThreshold, TP, FP, FN, null);
	}

	/**
	 * Calculate the match results for the given actual and predicted points.
	 * Points that are within the distance threshold are identified as a match.
	 * The number of true positives, false positives and false negatives are calculated.
	 * 
	 * @param actualPoints
	 * @param predictedPoints
	 * @param dThreshold
	 *            The distance threshold
	 * @param TP
	 *            True Positives
	 * @param FP
	 *            False Positives
	 * @param FN
	 *            False Negatives
	 * @param matches
	 *            The matched true positives (point1 = actual, point2 = predicted)
	 * @return The match results
	 */
	public static MatchResult analyseResults2D(Coordinate[] actualPoints, Coordinate[] predictedPoints,
			double dThreshold, List<Coordinate> TP, List<Coordinate> FP, List<Coordinate> FN, List<PointPair> matches)
	{
		dThreshold *= dThreshold; // We will use the squared distance

		final int predictedPointsLength = (predictedPoints != null) ? predictedPoints.length : 0;
		final int actualPointsLength = (actualPoints != null) ? actualPoints.length : 0;

		int tp = 0; // true positives (actual with matched predicted point)
		int fp = predictedPointsLength; // false positives (actual with no matched predicted point)
		int fn = actualPointsLength; // false negatives (predicted point with no actual point)
		double rmsd = 0;

		if (predictedPointsLength == 0 || actualPointsLength == 0)
		{
			if (FP != null)
				FP.addAll(asList(predictedPoints));
			if (FN != null)
				FN.addAll(asList(actualPoints));
			return new MatchResult(tp, fp, fn, rmsd);
		}

		// loop over the two arrays assigning the closest unassigned pair
		boolean[] resultAssignment = new boolean[predictedPointsLength];
		boolean[] roiAssignment = new boolean[fn];
		ArrayList<Assignment> assignments = new ArrayList<Assignment>(predictedPointsLength);

		int[] falsePositives = null, falseNegatives = null;
		if (FP != null)
		{
			falsePositives = ascendingArray(predictedPointsLength);
		}
		if (FN != null)
		{
			falseNegatives = ascendingArray(actualPointsLength);
		}

		// Pre-calculate all-vs-all distance matrix if it can fit in memory
		int size = predictedPointsLength * actualPointsLength;
		float[][] dMatrix = null;
		if (size < 200 * 200)
		{
			dMatrix = new float[predictedPointsLength][actualPointsLength];
			for (int predictedId = predictedPointsLength; predictedId-- > 0;)
			{
				final float x = predictedPoints[predictedId].getX();
				final float y = predictedPoints[predictedId].getY();
				for (int actualId = actualPointsLength; actualId-- > 0;)
				{
					dMatrix[predictedId][actualId] = (float) actualPoints[actualId].distance2(x, y);
				}
			}
		}

		do
		{
			assignments.clear();

			// Process each result
			for (int predictedId = predictedPointsLength; predictedId-- > 0;)
			{
				if (resultAssignment[predictedId])
					continue; // Already assigned

				final float x = predictedPoints[predictedId].getX();
				final float y = predictedPoints[predictedId].getY();

				// Find closest ROI point
				float d2Min = (float) dThreshold; //Float.MAX_VALUE;
				int targetId = -1;
				for (int actualId = actualPointsLength; actualId-- > 0;)
				{
					if (roiAssignment[actualId])
						continue; // Already assigned

					if (dMatrix != null)
					{
						if (dMatrix[predictedId][actualId] <= d2Min)
						{
							d2Min = dMatrix[predictedId][actualId];
							targetId = actualId;
						}
					}
					else
					{
						Coordinate actualPoint = actualPoints[actualId];

						// Calculate in steps for increased speed (allows early exit)
						float dx = actualPoint.getX() - x;
						dx *= dx;
						if (dx <= d2Min)
						{
							float dy = actualPoint.getY() - y;
							dy *= dy;
							if (dy <= d2Min)
							{
								final float d2 = dx + dy;
								if (d2 <= d2Min)
								{
									d2Min = d2;
									targetId = actualId;
								}
							}
						}
					}
				}

				// Store closest ROI point
				if (targetId > -1)
				{
					assignments.add(new Assignment(targetId, predictedId, d2Min));
				}
			}

			// If there are assignments
			if (!assignments.isEmpty())
			{
				// Pick the closest pair to be assigned
				Collections.sort(assignments);

				// Process in order
				for (Assignment closest : assignments)
				{
					// Skip those that have already been assigned since this will be a lower score.
					// Note at least one assignment should be processed as potential assignments are made 
					// using only unassigned points.
					if (resultAssignment[closest.getPredictedId()] || roiAssignment[closest.getTargetId()])
						continue;
					resultAssignment[closest.getPredictedId()] = true;
					roiAssignment[closest.getTargetId()] = true;

					// If within accuracy then classify as a match
					if (closest.getDistance() < dThreshold)
					{
						tp++;
						fn--;
						fp--;
						rmsd += closest.getDistance(); // Already a squared distance

						if (TP != null)
						{
							TP.add(predictedPoints[closest.getPredictedId()]);
						}
						if (FP != null)
						{
							falsePositives[closest.getPredictedId()] = -1;
						}
						if (FN != null)
						{
							falseNegatives[closest.getTargetId()] = -1;
						}
						if (matches != null)
							matches.add(new PointPair(actualPoints[closest.getTargetId()], predictedPoints[closest
									.getPredictedId()]));
					}
					else
					{
						// No more assignments within the distance threshold
						break;
					}
				}
			}

		} while (!assignments.isEmpty());

		// Add to lists
		if (FP != null)
		{
			for (int i = 0; i < predictedPointsLength; i++)
			{
				if (falsePositives[i] >= 0)
					FP.add(predictedPoints[i]);
			}
		}
		if (FN != null)
		{
			for (int i = 0; i < actualPointsLength; i++)
			{
				if (falseNegatives[i] >= 0)
					FN.add(actualPoints[i]);
			}
		}

		if (tp > 0)
			rmsd = Math.sqrt(rmsd / tp);
		return new MatchResult(tp, fp, fn, rmsd);
	}

	private static Collection<Coordinate> asList(Coordinate[] points)
	{
		if (points != null)
			return Arrays.asList(points);
		return new ArrayList<Coordinate>(0);
	}

	/**
	 * Calculate the match results for the given actual and predicted points.
	 * Points that are within the distance threshold are identified as a match.
	 * The number of true positives, false positives and false negatives are calculated.
	 * 
	 * @param actualPoints
	 * @param predictedPoints
	 * @param dThreshold
	 *            The distance threshold
	 * @return The match results
	 */
	public static MatchResult analyseResults3D(Coordinate[] actualPoints, Coordinate[] predictedPoints,
			double dThreshold)
	{
		return analyseResults3D(actualPoints, predictedPoints, dThreshold, null, null, null, null);
	}

	/**
	 * Calculate the match results for the given actual and predicted points.
	 * Points that are within the distance threshold are identified as a match.
	 * The number of true positives, false positives and false negatives are calculated.
	 * 
	 * @param actualPoints
	 * @param predictedPoints
	 * @param dThreshold
	 *            The distance threshold
	 * @param TP
	 *            True Positives
	 * @param FP
	 *            False Positives
	 * @param FN
	 *            False Negatives
	 * @return The match results
	 */
	public static MatchResult analyseResults3D(Coordinate[] actualPoints, Coordinate[] predictedPoints,
			double dThreshold, List<Coordinate> TP, List<Coordinate> FP, List<Coordinate> FN)
	{
		return analyseResults3D(actualPoints, predictedPoints, dThreshold, TP, FP, FN, null);
	}

	/**
	 * Calculate the match results for the given actual and predicted points.
	 * Points that are within the distance threshold are identified as a match.
	 * The number of true positives, false positives and false negatives are calculated.
	 * 
	 * @param actualPoints
	 * @param predictedPoints
	 * @param dThreshold
	 *            The distance threshold
	 * @param TP
	 *            True Positives
	 * @param FP
	 *            False Positives
	 * @param FN
	 *            False Negatives
	 * @param matches
	 *            The matched true positives (point1 = actual, point2 = predicted)
	 * @return The match results
	 */
	public static MatchResult analyseResults3D(Coordinate[] actualPoints, Coordinate[] predictedPoints,
			double dThreshold, List<Coordinate> TP, List<Coordinate> FP, List<Coordinate> FN, List<PointPair> matches)
	{
		dThreshold *= dThreshold; // We will use the squared distance

		final int predictedPointsLength = (predictedPoints != null) ? predictedPoints.length : 0;
		final int actualPointsLength = (actualPoints != null) ? actualPoints.length : 0;

		int tp = 0; // true positives (actual with matched predicted point)
		int fp = predictedPointsLength; // false positives (actual with no matched predicted point)
		int fn = actualPointsLength; // false negatives (predicted point with no actual point)
		double rmsd = 0;

		if (predictedPointsLength == 0 || actualPointsLength == 0)
		{
			if (FP != null)
				FP.addAll(asList(predictedPoints));
			if (FN != null)
				FN.addAll(asList(actualPoints));
			return new MatchResult(tp, fp, fn, rmsd);
		}

		// loop over the two arrays assigning the closest unassigned pair
		boolean[] resultAssignment = new boolean[predictedPointsLength];
		boolean[] roiAssignment = new boolean[fn];
		ArrayList<Assignment> assignments = new ArrayList<Assignment>(predictedPointsLength);

		int[] falsePositives = null, falseNegatives = null;
		if (FP != null)
		{
			falsePositives = ascendingArray(predictedPointsLength);
		}
		if (FN != null)
		{
			falseNegatives = ascendingArray(actualPointsLength);
		}

		// Pre-calculate all-vs-all distance matrix if it can fit in memory
		int size = predictedPointsLength * actualPointsLength;
		float[][] dMatrix = null;
		if (size < 200 * 200)
		{
			dMatrix = new float[predictedPointsLength][actualPointsLength];
			for (int predictedId = predictedPointsLength; predictedId-- > 0;)
			{
				float x = predictedPoints[predictedId].getX();
				float y = predictedPoints[predictedId].getY();
				float z = predictedPoints[predictedId].getZ();
				for (int actualId = actualPointsLength; actualId-- > 0;)
				{
					double d2 = actualPoints[actualId].distance2(x, y, z);
					dMatrix[predictedId][actualId] = (float) d2;
				}
			}
		}

		do
		{
			assignments.clear();

			// Process each result
			for (int predictedId = predictedPointsLength; predictedId-- > 0;)
			{
				if (resultAssignment[predictedId])
					continue; // Already assigned

				final float x = predictedPoints[predictedId].getX();
				final float y = predictedPoints[predictedId].getY();
				final float z = predictedPoints[predictedId].getZ();

				// Find closest ROI point
				float d2Min = (float) dThreshold; //Float.MAX_VALUE;
				int targetId = -1;
				for (int actualId = actualPointsLength; actualId-- > 0;)
				{
					if (roiAssignment[actualId])
						continue; // Already assigned

					if (dMatrix != null)
					{
						if (dMatrix[predictedId][actualId] <= d2Min)
						{
							d2Min = dMatrix[predictedId][actualId];
							targetId = actualId;
						}
					}
					else
					{
						Coordinate actualPoint = actualPoints[actualId];

						// Calculate in steps for increased speed (allows early exit)
						float dx = actualPoint.getX() - x;
						dx *= dx;
						if (dx <= d2Min)
						{
							float dy = actualPoint.getY() - y;
							dy *= dy;
							if (dy <= d2Min)
							{
								float dz = actualPoint.getZ() - z;
								dz *= dz;
								if (dz <= d2Min)
								{
									final float d2 = dx + dy + dz;
									if (d2 <= d2Min)
									{
										d2Min = d2;
										targetId = actualId;
									}
								}
							}
						}
					}
				}

				// Store closest ROI point
				if (targetId > -1)
				{
					assignments.add(new Assignment(targetId, predictedId, d2Min));
				}
			}

			// If there are assignments
			if (!assignments.isEmpty())
			{
				// Pick the closest pair to be assigned
				Collections.sort(assignments);

				// Process in order
				for (Assignment closest : assignments)
				{
					// Skip those that have already been assigned since this will be a lower score.
					// Note at least one assignment should be processed as potential assignments are made 
					// using only unassigned points.
					if (resultAssignment[closest.getPredictedId()] || roiAssignment[closest.getTargetId()])
						continue;

					resultAssignment[closest.getPredictedId()] = true;
					roiAssignment[closest.getTargetId()] = true;

					// If within accuracy then classify as a match
					if (closest.getDistance() < dThreshold)
					{
						tp++;
						fn--;
						fp--;
						rmsd += closest.getDistance();

						if (TP != null)
						{
							TP.add(predictedPoints[closest.getPredictedId()]);
						}
						if (FP != null)
						{
							falsePositives[closest.getPredictedId()] = -1;
						}
						if (FN != null)
						{
							falseNegatives[closest.getTargetId()] = -1;
						}
						if (matches != null)
							matches.add(new PointPair(actualPoints[closest.getTargetId()], predictedPoints[closest
									.getPredictedId()]));
					}
					else
					{
						// No more assignments within the distance threshold
						break;
					}
				}
			}

		} while (!assignments.isEmpty());

		// Add to lists
		if (FP != null)
		{
			for (int i = 0; i < predictedPointsLength; i++)
			{
				if (falsePositives[i] >= 0)
					FP.add(predictedPoints[i]);
			}
		}
		if (FN != null)
		{
			for (int i = 0; i < actualPointsLength; i++)
			{
				if (falseNegatives[i] >= 0)
					FN.add(actualPoints[i]);
			}
		}

		if (tp > 0)
			rmsd = Math.sqrt(rmsd / tp);
		return new MatchResult(tp, fp, fn, rmsd);
	}

	private static int[] ascendingArray(int length)
	{
		int[] array = new int[length];
		for (int i = 0; i < length; i++)
			array[i] = i;
		return array;
	}

	/**
	 * Calculate the match results for the given actual and predicted fluorophore pulses.
	 * Points that are within the distance threshold are identified as a match. The score is calculated using half the
	 * distance threshold and the overlap in time. Assignments are made using the highest scoring matches.
	 * <p>
	 * The total score is stored in the RMSD field of the MatchResult. The number of true positives, false positives and
	 * false negatives are calculated.
	 * 
	 * @param actualPoints
	 * @param predictedPoints
	 * @param dThreshold
	 *            The distance threshold
	 * @param TP
	 *            True Positives
	 * @param FP
	 *            False Positives
	 * @param FN
	 *            False Negatives
	 * @param matches
	 *            The matched true positives (point1 = actual, point2 = predicted)
	 * @return The match results
	 */
	public static MatchResult analyseResults2D(Pulse[] actualPoints, Pulse[] predictedPoints, final double dThreshold,
			List<Pulse> TP, List<Pulse> FP, List<Pulse> FN, List<PointPair> matches)
	{
		// We will use the squared distance for speed
		final double halfDThreshold2 = (dThreshold * 0.5) * (dThreshold * 0.5);
		final float floatDThreshold = (float) dThreshold;
		final double dThreshold2 = dThreshold * dThreshold;

		final int predictedPointsLength = (predictedPoints != null) ? predictedPoints.length : 0;
		final int actualPointsLength = (actualPoints != null) ? actualPoints.length : 0;

		int tp = 0; // true positives (actual with matched predicted point)
		int fp = predictedPointsLength; // false positives (actual with no matched predicted point)
		int fn = actualPointsLength; // false negatives (predicted point with no actual point)
		double score = 0;

		if (predictedPointsLength == 0 || actualPointsLength == 0)
		{
			if (FP != null)
				FP.addAll(asList(predictedPoints));
			if (FN != null)
				FN.addAll(asList(actualPoints));
			return new MatchResult(tp, fp, fn, score);
		}

		// loop over the two arrays assigning the closest unassigned pair
		boolean[] resultAssignment = new boolean[predictedPointsLength];
		boolean[] roiAssignment = new boolean[fn];
		ArrayList<Assignment> assignments = new ArrayList<Assignment>(predictedPointsLength);

		int[] falsePositives = null, falseNegatives = null;
		if (FP != null)
		{
			falsePositives = ascendingArray(predictedPointsLength);
		}
		if (FN != null)
		{
			falseNegatives = ascendingArray(actualPointsLength);
		}

		// Sort by time to allow efficient looping
		Arrays.sort(actualPoints);
		Arrays.sort(predictedPoints);

		// Pre-calculate all-vs-all distance matrix if it can fit in memory
		int size = predictedPointsLength * actualPointsLength;
		float[][] dMatrix = null;
		if (size < 200 * 200)
		{
			dMatrix = new float[predictedPointsLength][actualPointsLength];
			for (int predictedId = 0; predictedId < predictedPointsLength; predictedId++)
			{
				final float x = predictedPoints[predictedId].getX();
				final float y = predictedPoints[predictedId].getY();
				for (int actualId = 0; actualId < actualPointsLength; actualId++)
				{
					dMatrix[predictedId][actualId] = (float) actualPoints[actualId].distance2(x, y);
				}
			}
		}

		do
		{
			assignments.clear();

			// Process each result
			for (int predictedId = 0; predictedId < predictedPointsLength; predictedId++)
			{
				if (resultAssignment[predictedId])
					continue; // Already assigned

				final float x = predictedPoints[predictedId].getX();
				final float y = predictedPoints[predictedId].getY();
				final int start = predictedPoints[predictedId].getStart();
				final int end = predictedPoints[predictedId].getEnd();

				// Find first overlapping pulse
				int actualId = 0;
				while (actualId < actualPointsLength && actualPoints[actualId].getEnd() < start)
				{
					actualId++;
				}

				// Find highest scoring point within the distance limit
				double scoreMax = 0;
				int targetId = -1;
				for (; actualId < actualPointsLength; actualId++)
				{
					if (roiAssignment[actualId])
						continue; // Already assigned
					if (actualPoints[actualId].getStart() > end)
						break; // No more overlap in time

					double d2;
					if (dMatrix != null)
					{
						d2 = dMatrix[predictedId][actualId];
					}
					else
					{
						Coordinate actualPoint = actualPoints[actualId];

						// Calculate in steps for increased speed (allows early exit)
						float dx = Math.abs(actualPoint.getX() - x);
						if (dx > floatDThreshold)
							continue;
						float dy = Math.abs(actualPoint.getY() - y);
						if (dy > floatDThreshold)
							continue;
						d2 = dx * dx + dy * dy;
					}

					// Do we need to exclude using the distance threshold? This is useful for binary classification
					// but will truncate the continuous nature of the score.
					if (d2 > dThreshold2)
						continue;

					double s = predictedPoints[predictedId].score(actualPoints[actualId], d2, halfDThreshold2);
					if (scoreMax < s)
					{
						scoreMax = s;
						targetId = actualId;
					}
				}

				// Store highest scoring point
				if (targetId > -1)
				{
					assignments.add(new Assignment(targetId, predictedId, scoreMax));
				}
			}

			// If there are assignments
			if (!assignments.isEmpty())
			{
				// Process highest scoring first
				Collections.sort(assignments);
				Collections.reverse(assignments);

				// Process in order of score
				for (Assignment closest : assignments)
				{
					// Skip those that have already been assigned since this will be a lower score.
					// Note at least one assignment should be processed as potential assignments are made 
					// using only unassigned points.
					if (resultAssignment[closest.getPredictedId()] || roiAssignment[closest.getTargetId()])
						continue;

					resultAssignment[closest.getPredictedId()] = true;
					roiAssignment[closest.getTargetId()] = true;

					tp++;
					fn--;
					fp--;
					score += closest.getDistance(); // This is the scoreMax (not the distance)

					if (TP != null)
					{
						TP.add(predictedPoints[closest.getPredictedId()]);
					}
					if (FP != null)
					{
						falsePositives[closest.getPredictedId()] = -1;
					}
					if (FN != null)
					{
						falseNegatives[closest.getTargetId()] = -1;
					}
					if (matches != null)
						matches.add(new PointPair(actualPoints[closest.getTargetId()], predictedPoints[closest
								.getPredictedId()]));
				}
			}

		} while (!assignments.isEmpty());

		// Add to lists
		if (FP != null)
		{
			for (int i = 0; i < predictedPointsLength; i++)
			{
				if (falsePositives[i] >= 0)
					FP.add(predictedPoints[i]);
			}
		}
		if (FN != null)
		{
			for (int i = 0; i < actualPointsLength; i++)
			{
				if (falseNegatives[i] >= 0)
					FN.add(actualPoints[i]);
			}
		}

		// Every time-point has the chance to contribute to the score.
		// Normalise score by the maximum of the number of actual/predicted time points.
		// This penalises too few or too many predictions
		int p1 = countTimePoints(actualPoints);
		int p2 = countTimePoints(predictedPoints);
		score /= FastMath.max(p1, p2);

		return new MatchResult(tp, fp, fn, score);
	}

	private static int countTimePoints(Pulse[] actualPoints)
	{
		int p1 = 0;
		for (Pulse p : actualPoints)
			p1 += p.getEnd() - p.getStart() + 1;
		return p1;
	}

	private static Collection<Pulse> asList(Pulse[] points)
	{
		if (points != null)
			return Arrays.asList(points);
		return new ArrayList<Pulse>(0);
	}
}
