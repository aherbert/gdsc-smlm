package gdsc.smlm.ij.plugins;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.math3.util.FastMath;

import gdsc.core.ij.BufferedTextWindow;
import gdsc.core.ij.Utils;
import gdsc.core.match.AUCCalculator;
import gdsc.core.match.BasePoint;
import gdsc.core.match.Coordinate;
import gdsc.core.match.FractionClassificationResult;
import gdsc.core.match.FractionalAssignment;
import gdsc.core.match.ImmutableFractionalAssignment;
import gdsc.core.match.RankedScoreCalculator;
import gdsc.core.utils.FastCorrelator;
import gdsc.core.utils.Maths;
import gdsc.core.utils.RampedScore;
import gdsc.core.utils.Settings;
import gdsc.core.utils.Statistics;
import gdsc.core.utils.StoredData;
import gdsc.smlm.data.config.ConfigurationException;
import gdsc.smlm.data.config.UnitConfig.DistanceUnit;
import gdsc.smlm.data.config.UnitConfig.IntensityUnit;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2014 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.smlm.data.config.FitConfig.DataFilterMethod;
import gdsc.smlm.data.config.FitConfig.DataFilterType;
import gdsc.smlm.engine.FitConfiguration;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.filters.MaximaSpotFilter;
import gdsc.smlm.filters.Spot;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.GaussianOverlapAnalysis;
import gdsc.smlm.ij.plugins.ResultsMatchCalculator.PeakResultPoint;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.ImageConverter;
import gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.procedures.StandardResultProcedure;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.procedure.TIntObjectProcedure;
import gnu.trove.procedure.TIntProcedure;
import gnu.trove.procedure.TObjectProcedure;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.ExtendedGenericDialog;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.Plot2;
import ij.gui.PlotWindow;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;
import ij.text.TextWindow;

/**
 * Filters the benchmark spot image created by CreateData plugin to identify candidates and then assess the filter.
 */
public class BenchmarkSpotFilter implements PlugIn
{
	public static final String TITLE = "Filter Spot Data";

	private static FitConfiguration fitConfig;
	private static FitEngineConfiguration config;
	private static double search = 1;
	private static double minSearch = 1;
	private static double maxSearch = 1;
	private static double border = 1;
	private static boolean[] batchPlot;
	private static String[] batchPlotNames;
	private static String[] SELECTION;
	private static int selection = 2;
	static
	{
		config = new FitEngineConfiguration();
		fitConfig = config.getFitConfiguration();
		batchPlot = new boolean[BatchResult.getNumberOfScores()];
		batchPlotNames = new String[batchPlot.length];
		for (int i = 0; i < batchPlot.length; i++)
			batchPlotNames[i] = BatchResult.getScoreName(i);
		SELECTION = new String[3];
		SELECTION[0] = BatchResult.getScoreName(0);
		SELECTION[1] = BatchResult.getScoreName(1);
		SELECTION[2] = BatchResult.getScoreName(0) + "+" + BatchResult.getScoreName(1);
	}

	private static final String[] MATCHING_METHOD = { "Single", "Multi", "Greedy" };
	@SuppressWarnings("unused")
	private static final int METHOD_SINGLE = 0;
	private static final int METHOD_MULTI = 1;
	private static final int METHOD_GREEDY = 2;

	private static double sAnalysisBorder = 2;
	private static boolean hardBorder = true;
	private static int matchingMethod = METHOD_MULTI;
	/**
	 * The last border used in analysis
	 */
	static Rectangle lastAnalysisBorder;
	private static double upperDistance = 1.5;
	private static double lowerDistance = 0.5;
	private static double upperSignalFactor = 2;
	private static double lowerSignalFactor = 1;
	private static boolean filterRelativeDistances = true;
	private static boolean scoreRelativeDistances = true;
	private double matchDistance;
	private double lowerMatchDistance;
	private static double recallFraction = 100;
	private static boolean showPlot = true;
	private static boolean rankByIntensity = false;
	private static boolean showFailuresPlot = false;
	private static boolean showTP = false;
	private static boolean showFP = false;
	private static boolean showFN = false;
	private static boolean sDebug = false;

	private static boolean batchMean = true;
	private static boolean batchGaussian = true;
	private static boolean batchCircular = false;
	private static boolean batchMedian = false;

	private boolean extraOptions, debug = false, batchMode = false;
	private long time = 0;

	// Cache batch results
	private static Settings batchSettings = null;
	private static ArrayList<BatchResult[]> cachedBatchResults = new ArrayList<BatchResult[]>();

	private static int id = 1;

	private static BufferedTextWindow summaryTable = null, batchSummaryTable = null;

	private ImagePlus imp;
	private MemoryPeakResults results;
	private Gaussian2DPeakResultCalculator calculator;
	private CreateData.SimulationParameters simulationParameters;

	private static TIntObjectHashMap<PSFSpot[]> actualCoordinates = null;
	private static int lastId = -1;
	//private static boolean lastRelativeDistances = false;

	private WindowOrganiser windowOrganiser;

	public static String tablePrefix, resultPrefix;

	// Used by the Benchmark Spot Fit plugin
	private static int filterResultsId = 0;
	static BenchmarkFilterResult filterResult = null;

	public class BenchmarkFilterResult
	{
		public final int simulationId;
		public final int id = ++filterResultsId;
		public TIntObjectHashMap<FilterResult> filterResults;
		public FitEngineConfiguration config;
		public MaximaSpotFilter spotFilter;
		public double auc;
		double[] r, p, j, c;
		int maxIndex;
		int fractionIndex;
		public double[][] cumul;
		public StoredData stats;
		public double auc2;
		public double slope;
		public double[] i1;
		public double[] i2;
		public double[] intensity;
		public boolean relativeDistances;
		public long time;

		public BenchmarkFilterResult(TIntObjectHashMap<FilterResult> filterResults, FitEngineConfiguration config,
				MaximaSpotFilter spotFilter)
		{
			this.simulationId = simulationParameters.id;
			this.filterResults = filterResults;
			this.config = config;
			this.spotFilter = spotFilter;
		}
	}

	private static class BatchResult
	{
		public double auc, j, p, r;
		public long time;
		public DataFilterMethod dataFilter;
		public double param;
		public double search;

		public BatchResult(BenchmarkFilterResult filterResult, DataFilterMethod dataFilter, double param, double search)
		{
			if (filterResult != null)
			{
				this.auc = filterResult.auc;
				this.j = filterResult.j[filterResult.maxIndex];
				this.p = filterResult.p[filterResult.maxIndex];
				this.r = filterResult.r[filterResult.maxIndex];
				this.time = filterResult.time;
			}
			this.dataFilter = dataFilter;
			this.param = param;
			this.search = search;
		}

		public double getScore(int i)
		{
			if (i == 0)
				return auc;
			if (i == 1)
				return j;
			if (i == 2)
				return p;
			if (i == 3)
				return r;
			if (i == 4)
				return time / 1e6;
			return 0;
		}

		public static String getScoreName(int i)
		{
			if (i == 0)
				return "AUC";
			if (i == 1)
				return "Max Jaccard";
			if (i == 2)
				return "Precision (at Max Jaccard)";
			if (i == 3)
				return "Recall (at Max Jaccard)";
			if (i == 4)
				return "Time (ms)";
			return "";
		}

		public static int getNumberOfScores()
		{
			return 5;
		}

		public String getName()
		{
			return String.format("%s:%s", dataFilter.toString(), Utils.rounded(search));
		}
	}

	public class ScoredSpot implements Comparable<ScoredSpot>
	{
		final boolean match;
		double[] scores;
		double score;
		// Total intensity of spots we matched
		private double intensity;
		final float background;
		final Spot spot;
		int fails;

		public ScoredSpot(boolean match, double score, double intensity, Spot spot, float background)
		{
			this.match = match;
			this.spot = spot;
			this.background = background;
			this.fails = 0;
			add(score, intensity);
		}

		public ScoredSpot(boolean match, Spot spot, float background)
		{
			this.match = match;
			this.spot = spot;
			this.background = background;
			this.fails = 0;
		}

		public ScoredSpot(boolean match, Spot spot, float background, int fails)
		{
			this.match = match;
			this.spot = spot;
			this.background = background;
			this.fails = fails;
		}

		/**
		 * Adds the result to the scored spot
		 *
		 * @param d
		 *            the d
		 * @param score
		 *            the score
		 * @param intensity
		 *            the intensity
		 */
		void add(double score, double intensity)
		{
			if (scores == null)
			{
				scores = new double[] { score };
				this.score = score;
				this.intensity = intensity;
			}
			else
			{
				final int size = scores.length;
				scores = Arrays.copyOf(scores, size + 1);
				scores[size] = score;
				this.score += score;
				this.intensity += intensity;
			}
		}

		public double getScore(int i)
		{
			return (scores != null && i < scores.length) ? scores[i] : 0;
		}

		/**
		 * Get the score
		 * 
		 * @return The score
		 */
		double getScore()
		{
			return score;
		}

		/**
		 * Get the opposite of the score.
		 *
		 * @return the anti-score
		 */
		public double antiScore()
		{
			// The use of partial scoring mimicks using multiple distance 
			// thresholds and taking an average of scores. 
			// For a single experiment at a single distance threshold a spot can 
			// match 0, 1 or more actual results. This would be:
			// 0  = TP=0 ,FP=1
			// 1  = TP=1 ,FP=0
			// 2+ = TP=2+,FP=0 (because it doesn't 'not match' anything)
			// So for now I will use FP (anti-score) in the range 0-1.

			// Only ever allow an anti-score in the range 0-1
			return FastMath.max(0, 1 - score);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		public int compareTo(ScoredSpot o)
		{
			if (spot.intensity > o.spot.intensity)
				return -1;
			if (spot.intensity < o.spot.intensity)
				return 1;
			return 0;
		}

		/**
		 * Gets the intensity of the spot without the background.
		 *
		 * @return the intensity
		 */
		public float getIntensity()
		{
			return spot.intensity - background;
		}
	}

	public class FilterResult
	{
		final int frame;
		final RankedScoreCalculator calc;
		final FractionClassificationResult result;
		final ScoredSpot[] spots;
		final PSFSpot[] actual;
		final boolean[] actualAssignment;

		public FilterResult(int frame, RankedScoreCalculator calc, FractionClassificationResult result,
				ScoredSpot[] spots, PSFSpot[] actual, boolean[] actualAssignment)
		{
			this.frame = frame;
			this.calc = calc;
			this.result = result;
			this.spots = spots;
			this.actual = actual;
			this.actualAssignment = actualAssignment;
		}
	}

	/**
	 * Used to allow multi-threading of the PSF overlap computation.
	 */
	private class OverlapWorker implements Runnable
	{
		volatile boolean finished = false;
		final BlockingQueue<Integer> jobs;
		final TIntObjectHashMap<ArrayList<Coordinate>> originalCoordinates;
		final TIntObjectHashMap<PSFSpot[]> coordinates;

		public OverlapWorker(BlockingQueue<Integer> jobs, TIntObjectHashMap<ArrayList<Coordinate>> originalCoordinates)
		{
			this.jobs = jobs;
			this.originalCoordinates = originalCoordinates;
			this.coordinates = new TIntObjectHashMap<PSFSpot[]>();
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		public void run()
		{
			try
			{
				while (true)
				{
					Integer job = jobs.take();
					if (job == null || job.intValue() < 0)
						break;
					if (!finished)
						// Only run if not finished to allow queue to be emptied
						run(job.intValue());
				}
			}
			catch (InterruptedException e)
			{
				System.out.println(e.toString());
				throw new RuntimeException(e);
			}
			finally
			{
				finished = true;
			}
		}

		private void run(int frame)
		{
			if (Utils.isInterrupted())
			{
				finished = true;
				return;
			}

			showProgress();

			// Extract the data
			PSFSpot[] actual = getCoordinates(originalCoordinates, frame);
			coordinates.put(frame, actual);

			if (actual.length == 0)
				return;

			// Here we approximate the PSF as a Gaussian adjusted for square pixels.			

			// Determine spots (1) that have an overlap with other spots (2).
			// In practice this will be any spot within 2SD1 + 2SD2.

			// Pre-compute the adjusted pixel widths for each PSF
			double[] sa = new double[actual.length];
			double[] sa2 = new double[actual.length];
			for (int i = 0; i < actual.length; i++)
			{
				sa[i] = PSFCalculator.squarePixelAdjustment(
						calculator.getStandardDeviation(actual[i].peakResult.getParameters()) * simulationParameters.a,
						simulationParameters.a) / simulationParameters.a;
				sa2[i] = 2 * sa[i];
			}

			double[] allParams = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK * actual.length];
			final int flags = GaussianFunctionFactory.FIT_SIMPLE_NS_NB_FIXED;
			for (int i = 0; i < actual.length; i++)
			{
				// Bounding rectangle for the spot. This serves as the reference frame 0,0
				final float cx = actual[i].getX();
				final float cy = actual[i].getY();

				GaussianOverlapAnalysis overlapAnalysis = null;

				// Check for overlap
				int offset = 0;
				for (int j = 0; j < actual.length; j++)
				{
					if (i == j)
						continue;
					final double dx = (actual[j].getX() - cx);
					final double dy = (actual[j].getY() - cy);
					final double d2 = dx * dx + dy * dy;
					final double threshold = sa2[i] + sa2[j];
					if (d2 <= threshold * threshold)
					{
						// These overlap.

						// Initialise a Gaussian2D function for i
						if (overlapAnalysis == null)
						{
							double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
							params[Gaussian2DFunction.SIGNAL] = actual[i].peakResult.getSignal();
							params[Gaussian2DFunction.X_POSITION] = cx;
							params[Gaussian2DFunction.Y_POSITION] = cy;
							params[Gaussian2DFunction.X_SD] = params[Gaussian2DFunction.Y_SD] = sa[i];
							overlapAnalysis = new GaussianOverlapAnalysis(flags, null, params, 2);
						}

						// Accumulate the function for j
						allParams[offset + Gaussian2DFunction.SIGNAL] = actual[j].peakResult.getSignal();
						allParams[offset + Gaussian2DFunction.X_POSITION] = actual[j].peakResult.getXPosition();
						allParams[offset + Gaussian2DFunction.Y_POSITION] = actual[j].peakResult.getYPosition();
						allParams[offset +
								Gaussian2DFunction.X_SD] = allParams[offset + Gaussian2DFunction.Y_SD] = sa[j];
						offset += Gaussian2DFunction.PARAMETERS_PER_PEAK;
					}
				}

				if (offset != 0)
				{
					final double[] overlapParams = Arrays.copyOf(allParams, 1 + offset);
					overlapAnalysis.add(overlapParams, false);
					actual[i].backgroundOffset = (float) overlapAnalysis.getWeightedbackground();

					// This is not currently used.
					// Computation of this would depend on how a filter is estimating the signal. 
					// The offset should be computed with the same method to create a 'fair' signal offset.
					//actual[i].intensityOffset = ?
				}
			}
		}

		/**
		 * Return an array of PSF spots for the given time point. Returns an empty array if there are no coordinates.
		 * 
		 * @param coords
		 * @param t
		 * @return The array list
		 */
		public PSFSpot[] getCoordinates(TIntObjectHashMap<ArrayList<Coordinate>> coords, Integer t)
		{
			ArrayList<Coordinate> list1 = coords.get(t);
			if (list1 != null)
			{
				PSFSpot[] list2 = new PSFSpot[list1.size()];
				int i = 0;
				for (Coordinate c : list1)
				{
					PeakResultPoint p = (PeakResultPoint) c;
					list2[i++] = new PSFSpot(p.getTime(), p.getX(), p.getY(), p.peakResult);
				}
				return list2;
			}
			else
			{
				return new PSFSpot[0];
			}
		}
	}

	/**
	 * Used to allow multi-threading of the fitting method
	 */
	private class Worker implements Runnable
	{
		volatile boolean finished = false;
		final BlockingQueue<Integer> jobs;
		final ImageStack stack;
		final MaximaSpotFilter spotFilter;
		final float background;
		final TIntObjectHashMap<FilterResult> results;

		float[] data = null;
		long time = 0;

		public Worker(BlockingQueue<Integer> jobs, ImageStack stack, MaximaSpotFilter spotFilter, float background)
		{
			this.jobs = jobs;
			this.stack = stack;
			this.spotFilter = spotFilter.clone();
			this.results = new TIntObjectHashMap<FilterResult>();
			this.background = background;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		public void run()
		{
			try
			{
				while (true)
				{
					Integer job = jobs.take();
					if (job == null || job.intValue() < 0)
						break;
					if (!finished)
						// Only run if not finished to allow queue to be emptied
						run(job.intValue());
				}
			}
			catch (InterruptedException e)
			{
				System.out.println(e.toString());
				throw new RuntimeException(e);
			}
			finally
			{
				finished = true;
			}
		}

		private void run(int frame)
		{
			if (Utils.isInterrupted())
			{
				finished = true;
				return;
			}

			showProgress();

			// Extract the data
			data = ImageConverter.getData(stack.getPixels(frame), stack.getWidth(), stack.getHeight(), null, data);

			// Use this code to subtract the background before filtering. This produces different results
			// from using the raw data.
			//float background = this.background;
			//if (background != 0)
			//{
			//	for (int i = stack.getWidth() * stack.getHeight(); i-- > 0;)
			//		data[i] -= background;
			//}
			//background = 0;

			long start = System.nanoTime();
			Spot[] spots = spotFilter.rank(data, stack.getWidth(), stack.getHeight());
			time += System.nanoTime() - start;

			// Debug the candidates
			//			if (debug && frame == 5)
			//			{
			//				StringBuilder sb = new StringBuilder();
			//				sb.append(spotFilter.getDescription()).append("\n");
			//				for (int i = 0; i < spots.length; i++)
			//				{
			//					sb.append(String.format("Fit %d [%d,%d = %.1f]\n", i, spots[i].x, spots[i].y, spots[i].intensity));
			//				}
			//				System.out.print(sb.toString());
			//			}

			// Score the spots that are matches
			PSFSpot[] actual = actualCoordinates.get(frame);
			if (actual == null)
				actual = new PSFSpot[0];

			// We do not remove results at the border from analysis.
			// We can just mark them to contribute less to the score.
			// TODO - Create a Tukey window weighting from the border to the edge.
			// The type of weighting could be user configurable, e.g. Hard, Tukey, Linear, etc.
			final double[] actualWeight = new double[actual.length];
			final double[] spotsWeight = new double[spots.length];
			double actualLength = actual.length;
			double spotsLength = spots.length;
			if (lastAnalysisBorder.x > 0)
			{
				actualLength = spotsLength = 0;

				final int analysisBorder = lastAnalysisBorder.x;
				final int xlimit = lastAnalysisBorder.x + lastAnalysisBorder.width;
				final int ylimit = lastAnalysisBorder.y + lastAnalysisBorder.height;
				for (int i = 0; i < actual.length; i++)
				{
					final PSFSpot c = actual[i];
					actualWeight[i] = 1;
					if (c.getX() < analysisBorder || c.getX() > xlimit || c.getY() < analysisBorder ||
							c.getY() > ylimit)
						// TODO - better weighting
						actualWeight[i] = 0;
					actualLength += actualWeight[i];
				}

				for (int i = 0; i < spots.length; i++)
				{
					final Spot s = spots[i];
					spotsWeight[i] = 1;
					if (s.x < analysisBorder || s.x > xlimit || s.y < analysisBorder || s.y > ylimit)
						// TODO - better weighting
						spotsWeight[i] = 0;
					spotsLength += spotsWeight[i];
				}

				// TODO - option for hard border to match the old scoring.
				// Create smaller arrays using only those with a weighting of 1.
				if (hardBorder)
				{
					PSFSpot[] actual2 = new PSFSpot[actual.length];
					int j = 0;
					for (int i = 0; i < actual.length; i++)
						if (actualWeight[i] == 1)
						{
							actualWeight[j] = 1;
							actual2[j++] = actual[i];
						}
					actual = Arrays.copyOf(actual2, j);

					Spot[] spots2 = new Spot[spots.length];
					j = 0;
					for (int i = 0; i < spots.length; i++)
						if (spotsWeight[i] == 1)
						{
							spotsWeight[j] = 1;
							spots2[j++] = spots[i];
						}
					spots = Arrays.copyOf(spots2, j);

					// Update lengths
					actualLength = actual.length;
					spotsLength = spots.length;
				}
			}
			else
			{
				Arrays.fill(actualWeight, 1);
				Arrays.fill(spotsWeight, 1);
			}

			ScoredSpot[] scoredSpots = new ScoredSpot[spots.length];
			FractionClassificationResult result;
			RankedScoreCalculator calc = null;

			// Store the count of false positives since the last true positive
			int fails = 0;

			final int nActual = actual.length;
			final boolean[] actualAssignment = new boolean[nActual];
			if (actual.length > 0)
			{
				SpotCoordinate[] predicted = getCoordinates(spots);

				// Use the distance to the true location to score the candidate
				final RampedScore score = new RampedScore(lowerMatchDistance, matchDistance);
				final RampedScore signalScore = (upperSignalFactor > 0)
						? new RampedScore(lowerSignalFactor, upperSignalFactor) : null;

				// Candidates may be close to many localisations. In order to compute the signal
				// factor correctly we have computed the signal offset for each spot with overlapping PSFs.
				// This is used to raise the spot intensity when computing the signal factor.

				// Compute assignments
				ArrayList<FractionalAssignment> fractionalAssignments = new ArrayList<FractionalAssignment>(
						predicted.length * 3);

				final double dmin = matchDistance * matchDistance;
				final int nPredicted = predicted.length;
				for (int j = 0; j < nPredicted; j++)
				{
					final float x = predicted[j].getX();
					final float y = predicted[j].getY();
					// Any spots that match 
					for (int i = 0; i < nActual; i++)
					{
						final double dx = (x - actual[i].getX());
						final double dy = (y - actual[i].getY());
						final double d2 = dx * dx + dy * dy;
						if (d2 <= dmin)
						{
							final double d = Math.sqrt(d2);
							double s = score.score(d);
							final double intensity = getIntensity(actual[i]);
							if (signalScore != null)
							{
								// Adjust intensity using the surrounding PSF contributions
								//final double rsf = intensity / (predicted[j].spot.intensity - background);
								final double rsf = (actual[i].backgroundOffset + intensity) /
										(predicted[j].spot.intensity - background);
								// Normalise so perfect is zero
								final double sf = Math.abs((rsf < 1) ? 1 - 1 / rsf : rsf - 1);
								s *= signalScore.score(sf);
							}
							s = RampedScore.flatten(s, 256);

							double distance = 1 - s;
							if (distance == 0)
							{
								// In the case of a match below the distance and signal factor thresholds
								// the distance will be 0. To distinguish between candidates all below 
								// the thresholds just take the closest.
								// We know d2 is below dmin so we subtract the delta.
								distance -= (dmin - d2);
							}

							// Store the match
							fractionalAssignments.add(new ImmutableFractionalAssignment(i, j, distance, s));
						}
					}
				}

				FractionalAssignment[] assignments = fractionalAssignments
						.toArray(new FractionalAssignment[fractionalAssignments.size()]);
				calc = new RankedScoreCalculator(assignments);

				// Assign matches
				double tp = 0;
				double fp;
				final double[] predictedScore = new double[nPredicted];
				if (matchingMethod == METHOD_GREEDY)
				{
					// Spots can match as many actual results as they can, first match wins
					int nA = nActual;

					for (FractionalAssignment a : assignments)
					{
						final int i = a.getTargetId();
						if (!actualAssignment[i])
						{
							final int j = a.getPredictedId();
							actualAssignment[i] = true;
							final double intensity = getIntensity(actual[i]);
							final double tpScore = a.getScore() * actualWeight[i];
							if (scoredSpots[j] == null)
								scoredSpots[j] = new ScoredSpot(true, tpScore, intensity, spots[j], background);
							else
								scoredSpots[j].add(tpScore, intensity);
							tp += tpScore;
							predictedScore[j] += tpScore;
							if (--nA == 0)
								break;
						}
					}
				}
				else if (matchingMethod == METHOD_MULTI)
				{
					// Spots can match as many actual results as they can. Matching is iterative
					// so only the best match is computed for each spot per round.
					int nA = nActual;

					// Flag to indicate a match was made
					boolean processAgain = true;
					OUTER: while (processAgain)
					{
						processAgain = false;
						final boolean[] predictedAssignment = new boolean[nPredicted];
						for (FractionalAssignment a : assignments)
						{
							final int i = a.getTargetId();
							if (!actualAssignment[i])
							{
								final int j = a.getPredictedId();
								if (!predictedAssignment[j])
								{
									actualAssignment[i] = true;
									predictedAssignment[j] = true;
									processAgain = true;
									final double intensity = getIntensity(actual[i]);
									final double tpScore = a.getScore() * actualWeight[i];
									if (scoredSpots[j] == null)
										scoredSpots[j] = new ScoredSpot(true, tpScore, intensity, spots[j], background);
									else
										scoredSpots[j].add(tpScore, intensity);
									tp += tpScore;
									predictedScore[j] += tpScore;
									if (--nA == 0)
										break OUTER;
								}
							}
						}
					}
				}
				else
				{
					// matchingMethod == METHOD_SINGLE
					// Spots can match only one actual result

					final boolean[] predictedAssignment = new boolean[nPredicted];

					int nP = nPredicted;
					int nA = nActual;

					for (FractionalAssignment a : assignments)
					{
						final int i = a.getTargetId();
						if (!actualAssignment[i])
						{
							final int j = a.getPredictedId();
							if (!predictedAssignment[j])
							{
								actualAssignment[i] = true;
								predictedAssignment[j] = true;
								final double tpScore = a.getScore() * actualWeight[i];
								scoredSpots[j] = new ScoredSpot(true, tpScore, getIntensity(actual[i]), spots[j],
										background);
								tp += tpScore;
								predictedScore[j] = tpScore;
								if (--nA == 0 || --nP == 0)
									break;
							}
						}
					}
				}

				// Compute the FP. 
				// Note: The TP score is the match score multiplied by the weight of the actual spot.
				// Although a predicted point can accumulate more than its weight for TP matches (due 
				// to multiple matching and the fuzzy border), no predicted point can score less than 
				// its weight. This means: 
				// TP + FN = nActual
				// TP + FP >= nPredicted (due to fuzzy border)
				//
				// Note: This score scenario is the same as if using a hard border exclusion and 
				// multi-matching where TP can be above 1 but FP is not.
				//
				// This does mean that a spot in the border that matches a result in the image
				// can have a high TP score and very low FP score. Hopefully this has little effect
				// in practice.				

				fp = spotsLength;
				for (int j = 0; j < predictedScore.length; j++)
				{
					if (predictedScore[j] > spotsWeight[j])
						predictedScore[j] = spotsWeight[j];
					fp -= predictedScore[j];
				}

				result = new FractionClassificationResult(tp, fp, 0, actualLength - tp);

				// Store the number of fails (negatives) before each positive 
				for (int i = 0; i < spots.length; i++)
				{
					if (scoredSpots[i] == null)
					{
						scoredSpots[i] = new ScoredSpot(false, spots[i], fails++);
					}
					else
					{
						scoredSpots[i].fails = fails++;
						if (scoredSpots[i].match)
						{
							//if (fails > 60)
							//	System.out.printf("%d @ %d : %d,%d\n", fails, frame, scoredSpots[i].spot.x,
							//		scoredSpots[i].spot.y);
							fails = 0;
						}
					}
				}
			}
			else
			{
				// No results.
				// All spots are false positives
				result = new FractionClassificationResult(0, spotsLength, 0, 0);
				for (int i = 0; i < spots.length; i++)
				{
					scoredSpots[i] = new ScoredSpot(false, spots[i], fails++);
				}
			}

			if (debug)
			{
				System.out.printf("Frame %d : N = %d, TP = %.2f, FP = %.2f, R = %.2f, P = %.2f\n", frame, actualLength,
						result.getTP(), result.getFP(), result.getRecall(), result.getPrecision());
			}

			results.put(frame, new FilterResult(frame, calc, result, scoredSpots, actual, actualAssignment));
		}

		private double getIntensity(final PSFSpot p)
		{
			// Use the amplitude as all spot filters currently estimate the height, not the total signal
			final double intensity = calculator.getAmplitude(p.peakResult.getParameters());
			return intensity;
		}

		private SpotCoordinate[] getCoordinates(Spot[] spots)
		{
			SpotCoordinate[] coords = new SpotCoordinate[spots.length];
			for (int i = 0; i < spots.length; i++)
				coords[i] = new SpotCoordinate(i, spots[i]);
			return coords;
		}

		private class SpotCoordinate extends BasePoint
		{
			Spot spot;

			public SpotCoordinate(int id, Spot spot)
			{
				// Centre ion the middle of the pixel
				super(spot.x + 0.5f, spot.y + 0.5f);
				this.spot = spot;
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		extraOptions = Utils.isExtraOptions();
		batchMode = "batch".equals(arg);

		simulationParameters = CreateData.simulationParameters;
		if (simulationParameters == null)
		{
			IJ.error(TITLE, "No benchmark spot parameters in memory");
			return;
		}
		imp = CreateData.getImage();
		if (imp == null)
		{
			IJ.error(TITLE, "No benchmark image");
			return;
		}
		results = CreateData.getResults();
		if (results == null)
		{
			IJ.error(TITLE, "No benchmark results in memory");
			return;
		}

		// Set-up the converters
		try
		{
			if (results.getCalibration() == null ||
					results.getCalibrationReader().getDistanceUnit() != DistanceUnit.PIXEL)
				throw new ConfigurationException("Require results in pixel distance units");

			int flags = Gaussian2DPeakResultHelper.AMPLITUDE;
			calculator = Gaussian2DPeakResultHelper.create(results.getPSF(), results.getCalibration(), flags);
		}
		catch (ConfigurationException e)
		{
			IJ.error(TITLE, "Bad configuration: " + e.getMessage());
			return;
		}

		if (!showDialog())
			return;

		// Clear old results to free memory
		if (filterResult != null)
		{
			filterResult.filterResults.clear();
			filterResult.filterResults = null;
			filterResult = null;
		}

		// For graphs
		windowOrganiser = new WindowOrganiser();

		if (batchMode)
		{
			// Batch mode to test enumeration of filters
			final double sd = simulationParameters.s / simulationParameters.a;
			final int limit = (int) Math.floor(3 * sd);

			double[] searchParam = getRange(minSearch, maxSearch, 1);

			// Continuous parameters
			double[] pEmpty = new double[0];
			double[] mParam = (batchMean) ? getRange(limit, 0.05) : pEmpty;
			double[] gParam = (batchGaussian) ? getRange(limit, 0.05) : pEmpty;

			// Less continuous parameters
			double[] cParam = (batchCircular) ? getRange(limit, 0.5) : pEmpty;

			// Discrete parameters
			double[] medParam = (batchMedian) ? getRange(limit, 1) : pEmpty;

			setupProgress(imp.getImageStackSize() * searchParam.length *
					(mParam.length + gParam.length + cParam.length + medParam.length), "Frame");

			ArrayList<BatchResult[]> batchResults = new ArrayList<BatchResult[]>(cachedBatchResults.size());
			config.setDataFilterType(DataFilterType.SINGLE);
			for (double search : searchParam)
			{
				// Run all, store the results for plotting.
				// Allow re-use of these if they are cached to allow quick reanalysis of results.
				config.setSearch(search);
				if (batchMean)
					batchResults.add(addToCache(DataFilterMethod.MEAN, mParam, search));
				if (batchGaussian)
					batchResults.add(addToCache(DataFilterMethod.GAUSSIAN, gParam, search));
				if (batchCircular)
					batchResults.add(addToCache(DataFilterMethod.CIRCULAR_MEAN, cParam, search));
				if (batchMean)
					batchResults.add(addToCache(DataFilterMethod.MEDIAN, medParam, search));
			}

			IJ.showProgress(-1);
			IJ.showStatus("");
			if (Utils.isInterrupted())
				return;

			// Analysis options
			GenericDialog gd = new GenericDialog(TITLE);
			gd.addMessage("Choose performance plots:");
			for (int i = 0; i < batchPlot.length; i++)
				gd.addCheckbox(batchPlotNames[i], batchPlot[i]);

			gd.addChoice("Selection", SELECTION, SELECTION[selection]);
			gd.addCheckbox("Show_plots", showPlot);
			gd.addCheckbox("Plot_rank_by_intensity", rankByIntensity);
			gd.addCheckbox("Show_failures_plots", showFailuresPlot);
			gd.addCheckbox("Show_TP", showTP);
			gd.addCheckbox("Show_FP", showFP);
			gd.addCheckbox("Show_FN", showFN);

			gd.showDialog();
			if (gd.wasCanceled())
				return;

			for (int i = 0; i < batchPlot.length; i++)
				batchPlot[i] = gd.getNextBoolean();
			selection = gd.getNextChoiceIndex();
			showPlot = gd.getNextBoolean();
			rankByIntensity = gd.getNextBoolean();
			showFailuresPlot = gd.getNextBoolean();
			showTP = gd.getNextBoolean();
			showFP = gd.getNextBoolean();
			showFN = gd.getNextBoolean();

			// Plot charts			
			for (int i = 0; i < batchPlot.length; i++)
				plot(i, batchResults);

			// Store in global singleton
			filterResult = analyse(batchResults);
		}
		else
		{
			// Single filter mode
			setupProgress(imp.getImageStackSize(), "Frame");

			filterResult = run(config, filterRelativeDistances);
		}

		IJ.showProgress(-1);
		IJ.showStatus("");

		getTable(false).flush();

		if (filterResult == null)
			return;

		// Store a clone of the config
		filterResult.config = filterResult.config.clone();

		// Debugging the matches
		if (debug)
			addSpotsToMemory(filterResult.filterResults);

		if (showFailuresPlot)
			showFailuresPlot(filterResult);
		if (showPlot)
			showPlot(filterResult);
		if (isShowOverlay())
			showOverlay(imp, filterResult);

		windowOrganiser.tile();
	}

	private BatchResult[] addToCache(DataFilterMethod dataFilter, double[] param, double search)
	{
		for (BatchResult[] batchResult : cachedBatchResults)
		{
			if (batchResult == null || batchResult.length == 0)
				continue;
			if (batchResult[0].dataFilter == dataFilter && batchResult[0].search == search)
				return batchResult;
		}
		BatchResult[] batchResult = run(dataFilter, param, search);
		cachedBatchResults.add(batchResult);
		return batchResult;
	}

	private void plot(int i, ArrayList<BatchResult[]> batchResults)
	{
		if (!batchPlot[i])
			return;

		Color[] colors = new Color[] { Color.red, Color.gray, Color.green, Color.blue, Color.magenta };

		String name = batchPlotNames[i];
		String title = TITLE + " Performance " + name;
		Plot plot = new Plot(title, "Relative width", name);
		final double scale = 1.0 / config.getHWHMMin();
		for (BatchResult[] batchResult : batchResults)
		{
			if (batchResult == null || batchResult.length == 0)
				continue;
			float[][] data = extractData(batchResult, i, scale);
			int colorIndex = batchResult[0].dataFilter.ordinal();
			plot.setColor(colors[colorIndex]);
			colors[colorIndex] = colors[colorIndex].darker();
			plot.addPoints(data[0], data[1], null, (batchResult.length > 1) ? Plot.LINE : Plot.CIRCLE,
					batchResult[0].getName());
		}
		plot.setColor(Color.black);
		plot.addLegend(null);
		if (name.contains("Time"))
			plot.setAxisYLog(true);
		PlotWindow pw = Utils.display(title, plot);
		plot.setLimitsToFit(true); // Seems to only work after drawing
		if (Utils.isNewWindow())
			windowOrganiser.add(pw);
	}

	private float[][] extractData(BatchResult[] batchResult, int index, double scale)
	{
		float[][] data = new float[2][batchResult.length];
		for (int i = 0; i < batchResult.length; i++)
		{
			data[0][i] = (float) (batchResult[i].param * scale);
			data[1][i] = (float) batchResult[i].getScore(index);
		}
		return data;
	}

	private BenchmarkFilterResult analyse(ArrayList<BatchResult[]> batchResults)
	{
		// Support z-score of AUC and Max. Jaccard combined.
		// For this wee need the statistics of the population of scores. 
		double[][] stats = getStats(batchResults);

		double max = 0;
		DataFilterMethod dataFilter = null;
		double search = 0;
		double param = 0;
		for (BatchResult[] batchResult : batchResults)
		{
			if (batchResult == null || batchResult.length == 0)
				continue;
			double[][] data = extractData(batchResult, selection, stats);
			int maxi = 0;
			for (int i = 0; i < batchResult.length; i++)
			{
				if (data[1][maxi] < data[1][i])
					maxi = i;
			}
			if (max < data[1][maxi])
			{
				max = data[1][maxi];
				dataFilter = batchResult[0].dataFilter;
				search = batchResult[0].search;
				param = data[0][maxi];
			}
		}

		if (dataFilter != null)
		{
			// Convert the absolute distance to be relative to the PSF width
			param = Maths.round(param / config.getHWHMMin(), 0.001);

			final double hwhmMax = config.getHWHMMax();
			// Convert absolute search distance to relative
			config.setSearch(Maths.round(search / hwhmMax, 0.001));

			if (filterRelativeDistances)
			{
				// If relative distances were specified then we can use the input values
				config.setBorder(border);
			}
			else
			{
				// Otherwise we must adjust the input values to convert the absolute values to relative
				config.setBorder(Maths.round(border / hwhmMax, 0.001));
			}

			// Run the filter using relative distances
			config.setDataFilter(dataFilter, param, false, 0);
			BenchmarkFilterResult result = run(config, true, true);
			getTable(true).flush();
			return result;
		}

		return null;
	}

	private double[][] getStats(ArrayList<BatchResult[]> batchResults)
	{
		if (selection < 2)
			return null;

		double[][] stats = new double[2][2];
		for (int index = 0; index < stats.length; index++)
		{
			Statistics s = new Statistics();
			for (BatchResult[] batchResult : batchResults)
			{
				if (batchResult == null || batchResult.length == 0)
					continue;
				for (int i = 0; i < batchResult.length; i++)
				{
					s.add(batchResult[i].getScore(index));
				}
			}
			stats[index][0] = s.getMean();
			stats[index][1] = s.getStandardDeviation();
		}
		return stats;
	}

	private double[][] extractData(BatchResult[] batchResult, int index, double[][] stats)
	{
		double[][] data = new double[2][batchResult.length];
		for (int i = 0; i < batchResult.length; i++)
		{
			data[0][i] = batchResult[i].param;
			data[1][i] = getScore(batchResult[i], index, stats);
		}
		return data;
	}

	private double getScore(BatchResult batchResult, int index, double[][] stats)
	{
		if (stats == null)
			return batchResult.getScore(index);
		// Z-score of all metrics combined
		double z = 0;
		for (int i = 0; i < stats.length; i++)
			z += (batchResult.getScore(i) - stats[i][0]) / stats[i][1];
		return z;
	}

	private double[] getRange(final double min, final double max, final double interval)
	{
		double[] param;
		int c = (int) ((max - min) / interval) + 2;
		param = new double[c];
		param[0] = min;
		int i = 0;
		while (param[i++] < max)
		{
			param[i] = min + i * interval;
		}
		return (i < c) ? Arrays.copyOf(param, i) : param;
	}

	private double[] getRange(final int limit, final double interval)
	{
		double[] param;
		int c = (int) (limit / interval);
		param = new double[c];
		for (int i = 1; i <= c; i++)
		{
			param[i - 1] = i * interval;
		}
		return param;
	}

	private BatchResult[] run(DataFilterMethod dataFilter, double[] param, double search)
	{
		progressPrefix = new BatchResult(null, dataFilter, 0, search).getName();
		BatchResult[] result = new BatchResult[param.length];

		config.setSearch(search);

		for (int i = 0; i < param.length; i++)
		{
			config.setDataFilter(dataFilter, param[i], false, 0);
			BenchmarkFilterResult filterResult = run(config, false);
			if (filterResult == null)
				return null;
			result[i] = new BatchResult(filterResult, dataFilter, param[i], search);
		}
		return result;
	}

	private boolean showDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		StringBuilder sb = new StringBuilder();
		sb.append("Finds spots in the benchmark image created by CreateData plugin.\n");
		final double sa = getSa() / simulationParameters.a;
		sb.append("PSF width = ").append(Utils.rounded(simulationParameters.s / simulationParameters.a))
				.append(" px (sa = ").append(Utils.rounded(sa)).append(" px). HWHM = ")
				.append(Utils.rounded(sa * Gaussian2DFunction.SD_TO_HWHM_FACTOR)).append(" px\n");
		sb.append("Simulation depth = ").append(Utils.rounded(simulationParameters.depth)).append(" nm");
		if (simulationParameters.fixedDepth)
			sb.append(" (fixed)");
		sb.append("\n \nConfigure the spot filter:");
		gd.addMessage(sb.toString());

		if (batchMode)
		{
			// Support enumeration of single spot filters
			gd.addCheckbox("Mean", batchMean);
			gd.addCheckbox("Gaussian", batchGaussian);
			gd.addCheckbox("Circular", batchCircular);
			gd.addCheckbox("Median", batchMedian);
			gd.addSlider("Min_search_width", 1, 4, minSearch);
			gd.addSlider("Max_search_width", 1, 4, maxSearch);
			gd.addCheckbox("Filter_relative_distances (to HWHM)", filterRelativeDistances);
		}
		else
		{
			gd.addChoice("Spot_filter_type", SettingsManager.getDataFilterTypeNames(),
					config.getDataFilterType().ordinal());
			gd.addChoice("Spot_filter", SettingsManager.getDataFilterMethodNames(),
					config.getDataFilterMethod(0).ordinal());

			gd.addCheckbox("Filter_relative_distances (to HWHM)", config.getDataFilterAbsolute(0));
			gd.addSlider("Smoothing", 0, 2.5, config.getSmooth(0));
			gd.addSlider("Search_width", 1, 4, search);
		}
		gd.addSlider("Border", 0, 5, border);
		gd.addCheckbox("Hard_border", hardBorder);

		gd.addMessage("Scoring options:");
		gd.addCheckbox("Score_relative_distances (to HWHM)", scoreRelativeDistances);
		gd.addSlider("Analysis_border", 0, 5, sAnalysisBorder);
		gd.addChoice("Matching_method", MATCHING_METHOD, MATCHING_METHOD[matchingMethod]);
		gd.addSlider("Match_distance", 0.5, 3.5, upperDistance);
		gd.addSlider("Lower_distance", 0, 3.5, lowerDistance);
		gd.addSlider("Signal_factor", 0, 3.5, upperSignalFactor);
		gd.addSlider("Lower_factor", 0, 3.5, lowerSignalFactor);
		gd.addSlider("Recall_fraction", 50, 100, recallFraction);
		if (!batchMode)
		{
			gd.addCheckbox("Show_plots", showPlot);
			gd.addCheckbox("Plot_rank_by_intensity", rankByIntensity);
			gd.addCheckbox("Show_failures_plots", showFailuresPlot);
			gd.addCheckbox("Show_TP", showTP);
			gd.addCheckbox("Show_FP", showFP);
			gd.addCheckbox("Show_FN", showFN);
		}
		if (extraOptions)
			gd.addCheckbox("Debug", sDebug);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		fitConfig.setInitialPeakStdDev(Maths.round(sa));

		if (batchMode)
		{
			batchMean = gd.getNextBoolean();
			batchGaussian = gd.getNextBoolean();
			batchCircular = gd.getNextBoolean();
			batchMedian = gd.getNextBoolean();

			if (!(batchMean || batchGaussian || batchCircular || batchMedian))
				return false;

			minSearch = gd.getNextNumber();
			maxSearch = gd.getNextNumber();
			filterRelativeDistances = gd.getNextBoolean();
		}
		else
		{
			config.setDataFilterType(SettingsManager.getDataFilterTypeValues()[gd.getNextChoiceIndex()]);
			filterRelativeDistances = gd.getNextBoolean();
			config.setDataFilter(SettingsManager.getDataFilterMethodValues()[gd.getNextChoiceIndex()],
					Maths.round(Math.abs(gd.getNextNumber()), 0.001), !filterRelativeDistances, 0);
			search = gd.getNextNumber();
		}
		border = gd.getNextNumber();
		hardBorder = gd.getNextBoolean();
		scoreRelativeDistances = gd.getNextBoolean();
		sAnalysisBorder = Math.abs(gd.getNextNumber());
		matchingMethod = gd.getNextChoiceIndex();
		upperDistance = Math.abs(gd.getNextNumber());
		lowerDistance = Math.abs(gd.getNextNumber());
		upperSignalFactor = Math.abs(gd.getNextNumber());
		lowerSignalFactor = Math.abs(gd.getNextNumber());
		recallFraction = Math.abs(gd.getNextNumber());
		if (!batchMode)
		{
			showPlot = gd.getNextBoolean();
			rankByIntensity = gd.getNextBoolean();
			showFailuresPlot = gd.getNextBoolean();
			showTP = gd.getNextBoolean();
			showFP = gd.getNextBoolean();
			showFN = gd.getNextBoolean();
		}
		if (extraOptions)
			debug = sDebug = gd.getNextBoolean();

		if (gd.invalidNumber())
			return false;

		if (lowerDistance > upperDistance)
			lowerDistance = upperDistance;
		if (lowerSignalFactor > upperSignalFactor)
			lowerSignalFactor = upperSignalFactor;

		if (batchMode)
		{
			// Clear the cached results if the setting changed
			Settings settings = new Settings(simulationParameters.id, filterRelativeDistances,
					//search, maxSearch, // Ignore search distance for smart caching 
					border, scoreRelativeDistances, sAnalysisBorder, hardBorder, matchingMethod, upperDistance,
					lowerDistance, upperSignalFactor, lowerSignalFactor, recallFraction);
			if (!settings.equals(batchSettings))
			{
				cachedBatchResults.clear();
			}
			batchSettings = settings;

			// Analysis during batch mode will always be done with absolute distances. 
			// However we must ensure the border distance is 
			// relative (if requested) so that the results are consistent with single-filter mode.
			if (filterRelativeDistances)
			{
				final double hwhmMax = config.getHWHMMax();
				config.setBorder(Maths.round(border * hwhmMax, 0.001));
			}
			else
			{
				config.setBorder(Maths.round(border, 0.001));
			}
		}
		else
		{
			config.setSearch(Maths.round(search, 0.001));
			config.setBorder(Maths.round(border, 0.001));

			// Single filter ...
			// Allow more complicated filters to be configured
			if (!PeakFit.configureDataFilter(config, PeakFit.FLAG_NO_SAVE))
				return false;
		}

		int analysisBorder;
		if (scoreRelativeDistances)
		{
			// Convert distance to PSF standard deviation units
			final double hwhmMax = config.getHWHMMax();
			matchDistance = upperDistance * hwhmMax;
			lowerMatchDistance = lowerDistance * hwhmMax;
			analysisBorder = (int) (sAnalysisBorder * hwhmMax);
		}
		else
		{
			matchDistance = upperDistance;
			lowerMatchDistance = lowerDistance;
			analysisBorder = (int) (sAnalysisBorder);
		}

		if (analysisBorder > 0)
		{
			lastAnalysisBorder = new Rectangle(analysisBorder, analysisBorder, imp.getWidth() - 2 * analysisBorder,
					imp.getHeight() - 2 * analysisBorder);
		}
		else
		{
			lastAnalysisBorder = new Rectangle(imp.getWidth(), imp.getHeight());
		}

		return true;
	}

	private double getSa()
	{
		final double sa = PSFCalculator.squarePixelAdjustment(simulationParameters.s, simulationParameters.a);
		return sa;
	}

	/** The total progress. */
	private int progress, stepProgress, totalProgress;
	private String progressPrefix;

	private void setupProgress(int total, String prefix)
	{
		totalProgress = total;
		stepProgress = Utils.getProgressInterval(totalProgress);
		progressPrefix = prefix;
		progress = 0;
	}

	/**
	 * Show progress.
	 */
	private synchronized void showProgress()
	{
		if (progress % stepProgress == 0)
		{
			//if (Utils.showStatus(String.format("%s: %d / %d", progressPrefix, progress, totalProgress)))
			if (Utils.showStatus(progressPrefix))
				IJ.showProgress(progress, totalProgress);
		}
		progress++;
	}

	private BenchmarkFilterResult run(FitEngineConfiguration config, boolean relativeDistances)
	{
		return run(config, relativeDistances, false);
	}

	private BenchmarkFilterResult run(FitEngineConfiguration config, boolean relativeDistances, boolean batchSummary)
	{
		if (Utils.isInterrupted())
			return null;

		MaximaSpotFilter spotFilter = config.createSpotFilter(relativeDistances);

		// Extract all the results in memory into a list per frame. This can be cached
		if (lastId != simulationParameters.id) // || lastRelativeDistances != relativeDistances)
		{
			// Always use float coordinates.
			// The Worker adds a pixel offset for the spot coordinates.
			TIntObjectHashMap<ArrayList<Coordinate>> coordinates = ResultsMatchCalculator.getCoordinates(results,
					false);
			actualCoordinates = new TIntObjectHashMap<PSFSpot[]>();
			lastId = simulationParameters.id;
			//lastRelativeDistances = relativeDistances;

			// Store these so we can reset them
			final int total = totalProgress;
			final String prefix = progressPrefix;

			// Spot PSFs may overlap so we must determine the amount of signal overlap and amplitude effect 
			// for each spot...
			IJ.showStatus("Computing PSF overlap ...");

			final int nThreads = Prefs.getThreads();
			final BlockingQueue<Integer> jobs = new ArrayBlockingQueue<Integer>(nThreads * 2);
			List<OverlapWorker> workers = new LinkedList<OverlapWorker>();
			List<Thread> threads = new LinkedList<Thread>();
			for (int i = 0; i < nThreads; i++)
			{
				OverlapWorker worker = new OverlapWorker(jobs, coordinates);
				Thread t = new Thread(worker);
				workers.add(worker);
				threads.add(t);
				t.start();
			}

			// Process the frames
			totalProgress = coordinates.size();
			stepProgress = Utils.getProgressInterval(totalProgress);
			progress = 0;
			coordinates.forEachKey(new TIntProcedure()
			{
				public boolean execute(int value)
				{
					put(jobs, value);
					return true;
				}
			});
			// Finish all the worker threads by passing in a null job
			for (int i = 0; i < threads.size(); i++)
			{
				put(jobs, -1);
			}

			// Wait for all to finish
			for (int i = 0; i < threads.size(); i++)
			{
				try
				{
					threads.get(i).join();
				}
				catch (InterruptedException e)
				{
					e.printStackTrace();
				}
				actualCoordinates.putAll(workers.get(i).coordinates);
			}
			threads.clear();

			IJ.showProgress(-1);
			IJ.showStatus("");

			setupProgress(total, prefix);
		}

		if (!batchMode)
			IJ.showStatus("Computing results ...");
		final ImageStack stack = imp.getImageStack();

		float background = 0;
		if (spotFilter.isAbsoluteIntensity())
		{
			// To allow the signal factor to be computed we need to lower the image by the background so 
			// that the intensities correspond to the results amplitude.
			// Just assume the background is uniform.
			StandardResultProcedure s = new StandardResultProcedure(results, IntensityUnit.COUNT);
			s.getBIXY();
			background = (float) (Maths.sum(s.background) / results.size());
		}

		// Create a pool of workers
		final int nThreads = Prefs.getThreads();
		BlockingQueue<Integer> jobs = new ArrayBlockingQueue<Integer>(nThreads * 2);
		List<Worker> workers = new LinkedList<Worker>();
		List<Thread> threads = new LinkedList<Thread>();
		for (int i = 0; i < nThreads; i++)
		{
			Worker worker = new Worker(jobs, stack, spotFilter, background);
			Thread t = new Thread(worker);
			workers.add(worker);
			threads.add(t);
			t.start();
		}

		// Fit the frames
		for (int i = 1; i <= stack.getSize(); i++)
		{
			put(jobs, i);
		}
		// Finish all the worker threads by passing in a null job
		for (int i = 0; i < threads.size(); i++)
		{
			put(jobs, -1);
		}

		// Wait for all to finish
		for (int i = 0; i < threads.size(); i++)
		{
			try
			{
				threads.get(i).join();
			}
			catch (InterruptedException e)
			{
				e.printStackTrace();
			}
		}
		threads.clear();

		if (Utils.isInterrupted())
			return null;

		if (!batchMode)
		{
			IJ.showProgress(-1);
			IJ.showStatus("Collecting results ...");
		}

		TIntObjectHashMap<FilterResult> filterResults = new TIntObjectHashMap<FilterResult>();
		time = 0;
		for (Worker w : workers)
		{
			time += w.time;
			filterResults.putAll(w.results);
		}

		// Show a table of the results
		BenchmarkFilterResult filterResult = summariseResults(filterResults, config, spotFilter, relativeDistances,
				batchSummary);

		if (!batchMode)
			IJ.showStatus("");

		return filterResult;
	}

	/**
	 * Add all the true-positives to memory as a new results set
	 * 
	 * @param filterResults
	 */
	void addSpotsToMemory(TIntObjectHashMap<FilterResult> filterResults)
	{
		final MemoryPeakResults results = new MemoryPeakResults();
		results.setName(TITLE + " TP " + id++);
		filterResults.forEachEntry(new TIntObjectProcedure<FilterResult>()
		{
			public boolean execute(int peak, FilterResult filterResult)
			{
				for (ScoredSpot spot : filterResult.spots)
				{
					if (spot.match)
					{
						final float[] params = new float[] { 0, spot.getIntensity(), 0, spot.spot.x, spot.spot.y, 0,
								0 };
						results.add(peak, spot.spot.x, spot.spot.y, spot.getIntensity(), 0d, 0f, params, null);
					}
				}
				return true;
			}
		});
		MemoryPeakResults.addResults(results);
	}

	/**
	 * Histogram the number of negatives preceeding each positive.
	 *
	 * @param filterResult
	 *            the filter result
	 * @return
	 */
	private double[][] histogramFailures(BenchmarkFilterResult filterResult)
	{
		final StoredData data = new StoredData();
		filterResult.filterResults.forEachEntry(new TIntObjectProcedure<FilterResult>()
		{
			public boolean execute(int peak, FilterResult filterResult)
			{
				for (ScoredSpot spot : filterResult.spots)
				{
					if (spot.match)
					{
						data.add(spot.fails);
					}
				}
				return true;
			}
		});

		double[][] h = Maths.cumulativeHistogram(data.getValues(), true);

		filterResult.cumul = h;
		filterResult.stats = data;

		return h;
	}

	private void showFailuresPlot(BenchmarkFilterResult filterResult)
	{
		double[][] h = filterResult.cumul;
		StoredData data = filterResult.stats;

		String xTitle = "Failures";
		final int id = Utils.showHistogram(TITLE, data, xTitle, 1, 0, 0);
		if (Utils.isNewWindow())
			windowOrganiser.add(id);

		String title = TITLE + " " + xTitle + " Cumulative";
		Plot2 plot = new Plot2(title, xTitle, "Frequency");
		double xMin = (data.size() == 0) ? 1 : h[0][0];
		double xMax = (data.size() == 0) ? 1 : h[0][h[0].length - 1] + 1;
		double xPadding = 0.05 * (xMax - xMin);
		plot.setLimits(xMin - xPadding, xMax, 0, 1.05);
		plot.setColor(Color.blue);
		plot.addPoints(h[0], h[1], Plot2.BAR);
		PlotWindow pw = Utils.display(title, plot);
		if (Utils.isNewWindow())
			windowOrganiser.add(pw);
	}

	private void put(BlockingQueue<Integer> jobs, int i)
	{
		try
		{
			jobs.put(i);
		}
		catch (InterruptedException e)
		{
			throw new RuntimeException("Unexpected interruption", e);
		}
	}

	private BenchmarkFilterResult summariseResults(TIntObjectHashMap<FilterResult> filterResults,
			FitEngineConfiguration config, MaximaSpotFilter spotFilter, boolean relativeDistances, boolean batchSummary)
	{
		BenchmarkFilterResult filterResult = new BenchmarkFilterResult(filterResults, config, spotFilter);

		// Note: 
		// Although we can compute the TP/FP score as each additional spot is added
		// using the RankedScoreCalculator this is not applicable to the PeakFit method.
		// The method relies on all spot candidates being present in order to make a
		// decision to fit the candidate as a multiple. So scoring the filter candidates using
		// for example the top 10 may get a better score than if all candidates were scored
		// and the scores accumulated for the top 10, it is not how the algorithm will use the 
		// candidate set. I.e. It does not use the top 10, then top 20 to refine the fit, etc. 
		// (the method is not iterative) .
		// We require an assessment of how a subset of the scored candidates
		// in ranked order contributes to the overall score, i.e. are the candidates ranked
		// in the correct order, those most contributing to the match to the underlying data 
		// should be higher up and those least contributing will be at the end.

		// TODO We could add some smart filtering of candidates before ranking. This would
		// allow assessment of the candidate set handed to PeakFit. E.g. Threshold the image
		// and only use candidates that are in the foreground region.

		double[][] cumul = histogramFailures(filterResult);

		// Create the overall match score
		final double[] total = new double[3];
		final ArrayList<ScoredSpot> allSpots = new ArrayList<BenchmarkSpotFilter.ScoredSpot>();
		filterResults.forEachValue(new TObjectProcedure<FilterResult>()
		{
			public boolean execute(FilterResult result)
			{
				total[0] += result.result.getTP();
				total[1] += result.result.getFP();
				total[2] += result.result.getFN();
				allSpots.addAll(Arrays.asList(result.spots));
				return true;
			}
		});
		double tp = total[0], fp = total[1], fn = total[2];
		FractionClassificationResult allResult = new FractionClassificationResult(tp, fp, 0, fn);
		// The number of actual results
		final double n = (tp + fn);

		StringBuilder sb = new StringBuilder();

		double signal = (simulationParameters.minSignal + simulationParameters.maxSignal) * 0.5;

		// Create the benchmark settings and the fitting settings
		sb.append(imp.getStackSize()).append('\t');
		final int w = lastAnalysisBorder.width;
		final int h = lastAnalysisBorder.height;
		sb.append(w).append('\t');
		sb.append(h).append('\t');
		sb.append(Utils.rounded(n)).append('\t');
		double density = (n / imp.getStackSize()) / (w * h) / (simulationParameters.a * simulationParameters.a / 1e6);
		sb.append(Utils.rounded(density)).append('\t');
		sb.append(Utils.rounded(signal)).append('\t');
		sb.append(Utils.rounded(simulationParameters.s)).append('\t');
		sb.append(Utils.rounded(simulationParameters.a)).append('\t');
		sb.append(Utils.rounded(simulationParameters.depth)).append('\t');
		sb.append(simulationParameters.fixedDepth).append('\t');
		sb.append(Utils.rounded(simulationParameters.gain)).append('\t');
		sb.append(Utils.rounded(simulationParameters.readNoise)).append('\t');
		sb.append(Utils.rounded(simulationParameters.b)).append('\t');
		sb.append(Utils.rounded(simulationParameters.b2)).append('\t');

		// Compute the noise
		double noise = simulationParameters.b2;
		if (simulationParameters.emCCD)
		{
			// The b2 parameter was computed without application of the EM-CCD noise factor of 2.
			//final double b2 = backgroundVariance + readVariance
			//                = simulationParameters.b + readVariance
			// This should be applied only to the background variance.
			final double readVariance = noise - simulationParameters.b;
			noise = simulationParameters.b * 2 + readVariance;
		}

		sb.append(Utils.rounded(signal / Math.sqrt(noise))).append('\t');
		sb.append(Utils.rounded(simulationParameters.s / simulationParameters.a)).append('\t');
		sb.append(config.getDataFilterType()).append('\t');
		//sb.append(spotFilter.getName()).append('\t');
		sb.append(spotFilter.getSearch()).append('\t');
		sb.append(spotFilter.getBorder()).append('\t');
		sb.append(Utils.rounded(spotFilter.getSpread())).append('\t');
		sb.append(config.getDataFilterMethod(0)).append('\t');
		final double param = config.getSmooth(0);
		final double hwhmMin = config.getHWHMMin();
		if (relativeDistances)
		{
			sb.append(Utils.rounded(param * hwhmMin)).append('\t');
			sb.append(Utils.rounded(param)).append('\t');
		}
		else
		{
			sb.append(Utils.rounded(param)).append('\t');
			sb.append(Utils.rounded(param / hwhmMin)).append('\t');
		}
		sb.append(spotFilter.getDescription()).append('\t');
		sb.append(lastAnalysisBorder.x).append('\t');
		sb.append(MATCHING_METHOD[matchingMethod]).append('\t');
		sb.append(Utils.rounded(lowerMatchDistance)).append('\t');
		sb.append(Utils.rounded(matchDistance)).append('\t');
		sb.append(Utils.rounded(lowerSignalFactor)).append('\t');
		sb.append(Utils.rounded(upperSignalFactor));

		resultPrefix = sb.toString();

		// Add the results
		sb.append('\t');

		// Rank the scored spots by intensity
		Collections.sort(allSpots);

		// Produce Recall, Precision, Jaccard for each cut of the spot candidates
		double[] r = new double[allSpots.size() + 1];
		double[] p = new double[r.length];
		double[] j = new double[r.length];
		double[] c = new double[r.length];
		double[] truePositives = new double[r.length];
		double[] falsePositives = new double[r.length];
		double[] intensity = new double[r.length];
		// Note: fn = n - tp
		tp = fp = 0;
		int i = 1;
		p[0] = 1;
		FastCorrelator corr = new FastCorrelator();
		double lastC = 0;
		double[] i1 = new double[r.length];
		double[] i2 = new double[r.length];
		int ci = 0;
		SimpleRegression regression = new SimpleRegression(false);
		for (ScoredSpot s : allSpots)
		{
			if (s.match)
			{
				// Score partial matches as part true-positive and part false-positive.
				// TP can be above 1 if we are allowing multiple matches.
				tp += s.getScore();
				fp += s.antiScore();
				// Just use a rounded intensity for now
				final double spotIntensity = s.getIntensity();
				final long v1 = (long) Math.round(spotIntensity);
				final long v2 = (long) Math.round(s.intensity);
				regression.addData(spotIntensity, s.intensity);
				i1[ci] = spotIntensity;
				i2[ci] = s.intensity;
				ci++;
				corr.add(v1, v2);
				lastC = corr.getCorrelation();
			}
			else
				fp++;
			r[i] = (double) tp / n;
			p[i] = (double) tp / (tp + fp);
			j[i] = (double) tp / (fp + n); // (tp+fp+fn) == (fp+n) since tp+fn=n;
			c[i] = lastC;
			truePositives[i] = tp;
			falsePositives[i] = fp;
			intensity[i] = s.getIntensity();
			i++;
		}
		i1 = Arrays.copyOf(i1, ci);
		i2 = Arrays.copyOf(i2, ci);

		final double slope = regression.getSlope();
		sb.append(Utils.rounded(slope)).append('\t');
		addResult(sb, allResult, c[c.length - 1]);

		// Output the match results when the recall achieves the fraction of the maximum.
		double target = r[r.length - 1];
		if (recallFraction < 100)
			target *= recallFraction / 100.0;
		int fractionIndex = 0;
		while (fractionIndex < r.length && r[fractionIndex] < target)
		{
			fractionIndex++;
		}
		if (fractionIndex == r.length)
			fractionIndex--;
		addResult(sb, new FractionClassificationResult(truePositives[fractionIndex], falsePositives[fractionIndex], 0,
				n - truePositives[fractionIndex]), c[fractionIndex]);

		// Output the match results at the maximum jaccard score
		int maxIndex = 0;
		for (int ii = 1; ii < r.length; ii++)
		{
			if (j[maxIndex] < j[ii])
				maxIndex = ii;
		}
		addResult(sb, new FractionClassificationResult(truePositives[maxIndex], falsePositives[maxIndex], 0,
				n - truePositives[maxIndex]), c[maxIndex]);

		sb.append(Utils.rounded(time / 1e6));

		// Calculate AUC (Average precision == Area Under Precision-Recall curve)
		final double auc = AUCCalculator.auc(p, r);
		// Compute the AUC using the adjusted precision curve
		// which uses the maximum precision for recall >= r
		final double[] maxp = new double[p.length];
		double max = 0;
		for (int k = maxp.length; k-- > 0;)
		{
			if (max < p[k])
				max = p[k];
			maxp[k] = max;
		}
		final double auc2 = AUCCalculator.auc(maxp, r);

		sb.append('\t').append(Utils.rounded(auc));
		sb.append('\t').append(Utils.rounded(auc2));

		// Output the number of fit failures that must be processed to capture fractions of the true positives
		if (cumul[0].length != 0)
		{
			sb.append('\t').append(Utils.rounded(getFailures(cumul, 0.80)));
			sb.append('\t').append(Utils.rounded(getFailures(cumul, 0.90)));
			sb.append('\t').append(Utils.rounded(getFailures(cumul, 0.95)));
			sb.append('\t').append(Utils.rounded(getFailures(cumul, 0.99)));
			sb.append('\t').append(Utils.rounded(cumul[0][cumul[0].length - 1]));
		}
		else
			sb.append("\t\t\t\t\t");

		BufferedTextWindow resultsTable = getTable(batchSummary);
		resultsTable.append(sb.toString());

		// Store results
		filterResult.auc = auc;
		filterResult.auc2 = auc2;
		filterResult.r = r;
		filterResult.p = p;
		filterResult.j = j;
		filterResult.c = c;
		filterResult.maxIndex = maxIndex;
		filterResult.fractionIndex = fractionIndex;
		filterResult.cumul = cumul;
		filterResult.slope = slope;
		filterResult.i1 = i1;
		filterResult.i2 = i2;
		filterResult.intensity = intensity;
		filterResult.relativeDistances = relativeDistances;
		filterResult.time = time;
		return filterResult;
	}

	private boolean isShowOverlay()
	{
		return (showTP || showFP || showFN);
	}

	private void showOverlay(ImagePlus imp, BenchmarkFilterResult filterResult)
	{
		final Overlay o = new Overlay();
		//int tp = 0, fp = 0, fn = 0, nn = 0;
		filterResult.filterResults.forEachValue(new TObjectProcedure<FilterResult>()
		{
			public boolean execute(FilterResult result)
			{
				final int size = result.spots.length;

				float[] tx = null, ty = null, fx = null, fy = null;
				if (showTP)
				{
					tx = new float[size];
					ty = new float[size];
				}
				if (showFP)
				{
					fx = new float[size];
					fy = new float[size];
				}
				int t = 0, f = 0;
				for (ScoredSpot s : result.spots)
				{
					if (s.match)
					{
						if (showTP)
						{
							tx[t] = s.spot.x + 0.5f;
							ty[t++] = s.spot.y + 0.5f;
						}
					}
					else
					{
						if (showFP)
						{
							fx[f] = s.spot.x + 0.5f;
							fy[f++] = s.spot.y + 0.5f;
						}
					}
				}
				//tp += t;
				//fp += f;
				if (showTP)
					SpotFinderPreview.addRoi(result.frame, o, tx, ty, t, Color.green);
				if (showFP)
					SpotFinderPreview.addRoi(result.frame, o, fx, fy, f, Color.red);
				if (showFN)
				{
					// We need the FN ...
					final PSFSpot[] actual = result.actual;
					final boolean[] actualAssignment = result.actualAssignment;
					//nn += actual.length;
					final float[] nx = new float[actual.length];
					final float[] ny = new float[actual.length];
					int n = 0;
					for (int i = 0; i < actual.length; i++)
					{
						if (!actualAssignment[i])
						{
							nx[n] = actual[i].getX();
							ny[n++] = actual[i].getY();
						}
					}
					//fn += n;
					SpotFinderPreview.addRoi(result.frame, o, nx, ny, n, Color.yellow);
				}
				return true;
			}
		});

		//System.out.printf("TP=%d, FP=%d, FN=%d, N=%d (%d) %d\n", tp, fp, fn, tp + fn, results.size(), nn);

		imp.setOverlay(o);
	}

	private void showPlot(BenchmarkFilterResult filterResult)
	{
		double[] r = filterResult.r;
		double[] p = filterResult.p;
		double[] j = filterResult.j;
		double[] c = filterResult.c;
		int fractionIndex = filterResult.fractionIndex;
		int maxIndex = filterResult.maxIndex;
		double auc = filterResult.auc;
		double auc2 = filterResult.auc2;
		double slope = filterResult.slope;
		double[] i1 = filterResult.i1;
		double[] i2 = filterResult.i2;
		double[] intensity = filterResult.intensity;

		double[] rank = new double[intensity.length];
		final double topIntensity = (intensity.length == 1) ? 0 : intensity[1];
		for (int i = 1; i < rank.length; i++)
		{
			if (rankByIntensity)
				rank[i] = topIntensity - intensity[i];
			else
				rank[i] = i;
		}

		String title = TITLE + " Performance";
		Plot2 plot = new Plot2(title, (rankByIntensity) ? "Relative Intensity" : "Spot Rank", "");
		double[] limits = Maths.limits(rank);
		plot.setLimits(limits[0], limits[1], 0, 1.05);
		plot.setColor(Color.blue);
		plot.addPoints(rank, p, Plot2.LINE);
		//plot.addPoints(rank, maxp, Plot2.DOT);
		plot.setColor(Color.red);
		plot.addPoints(rank, r, Plot2.LINE);
		plot.setColor(Color.black);
		plot.addPoints(rank, j, Plot2.LINE);
		// Plot correlation - update the scale to be 0-1?
		plot.setColor(Color.yellow);
		plot.addPoints(rank, c, Plot2.LINE);
		plot.setColor(Color.magenta);
		plot.drawLine(rank[fractionIndex], 0, rank[fractionIndex],
				Maths.max(p[fractionIndex], r[fractionIndex], j[fractionIndex], c[fractionIndex]));
		plot.setColor(Color.pink);
		plot.drawLine(rank[maxIndex], 0, rank[maxIndex], Maths.max(p[maxIndex], r[maxIndex], j[maxIndex], c[maxIndex]));
		plot.setColor(Color.black);
		plot.addLabel(0, 0, "Precision=Blue, Recall=Red, Jaccard=Black, Correlation=Yellow");

		PlotWindow pw = Utils.display(title, plot);
		if (Utils.isNewWindow())
			windowOrganiser.add(pw);

		title = TITLE + " Precision-Recall";
		plot = new Plot2(title, "Recall", "Precision");
		plot.setLimits(0, 1, 0, 1.05);
		plot.setColor(Color.red);
		plot.addPoints(r, p, Plot2.LINE);
		//plot.setColor(Color.magenta);
		//plot.addPoints(r, maxp, Plot2.LINE);
		plot.drawLine(r[r.length - 1], p[r.length - 1], r[r.length - 1], 0);
		plot.setColor(Color.black);
		plot.addLabel(0, 0, "AUC = " + Utils.rounded(auc) + ", AUC2 = " + Utils.rounded(auc2));
		PlotWindow pw2 = Utils.display(title, plot);
		if (Utils.isNewWindow())
			windowOrganiser.add(pw2);

		title = TITLE + " Intensity";
		plot = new Plot2(title, "Candidate", "Spot");
		double[] limits1 = Maths.limits(i1);
		double[] limits2 = Maths.limits(i2);
		plot.setLimits(limits1[0], limits1[1], limits2[0], limits2[1]);
		plot.addLabel(0, 0,
				String.format("Correlation=%s; Slope=%s", Utils.rounded(c[c.length - 1]), Utils.rounded(slope)));
		plot.setColor(Color.red);
		plot.addPoints(i1, i2, Plot.DOT);
		if (slope > 1)
			plot.drawLine(limits1[0], limits1[0] * slope, limits1[1], limits1[1] * slope);
		else
			plot.drawLine(limits2[0] / slope, limits2[0], limits2[1] / slope, limits2[1]);
		PlotWindow pw3 = Utils.display(title, plot);
		if (Utils.isNewWindow())
			windowOrganiser.add(pw3);
	}

	private double getFailures(double[][] cumul, double fraction)
	{
		int i = 0;
		final double[] sum = cumul[1];
		while (sum[i] < fraction && i < sum.length)
			i++;
		return i;
	}

	private void addResult(StringBuilder sb, FractionClassificationResult matchResult, double c)
	{
		addCount(sb, matchResult.getTP());
		addCount(sb, matchResult.getFP());
		sb.append(Utils.rounded(matchResult.getRecall())).append('\t');
		sb.append(Utils.rounded(matchResult.getPrecision())).append('\t');
		sb.append(Utils.rounded(matchResult.getJaccard())).append('\t');
		sb.append(Utils.rounded(c)).append('\t');
	}

	private static void addCount(StringBuilder sb, double value)
	{
		// Check if the double holds an integer count
		if ((int) value == value)
		{
			sb.append((int) value);
		}
		else
		{
			// Otherwise add the counts using at least 2 dp
			if (value > 100)
				sb.append(IJ.d2s(value));
			else
				sb.append(Utils.rounded(value));
		}
		sb.append('\t');
	}

	private BufferedTextWindow getTable(boolean batchSummary)
	{
		if (batchSummary)
		{
			if (batchSummaryTable == null || !batchSummaryTable.isVisible())
			{
				TextWindow table = new TextWindow(TITLE + " Batch", createHeader(), "", 1000, 300);
				table.setVisible(true);
				batchSummaryTable = new BufferedTextWindow(table);
			}
			return batchSummaryTable;
		}
		else
		{
			if (summaryTable == null || !summaryTable.isVisible())
			{
				TextWindow table = new TextWindow(TITLE, createHeader(), "", 1000, 300);
				table.setVisible(true);
				summaryTable = new BufferedTextWindow(table);
			}
			return summaryTable;
		}
	}

	private String createHeader()
	{
		StringBuilder sb = new StringBuilder(
				"Frames\tW\tH\tMolecules\tDensity (um^-2)\tN\ts (nm)\ta (nm)\tDepth (nm)\tFixed\tGain\tReadNoise (ADUs)\tB (photons)\tb2 (photons)\tSNR\ts (px)\t");
		sb.append(
				"Type\tSearch\tBorder\tWidth\tFilter\tAbs.Param\tRel.Param\tDescription\tA.Border\tMatching\tlower d\td\tlower sf\tsf");
		tablePrefix = sb.toString();
		sb.append("\tSlope\t");
		sb.append("TP\tFP\tRecall\tPrecision\tJaccard\tR\t");
		sb.append("TP\tFP\tRecall\tPrecision\tJaccard\tR\t");
		sb.append("TP\tFP\tRecall\tPrecision\tJaccard\tR\t");
		sb.append("Time (ms)\t");
		sb.append("AUC\tAUC2\t");
		sb.append("Fail80\tFail90\tFail95\tFail99\tFail100");
		return sb.toString();
	}

	/**
	 * Updates the given configuration using the latest settings used in benchmarking.
	 *
	 * @param pConfig
	 *            the configuration
	 * @return true, if successful
	 */
	public static boolean updateConfiguration(FitEngineConfiguration pConfig)
	{
		if (filterResult == null)
			return false;

		double scaleSearch = 1;
		double scaleSmooth = 1;
		if (!filterResult.relativeDistances)
		{
			// Distance were absolute. Convert using the HWHM so they are relative
			// to the configured fitting width.
			scaleSearch = 1 / pConfig.getHWHMMax();
			scaleSmooth = 1 / pConfig.getHWHMMin();
		}

		FitEngineConfiguration config = filterResult.config;

		pConfig.setDataFilterType(config.getDataFilterType());
		final int nFilters = config.getNumberOfFilters();
		for (int n = 0; n < nFilters; n++)
		{
			pConfig.setDataFilter(config.getDataFilterMethod(n), config.getSmooth(n) * scaleSmooth, false, n);
		}
		pConfig.setSearch(config.getSearch() * scaleSearch);
		pConfig.setBorder(config.getBorder() * scaleSearch);

		return true;
	}
}
