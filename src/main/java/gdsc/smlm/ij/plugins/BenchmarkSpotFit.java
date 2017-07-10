package gdsc.smlm.ij.plugins;

import java.awt.Checkbox;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import java.awt.Color;
import java.awt.Rectangle;
import java.awt.TextArea;
import java.awt.TextField;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import gdsc.core.ij.Utils;
import gdsc.core.match.Assignment;
import gdsc.core.match.AssignmentComparator;
import gdsc.core.match.BasePoint;
import gdsc.core.match.Coordinate;
import gdsc.core.match.FractionClassificationResult;
import gdsc.core.match.FractionalAssignment;
import gdsc.core.match.ImmutableFractionalAssignment;
import gdsc.core.match.PointPair;
import gdsc.core.utils.Correlator;
import gdsc.core.utils.FastCorrelator;
import gdsc.core.utils.Maths;
import gdsc.core.utils.RampedScore;
import gdsc.core.utils.Settings;
import gdsc.core.utils.Sort;
import gdsc.core.utils.StoredDataStatistics;
import gdsc.core.utils.XmlUtils;
import gdsc.smlm.data.config.FitProtos.FitSolver;
import gdsc.smlm.data.config.FitProtos.NoiseEstimatorMethod;
import gdsc.smlm.data.config.GUIProtos.GUIFilterSettings;
import gdsc.smlm.engine.FitConfiguration;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.engine.FitParameters;
import gdsc.smlm.engine.FitParameters.FitTask;
import gdsc.smlm.engine.FitWorker;
import gdsc.smlm.engine.ParameterisedFitJob;
import gdsc.smlm.filters.MaximaSpotFilter;
import gdsc.smlm.filters.Spot;
import gdsc.smlm.fitting.FitResult;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.ij.plugins.BenchmarkSpotFilter.FilterResult;
import gdsc.smlm.ij.plugins.BenchmarkSpotFilter.ScoredSpot;
import gdsc.smlm.ij.plugins.ResultsMatchCalculator.PeakResultPoint;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.ImageConverter;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResults;
import gdsc.smlm.results.filter.BasePreprocessedPeakResult;
import gdsc.smlm.results.filter.CoordinateStoreFactory;
import gdsc.smlm.results.filter.DirectFilter;
import gdsc.smlm.results.filter.EShiftFilter;
import gdsc.smlm.results.filter.Filter;
import gdsc.smlm.results.filter.FilterSet;
import gdsc.smlm.results.filter.MultiFilter2;
import gdsc.smlm.results.filter.MultiPathFilter;
import gdsc.smlm.results.filter.MultiPathFilter.FractionScoreStore;
import gdsc.smlm.results.filter.MultiPathFitResult;
import gdsc.smlm.results.filter.MultiPathFitResults;
import gdsc.smlm.results.filter.ParameterType;
import gdsc.smlm.results.filter.PeakFractionalAssignment;
import gdsc.smlm.results.filter.PrecisionFilter;
import gdsc.smlm.results.filter.PreprocessedPeakResult;
import gdsc.smlm.results.filter.ResultAssignment;
import gdsc.smlm.results.filter.SNRFilter;
import gdsc.smlm.results.filter.ShiftFilter;
import gdsc.smlm.results.filter.SignalFilter;
import gdsc.smlm.results.filter.WidthFilter;
import gdsc.smlm.results.filter.WidthFilter2;
import gdsc.smlm.results.filter.XStreamWrapper;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.procedure.TIntObjectProcedure;
import gnu.trove.procedure.TIntProcedure;
import gnu.trove.set.hash.TIntHashSet;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.ExtendedGenericDialog;
import ij.gui.Plot;
import ij.gui.Plot2;
import ij.gui.PlotWindow;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;
import ij.text.TextWindow;

/**
 * Fits all the candidate spots identified by the benchmark spot filter plugin.
 */
public class BenchmarkSpotFit implements PlugIn, ItemListener
{
	static final String TITLE = "Fit Spot Data";

	// Used to try and guess the range for filtering the results
	private enum LowerLimit
	{
		ZERO(false), ONE_PERCENT(false), MAX_NEGATIVE_CUMUL_DELTA(true), HALF_MAX_JACCARD_VALUE(false, true);

		final boolean requiresDelta;
		final boolean requiresJaccard;

		LowerLimit(boolean requiresDelta)
		{
			this(requiresDelta, false);
		}

		LowerLimit(boolean requiresDelta, boolean requiresJaccard)
		{
			this.requiresDelta = requiresDelta;
			this.requiresJaccard = requiresJaccard;
		}

		public boolean requiresDeltaHistogram()
		{
			return requiresDelta;
		}
	}

	private enum UpperLimit
	{
		ZERO(false), MAX_POSITIVE_CUMUL_DELTA(true), NINETY_NINE_PERCENT(false), NINETY_NINE_NINE_PERCENT(
				false), MAX_JACCARD2(false, true);

		final boolean requiresDelta;
		final boolean requiresJaccard;

		UpperLimit(boolean requiresDelta)
		{
			this(requiresDelta, false);
		}

		UpperLimit(boolean requiresDelta, boolean requiresJaccard)
		{
			this.requiresDelta = requiresDelta;
			this.requiresJaccard = requiresJaccard;
		}

		public boolean requiresDeltaHistogram()
		{
			return requiresDelta;
		}
	}

	private class FilterCriteria
	{
		final ParameterType type;
		final String name;
		final LowerLimit lower;
		final UpperLimit upper;
		final int minBinWidth;
		final boolean restrictRange;
		final boolean requireLabel;

		public FilterCriteria(ParameterType type, LowerLimit lower, UpperLimit upper)
		{
			this(type, type.toString(), lower, upper, 0, true, true);
		}

		public FilterCriteria(ParameterType type, String name, LowerLimit lower, UpperLimit upper, int minBinWidth,
				boolean restrictRange, boolean requireLabel)
		{
			this.type = type;
			this.name = name;
			this.lower = lower;
			this.upper = upper;
			this.minBinWidth = minBinWidth;
			this.restrictRange = restrictRange;
			this.requireLabel = requireLabel;
		}
	}

	private static FilterCriteria[] filterCriteria = null;
	private static final int FILTER_SIGNAL = 0;
	private static final int FILTER_SNR = 1;
	private static final int FILTER_MIN_WIDTH = 2;
	private static final int FILTER_MAX_WIDTH = 3;
	private static final int FILTER_SHIFT = 4;
	private static final int FILTER_ESHIFT = 5;
	private static final int FILTER_PRECISION = 6;
	private static final int FILTER_ITERATIONS = 7;
	private static final int FILTER_EVALUATIONS = 8;

	private FilterCriteria[] createFilterCriteria()
	{
		if (filterCriteria == null)
		{
			filterCriteria = new FilterCriteria[9];
			int i = 0;
			//@formatter:off
			filterCriteria[i++] = new FilterCriteria(ParameterType.SIGNAL,     LowerLimit.ONE_PERCENT, UpperLimit.MAX_POSITIVE_CUMUL_DELTA);
			filterCriteria[i++] = new FilterCriteria(ParameterType.SNR,        LowerLimit.ONE_PERCENT, UpperLimit.MAX_POSITIVE_CUMUL_DELTA);
			filterCriteria[i++] = new FilterCriteria(ParameterType.MIN_WIDTH,  LowerLimit.ONE_PERCENT, UpperLimit.ZERO);
			filterCriteria[i++] = new FilterCriteria(ParameterType.MAX_WIDTH,  LowerLimit.ZERO,        UpperLimit.NINETY_NINE_PERCENT);
			filterCriteria[i++] = new FilterCriteria(ParameterType.SHIFT,      LowerLimit.MAX_NEGATIVE_CUMUL_DELTA, UpperLimit.NINETY_NINE_PERCENT);
			filterCriteria[i++] = new FilterCriteria(ParameterType.ESHIFT,     LowerLimit.MAX_NEGATIVE_CUMUL_DELTA, UpperLimit.NINETY_NINE_PERCENT);
			// Precision has enough discrimination power to be able to use the jaccard score
			filterCriteria[i++] = new FilterCriteria(ParameterType.PRECISION,  LowerLimit.HALF_MAX_JACCARD_VALUE, UpperLimit.MAX_JACCARD2);
			// These are not filters but are used for stats analysis
			filterCriteria[i++] = new FilterCriteria(null, "Iterations", LowerLimit.ONE_PERCENT, UpperLimit.NINETY_NINE_NINE_PERCENT, 1, false, false);
			filterCriteria[i++] = new FilterCriteria(null, "Evaluations",LowerLimit.ONE_PERCENT, UpperLimit.NINETY_NINE_NINE_PERCENT, 1, false, false);
			//@formatter:on

			// Some parameter types may not be for DirectFilters so ignore this check...
			// We just have to be sure that we support all the types produced by any DirectFilter.

			//			// Do a check to ensure we have all the parameter types in the correct order.
			//			// This is needed so that all possible filters can be processed.
			//			ParameterType[] types = ParameterType.values();
			//			NEXT_TYPE: for (int k = 0; k < types.length; k++)
			//			{
			//				for (int j = 0; j < filterCriteria.length; j++)
			//					if (filterCriteria[j].type == types[k])
			//						continue NEXT_TYPE;
			//				filterCriteria = null;
			//				throw new RuntimeException("Missing parameter type: " + types[k].toString());
			//			}
		}
		return filterCriteria;
	}

	private static double[] min;
	private static double[] max;

	/**
	 * Gets the min value of the most recent fit data for the given parameter name.
	 *
	 * @param name
	 *            the name
	 * @return the min
	 */
	public static double getMin(ParameterType type)
	{
		return getValue(type, min, 0);
	}

	/**
	 * Gets the max value of the most recent fit data for the given parameter name.
	 *
	 * @param name
	 *            the name
	 * @return the max
	 */
	public static double getMax(ParameterType type)
	{
		return getValue(type, max, Double.MAX_VALUE);
	}

	private static double getValue(ParameterType type, double[] array, double defaultValue)
	{
		if (type == null || array == null || filterCriteria == null)
			return defaultValue;

		// Assume these are roughly the same
		if (type == ParameterType.PRECISION2)
			type = ParameterType.PRECISION;

		for (int j = 0; j < filterCriteria.length; j++)
			if (filterCriteria[j].type == type)
				return array[j];

		// All other types will have a default value
		return defaultValue;
	}

	static FitConfiguration fitConfig;
	static FitEngineConfiguration config;
	static MultiPathFilter multiFilter;
	private static final MultiPathFilter defaultMultiFilter;
	private static final double[] defaultParameters;
	private static MultiFilter2 minimalFilter;;

	static
	{
		config = new FitEngineConfiguration();
		fitConfig = config.getFitConfiguration();
		// Set some default fit settings here ...
		fitConfig.setDisableSimpleFilter(false);
		fitConfig.setMinPhotons(1); // Do not allow negative photons 
		fitConfig.setCoordinateShiftFactor(0);
		fitConfig.setPrecisionThreshold(0);
		fitConfig.setMinWidthFactor(0);
		fitConfig.setWidthFactor(0);

		fitConfig.setBackgroundFitting(true);
		fitConfig.setNoise(0);
		config.setNoiseMethod(NoiseEstimatorMethod.QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES);

		// Use bounded so that we can fit the neighbours
		config.setIncludeNeighbours(true);
		fitConfig.setFitSolver(FitSolver.LVM_LSE);

		// Add a filter to use for storing the slice results:
		// Use the standard configuration to ensure sensible fits are stored as the current slice results.
		FitConfiguration tmp = new FitConfiguration();
		tmp.setPrecisionUsingBackground(true); // So we get a MultiFilter2 to match the minimal filter

		// Add a minimum filter to use for storing estimates
		minimalFilter = FitWorker.createMinimalFilter();

		final DirectFilter primaryFilter = tmp.getDefaultSmartFilter();

		// We might as well use the doublet fits given we will compute them.
		final double residualsThreshold = 0.4;
		multiFilter = new MultiPathFilter(primaryFilter, minimalFilter, residualsThreshold);
		defaultMultiFilter = multiFilter;
		defaultParameters = createParameters();
	}

	private static double fractionPositives = 100;
	private static double fractionNegativesAfterAllPositives = 50;
	private static int negativesAfterAllPositives = 10;
	private static double distance = 1.5;
	private static double lowerDistance = 1.5;
	// Allow other plugins to access these
	static double signalFactor = 2;
	static double lowerSignalFactor = 1;

	private static boolean useBenchmarkSettings = false;
	static boolean computeDoublets = true; //config.getResidualsThreshold() < 1;
	private static boolean showFilterScoreHistograms = false;
	private static boolean saveFilterRange = true;
	private static boolean showCorrelation = false;
	private static boolean rankByIntensity = false;

	private TextArea taFilterXml;
	private TextField textFailLimit;
	private Checkbox cbIncludeNeighbours;
	private TextField textNeighbourHeight;
	private Checkbox cbComputeDoublets;

	private boolean extraOptions = false;
	// Flag used when being called by another plugin to suppress dialogs
	private boolean silent = false;
	// Flag used when being called by another plugin to idicate success
	boolean finished;

	private static TextWindow summaryTable = null;

	private ImagePlus imp;
	private MemoryPeakResults results;
	private CreateData.SimulationParameters simulationParameters;
	private MaximaSpotFilter spotFilter;

	private static TIntObjectHashMap<ArrayList<Coordinate>> actualCoordinates = null;
	private static TIntObjectHashMap<FilterCandidates> filterCandidates;
	private static double fP, fN;
	private static int nP, nN;

	static int lastId = -1, lastFilterId = -1;
	private static Settings lastSettings = null;

	// Allow other plugins to access the results
	static int fitResultsId = 0;
	static TIntObjectHashMap<FilterCandidates> fitResults;
	static double distanceInPixels;
	static double lowerDistanceInPixels;
	static double candidateTN, candidateFN;
	// Allow access to the time
	static StopWatch stopWatch;

	public static String tablePrefix, resultPrefix;

	public class MultiPathPoint extends BasePoint
	{
		final static int SPOT = -1;
		final static int SINGLE = 0;
		final static int MULTI = 1;
		final static int DOUBLET = 2;
		final static int MULTIDOUBLET = 3;

		final PreprocessedPeakResult result;
		final int id, type, i;

		public MultiPathPoint(PreprocessedPeakResult result, int id, int type, int i)
		{
			super(result.getX(), result.getY());
			this.result = result;
			this.id = id;
			this.type = type;
			this.i = i;
		}

		public MultiPathPoint(float x, float y, int id, int type, int i)
		{
			super(x, y);
			this.result = null;
			this.id = id;
			this.type = type;
			this.i = i;
		}
	}

	/**
	 * Store details of spot candidates that match actual spots
	 */
	public abstract class SpotMatch
	{
		/**
		 * The index for the spot candidate
		 */
		final int i;
		/**
		 * The distance to the spot
		 */
		final double d;
		/**
		 * The depth of the actual spot
		 */
		final double z;
		/**
		 * The score
		 */
		double score;

		public SpotMatch(int i, double d, double z)
		{
			this.i = i;
			this.d = d;
			this.z = z;
		}

		/**
		 * @return True if the spot candidate was successfully fitted
		 */
		public abstract boolean isFitResult();

		/**
		 * Return a score for the difference between the fitted and actual signal. Zero is no difference. Negative is
		 * the fitted is below the actual. Positive means the fitted is above the actual.
		 * 
		 * @return The factor difference between the successfully fitted signal and the actual signal.
		 */
		public abstract double getSignalFactor();

		/**
		 * Return a score for the difference between the fitted and actual signal. Zero is no difference.
		 * Positive means the fitted is difference from the actual.
		 * 
		 * @return The factor difference between the successfully fitted signal and the actual signal.
		 */
		public double getAbsoluteSignalFactor()
		{
			return Math.abs(getSignalFactor());
		}
	}

	/**
	 * Store details of a fitted spot candidate that matches an actual spot
	 */
	public class FitMatch extends SpotMatch
	{
		final MultiPathPoint point;
		final double predictedSignal, actualSignal;
		final double sf;

		public FitMatch(MultiPathPoint point, double d, double z, double predictedSignal, double actualSignal)
		{
			super(point.i, d, z);
			this.point = point;
			this.predictedSignal = predictedSignal;
			this.actualSignal = actualSignal;
			this.sf = BenchmarkSpotFit.getSignalFactor(predictedSignal, actualSignal);
		}

		@Override
		public boolean isFitResult()
		{
			return true;
		}

		@Override
		public double getSignalFactor()
		{
			return sf;
		}
	}

	static double getSignalFactor(double predictedSignal, double actualSignal)
	{
		final double rsf = predictedSignal / actualSignal;
		// The relative signal factor is 1 for a perfect fit, less than 1 for below and above 1 for above.
		// Reset the signal factor from 0
		double sf = (rsf < 1) ? 1 - 1 / rsf : rsf - 1;
		return sf;
	}

	/**
	 * Store details of a spot candidate that matches an actual spot
	 */
	public class CandidateMatch extends SpotMatch
	{
		public CandidateMatch(int i, double d, double z)
		{
			super(i, d, z);
		}

		@Override
		public boolean isFitResult()
		{
			return false;
		}

		@Override
		public double getSignalFactor()
		{
			return 0;
		}
	}

	public class FilterCandidates implements Cloneable
	{
		// Integer counts of positives (matches) and negatives
		final int p, n;
		// Double sums of the fractions match score and antiscore 
		final double np, nn;
		final ScoredSpot[] spots;
		final int maxCandidate;
		double tp, fp, tn, fn;
		MultiPathFitResult[] fitResult;
		float noise;

		/** Store the z-position of the actual spots for later analysis. Size is the number of actual spots */
		double[] zPosition;

		/**
		 * Store details about the actual spots that were matched by spot candidates or fitted spot candidates
		 */
		SpotMatch[] match;

		public FilterCandidates(int p, int n, double np, double nn, ScoredSpot[] spots, int maxCandidate)
		{
			this.p = p;
			this.n = n;
			this.np = np;
			this.nn = nn;
			this.spots = spots;
			this.maxCandidate = maxCandidate;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Object#clone()
		 */
		public FilterCandidates clone()
		{
			try
			{
				return (FilterCandidates) super.clone();
			}
			catch (CloneNotSupportedException e)
			{
				return null;
			}
		}
	}

	private class Ranking implements Comparable<Ranking>
	{
		final double value;
		final int index;

		Ranking(double value, int index)
		{
			this.value = value;
			this.index = index;
		}

		public int compareTo(Ranking that)
		{
			return Double.compare(this.value, that.value);
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
		final FitWorker fitWorker;
		final TIntObjectHashMap<ArrayList<Coordinate>> actualCoordinates;
		final TIntObjectHashMap<FilterCandidates> filterCandidates;
		final TIntObjectHashMap<FilterCandidates> results;
		final Rectangle bounds;
		final MultiPathFilter multiFilter;

		float[] data = null;
		List<PointPair> matches = new ArrayList<PointPair>();

		public Worker(BlockingQueue<Integer> jobs, ImageStack stack,
				TIntObjectHashMap<ArrayList<Coordinate>> actualCoordinates,
				TIntObjectHashMap<FilterCandidates> filterCandidates, PeakResults peakResults)
		{
			this.jobs = jobs;
			this.stack = stack;
			this.fitWorker = new FitWorker(config.clone(), peakResults, null);

			final int fitting = config.getRelativeFitting();
			fitWorker.setSearchParameters(spotFilter.clone(), fitting);

			this.actualCoordinates = actualCoordinates;
			this.filterCandidates = filterCandidates;
			this.results = new TIntObjectHashMap<FilterCandidates>();
			bounds = new Rectangle(0, 0, stack.getWidth(), stack.getHeight());
			// Instance copy
			multiFilter = BenchmarkSpotFit.multiFilter.clone();
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
						// Only run jobs when not finished. This allows the queue to be emptied.
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

			FilterCandidates candidates = filterCandidates.get(frame);
			int totalCandidates = candidates.spots.length;
			final MultiPathFitResult[] fitResult = new MultiPathFitResult[totalCandidates];

			// Fit the candidates and store the results
			final FitParameters parameters = new FitParameters();
			final Spot[] spots = new Spot[candidates.spots.length];
			for (int i = 0; i < spots.length; i++)
			{
				spots[i] = candidates.spots[i].spot;
				// Debug candidates...
				//if (frame == 5)
				//	System.out.printf("Fit %d [%d,%d = %.1f]\n", i, spots[i].x, spots[i].y, spots[i].intensity);
			}
			parameters.spots = spots;
			parameters.maxCandidate = candidates.maxCandidate;
			parameters.fitTask = FitTask.BENCHMARKING;
			parameters.benchmarkFilter = multiFilter;

			final ParameterisedFitJob job = new ParameterisedFitJob(parameters, frame, data, bounds);
			fitWorker.run(job); // Results will be stored in the fit job 

			for (int i = 0; i < totalCandidates; i++)
			{
				fitResult[i] = job.getMultiPathFitResult(i);
			}

			// Compute the matches of the fitted spots to the simulated positions.
			// We will match all fitting results so providing the upper limit for the match score after filtering.
			final Coordinate[] actual = ResultsMatchCalculator.getCoordinates(actualCoordinates, frame);
			final double[] zPosition = new double[actual.length];
			SpotMatch[] match = new SpotMatch[actual.length];
			int matchCount = 0;
			final RampedScore rampedScore = new RampedScore(lowerDistanceInPixels, distanceInPixels);
			final RampedScore signalScore = (signalFactor > 0) ? new RampedScore(lowerSignalFactor, signalFactor)
					: null;
			if (actual.length > 0)
			{
				// Build a list of the coordinates z-depth using the PeakResultPoint
				for (int i = 0; i < actual.length; i++)
				{
					PeakResultPoint p = (PeakResultPoint) actual[i];
					zPosition[i] = p.peakResult.error;
				}

				// Allow for doublets the predicted array
				final ArrayList<MultiPathPoint> predicted = new ArrayList<MultiPathPoint>(totalCandidates * 2);
				matches.clear();

				for (int i = 0; i < totalCandidates; i++)
				{
					// Use all the results. 
					// Store results using the candidate Id. The best match for each Id will be chosen.
					final int size = predicted.size();

					// Allow all new results to be processed as multi-fit can return more than 1 new result.
					add(predicted, fitResult[i].getSingleFitResult(), MultiPathPoint.SINGLE, i);
					add(predicted, fitResult[i].getMultiFitResult(), MultiPathPoint.MULTI, i);
					add(predicted, fitResult[i].getDoubletFitResult(), MultiPathPoint.DOUBLET, i);
					add(predicted, fitResult[i].getMultiDoubletFitResult(), MultiPathPoint.MULTI, i);
					if (size == predicted.size())
					{
						// Use the candidate position instead
						predicted.add(
								new MultiPathPoint(spots[i].x + 0.5f, spots[i].y + 0.5f, i, MultiPathPoint.SPOT, i));
					}
				}
				// If we made any fits then score them
				if (!predicted.isEmpty())
				{
					// Match fit results/candidates with their closest actual spot

					// TODO - Is using the closest match the best way to do this for high density data?
					// Perhaps we should pair up the closest matches using the signal factor as well.

					final double matchDistance = distanceInPixels * distanceInPixels;
					ArrayList<Assignment> assignments = new ArrayList<Assignment>();

					// Match all the fit results to spots. We want to match all fit results to actual spots.
					// All remaining candidate spots can then be matched to any remaining actual spots.
					ResultAssignment[][] resultAssignments = new ResultAssignment[predicted.size()][];
					int[] resultAssignmentsSize = new int[resultAssignments.length];
					for (int ii = 0; ii < actual.length; ii++)
					{
						final float x = actual[ii].getX();
						final float y = actual[ii].getY();
						for (int jj = 0; jj < predicted.size(); jj++)
						{
							final double d2 = predicted.get(jj).distance2(x, y);
							if (d2 <= matchDistance)
							{
								// Get the score
								double distance = d2;
								double score = 0;

								if (predicted.get(jj).type != MultiPathPoint.SPOT)
								{
									// Use the signal and ramped distance scoring
									score = rampedScore.score(Math.sqrt(d2));
									if (signalScore != null)
									{
										final PeakResultPoint p3 = (PeakResultPoint) actual[ii];
										double sf = getSignalFactor(predicted.get(jj).result.getSignal(),
												p3.peakResult.getSignal());
										score *= signalScore.score(Math.abs(sf));

										if (score == 0)
											// This doesn't match
											continue;
									}

									// Invert for the ranking (i.e. low is best)
									distance = 1 - score;

									// Ensure a perfect match can still be ranked ... closest first
									if (distance == 0)
										distance -= (matchDistance - d2);

									// Store the assignments for this result for later filter analysis
									if (resultAssignments[jj] == null)
										resultAssignments[jj] = new ResultAssignment[5];
									else if (resultAssignments[jj].length == resultAssignmentsSize[jj])
										resultAssignments[jj] = Arrays.copyOf(resultAssignments[jj],
												resultAssignmentsSize[jj] * 2);
									resultAssignments[jj][resultAssignmentsSize[jj]++] = new ResultAssignment(ii,
											distance, score);
								}
								else
								{
									// This is not a fit. Ensure that remaining candidate spots are assigned
									// after any fit results.
									distance += matchDistance + 1;
								}

								// Store the index of the predicted point in the score
								assignments
										.add(new ImmutableFractionalAssignment(ii, predicted.get(jj).id, distance, jj));
							}
						}
					}

					// Set assignments
					for (int jj = 0; jj < predicted.size(); jj++)
					{
						if (predicted.get(jj).type == MultiPathPoint.SPOT)
							continue;
						if (resultAssignmentsSize[jj] != 0)
						{
							BasePreprocessedPeakResult result = (BasePreprocessedPeakResult) predicted.get(jj).result;
							result.setAssignments(Arrays.copyOf(resultAssignments[jj], resultAssignmentsSize[jj]));
						}
					}

					AssignmentComparator.sort(assignments);

					final boolean[] actualAssignment = new boolean[actual.length];
					// We use the candidate Id as the id
					final boolean[] predictedAssignment = new boolean[fitResult.length];

					for (Assignment assignment : assignments)
					{
						if (!actualAssignment[assignment.getTargetId()])
						{
							if (!predictedAssignment[assignment.getPredictedId()])
							{
								actualAssignment[assignment.getTargetId()] = true;
								predictedAssignment[assignment.getPredictedId()] = true;

								final PeakResultPoint p3 = (PeakResultPoint) actual[assignment.getTargetId()];
								int jj = (int) ((ImmutableFractionalAssignment) assignment).getScore();
								final MultiPathPoint point = predicted.get(jj);

								final double d = point.distanceXY(p3);

								if (point.type != MultiPathPoint.SPOT)
								{
									// This is a fitted candidate

									final double a = p3.peakResult.getSignal(); // Should be in photons
									final double p = point.result.getPhotons();

									match[matchCount++] = new FitMatch(point, d, p3.peakResult.error, p, a);
								}
								else
								{
									// This is a candidate that could not be fitted
									match[matchCount++] = new CandidateMatch(point.i, d, p3.peakResult.error);
								}
							}
						}
					}
				}
			}
			match = Arrays.copyOf(match, matchCount);

			// Mark the results 
			double tp = 0;
			double fp = 0;
			double tn = 0;
			double fn = 0;

			for (int i = 0; i < match.length; i++)
			{
				// Score using just the distance.
				// The signal has been used only to compute the best match when two spots are close. 
				final double s = rampedScore.scoreAndFlatten(match[i].d, 256);
				match[i].score = s;

				if (match[i].isFitResult())
				{
					// This is a fitted result so is a positive
					tp += s;
					fp += 1 - s;
				}
				else
				{
					// This is an unfitted result that matches so is a negative
					fn += s;
					tn += 1 - s;
				}
			}

			// Store the results using a copy of the original (to preserve the candidates for repeat analysis)
			candidates = candidates.clone();
			candidates.tp = tp;
			candidates.fp = fp;
			candidates.tn = tn;
			candidates.fn = fn;
			candidates.fitResult = fitResult;
			candidates.match = match;
			candidates.zPosition = zPosition;
			candidates.noise = fitWorker.getNoise();
			results.put(frame, candidates);
		}

		private void add(ArrayList<MultiPathPoint> predicted,
				gdsc.smlm.results.filter.MultiPathFitResult.FitResult fitResult, int type, int spotId)
		{
			if (fitResult == null || fitResult.status != 0 || fitResult.results == null)
				return;

			for (int i = 0; i < fitResult.results.length; i++)
			{
				if (fitResult.results[i].isNewResult())
					predicted.add(new MultiPathPoint(fitResult.results[i], fitResult.results[i].getCandidateId(), type,
							spotId));
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

		silent = false;
		finished = false;
		if (!initialise())
			return;

		if (!showDialog())
			return;

		run();
		finished = true;
	}

	private boolean initialise()
	{
		simulationParameters = CreateData.simulationParameters;
		if (simulationParameters == null)
		{
			if (!silent)
				IJ.error(TITLE, "No benchmark spot parameters in memory");
			return false;
		}
		imp = CreateData.getImage();
		if (imp == null)
		{
			if (!silent)
				IJ.error(TITLE, "No benchmark image");
			return false;
		}
		results = CreateData.getResults();
		if (results == null)
		{
			if (!silent)
				IJ.error(TITLE, "No benchmark results in memory");
			return false;
		}
		if (BenchmarkSpotFilter.filterResult == null)
		{
			if (!silent)
				IJ.error(TITLE, "No benchmark spot candidates in memory");
			return false;
		}
		if (BenchmarkSpotFilter.filterResult.simulationId != simulationParameters.id)
		{
			if (!silent)
				IJ.error(TITLE, "Update the benchmark spot candidates for the latest simulation");
			return false;
		}
		// This is required to initialise the FitWorker
		spotFilter = BenchmarkSpotFilter.filterResult.spotFilter;
		return true;
	}

	private boolean showDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage(String.format(
				"Fit candidate spots in the benchmark image created by " + CreateData.TITLE +
						" plugin\nand identified by the " + BenchmarkSpotFilter.TITLE +
						" plugin.\nPSF width = %s nm (Square pixel adjustment = %s nm)\n \nConfigure the fitting:",
				Utils.rounded(simulationParameters.s), Utils.rounded(getSa())));

		gd.addSlider("Fraction_positives", 50, 100, fractionPositives);
		gd.addSlider("Fraction_negatives_after_positives", 0, 100, fractionNegativesAfterAllPositives);
		gd.addSlider("Min_negatives_after_positives", 0, 10, negativesAfterAllPositives);
		gd.addSlider("Match_distance", 0.5, 3.5, distance);
		gd.addSlider("Lower_distance", 0, 3.5, lowerDistance);
		gd.addSlider("Match_signal", 0, 3.5, signalFactor);
		gd.addSlider("Lower_signal", 0, 3.5, lowerSignalFactor);

		// Collect options for fitting
		final double sa = getSa();
		fitConfig.setInitialPeakStdDev(Maths.round(sa / simulationParameters.a));
		PeakFit.addPSFOptions(gd, fitConfig);
		gd.addSlider("Fitting_width", 2, 4.5, config.getFitting());
		gd.addChoice("Fit_solver", SettingsManager.getFitSolverNames(), fitConfig.getFitSolver().ordinal());

		gd.addMessage("Multi-path filter (used to pick optimum results during fitting)");

		// Allow loading the best filter for these results
		boolean benchmarkSettingsCheckbox = fitResultsId == BenchmarkFilterAnalysis.lastId;

		// This should always be an opt-in decision. Otherwise the user cannot use the previous settings
		useBenchmarkSettings = false;
		Checkbox cbBenchmark = null;
		if (benchmarkSettingsCheckbox)
			cbBenchmark = gd.addAndGetCheckbox("Benchmark_settings", useBenchmarkSettings);

		gd.addTextAreas(XmlUtils.convertQuotes(multiFilter.toXML()), null, 6, 60);

		textFailLimit = gd.addAndGetNumericField("Fail_limit", config.getFailuresLimit(), 0);
		cbIncludeNeighbours = gd.addAndGetCheckbox("Include_neighbours", config.isIncludeNeighbours());
		gd.addAndGetSlider("Neighbour_height", 0.01, 1, config.getNeighbourHeightThreshold());
		textNeighbourHeight = gd.getLastTextField();
		//gd.addSlider("Residuals_threshold", 0.01, 1, config.getResidualsThreshold());
		cbComputeDoublets = gd.addAndGetCheckbox("Compute_doublets", computeDoublets);
		gd.addNumericField("Duplicate_distance", config.getDuplicateDistance(), 2);
		gd.addCheckbox("Show_score_histograms", showFilterScoreHistograms);
		gd.addCheckbox("Show_correlation", showCorrelation);
		gd.addCheckbox("Plot_rank_by_intensity", rankByIntensity);
		gd.addCheckbox("Save_filter_range", saveFilterRange);

		if (extraOptions)
		{
		}

		// Add a mouse listener to the config file field
		if (benchmarkSettingsCheckbox && Utils.isShowGenericDialog())
		{
			taFilterXml = gd.getTextArea1();
			cbBenchmark.addItemListener(this);

			if (useBenchmarkSettings)
			{
				FitEngineConfiguration tmp = new FitEngineConfiguration();
				FitConfiguration tmpFitConfig = tmp.getFitConfiguration();
				tmpFitConfig.setComputeResiduals(true); // Collect the residuals threshold
				if (BenchmarkFilterAnalysis.updateConfiguration(tmp, false))
				{
					textFailLimit.setText("" + tmp.getFailuresLimit());
					cbIncludeNeighbours.setState(tmp.isIncludeNeighbours());
					textNeighbourHeight.setText(Utils.rounded(tmp.getNeighbourHeightThreshold()));
					cbComputeDoublets.setState(tmp.getResidualsThreshold() < 1);

					final DirectFilter primaryFilter = tmpFitConfig.getSmartFilter();
					final double residualsThreshold = tmp.getResidualsThreshold();
					taFilterXml.setText(new MultiPathFilter(primaryFilter, minimalFilter, residualsThreshold).toXML());
				}
			}
		}

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		fractionPositives = Math.abs(gd.getNextNumber());
		fractionNegativesAfterAllPositives = Math.abs(gd.getNextNumber());
		negativesAfterAllPositives = (int) Math.abs(gd.getNextNumber());
		distance = Math.abs(gd.getNextNumber());
		lowerDistance = Math.abs(gd.getNextNumber());
		signalFactor = Math.abs(gd.getNextNumber());
		lowerSignalFactor = Math.abs(gd.getNextNumber());

		fitConfig.setPSFType(PeakFit.getPSFTypeValues()[gd.getNextChoiceIndex()]);
		config.setFitting(gd.getNextNumber());
		fitConfig.setFitSolver(gd.getNextChoiceIndex());

		boolean myUseBenchmarkSettings = false;
		if (benchmarkSettingsCheckbox)
			//useBenchmarkSettings = 
			myUseBenchmarkSettings = gd.getNextBoolean();

		// Read dialog settings
		String xml = gd.getNextText();
		int failLimit = (int) gd.getNextNumber();
		boolean includeNeighbours = gd.getNextBoolean();
		double neighbourHeightThreshold = gd.getNextNumber();
		boolean myComputeDoublets = gd.getNextBoolean();
		double myDuplicateDistance = gd.getNextNumber();

		gd.collectOptions();

		MultiPathFilter myMultiFilter = null;
		if (myUseBenchmarkSettings && !Utils.isShowGenericDialog())
		{
			// Only copy the benchmark settings if not interactive
			FitEngineConfiguration tmp = new FitEngineConfiguration();
			FitConfiguration tmpFitConfig = tmp.getFitConfiguration();
			tmpFitConfig.setComputeResiduals(true); // Collect the residuals threshold
			if (BenchmarkFilterAnalysis.updateConfiguration(tmp, false))
			{
				config.setFailuresLimit(tmp.getFailuresLimit());
				config.setIncludeNeighbours(tmp.isIncludeNeighbours());
				config.setNeighbourHeightThreshold(tmp.getNeighbourHeightThreshold());
				computeDoublets = (tmp.getResidualsThreshold() < 1);
				config.setDuplicateDistance(tmp.getDuplicateDistance());

				final DirectFilter primaryFilter = tmpFitConfig.getSmartFilter();
				final double residualsThreshold = tmp.getResidualsThreshold();
				myMultiFilter = new MultiPathFilter(primaryFilter, minimalFilter, residualsThreshold);
			}
		}
		else
		{
			myMultiFilter = MultiPathFilter.fromXML(xml);

			config.setFailuresLimit(failLimit);
			config.setIncludeNeighbours(includeNeighbours);
			config.setNeighbourHeightThreshold(neighbourHeightThreshold);
			computeDoublets = myComputeDoublets;
			config.setDuplicateDistance(myDuplicateDistance);
		}

		if (myMultiFilter == null)
		{
			gd = new ExtendedGenericDialog(TITLE);
			gd.addMessage("The multi-path filter was invalid.\n \nContinue with a default filter?");
			gd.enableYesNoCancel();
			gd.hideCancelButton();
			gd.showDialog();
			if (!gd.wasOKed())
				return false;
		}
		else
		{
			multiFilter = myMultiFilter;
		}

		if (computeDoublets)
		{
			//config.setComputeResiduals(true);
			config.setResidualsThreshold(0);
			fitConfig.setComputeResiduals(true);
		}
		else
		{
			config.setResidualsThreshold(1);
			fitConfig.setComputeResiduals(false);
		}
		showFilterScoreHistograms = gd.getNextBoolean();
		showCorrelation = gd.getNextBoolean();
		rankByIntensity = gd.getNextBoolean();
		saveFilterRange = gd.getNextBoolean();

		// Avoid stupidness, i.e. things that move outside the fit window and are bad widths

		// TODO - Fix this for simple or smart filter...
		fitConfig.setDisableSimpleFilter(false);

		fitConfig.setMinPhotons(15); // Realistically we cannot fit lower than this
		// Disable shift as candidates may be re-mapped to alternative candidates so the initial position is wrong.
		fitConfig.setCoordinateShiftFactor(0);
		fitConfig.setMinWidthFactor(1.0 / 5);
		fitConfig.setWidthFactor(5);
		// Disable the direct filter
		fitConfig.setDirectFilter(null);

		if (extraOptions)
		{
		}

		if (gd.invalidNumber())
			return false;

		if (lowerDistance > distance)
			lowerDistance = distance;
		if (lowerSignalFactor > signalFactor)
			lowerSignalFactor = signalFactor;

		// Distances relative to sa (not s) as this is the same as the BenchmarkSpotFilter plugin 
		distanceInPixels = distance * sa / simulationParameters.a;
		lowerDistanceInPixels = lowerDistance * sa / simulationParameters.a;

		// Copy simulation defaults if a new simulation
		if (lastId != simulationParameters.id)
		{
			// This is needed to configure the fit solver
			fitConfig.setNmPerPixel(simulationParameters.a);
			fitConfig.setGain(simulationParameters.gain);
			fitConfig.setAmplification(simulationParameters.amplification);
			fitConfig.setReadNoise(simulationParameters.readNoise);
			fitConfig.setBias(simulationParameters.bias);
			fitConfig.setEmCCD(simulationParameters.emCCD);
		}
		if (!PeakFit.configureFitSolver(config, (extraOptions) ? PeakFit.FLAG_EXTRA_OPTIONS : 0))
			return false;

		return true;
	}

	/** The total progress. */
	int progress, stepProgress, totalProgress;

	/**
	 * Show progress.
	 */
	private synchronized void showProgress()
	{
		if (progress % stepProgress == 0)
		{
			if (Utils.showStatus("Fitting frame: " + progress + " / " + totalProgress))
				IJ.showProgress(progress, totalProgress);
		}
		progress++;
	}

	private void run()
	{
		// Extract all the results in memory into a list per frame. This can be cached
		boolean refresh = false;
		if (lastId != simulationParameters.id)
		{
			// Do not get integer coordinates
			// The Coordinate objects will be PeakResultPoint objects that store the original PeakResult
			// from the MemoryPeakResults
			actualCoordinates = ResultsMatchCalculator.getCoordinates(results, false);
			lastId = simulationParameters.id;
			refresh = true;
		}

		// Extract all the candidates into a list per frame. This can be cached if the settings have not changed
		final int width = (config.isIncludeNeighbours()) ? config.getRelativeFitting() : 0;
		final Settings settings = new Settings(BenchmarkSpotFilter.filterResult.id, fractionPositives,
				fractionNegativesAfterAllPositives, negativesAfterAllPositives, width);
		if (refresh || !settings.equals(lastSettings))
		{
			filterCandidates = subsetFilterResults(BenchmarkSpotFilter.filterResult.filterResults, width);
			lastSettings = settings;
			lastFilterId = BenchmarkSpotFilter.filterResult.id;
		}

		stopWatch = StopWatch.createStarted();
		final ImageStack stack = imp.getImageStack();

		clearFitResults();

		// Save results to memory
		MemoryPeakResults peakResults = new MemoryPeakResults();
		peakResults.copySettings(this.results);
		peakResults.setName(TITLE);
		MemoryPeakResults.addResults(peakResults);

		// Create a pool of workers
		final int nThreads = Prefs.getThreads();
		BlockingQueue<Integer> jobs = new ArrayBlockingQueue<Integer>(nThreads * 2);
		List<Worker> workers = new LinkedList<Worker>();
		List<Thread> threads = new LinkedList<Thread>();
		for (int i = 0; i < nThreads; i++)
		{
			Worker worker = new Worker(jobs, stack, actualCoordinates, filterCandidates, peakResults);
			Thread t = new Thread(worker);
			workers.add(worker);
			threads.add(t);
			t.start();
		}

		// Fit the frames
		long runTime = System.nanoTime();
		totalProgress = stack.getSize();
		stepProgress = Utils.getProgressInterval(totalProgress);
		progress = 0;
		for (int i = 1; i <= totalProgress; i++)
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

		IJ.showProgress(1);
		runTime = System.nanoTime() - runTime;

		if (Utils.isInterrupted())
		{
			return;
		}

		stopWatch.stop();
		final String timeString = stopWatch.toString();
		IJ.log("Spot fit time : " + timeString);

		IJ.showStatus("Collecting results ...");

		fitResultsId++;
		fitResults = new TIntObjectHashMap<FilterCandidates>();
		for (Worker w : workers)
		{
			fitResults.putAll(w.results);
		}

		// Assign a unique ID to each result
		int count = 0;
		// Materialise into an array since we use it twice
		FilterCandidates[] candidates = fitResults.values(new FilterCandidates[fitResults.size()]);
		for (FilterCandidates result : candidates)
		{
			for (int i = 0; i < result.fitResult.length; i++)
			{
				final MultiPathFitResult fitResult = result.fitResult[i];
				count += count(fitResult.getSingleFitResult());
				count += count(fitResult.getMultiFitResult());
				count += count(fitResult.getDoubletFitResult());
				count += count(fitResult.getMultiDoubletFitResult());
			}
		}
		PreprocessedPeakResult[] preprocessedPeakResults = new PreprocessedPeakResult[count];
		count = 0;
		for (FilterCandidates result : candidates)
		{
			for (int i = 0; i < result.fitResult.length; i++)
			{
				final MultiPathFitResult fitResult = result.fitResult[i];
				count = store(fitResult.getSingleFitResult(), count, preprocessedPeakResults);
				count = store(fitResult.getMultiFitResult(), count, preprocessedPeakResults);
				count = store(fitResult.getDoubletFitResult(), count, preprocessedPeakResults);
				count = store(fitResult.getMultiDoubletFitResult(), count, preprocessedPeakResults);
			}
		}

		summariseResults(fitResults, runTime, preprocessedPeakResults, count);

		IJ.showStatus("");
	}

	/**
	 * Clear old results to free memory
	 */
	private void clearFitResults()
	{
		if (fitResults != null)
			fitResults.clear();
		fitResults = null;
	}

	private int count(gdsc.smlm.results.filter.MultiPathFitResult.FitResult fitResult)
	{
		if (fitResult == null || fitResult.results == null)
			return 0;
		int count = 0;
		for (int i = 0; i < fitResult.results.length; i++)
		{
			PreprocessedPeakResult result = fitResult.results[i];
			if (result.isNewResult())
				count++;
		}
		return count;
	}

	private int store(gdsc.smlm.results.filter.MultiPathFitResult.FitResult fitResult, int count,
			PreprocessedPeakResult[] preprocessedPeakResults)
	{
		if (fitResult == null || fitResult.results == null)
			return count;
		for (int i = 0; i < fitResult.results.length; i++)
		{
			BasePreprocessedPeakResult result = (BasePreprocessedPeakResult) fitResult.results[i];
			if (result.isNewResult())
			{
				result.uniqueId = count++;
				preprocessedPeakResults[result.uniqueId] = result;
			}
		}
		return count;
	}

	/**
	 * Extract all the filter candidates in order until the desired number of positives have been reached and the number
	 * of negatives matches the configured parameters.
	 * 
	 * @param filterResults
	 * @return The filter candidates
	 */
	private TIntObjectHashMap<FilterCandidates> subsetFilterResults(TIntObjectHashMap<FilterResult> filterResults,
			int fitting)
	{
		// Convert fractions from percent 
		final double f1 = Math.min(1, fractionPositives / 100.0);
		final double f2 = fractionNegativesAfterAllPositives / 100.0;

		final int[] counter = new int[2];

		final TIntObjectHashMap<FilterCandidates> subset = new TIntObjectHashMap<FilterCandidates>();
		fP = fN = 0;
		nP = nN = 0;
		final double[] fX = new double[2];
		final int[] nX = new int[2];
		filterResults.forEachEntry(new TIntObjectProcedure<FilterResult>()
		{
			public boolean execute(int frame, FilterResult r)
			{
				// Determine the number of positives to find. This score may be fractional.
				fX[0] += r.result.getTP();
				fX[1] += r.result.getFP();

				// Q. Is r.result.getTP() not the same as the total of r.spots[i].match?
				// A. Not if we used fractional scoring.
				int c = 0;
				for (int i = r.spots.length; i-- > 0;)
				{
					if (r.spots[i].match)
						c++;
				}
				nX[0] += c;
				nX[1] += (r.spots.length - c);

				// Make the target use the fractional score
				final double np2 = r.result.getTP() * f1;
				double targetP = np2;

				// Set the target using the closest
				if (f1 < 1)
				{
					double np = 0;
					double min = r.result.getTP();
					for (ScoredSpot spot : r.spots)
					{
						if (spot.match)
						{
							np += spot.getScore();
							double d = np2 - np;
							if (d < min)
							{
								min = d;
								targetP = np;
							}
							else
							{
								break;
							}
						}
					}

					//if (targetP < np2)
					//	System.out.printf("np2 = %.2f, targetP = %.2f\n", np2, targetP);
				}

				// Count the number of positive & negatives
				int p = 0, n = 0;
				double np = 0, nn = 0;

				boolean reachedTarget = false;
				int nAfter = 0;

				int count = 0;
				for (ScoredSpot spot : r.spots)
				{
					count++;
					nn += spot.antiScore();
					if (spot.match)
					{
						np += spot.getScore();
						p++;
						if (!reachedTarget)
						{
							reachedTarget = np >= targetP;
						}
					}
					else
					{
						n++;
						if (reachedTarget)
						{
							nAfter++;
						}
					}

					if (reachedTarget)
					{
						// Check if we have reached both the limits
						if (nAfter >= negativesAfterAllPositives && (double) n / (n + p) >= f2)
							break;
					}
				}

				counter[0] += count;
				counter[1] += r.spots.length;

				// Debug
				//System.out.printf("Frame %d : %.1f / (%.1f + %.1f). p=%d, n=%d, after=%d, f=%.1f\n", result.getKey().intValue(),
				//		r.result.getTP(), r.result.getTP(), r.result.getFP(), p, n,
				//		nAfter, (double) n / (n + p));

				// We can use all the candidates but only fit up to count
				subset.put(frame, new FilterCandidates(p, n, np, nn, r.spots, count));
				return true;
			}
		});

		fP = fX[0];
		fN = fX[1];
		nP = nX[0];
		nN = nX[1];

		// We now add all the candidates but only fit the first N
		int target = counter[0];
		int total = counter[1];
		int added = total - target;

		if (extraOptions && added > target)
			Utils.log("Added %s to %s (total = %d)", Utils.pleural(added, "neighbour"),
					Utils.pleural(target, "candidate"), total);

		return subset;
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

	/**
	 * Create an abstract class to allow a count to be passed to the constructor. The procedure can be coded inline
	 * using final object references.
	 */
	private abstract class CustomTIntProcedure implements TIntProcedure
	{
		int c;

		CustomTIntProcedure(int count)
		{
			c = count;
		}
	}

	private void summariseResults(TIntObjectHashMap<FilterCandidates> filterCandidates, long runTime,
			final PreprocessedPeakResult[] preprocessedPeakResults, int nUniqueIDs)
	{
		createTable();

		// Summarise the fitting results. N fits, N failures. 
		// Optimal match statistics if filtering is perfect (since fitting is not perfect).
		StoredDataStatistics distanceStats = new StoredDataStatistics();
		StoredDataStatistics depthStats = new StoredDataStatistics();

		// Get stats for all fitted results and those that match 
		// Signal, SNR, Width, xShift, yShift, Precision
		createFilterCriteria();
		StoredDataStatistics[][] stats = new StoredDataStatistics[3][filterCriteria.length];
		for (int i = 0; i < stats.length; i++)
			for (int j = 0; j < stats[i].length; j++)
				stats[i][j] = new StoredDataStatistics();

		final double nmPerPixel = simulationParameters.a;
		double tp = 0, fp = 0;
		int failcTP = 0, failcFP = 0;
		int cTP = 0, cFP = 0;
		int[] singleStatus = null, multiStatus = null, doubletStatus = null, multiDoubletStatus = null;
		singleStatus = new int[FitStatus.values().length];
		multiStatus = new int[singleStatus.length];
		doubletStatus = new int[singleStatus.length];
		multiDoubletStatus = new int[singleStatus.length];

		// Easier to materialise the values since we have a lot of non final variables to manipulate
		final int[] frames = new int[filterCandidates.size()];
		final FilterCandidates[] candidates = new FilterCandidates[filterCandidates.size()];
		final int[] counter = new int[1];
		filterCandidates.forEachEntry(new TIntObjectProcedure<FilterCandidates>()
		{
			public boolean execute(int a, FilterCandidates b)
			{
				frames[counter[0]] = a;
				candidates[counter[0]] = b;
				counter[0]++;
				return true;
			}
		});

		for (FilterCandidates result : candidates)
		{
			// Count the number of fit results that matched (tp) and did not match (fp)
			tp += result.tp;
			fp += result.fp;

			for (int i = 0; i < result.fitResult.length; i++)
			{
				if (result.spots[i].match)
					cTP++;
				else
					cFP++;
				final MultiPathFitResult fitResult = result.fitResult[i];

				if (singleStatus != null && result.spots[i].match)
				{
					// Debugging reasons for fit failure
					addStatus(singleStatus, fitResult.getSingleFitResult());
					addStatus(multiStatus, fitResult.getMultiFitResult());
					addStatus(doubletStatus, fitResult.getDoubletFitResult());
					addStatus(multiDoubletStatus, fitResult.getMultiDoubletFitResult());
				}

				if (noMatch(fitResult))
				{
					if (result.spots[i].match)
						failcTP++;
					else
						failcFP++;
				}

				// We have multi-path results.
				// We want statistics for:
				// [0] all fitted spots
				// [1] fitted spots that match a result
				// [2] fitted spots that do not match a result
				addToStats(fitResult.getSingleFitResult(), stats);
				addToStats(fitResult.getMultiFitResult(), stats);
				addToStats(fitResult.getDoubletFitResult(), stats);
				addToStats(fitResult.getMultiDoubletFitResult(), stats);
			}

			// Statistics on spots that fit an actual result
			for (int i = 0; i < result.match.length; i++)
			{
				if (!result.match[i].isFitResult())
					// For now just ignore the candidates that matched
					continue;

				FitMatch fitMatch = (FitMatch) result.match[i];
				distanceStats.add(fitMatch.d * nmPerPixel);
				depthStats.add(fitMatch.z * nmPerPixel);
			}
		}

		// Store data for computing correlation
		double[] i1 = new double[depthStats.getN()];
		double[] i2 = new double[i1.length];
		double[] is = new double[i1.length];
		int ci = 0;
		for (FilterCandidates result : candidates)
		{
			for (int i = 0; i < result.match.length; i++)
			{
				if (!result.match[i].isFitResult())
					// For now just ignore the candidates that matched
					continue;

				FitMatch fitMatch = (FitMatch) result.match[i];
				ScoredSpot spot = result.spots[fitMatch.i];
				i1[ci] = fitMatch.predictedSignal;
				i2[ci] = fitMatch.actualSignal;
				is[ci] = spot.spot.intensity;
				ci++;
			}
		}

		// We want to compute the Jaccard against the spot metric

		// Filter the results using the multi-path filter
		ArrayList<MultiPathFitResults> multiPathResults = new ArrayList<MultiPathFitResults>(filterCandidates.size());
		for (int i = 0; i < frames.length; i++)
		{
			int frame = frames[i];
			MultiPathFitResult[] multiPathFitResults = candidates[i].fitResult;
			int totalCandidates = candidates[i].spots.length;
			int nActual = actualCoordinates.get(frame).size();
			multiPathResults.add(new MultiPathFitResults(frame, multiPathFitResults, totalCandidates, nActual));
		}
		// Score the results and count the number returned
		List<FractionalAssignment[]> assignments = new ArrayList<FractionalAssignment[]>();
		final TIntHashSet set = new TIntHashSet(nUniqueIDs);
		FractionScoreStore scoreStore = new FractionScoreStore()
		{
			public void add(int uniqueId)
			{
				set.add(uniqueId);
			}
		};
		MultiPathFitResults[] multiResults = multiPathResults.toArray(new MultiPathFitResults[multiPathResults.size()]);
		// Filter with no filter
		MultiPathFilter mpf = new MultiPathFilter(new SignalFilter(0), null, multiFilter.residualsThreshold);
		FractionClassificationResult fractionResult = mpf.fractionScoreSubset(multiResults, Integer.MAX_VALUE,
				this.results.size(), assignments, scoreStore,
				CoordinateStoreFactory.create(imp.getWidth(), imp.getHeight(), config.getDuplicateDistance()));
		double nPredicted = fractionResult.getTP() + fractionResult.getFP();

		final double[][] matchScores = new double[set.size()][];
		int count = 0;
		for (int i = 0; i < assignments.size(); i++)
		{
			FractionalAssignment[] a = assignments.get(i);
			if (a == null)
				continue;
			for (int j = 0; j < a.length; j++)
			{
				final PreprocessedPeakResult r = ((PeakFractionalAssignment) a[j]).peakResult;
				set.remove(r.getUniqueId());

				final double precision = Math.sqrt(r.getLocationVariance());
				final double signal = r.getSignal();
				final double snr = r.getSNR();
				final double width = r.getXSDFactor();
				final double xShift = r.getXRelativeShift2();
				final double yShift = r.getYRelativeShift2();
				// Since these two are combined for filtering and the max is what matters.
				final double shift = (xShift > yShift) ? Math.sqrt(xShift) : Math.sqrt(yShift);
				final double eshift = Math.sqrt(xShift + yShift);

				final double[] score = new double[8];
				score[FILTER_SIGNAL] = signal;
				score[FILTER_SNR] = snr;
				score[FILTER_MIN_WIDTH] = width;
				score[FILTER_MAX_WIDTH] = width;
				score[FILTER_SHIFT] = shift;
				score[FILTER_ESHIFT] = eshift;
				score[FILTER_PRECISION] = precision;
				score[FILTER_PRECISION + 1] = a[j].getScore();
				matchScores[count++] = score;
			}
		}
		// Add the rest
		set.forEach(new CustomTIntProcedure(count)
		{
			public boolean execute(int uniqueId)
			{
				// This should not be null or something has gone wrong
				PreprocessedPeakResult r = preprocessedPeakResults[uniqueId];
				if (r == null)
					throw new RuntimeException("Missing result: " + uniqueId);
				final double precision = Math.sqrt(r.getLocationVariance());
				final double signal = r.getSignal();
				final double snr = r.getSNR();
				final double width = r.getXSDFactor();
				final double xShift = r.getXRelativeShift2();
				final double yShift = r.getYRelativeShift2();
				// Since these two are combined for filtering and the max is what matters.
				final double shift = (xShift > yShift) ? Math.sqrt(xShift) : Math.sqrt(yShift);
				final double eshift = Math.sqrt(xShift + yShift);

				final double[] score = new double[8];
				score[FILTER_SIGNAL] = signal;
				score[FILTER_SNR] = snr;
				score[FILTER_MIN_WIDTH] = width;
				score[FILTER_MAX_WIDTH] = width;
				score[FILTER_SHIFT] = shift;
				score[FILTER_ESHIFT] = eshift;
				score[FILTER_PRECISION] = precision;
				matchScores[c++] = score;
				return true;
			}
		});

		// Debug the reasons the fit failed
		if (singleStatus != null)
		{
			String name = PeakFit.getSolverName(fitConfig);
			if (fitConfig.getFitSolver() == FitSolver.MLE && fitConfig.isModelCamera())
				name += " Camera";
			System.out.println("Failure counts: " + name);
			printFailures("Single", singleStatus);
			printFailures("Multi", multiStatus);
			printFailures("Doublet", doubletStatus);
			printFailures("Multi doublet", multiDoubletStatus);
		}

		StringBuilder sb = new StringBuilder(300);

		// Add information about the simulation
		final double signal = simulationParameters.signalPerFrame; //(simulationParameters.minSignal + simulationParameters.maxSignal) * 0.5;
		final int n = results.size();
		sb.append(imp.getStackSize()).append('\t');
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		sb.append(w).append('\t');
		sb.append(h).append('\t');
		sb.append(n).append('\t');
		double density = ((double) n / imp.getStackSize()) / (w * h) /
				(simulationParameters.a * simulationParameters.a / 1e6);
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

		if (simulationParameters.fullSimulation)
		{
			// The total signal is spread over frames
		}

		sb.append(Utils.rounded(signal / Math.sqrt(noise))).append('\t');
		sb.append(Utils.rounded(simulationParameters.s / simulationParameters.a)).append('\t');

		sb.append(spotFilter.getDescription());

		// nP and nN is the fractional score of the spot candidates 
		addCount(sb, nP + nN);
		addCount(sb, nP);
		addCount(sb, nN);
		addCount(sb, fP);
		addCount(sb, fN);
		String name = PeakFit.getSolverName(fitConfig);
		if (fitConfig.getFitSolver() == FitSolver.MLE && fitConfig.isModelCamera())
			name += " Camera";
		add(sb, name);
		add(sb, config.getFitting());

		resultPrefix = sb.toString();

		// Q. Should I add other fit configuration here?

		// The fraction of positive and negative candidates that were included
		add(sb, (100.0 * cTP) / nP);
		add(sb, (100.0 * cFP) / nN);

		// Score the fitting results compared to the original simulation.

		// Score the candidate selection:
		add(sb, cTP + cFP);
		add(sb, cTP);
		add(sb, cFP);
		// TP are all candidates that can be matched to a spot
		// FP are all candidates that cannot be matched to a spot
		// FN = The number of missed spots
		FractionClassificationResult m = new FractionClassificationResult(cTP, cFP, 0,
				simulationParameters.molecules - cTP);
		add(sb, m.getRecall());
		add(sb, m.getPrecision());
		add(sb, m.getF1Score());
		add(sb, m.getJaccard());

		// Score the fitting results:
		add(sb, failcTP);
		add(sb, failcFP);

		// TP are all fit results that can be matched to a spot
		// FP are all fit results that cannot be matched to a spot
		// FN = The number of missed spots
		add(sb, tp);
		add(sb, fp);
		m = new FractionClassificationResult(tp, fp, 0, simulationParameters.molecules - tp);
		add(sb, m.getRecall());
		add(sb, m.getPrecision());
		add(sb, m.getF1Score());
		add(sb, m.getJaccard());

		// Do it again but pretend we can perfectly filter all the false positives
		//add(sb, tp);
		m = new FractionClassificationResult(tp, 0, 0, simulationParameters.molecules - tp);
		// Recall is unchanged
		// Precision will be 100%
		add(sb, m.getF1Score());
		add(sb, m.getJaccard());

		// The mean may be subject to extreme outliers so use the median
		double median = distanceStats.getMedian();
		add(sb, median);

		WindowOrganiser wo = new WindowOrganiser();

		String label = String.format("Recall = %s. n = %d. Median = %s nm. SD = %s nm", Utils.rounded(m.getRecall()),
				distanceStats.getN(), Utils.rounded(median), Utils.rounded(distanceStats.getStandardDeviation()));
		int id = Utils.showHistogram(TITLE, distanceStats, "Match Distance (nm)", 0, 0, 0, label);
		if (Utils.isNewWindow())
			wo.add(id);

		median = depthStats.getMedian();
		add(sb, median);

		// Sort by spot intensity and produce correlation
		int[] indices = Utils.newArray(i1.length, 0, 1);
		if (showCorrelation)
			Sort.sort(indices, is, rankByIntensity);
		double[] r = (showCorrelation) ? new double[i1.length] : null;
		double[] sr = (showCorrelation) ? new double[i1.length] : null;
		double[] rank = (showCorrelation) ? new double[i1.length] : null;
		ci = 0;
		FastCorrelator fastCorrelator = new FastCorrelator();
		ArrayList<Ranking> pc1 = new ArrayList<Ranking>();
		ArrayList<Ranking> pc2 = new ArrayList<Ranking>();
		for (int ci2 : indices)
		{
			fastCorrelator.add((long) Math.round(i1[ci2]), (long) Math.round(i2[ci2]));
			pc1.add(new Ranking(i1[ci2], ci));
			pc2.add(new Ranking(i2[ci2], ci));
			if (showCorrelation)
			{
				r[ci] = fastCorrelator.getCorrelation();
				sr[ci] = Correlator.correlation(rank(pc1), rank(pc2));
				if (rankByIntensity)
					rank[ci] = is[0] - is[ci];
				else
					rank[ci] = ci;
			}
			ci++;
		}

		final double pearsonCorr = fastCorrelator.getCorrelation();
		final double rankedCorr = Correlator.correlation(rank(pc1), rank(pc2));

		// Get the regression
		SimpleRegression regression = new SimpleRegression(false);
		for (int i = 0; i < pc1.size(); i++)
			regression.addData(pc1.get(i).value, pc2.get(i).value);
		//final double intercept = regression.getIntercept();
		final double slope = regression.getSlope();

		if (showCorrelation)
		{
			String title = TITLE + " Intensity";
			Plot plot = new Plot(title, "Candidate", "Spot");
			double[] limits1 = Maths.limits(i1);
			double[] limits2 = Maths.limits(i2);
			plot.setLimits(limits1[0], limits1[1], limits2[0], limits2[1]);
			label = String.format("Correlation=%s; Ranked=%s; Slope=%s", Utils.rounded(pearsonCorr),
					Utils.rounded(rankedCorr), Utils.rounded(slope));
			plot.addLabel(0, 0, label);
			plot.setColor(Color.red);
			plot.addPoints(i1, i2, Plot.DOT);
			if (slope > 1)
				plot.drawLine(limits1[0], limits1[0] * slope, limits1[1], limits1[1] * slope);
			else
				plot.drawLine(limits2[0] / slope, limits2[0], limits2[1] / slope, limits2[1]);
			PlotWindow pw = Utils.display(title, plot);
			if (Utils.isNewWindow())
				wo.add(pw);

			title = TITLE + " Correlation";
			plot = new Plot(title, "Spot Rank", "Correlation");
			double[] xlimits = Maths.limits(rank);
			double[] ylimits = Maths.limits(r);
			ylimits = Maths.limits(ylimits, sr);
			plot.setLimits(xlimits[0], xlimits[1], ylimits[0], ylimits[1]);
			plot.setColor(Color.red);
			plot.addPoints(rank, r, Plot.LINE);
			plot.setColor(Color.blue);
			plot.addPoints(rank, sr, Plot.LINE);
			plot.setColor(Color.black);
			plot.addLabel(0, 0, label);
			pw = Utils.display(title, plot);
			if (Utils.isNewWindow())
				wo.add(pw);
		}

		add(sb, pearsonCorr);
		add(sb, rankedCorr);
		add(sb, slope);

		label = String.format("n = %d. Median = %s nm", depthStats.getN(), Utils.rounded(median));
		id = Utils.showHistogram(TITLE, depthStats, "Match Depth (nm)", 0, 1, 0, label);
		if (Utils.isNewWindow())
			wo.add(id);

		// Plot histograms of the stats on the same window
		double[] lower = new double[filterCriteria.length];
		double[] upper = new double[lower.length];
		min = new double[lower.length];
		max = new double[lower.length];
		for (int i = 0; i < stats[0].length; i++)
		{
			double[] limits = showDoubleHistogram(stats, i, wo, matchScores, nPredicted);
			lower[i] = limits[0];
			upper[i] = limits[1];
			min[i] = limits[2];
			max[i] = limits[3];
		}

		// Reconfigure some of the range limits
		upper[FILTER_SIGNAL] *= 2; // Make this a bit bigger
		upper[FILTER_SNR] *= 2; // Make this a bit bigger
		double factor = 0.25;
		if (lower[FILTER_MIN_WIDTH] != 0)
			upper[FILTER_MIN_WIDTH] = 1 - Math.max(0, factor * (1 - lower[FILTER_MIN_WIDTH])); // (assuming lower is less than 1)
		if (upper[FILTER_MIN_WIDTH] != 0)
			lower[FILTER_MAX_WIDTH] = 1 + Math.max(0, factor * (upper[FILTER_MAX_WIDTH] - 1)); // (assuming upper is more than 1)

		// Round the ranges
		final double[] interval = new double[stats[0].length];
		interval[FILTER_SIGNAL] = SignalFilter.DEFAULT_INCREMENT;
		interval[FILTER_SNR] = SNRFilter.DEFAULT_INCREMENT;
		interval[FILTER_MIN_WIDTH] = WidthFilter2.DEFAULT_MIN_INCREMENT;
		interval[FILTER_MAX_WIDTH] = WidthFilter.DEFAULT_INCREMENT;
		interval[FILTER_SHIFT] = ShiftFilter.DEFAULT_INCREMENT;
		interval[FILTER_ESHIFT] = EShiftFilter.DEFAULT_INCREMENT;
		interval[FILTER_PRECISION] = PrecisionFilter.DEFAULT_INCREMENT;
		interval[FILTER_ITERATIONS] = 0.1;
		interval[FILTER_EVALUATIONS] = 0.1;

		// Create a range increment
		double[] increment = new double[lower.length];
		for (int i = 0; i < increment.length; i++)
		{
			lower[i] = Maths.floor(lower[i], interval[i]);
			upper[i] = Maths.ceil(upper[i], interval[i]);
			double range = upper[i] - lower[i];
			// Allow clipping if the range is small compared to the min increment
			double multiples = range / interval[i];
			// Use 8 multiples for the equivalent of +/- 4 steps around the centre
			if (multiples < 8)
			{
				multiples = Math.ceil(multiples);
			}
			else
				multiples = 8;
			increment[i] = Maths.ceil(range / multiples, interval[i]);

			if (i == FILTER_MIN_WIDTH)
				// Requires clipping based on the upper limit
				lower[i] = upper[i] - increment[i] * multiples;
			else
				upper[i] = lower[i] + increment[i] * multiples;
		}

		for (int i = 0; i < stats[0].length; i++)
		{
			lower[i] = Maths.round(lower[i]);
			upper[i] = Maths.round(upper[i]);
			min[i] = Maths.round(min[i]);
			max[i] = Maths.round(max[i]);
			increment[i] = Maths.round(increment[i]);
			sb.append('\t').append(min[i]).append(':').append(lower[i]).append('-').append(upper[i]).append(':')
					.append(max[i]);
		}

		// Disable some filters
		increment[FILTER_SIGNAL] = Double.POSITIVE_INFINITY;
		//increment[FILTER_SHIFT] = Double.POSITIVE_INFINITY;
		increment[FILTER_ESHIFT] = Double.POSITIVE_INFINITY;

		wo.tile();

		sb.append('\t').append(Utils.timeToString(runTime / 1000000.0));

		summaryTable.append(sb.toString());

		if (saveFilterRange)
		{
			GUIFilterSettings filterSettings = SettingsManager.readGUIFilterSettings(0);

			String filename = (silent) ? filterSettings.getFilterSetFilename()
					: Utils.getFilename("Filter_range_file", filterSettings.getFilterSetFilename());
			if (filename == null)
				return;
			// Remove extension to store the filename
			filename = Utils.replaceExtension(filename, ".xml");
			filterSettings = filterSettings.toBuilder().setFilterSetFilename(filename).build();

			// Create a filter set using the ranges
			ArrayList<Filter> filters = new ArrayList<Filter>(3);
			filters.add(new MultiFilter2(lower[0], (float) lower[1], lower[2], lower[3], lower[4], lower[5], lower[6]));
			filters.add(new MultiFilter2(upper[0], (float) upper[1], upper[2], upper[3], upper[4], upper[5], upper[6]));
			filters.add(new MultiFilter2(increment[0], (float) increment[1], increment[2], increment[3], increment[4],
					increment[5], increment[6]));
			if (saveFilters(filename, filters))
				SettingsManager.writeSettings(filterSettings);

			// Create a filter set using the min/max and the initial bounds.
			// Set sensible limits
			min[FILTER_SIGNAL] = Math.max(min[FILTER_SIGNAL], 30);
			max[FILTER_PRECISION] = Math.min(max[FILTER_PRECISION], 100);

			// Commented this out so that the 4-set filters are the same as the 3-set filters.
			// The difference leads to differences when optimising.
			//			// Use half the initial bounds (hoping this is a good starting guess for the optimum)
			//			final boolean[] limitToLower = new boolean[min.length];
			//			limitToLower[FILTER_SIGNAL] = true;
			//			limitToLower[FILTER_SNR] = true;
			//			limitToLower[FILTER_MIN_WIDTH] = true;
			//			limitToLower[FILTER_MAX_WIDTH] = false;
			//			limitToLower[FILTER_SHIFT] = false;
			//			limitToLower[FILTER_ESHIFT] = false;
			//			limitToLower[FILTER_PRECISION] = true;
			//			for (int i = 0; i < limitToLower.length; i++)
			//			{
			//				final double range = (upper[i] - lower[i]) / 2;
			//				if (limitToLower[i])
			//					upper[i] = lower[i] + range;
			//				else
			//					lower[i] = upper[i] - range;
			//			}

			filters = new ArrayList<Filter>(4);
			filters.add(new MultiFilter2(min[0], (float) min[1], min[2], min[3], min[4], min[5], min[6]));
			filters.add(new MultiFilter2(lower[0], (float) lower[1], lower[2], lower[3], lower[4], lower[5], lower[6]));
			filters.add(new MultiFilter2(upper[0], (float) upper[1], upper[2], upper[3], upper[4], upper[5], upper[6]));
			filters.add(new MultiFilter2(max[0], (float) max[1], max[2], max[3], max[4], max[5], max[6]));
			saveFilters(Utils.replaceExtension(filename, ".4.xml"), filters);
		}
	}

	private boolean saveFilters(String filename, ArrayList<Filter> filters)
	{
		ArrayList<FilterSet> filterList = new ArrayList<FilterSet>(1);
		// Add Range keyword to identify as a range filter set 
		filterList.add(new FilterSet("Range", filters));
		FileOutputStream fos = null;
		try
		{
			fos = new FileOutputStream(filename);
			// Use the instance (not .toXML() method) to allow the exception to be caught
			XStreamWrapper.getInstance().toXML(filterList, fos);
			return true;
		}
		catch (Exception e)
		{
			IJ.log("Unable to save the filter set to file: " + e.getMessage());
		}
		finally
		{
			if (fos != null)
			{
				try
				{
					fos.close();
				}
				catch (IOException e)
				{
					// Ignore
				}
			}
		}
		return false;
	}

	private void printFailures(String title, int[] status)
	{
		int total = 0;
		// Count failures
		for (int i = 1; i < status.length; i++)
		{
			if (status[i] != 0)
			{
				total += status[i];
			}
		}
		// Print failures
		if (total != 0)
			for (int i = 1; i < status.length; i++)
			{
				if (status[i] != 0)
				{
					System.out.printf("%s %s = %d / %d  (%.2f)\n", title, FitStatus.values()[i].toString(), status[i],
							total, 100.0 * status[i] / total);
				}
			}
		// Print total
		int all = total + status[0];
		if (all != 0)
			System.out.printf("%s %s = %d / %d  (%.2f)\n", title, "Total", total, all, 100.0 * total / all);
	}

	private void addStatus(int[] status, gdsc.smlm.results.filter.MultiPathFitResult.FitResult fitResult)
	{
		if (fitResult != null)
			status[fitResult.status]++;
	}

	private void addToStats(gdsc.smlm.results.filter.MultiPathFitResult.FitResult fitResult,
			StoredDataStatistics[][] stats)
	{
		if (fitResult == null)
			return;

		final FitResult actualFitResult = (FitResult) fitResult.getData();

		if (fitResult.status != 0)
		{
			// Add the evaluations for spots that were not OK
			stats[0][FILTER_ITERATIONS].add(actualFitResult.getIterations());
			stats[0][FILTER_EVALUATIONS].add(actualFitResult.getEvaluations());
			return;
		}

		if (fitResult.results == null)
			return;

		boolean isMatch = false;

		for (int resultIndex = 0; resultIndex < fitResult.results.length; resultIndex++)
		{
			BasePreprocessedPeakResult result = (BasePreprocessedPeakResult) fitResult.results[resultIndex];

			// Q. Only build stats on new results?
			if (!result.isNewResult())
				continue;

			// This was fit - Get statistics
			final double precision = Math.sqrt(result.getLocationVariance());

			final double signal = result.getSignal();
			final double snr = result.getSNR();
			final double width = result.getXSDFactor();
			final double xShift = result.getXRelativeShift2();
			final double yShift = result.getYRelativeShift2();
			// Since these two are combined for filtering and the max is what matters.
			final double shift = (xShift > yShift) ? Math.sqrt(xShift) : Math.sqrt(yShift);
			final double eshift = Math.sqrt(xShift + yShift);

			stats[0][FILTER_SIGNAL].add(signal);
			stats[0][FILTER_SNR].add(snr);
			if (width < 1)
				stats[0][FILTER_MIN_WIDTH].add(width);
			else
				stats[0][FILTER_MAX_WIDTH].add(width);
			stats[0][FILTER_SHIFT].add(shift);
			stats[0][FILTER_ESHIFT].add(eshift);
			stats[0][FILTER_PRECISION].add(precision);
			if (resultIndex == 0)
			{
				stats[0][FILTER_ITERATIONS].add(actualFitResult.getIterations());
				stats[0][FILTER_EVALUATIONS].add(actualFitResult.getEvaluations());
			}

			// Add to the TP or FP stats 
			// If it has assignments then it was a match to something
			isMatch |= result.hasAssignments();
			final int index = (result.hasAssignments()) ? 1 : 2;
			stats[index][FILTER_SIGNAL].add(signal);
			stats[index][FILTER_SNR].add(snr);
			if (width < 1)
				stats[index][FILTER_MIN_WIDTH].add(width);
			else
				stats[index][FILTER_MAX_WIDTH].add(width);
			stats[index][FILTER_SHIFT].add(shift);
			stats[index][FILTER_ESHIFT].add(eshift);
			stats[index][FILTER_PRECISION].add(precision);
			if (resultIndex == 0)
			{
			}
		}

		final int index = (isMatch) ? 1 : 2;
		stats[index][FILTER_ITERATIONS].add(actualFitResult.getIterations());
		stats[index][FILTER_EVALUATIONS].add(actualFitResult.getEvaluations());
	}

	private boolean noMatch(MultiPathFitResult fitResult)
	{
		if (isMatch(fitResult.getSingleFitResult()))
			return false;
		if (isMatch(fitResult.getMultiFitResult()))
			return false;
		if (isMatch(fitResult.getDoubletFitResult()))
			return false;
		if (isMatch(fitResult.getMultiDoubletFitResult()))
			return false;

		return true;
	}

	private boolean isMatch(gdsc.smlm.results.filter.MultiPathFitResult.FitResult fitResult)
	{
		if (fitResult == null)
			return false;
		if (fitResult.status != 0)
			return false;
		if (fitResult.results == null)
			return false;
		for (int resultIndex = 0; resultIndex < fitResult.results.length; resultIndex++)
		{
			BasePreprocessedPeakResult result = (BasePreprocessedPeakResult) fitResult.results[resultIndex];
			if (!result.isNewResult())
				continue;
			if (result.hasAssignments())
				return true;
		}
		return false;
	}

	private int[] rank(ArrayList<Ranking> pc)
	{
		Collections.sort(pc);
		int[] ranking = new int[pc.size()];
		int rank = 1;
		for (Ranking r : pc)
		{
			ranking[r.index] = rank++;
		}
		return ranking;
	}

	private double[] showDoubleHistogram(StoredDataStatistics[][] stats, final int i, WindowOrganiser wo,
			double[][] matchScores, double nPredicted)
	{
		String xLabel = filterCriteria[i].name;
		LowerLimit lower = filterCriteria[i].lower;
		UpperLimit upper = filterCriteria[i].upper;

		double[] j = null;
		double[] metric = null;
		double maxJ = 0;
		if (i <= FILTER_PRECISION && (showFilterScoreHistograms || upper.requiresJaccard || lower.requiresJaccard))
		{
			// Jaccard score verses the range of the metric
			Arrays.sort(matchScores, new Comparator<double[]>()
			{
				public int compare(double[] o1, double[] o2)
				{
					if (o1[i] < o2[i])
						return -1;
					if (o1[i] > o2[i])
						return 1;
					return 0;
				}
			});

			final int scoreIndex = FILTER_PRECISION + 1;
			int n = results.size();
			double tp = 0;
			double fp = 0;
			j = new double[matchScores.length + 1];
			metric = new double[j.length];
			for (int k = 0; k < matchScores.length; k++)
			{
				final double score = matchScores[k][scoreIndex];
				tp += score;
				fp += (1 - score);
				j[k + 1] = tp / (fp + n);
				metric[k + 1] = matchScores[k][i];
			}
			metric[0] = metric[1];
			maxJ = Maths.max(j);

			if (showFilterScoreHistograms)
			{
				String title = TITLE + " Jaccard " + xLabel;
				Plot plot = new Plot(title, xLabel, "Jaccard", metric, j);
				// Remove outliers
				double[] limitsx = Maths.limits(metric);
				Percentile p = new Percentile();
				double l = p.evaluate(metric, 25);
				double u = p.evaluate(metric, 75);
				double iqr = 1.5 * (u - l);
				limitsx[1] = Math.min(limitsx[1], u + iqr);
				plot.setLimits(limitsx[0], limitsx[1], 0, Maths.max(j));
				PlotWindow pw = Utils.display(title, plot);
				if (Utils.isNewWindow())
					wo.add(pw);
			}
		}

		// [0] is all
		// [1] is matches
		// [2] is no match
		StoredDataStatistics s1 = stats[0][i];
		StoredDataStatistics s2 = stats[1][i];
		StoredDataStatistics s3 = stats[2][i];

		if (s1.getN() == 0)
			return new double[4];

		DescriptiveStatistics d = s1.getStatistics();
		double median = 0;
		Plot2 plot = null;
		String title = null;

		if (showFilterScoreHistograms)
		{
			median = d.getPercentile(50);
			String label = String.format("n = %d. Median = %s nm", s1.getN(), Utils.rounded(median));
			int id = Utils.showHistogram(TITLE, s1, xLabel, filterCriteria[i].minBinWidth,
					(filterCriteria[i].restrictRange) ? 1 : 0, 0, label);
			if (id == 0)
			{
				IJ.log("Failed to show the histogram: " + xLabel);
				return new double[4];
			}

			if (Utils.isNewWindow())
				wo.add(id);

			title = WindowManager.getImage(id).getTitle();

			// Reverse engineer the histogram settings
			plot = Utils.plot;
			double[] xValues = Utils.xValues;
			int bins = xValues.length;
			double yMin = xValues[0];
			double binSize = xValues[1] - xValues[0];
			double yMax = xValues[0] + (bins - 1) * binSize;

			if (s2.getN() > 0)
			{
				double[] values = s2.getValues();
				double[][] hist = Utils.calcHistogram(values, yMin, yMax, bins);

				if (hist[0].length > 0)
				{
					plot.setColor(Color.red);
					plot.addPoints(hist[0], hist[1], Plot2.BAR);
					Utils.display(title, plot);
				}
			}

			if (s3.getN() > 0)
			{
				double[] values = s3.getValues();
				double[][] hist = Utils.calcHistogram(values, yMin, yMax, bins);

				if (hist[0].length > 0)
				{
					plot.setColor(Color.blue);
					plot.addPoints(hist[0], hist[1], Plot2.BAR);
					Utils.display(title, plot);
				}
			}
		}

		// Do cumulative histogram
		double[][] h1 = Maths.cumulativeHistogram(s1.getValues(), true);
		double[][] h2 = Maths.cumulativeHistogram(s2.getValues(), true);
		double[][] h3 = Maths.cumulativeHistogram(s3.getValues(), true);

		if (showFilterScoreHistograms)
		{
			title = TITLE + " Cumul " + xLabel;
			plot = new Plot2(title, xLabel, "Frequency");
			// Find limits
			double[] xlimit = Maths.limits(h1[0]);
			xlimit = Maths.limits(xlimit, h2[0]);
			xlimit = Maths.limits(xlimit, h3[0]);
			// Restrict using the inter-quartile range 
			if (filterCriteria[i].restrictRange)
			{
				double q1 = d.getPercentile(25);
				double q2 = d.getPercentile(75);
				double iqr = (q2 - q1) * 2.5;
				xlimit[0] = Maths.max(xlimit[0], median - iqr);
				xlimit[1] = Maths.min(xlimit[1], median + iqr);
			}
			plot.setLimits(xlimit[0], xlimit[1], 0, 1.05);
			plot.addPoints(h1[0], h1[1], Plot.LINE);
			plot.setColor(Color.red);
			plot.addPoints(h2[0], h2[1], Plot.LINE);
			plot.setColor(Color.blue);
			plot.addPoints(h3[0], h3[1], Plot.LINE);
		}

		// Determine the maximum difference between the TP and FP
		double maxx1 = 0;
		double maxx2 = 0;
		double max1 = 0;
		double max2 = 0;

		// We cannot compute the delta histogram, or use percentiles
		if (s2.getN() == 0)
		{
			upper = UpperLimit.ZERO;
			lower = LowerLimit.ZERO;
		}

		final boolean requireLabel = (showFilterScoreHistograms && filterCriteria[i].requireLabel);
		if (requireLabel || upper.requiresDeltaHistogram() || lower.requiresDeltaHistogram())
		{
			if (s2.getN() != 0 && s3.getN() != 0)
			{
				LinearInterpolator li = new LinearInterpolator();
				PolynomialSplineFunction f1 = li.interpolate(h2[0], h2[1]);
				PolynomialSplineFunction f2 = li.interpolate(h3[0], h3[1]);
				for (double x : h1[0])
				{
					if (x < h2[0][0] || x < h3[0][0])
						continue;
					try
					{
						double v1 = f1.value(x);
						double v2 = f2.value(x);
						double diff = v2 - v1;
						if (diff > 0)
						{
							if (max1 < diff)
							{
								max1 = diff;
								maxx1 = x;
							}
						}
						else
						{
							if (max2 > diff)
							{
								max2 = diff;
								maxx2 = x;
							}
						}
					}
					catch (OutOfRangeException e)
					{
						// Because we reached the end
						break;
					}
				}
			}
			else
			{
				// Switch to percentiles if we have no delta histogram
				if (upper.requiresDeltaHistogram())
					upper = UpperLimit.NINETY_NINE_PERCENT;
				if (lower.requiresDeltaHistogram())
					lower = LowerLimit.ONE_PERCENT;
			}

			//			System.out.printf("Bounds %s : %s, pos %s, neg %s, %s\n", xLabel, Utils.rounded(getPercentile(h2, 0.01)),
			//					Utils.rounded(maxx1), Utils.rounded(maxx2), Utils.rounded(getPercentile(h1, 0.99)));
		}

		if (showFilterScoreHistograms)
		{
			// We use bins=1 on charts where we do not need a label
			if (requireLabel)
			{
				String label = String.format("Max+ %s @ %s, Max- %s @ %s", Utils.rounded(max1), Utils.rounded(maxx1),
						Utils.rounded(max2), Utils.rounded(maxx2));
				plot.setColor(Color.black);
				plot.addLabel(0, 0, label);
			}
			PlotWindow pw = Utils.display(title, plot);
			if (Utils.isNewWindow())
				wo.add(pw.getImagePlus().getID());
		}

		// Now compute the bounds using the desired limit
		double l, u;
		switch (lower)
		{
			case ONE_PERCENT:
				l = getPercentile(h2, 0.01);
				break;
			case MAX_NEGATIVE_CUMUL_DELTA:
				l = maxx2;
				break;
			case ZERO:
				l = 0;
				break;
			case HALF_MAX_JACCARD_VALUE:
				l = getValue(metric, j, maxJ * 0.5);
				break;
			default:
				throw new RuntimeException("Missing lower limit method");
		}
		switch (upper)
		{
			case MAX_POSITIVE_CUMUL_DELTA:
				u = maxx1;
				break;
			case NINETY_NINE_PERCENT:
				u = getPercentile(h2, 0.99);
				break;
			case NINETY_NINE_NINE_PERCENT:
				u = getPercentile(h2, 0.999);
				break;
			case ZERO:
				u = 0;
				break;
			case MAX_JACCARD2:
				u = getValue(metric, j, maxJ) * 2;
				//System.out.printf("MaxJ = %.4f @ %.3f\n", maxJ, u / 2);
				break;
			default:
				throw new RuntimeException("Missing upper limit method");
		}
		double min = getPercentile(h1, 0);
		double max = getPercentile(h1, 1);
		return new double[] { l, u, min, max };
	}

	/**
	 * @param h
	 *            The cumulative histogram
	 * @param p
	 *            The fraction
	 * @return The value for the given fraction
	 */
	private double getPercentile(double[][] h, double p)
	{
		double[] x = h[0];
		double[] y = h[1];
		return getValue(x, y, p);
	}

	/**
	 * Gets the value from x corresponding to the value p in the y values.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param p
	 *            the p
	 * @return the value
	 */
	private double getValue(double[] x, double[] y, double p)
	{
		for (int i = 0; i < x.length; i++)
		{
			if (y[i] >= p)
			{
				if (i == 0 || y[i] == p)
					return x[i];
				// Interpolation
				double delta = (p - y[i - 1]) / (y[i] - y[i - 1]);
				return x[i - 1] + delta * (x[i] - x[i - 1]);
			}
		}
		return x[x.length - 1];
	}

	private static void add(StringBuilder sb, String value)
	{
		sb.append('\t').append(value);
	}

	private static void add(StringBuilder sb, int value)
	{
		sb.append('\t').append(value);
	}

	private static void add(StringBuilder sb, double value)
	{
		add(sb, Utils.rounded(value));
	}

	private static void addCount(StringBuilder sb, double value)
	{
		// Check if the double holds an integer count
		if ((int) value == value)
		{
			sb.append('\t').append((int) value);
		}
		else
		{
			// Otherwise add the counts using at least 2 dp
			if (value > 100)
				sb.append('\t').append(IJ.d2s(value));
			else
				add(sb, Utils.rounded(value));
		}
	}

	private void createTable()
	{
		if (summaryTable == null || !summaryTable.isVisible())
		{
			summaryTable = new TextWindow(TITLE, createHeader(false), "", 1000, 300);
			summaryTable.setVisible(true);
		}
	}

	private String createHeader(boolean extraRecall)
	{
		StringBuilder sb = new StringBuilder(
				"Frames\tW\tH\tMolecules\tDensity (um^-2)\tN\ts (nm)\ta (nm)\tDepth (nm)\tFixed\tGain\tReadNoise (ADUs)\tB (photons)\tb2 (photons)\tSNR\ts (px)\t");
		sb.append("Filter\t");
		sb.append("Spots\t");
		sb.append("nP\t");
		sb.append("nN\t");
		sb.append("fP\t");
		sb.append("fN\t");
		sb.append("Solver\t");
		sb.append("Fitting");

		tablePrefix = sb.toString();

		sb.append('\t');
		sb.append("% nP\t");
		sb.append("% nN\t");
		sb.append("Total\t");
		sb.append("cTP\t");
		sb.append("cFP\t");

		sb.append("cRecall\t");
		sb.append("cPrecision\t");
		sb.append("cF1\t");
		sb.append("cJaccard\t");

		sb.append("Fail cTP\t");
		sb.append("Fail cFP\t");
		sb.append("TP\t");
		sb.append("FP\t");
		sb.append("Recall\t");
		sb.append("Precision\t");
		sb.append("F1\t");
		sb.append("Jaccard\t");

		sb.append("pF1\t");
		sb.append("pJaccard\t");

		sb.append("Med.Distance (nm)\t");
		sb.append("Med.Depth (nm)\t");
		sb.append("Correlation\t");
		sb.append("Ranked\t");
		sb.append("Slope\t");

		createFilterCriteria();
		for (FilterCriteria f : filterCriteria)
			sb.append(f.name).append('\t');

		sb.append("Run time");
		return sb.toString();
	}

	private double getSa()
	{
		final double sa = PSFCalculator.squarePixelAdjustment(simulationParameters.s, simulationParameters.a);
		return sa;
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
		final FitConfiguration pFitConfig = pConfig.getFitConfiguration();

		pFitConfig.setPSF(pFitConfig.getPSF());
		pFitConfig.setFitSolverSettings(pFitConfig.getFitSolverSettings());
		pFitConfig.setFilterSettings(pFitConfig.getFilterSettings());

		// Set the fit engine settings manually to avoid merging all child settings
		pConfig.setFitting(config.getFitting());
		pConfig.setIncludeNeighbours(config.isIncludeNeighbours());
		pConfig.setNeighbourHeightThreshold(config.getNeighbourHeightThreshold());
		pConfig.setDuplicateDistance(config.getDuplicateDistance());

		if (computeDoublets)
		{
			//config.setComputeResiduals(true);
			pConfig.setResidualsThreshold(0);
			pFitConfig.setComputeResiduals(true);
		}
		else
		{
			pConfig.setResidualsThreshold(1);
			pFitConfig.setComputeResiduals(false);
		}

		// We used simple filtering. 
		pFitConfig.setSmartFilter(false);

		return true;
	}

	public void itemStateChanged(ItemEvent e)
	{
		if (e.getSource() instanceof Checkbox)
		{
			Checkbox checkbox = (Checkbox) e.getSource();

			int failLimit;
			boolean includeNeighbours;
			double neighbourHeightThrehsold;
			boolean computeDoublets;
			MultiPathFilter myMultiFilter;

			if (checkbox.getState())
			{
				FitEngineConfiguration tmp = new FitEngineConfiguration();
				FitConfiguration tmpFitConfig = tmp.getFitConfiguration();
				tmpFitConfig.setComputeResiduals(true); // Collect residuals threshold
				if (BenchmarkFilterAnalysis.updateConfiguration(tmp, false))
				{
					failLimit = tmp.getFailuresLimit();
					includeNeighbours = tmp.isIncludeNeighbours();
					neighbourHeightThrehsold = tmp.getNeighbourHeightThreshold();
					computeDoublets = tmp.getResidualsThreshold() < 1;

					final DirectFilter primaryFilter = tmpFitConfig.getSmartFilter();
					final double residualsThreshold = tmp.getResidualsThreshold();
					myMultiFilter = new MultiPathFilter(primaryFilter, minimalFilter, residualsThreshold);
				}
				else
				{
					IJ.log("Failed to update settings using the filter analysis");
					checkbox.setState(false);
					return;
				}
			}
			else
			{
				failLimit = config.getFailuresLimit();
				includeNeighbours = config.isIncludeNeighbours();
				neighbourHeightThrehsold = config.getNeighbourHeightThreshold();
				computeDoublets = BenchmarkSpotFit.computeDoublets;
				myMultiFilter = multiFilter;
			}

			// Update the dialog
			taFilterXml.setText(myMultiFilter.toXML());
			textFailLimit.setText("" + failLimit);
			cbIncludeNeighbours.setState(includeNeighbours);
			textNeighbourHeight.setText(Utils.rounded(neighbourHeightThrehsold));
			cbComputeDoublets.setState(computeDoublets);
		}
	}

	/**
	 * Run the analysis non-interactively using the given filter settings.
	 *
	 * @param filter
	 *            the filter
	 * @param residualsThreshold
	 *            the residuals threshold
	 * @param failuresLimit
	 *            the failures limit
	 * @param duplicateDistance
	 *            the duplicate distance
	 */
	public void run(DirectFilter filter, double residualsThreshold, int failuresLimit, double duplicateDistance)
	{
		multiFilter = new MultiPathFilter(filter, minimalFilter, residualsThreshold);
		config.setFailuresLimit(failuresLimit);
		config.setDuplicateDistance(duplicateDistance);

		clearFitResults();

		silent = true;
		if (!initialise())
			return;

		run();
	}

	/**
	 * Reset the multi path filter and non-filter parameters. This may have been updated when copying benchmark filter
	 * settings or by the user within the dialog.
	 *
	 * @return true, if a reset was required
	 */
	public boolean resetMultiPathFilter()
	{
		if (defaultMultiFilter.equals(multiFilter))
		{
			if (equals(createParameters(), defaultParameters))
				return false;
		}
		multiFilter = defaultMultiFilter;
		config.setFailuresLimit((int) defaultParameters[0]);
		// Note we are not resetting the residuals threshold in the config.
		// Only the threshold in the multi-filter matters.
		config.setDuplicateDistance(defaultParameters[1]);
		return true;
	}

	private boolean equals(double[] currentParameters, double[] previousParameters)
	{
		for (int i = 0; i < previousParameters.length; i++)
			if (previousParameters[i] != currentParameters[i])
				return false;
		return true;
	}

	private static double[] createParameters()
	{
		return new double[] { config.getFailuresLimit(), config.getDuplicateDistance() };
	}
}
