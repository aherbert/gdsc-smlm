package gdsc.smlm.ij.plugins;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import org.apache.commons.math3.util.FastMath;

import gdsc.core.ij.Utils;
import gdsc.core.match.BasePoint;
import gdsc.core.match.Coordinate;
import gdsc.core.match.FractionClassificationResult;
import gdsc.core.match.MatchCalculator;
import gdsc.core.match.PointPair;
import gdsc.core.utils.Maths;
import gdsc.core.utils.RampedScore;
import gdsc.core.utils.StoredDataStatistics;

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

import gdsc.smlm.engine.DataFilter;
import gdsc.smlm.engine.DataFilterType;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.filters.MaximaSpotFilter;
import gdsc.smlm.filters.Spot;
import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.ImageConverter;
import gdsc.smlm.results.MemoryPeakResults;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
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
	static
	{
		fitConfig = new FitConfiguration();
		config = new FitEngineConfiguration(fitConfig);
		config.setSearch(1);
	}

	private static double sAnalysisBorder = 0;
	private static boolean multipleMatches = true;
	private int analysisBorder;
	private static double upperDistance = 1.5;
	private static double lowerDistance = 1.5;
	private static boolean relativeDistances = true;
	private double matchDistance;
	private double lowerMatchDistance;
	private float pixelOffset;
	private static double recallFraction = 100;
	private static boolean showPlot = true;
	private static boolean showFailuresPlot = true;
	private static boolean sDebug = false;
	private boolean extraOptions, debug = false;
	private long time = 0;

	private static int id = 1;

	private static TextWindow summaryTable = null;

	private ImagePlus imp;
	private MemoryPeakResults results;
	private CreateData.SimulationParameters simulationParameters;

	private static HashMap<Integer, ArrayList<Coordinate>> actualCoordinates = null;
	private static int lastId = -1;
	private static boolean lastRelativeDistances = false;

	// Used by the Benchmark Spot Fit plugin
	static int simulationId = 0;
	static int filterResultsId = 0;
	static HashMap<Integer, FilterResult> filterResults = null;
	static MaximaSpotFilter spotFilter = null;

	private int idCount = 0;
	private int[] idList = new int[4];

	public class ScoredSpot implements Comparable<ScoredSpot>
	{
		final boolean match;
		double[] distances, scores;
		double score;
		final Spot spot;
		int fails;

		public ScoredSpot(boolean match, double d, double score, Spot spot)
		{
			this.match = match;
			this.spot = spot;
			this.fails = 0;
			add(d, score);
		}

		public ScoredSpot(boolean match, Spot spot)
		{
			this.match = match;
			this.spot = spot;
			this.fails = 0;
		}

		public ScoredSpot(boolean match, Spot spot, int fails)
		{
			this.match = match;
			this.spot = spot;
			this.fails = fails;
		}

		/**
		 * Adds the result to the scored spot
		 *
		 * @param d
		 *            the d
		 * @param score
		 *            the score
		 */
		void add(double d, double score)
		{
			if (distances == null)
			{
				distances = new double[] { d };
				scores = new double[] { score };
				this.score = score;
			}
			else
			{
				final int size = distances.length;
				distances = Arrays.copyOf(distances, size + 1);
				scores = Arrays.copyOf(scores, size + 1);
				distances[size] = d;
				scores[size] = score;
				this.score += score;
			}
		}

		public double getDistance(int i)
		{
			return (distances != null && i < distances.length) ? distances[i] : 0;
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
			// XXX - I am not sure what result I would like here.
			// One method allows the anti-score to be restricted to the range 0-1.
			// i.e. a spot not matching anything should get a penalty of 1.
			// The other method says that a spot not matching anything can have a 
			// penalty above 1, if it is further away from multiple actual results
			// than the mid point between the lower and upper match distance.
			
			// The use of partial scoring mimicks using multiple distance 
			// thresholds and taking an average of scores. 
			// For a single experiment at a single distance threshold a spot can 
			// match 0, 1 or more actual results. This would be:
			// 0  = TP=0 ,FP=1
			// 1  = TP=1 ,FP=0
			// 2+ = TP=2+,FP=0 (because it doesn't 'not match' anything)
			// So for now I will use FP (anti-score) in the range 0-1.
			
			// Ensure score + anti-score add up to number number of scores
			//return (distances == null) ? 1 : distances.length - score;
			
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
	}

	public class FilterResult
	{
		final FractionClassificationResult result;
		final ScoredSpot[] spots;

		public FilterResult(FractionClassificationResult result, ScoredSpot[] spots)
		{
			this.result = result;
			this.spots = spots;
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
		final HashMap<Integer, ArrayList<Coordinate>> actualCoordinates;
		List<Coordinate> FP = new ArrayList<Coordinate>();
		final HashMap<Integer, FilterResult> results;
		List<PointPair> matches = new ArrayList<PointPair>();

		float[] data = null;
		long time = 0;

		public Worker(BlockingQueue<Integer> jobs, ImageStack stack, MaximaSpotFilter spotFilter,
				HashMap<Integer, ArrayList<Coordinate>> actualCoordinates)
		{
			this.jobs = jobs;
			this.stack = stack;
			this.spotFilter = spotFilter.clone();
			this.actualCoordinates = actualCoordinates;
			this.results = new HashMap<Integer, FilterResult>();
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
				while (!finished)
				{
					Integer job = jobs.take();
					if (job == null || job.intValue() < 0 || finished)
						break;
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
			// Extract the data
			data = ImageConverter.getData(stack.getPixels(frame), stack.getWidth(), stack.getHeight(), null, data);

			long start = System.nanoTime();
			Spot[] spots = spotFilter.rank(data, stack.getWidth(), stack.getHeight());
			time += System.nanoTime() - start;

			// Score the spots that are matches
			Coordinate[] actual = ResultsMatchCalculator.getCoordinates(actualCoordinates, frame);

			// Remove spots at the border from analysis
			if (analysisBorder > 0)
			{
				final int xlimit = stack.getWidth() - analysisBorder;
				final int ylimit = stack.getHeight() - analysisBorder;
				Coordinate[] actual2 = new Coordinate[actual.length];
				int count = 0;
				for (Coordinate c : actual)
				{
					if (c.getX() < analysisBorder || c.getX() > xlimit || c.getY() < analysisBorder ||
							c.getY() > ylimit)
						continue;
					actual2[count++] = c;
				}
				actual = Arrays.copyOf(actual2, count);
				count = 0;
				Spot[] spots2 = new Spot[spots.length];
				for (Spot s : spots)
				{
					if (s.x < analysisBorder || s.x > xlimit || s.y < analysisBorder || s.y > ylimit)
						continue;
					spots2[count++] = s;
				}
				spots = Arrays.copyOf(spots2, count);
			}

			ScoredSpot[] scoredSpots = new ScoredSpot[spots.length];
			FractionClassificationResult result;

			// Store the count of false positives since the last true positive
			int fails = 0;

			// Use the distance to the true location to score the candidate
			matches.clear();
			RampedScore score = new RampedScore(lowerMatchDistance, matchDistance);

			if (actual.length > 0)
			{
				Coordinate[] predicted = getCoordinates(spots);

				double tp = 0, fp = 0;

				if (!multipleMatches)
				{
					MatchCalculator.analyseResults2D(actual, predicted, matchDistance, null, FP, null, matches);

					// Store the true and false positives. Maintain the original ranked order.
					fp = FP.size();
					for (PointPair pair : matches)
					{
						SpotCoordinate sc = (SpotCoordinate) pair.getPoint2();
						final double d = pair.getXYDistance();
						final double s = RampedScore.flatten(score.score(d), 256);
						scoredSpots[sc.id] = new ScoredSpot(true, d, s, sc.spot);
						// Score partial matches as part true-positive and part false-positive
						tp += s;
						fp += (1 - s);
					}
					for (Coordinate c : FP)
					{
						SpotCoordinate sc = (SpotCoordinate) c;
						scoredSpots[sc.id] = new ScoredSpot(false, sc.spot);
					}
				}
				else
				{
					// Assign each actual result to the closest spot. 
					// This allows a spot to have more than one match (i.e. a doublet) 
					for (int i = 0; i < actual.length; i++)
					{
						double dmin = matchDistance * matchDistance;
						int match = -1;
						final double x = actual[i].getX();
						final double y = actual[i].getY();
						for (int j = 0; j < predicted.length; j++)
						{
							final double dx = (x - predicted[j].getX());
							final double dy = (y - predicted[j].getY());
							final double d2 = dx * dx + dy * dy;
							if (d2 <= dmin)
							{
								dmin = d2;
								match = j;
							}
						}

						if (match != -1)
						{
							final double d = Math.sqrt(dmin);
							final double s = RampedScore.flatten(score.score(d), 256);
							final ScoredSpot ss = scoredSpots[match];
							if (ss == null)
							{
								scoredSpots[match] = new ScoredSpot(true, d, s, spots[match]);
							}
							else
							{
								// This is a second match to the same spot
								ss.add(d, s);
							}
							tp += s;
						}
					}
					for (int j = 0; j < spots.length; j++)
					{
						if (scoredSpots[j] == null)
						{
							scoredSpots[j] = new ScoredSpot(false, spots[j]);
							fp++;
						}
						else
						{
							fp += scoredSpots[j].antiScore();
						}
					}
				}

				result = new FractionClassificationResult(tp, fp, 0, actual.length - tp);

				// Store the number of fails (negatives) before each positive 
				for (int i = 0; i < spots.length; i++)
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
			else
			{
				// All spots are false positives
				result = new FractionClassificationResult(0, spots.length, 0, 0);

				for (int i = 0; i < spots.length; i++)
				{
					scoredSpots[i] = new ScoredSpot(false, spots[i], fails++);
				}
			}

			if (debug)
			{
				System.out.printf("Frame %d : N = %d, TP = %.2f, FP = %.2f, R = %.2f, P = %.2f\n", frame, actual.length,
						result.getTP(), result.getFP(), result.getRecall(), result.getPrecision());
			}

			results.put(frame, new FilterResult(result, scoredSpots));
		}

		private Coordinate[] getCoordinates(Spot[] spots)
		{
			Coordinate[] coords = new Coordinate[spots.length];
			for (int i = 0; i < spots.length; i++)
				coords[i] = new SpotCoordinate(i, spots[i]);
			return coords;
		}

		private class SpotCoordinate extends BasePoint
		{
			int id;
			Spot spot;

			public SpotCoordinate(int id, Spot spot)
			{
				super(spot.x + pixelOffset, spot.y + pixelOffset);
				this.id = id;
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
		results = MemoryPeakResults.getResults(CreateData.CREATE_DATA_IMAGE_TITLE + " (Create Data)");
		if (results == null)
		{
			IJ.error(TITLE, "No benchmark results in memory");
			return;
		}

		if (!showDialog())
			return;

		run();
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		StringBuilder sb = new StringBuilder();
		sb.append("Finds spots in the benchmark image created by CreateData plugin.\n");
		sb.append("PSF width = ").append(Utils.rounded(simulationParameters.s / simulationParameters.a))
				.append(" px\n");
		sb.append("Simulation depth = ").append(Utils.rounded(simulationParameters.depth)).append(" nm");
		if (simulationParameters.fixedDepth)
			sb.append(" (fixed)");
		sb.append("\n \nConfigure the spot filter:");
		gd.addMessage(sb.toString());

		String[] filterTypes = SettingsManager.getNames((Object[]) DataFilterType.values());
		gd.addChoice("Spot_filter_type", filterTypes, filterTypes[config.getDataFilterType().ordinal()]);
		String[] filterNames = SettingsManager.getNames((Object[]) DataFilter.values());
		gd.addChoice("Spot_filter", filterNames, filterNames[config.getDataFilter(0).ordinal()]);

		gd.addCheckbox("Relative_distances", relativeDistances);

		gd.addSlider("Smoothing", 0, 2.5, config.getSmooth(0));
		gd.addSlider("Search_width", 1, 4, config.getSearch());
		gd.addSlider("Border", 0, 5, config.getBorder());

		gd.addMessage("Scoring options:");
		gd.addSlider("Analysis_border", 0, 5, sAnalysisBorder);
		gd.addCheckbox("Multiple_matches", multipleMatches);
		gd.addSlider("Match_distance", 1, 3, upperDistance);
		gd.addSlider("Lower_distance", 1, 3, lowerDistance);
		gd.addSlider("Recall_fraction", 50, 100, recallFraction);
		gd.addCheckbox("Show_plots", showPlot);
		gd.addCheckbox("Show_failures_plots", showFailuresPlot);
		if (extraOptions)
			gd.addCheckbox("Debug", sDebug);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		config.setDataFilterType(gd.getNextChoiceIndex());
		config.setDataFilter(gd.getNextChoiceIndex(), Math.abs(gd.getNextNumber()), 0);
		relativeDistances = gd.getNextBoolean();
		config.setSearch(gd.getNextNumber());
		config.setBorder(gd.getNextNumber());
		sAnalysisBorder = Math.abs(gd.getNextNumber());
		multipleMatches = gd.getNextBoolean();
		upperDistance = Math.abs(gd.getNextNumber());
		lowerDistance = Math.abs(gd.getNextNumber());
		recallFraction = Math.abs(gd.getNextNumber());
		showPlot = gd.getNextBoolean();
		showFailuresPlot = gd.getNextBoolean();
		if (extraOptions)
			debug = sDebug = gd.getNextBoolean();

		if (gd.invalidNumber())
			return false;

		if (lowerDistance > upperDistance)
			lowerDistance = upperDistance;

		GlobalSettings settings = new GlobalSettings();
		settings.setFitEngineConfiguration(config);
		if (!PeakFit.configureDataFilter(settings, null, false))
			return false;

		if (relativeDistances)
		{
			// Convert distance to PSF standard deviation units
			final double sd = simulationParameters.s / simulationParameters.a;
			matchDistance = upperDistance * sd;
			lowerMatchDistance = lowerDistance * sd;
			analysisBorder = (int) (analysisBorder * sd);
			// Add 0.5 offset to centre the spot in the pixel
			pixelOffset = 0.5f;
		}
		else
		{
			matchDistance = upperDistance;
			lowerMatchDistance = lowerDistance;
			analysisBorder = (int) (analysisBorder);
			// Absolute distances in pixels will use integer coordinates so no offset
			pixelOffset = 0;
		}

		return true;
	}

	private void run()
	{
		spotFilter = config.createSpotFilter(relativeDistances);

		// Extract all the results in memory into a list per frame. This can be cached
		if (lastId != simulationParameters.id || lastRelativeDistances != relativeDistances)
		{
			// When using pixel distances then get the coordinates in integers.
			// The Worker adds a pixel offset if appropriate.
			actualCoordinates = ResultsMatchCalculator.getCoordinates(results.getResults(), !relativeDistances);
			lastId = simulationParameters.id;
			lastRelativeDistances = relativeDistances;
		}

		final ImageStack stack = imp.getImageStack();

		// Clear old results to free memory
		if (filterResults != null)
			filterResults.clear();
		filterResults = null;

		// Create a pool of workers
		final int nThreads = Prefs.getThreads();
		BlockingQueue<Integer> jobs = new ArrayBlockingQueue<Integer>(nThreads * 2);
		List<Worker> workers = new LinkedList<Worker>();
		List<Thread> threads = new LinkedList<Thread>();
		for (int i = 0; i < nThreads; i++)
		{
			Worker worker = new Worker(jobs, stack, spotFilter, actualCoordinates);
			Thread t = new Thread(worker);
			workers.add(worker);
			threads.add(t);
			t.start();
		}

		final int totalFrames = stack.getSize();

		// Fit the frames
		final int step = Utils.getProgressInterval(totalFrames);
		for (int i = 1; i <= totalFrames; i++)
		{
			// TODO : Should we only process the frame if there were simulated spots?

			put(jobs, i);
			if (i % step == 0)
			{
				if (Utils.showStatus("Frame: " + i + " / " + totalFrames))
					IJ.showProgress(i, totalFrames);
			}
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
		IJ.showStatus("Collecting results ...");

		filterResultsId++;
		simulationId = simulationParameters.id;
		filterResults = new HashMap<Integer, FilterResult>();
		for (Worker w : workers)
		{
			time += w.time;
			filterResults.putAll(w.results);
		}

		double[][] h = histogramFailures(filterResults);

		// Show a table of the results
		summariseResults(filterResults, spotFilter, h);

		if (idCount > 0)
		{
			idList = Arrays.copyOf(idList, idCount);
			new WindowOrganiser().tileWindows(idList);
		}

		// Debugging the matches
		if (debug)
			addSpotsToMemory(filterResults);

		IJ.showStatus("");
	}

	/**
	 * Add all the true-positives to memory as a new results set
	 * 
	 * @param filterResults
	 */
	void addSpotsToMemory(HashMap<Integer, FilterResult> filterResults)
	{
		MemoryPeakResults results = new MemoryPeakResults();
		results.setName(TITLE + " TP " + id++);
		for (Entry<Integer, FilterResult> result : filterResults.entrySet())
		{
			int peak = result.getKey();
			for (ScoredSpot spot : result.getValue().spots)
			{
				if (spot.match)
				{
					final float[] params = new float[] { 0, spot.spot.intensity, 0, spot.spot.x, spot.spot.y, 0, 0 };
					results.add(peak, spot.spot.x, spot.spot.y, spot.spot.intensity, 0d, 0f, params, null);
				}
			}
		}
		MemoryPeakResults.addResults(results);
	}

	/**
	 * Histogram the number of negatives preceeding each positive
	 * 
	 * @param filterResults
	 */
	private double[][] histogramFailures(HashMap<Integer, FilterResult> filterResults)
	{
		StoredDataStatistics stats = new StoredDataStatistics();
		for (Entry<Integer, FilterResult> result : filterResults.entrySet())
		{
			for (ScoredSpot spot : result.getValue().spots)
			{
				if (spot.match)
				{
					stats.add(spot.fails);
				}
			}
		}
		String xTitle = "Failures";

		double[][] h = Maths.cumulativeHistogram(stats.getValues(), true);

		if (showFailuresPlot)
		{
			final int id = Utils.showHistogram(TITLE, stats, xTitle, 1, 0, 0);
			if (Utils.isNewWindow())
				idList[idCount++] = id;

			String title = TITLE + " " + xTitle + " Cumulative";
			Plot2 plot = new Plot2(title, xTitle, "Frequency");
			double xMax = h[0][h[0].length - 1] + 1;
			double xPadding = 0.05 * (xMax - h[0][0]);
			plot.setLimits(h[0][0] - xPadding, xMax, 0, 1.05);
			plot.setColor(Color.blue);
			plot.addPoints(h[0], h[1], Plot2.BAR);
			PlotWindow pw = Utils.display(title, plot);
			if (Utils.isNewWindow())
				idList[idCount++] = pw.getImagePlus().getID();
		}

		return h;
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

	private void summariseResults(HashMap<Integer, FilterResult> filterResults, MaximaSpotFilter spotFilter,
			double[][] cumul)
	{
		createTable();

		// Create the overall match score
		double tp = 0, fp = 0, fn = 0;
		ArrayList<ScoredSpot> allSpots = new ArrayList<BenchmarkSpotFilter.ScoredSpot>();
		for (FilterResult result : filterResults.values())
		{
			tp += result.result.getTP();
			fp += result.result.getFP();
			fn += result.result.getFN();
			allSpots.addAll(Arrays.asList(result.spots));
		}
		FractionClassificationResult allResult = new FractionClassificationResult(tp, fp, 0, fn);
		// The number of actual results
		final int n = (int) (tp + fn);

		StringBuilder sb = new StringBuilder();

		double signal = (simulationParameters.minSignal + simulationParameters.maxSignal) * 0.5;

		// Create the benchmark settings and the fitting settings
		sb.append(imp.getStackSize()).append("\t");
		final int w = imp.getWidth() - 2 * analysisBorder;
		final int h = imp.getHeight() - 2 * analysisBorder;
		sb.append(w).append("\t");
		sb.append(h).append("\t");
		sb.append(n).append("\t");
		double density = ((double) n / imp.getStackSize()) / (w * h) /
				(simulationParameters.a * simulationParameters.a / 1e6);
		sb.append(Utils.rounded(density)).append("\t");
		sb.append(Utils.rounded(signal)).append("\t");
		sb.append(Utils.rounded(simulationParameters.s)).append("\t");
		sb.append(Utils.rounded(simulationParameters.a)).append("\t");
		sb.append(Utils.rounded(simulationParameters.depth)).append("\t");
		sb.append(simulationParameters.fixedDepth).append("\t");
		sb.append(Utils.rounded(simulationParameters.gain)).append("\t");
		sb.append(Utils.rounded(simulationParameters.readNoise)).append("\t");
		sb.append(Utils.rounded(simulationParameters.b)).append("\t");
		sb.append(Utils.rounded(simulationParameters.b2)).append("\t");

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

		sb.append(Utils.rounded(signal / Math.sqrt(noise))).append("\t");
		sb.append(Utils.rounded(simulationParameters.s / simulationParameters.a)).append("\t");
		sb.append(config.getDataFilterType()).append("\t");
		//sb.append(spotFilter.getName()).append("\t");
		sb.append(spotFilter.getSearch()).append("\t");
		sb.append(spotFilter.getBorder()).append("\t");
		sb.append(Utils.rounded(spotFilter.getSpread())).append("\t");
		sb.append(config.getDataFilter(0)).append("\t");
		sb.append(Utils.rounded(config.getSmooth(0))).append("\t");
		sb.append(spotFilter.getDescription()).append("\t");
		sb.append(analysisBorder).append("\t");
		sb.append(multipleMatches).append("\t");
		sb.append(Utils.rounded(lowerMatchDistance)).append("\t");
		sb.append(Utils.rounded(matchDistance)).append("\t");

		// Add the results

		addResult(sb, allResult);

		// Rank the scored spots by intensity
		Collections.sort(allSpots);

		// Produce Recall, Precision, Jaccard for each cut of the spot candidates
		double[] r = new double[allSpots.size() + 1];
		double[] p = new double[r.length];
		double[] j = new double[r.length];
		double[] truePositives = new double[r.length];
		double[] falsePositives = new double[r.length];
		double[] rank = new double[r.length];
		// Note: fn = n - tp
		tp = fp = 0;
		int i = 1;
		p[0] = 1;
		for (ScoredSpot s : allSpots)
		{
			if (s.match)
			{
				// Score partial matches as part true-positive and part false-positive.
				// TP can be above 1 if we are allowing multiple matches.
				tp += s.getScore();
				fp += s.antiScore();
			}
			else
				fp++;
			r[i] = (double) tp / n;
			p[i] = (double) tp / (tp + fp);
			j[i] = (double) tp / (fp + n); // (tp+fp+fn) == (fp+n) since tp+fn=n;
			truePositives[i] = tp;
			falsePositives[i] = fp;
			rank[i] = i;
			i++;
		}

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
				n - truePositives[fractionIndex]));

		// Output the match results at the maximum jaccard score
		int maxIndex = 0;
		for (int ii = 1; ii < r.length; ii++)
		{
			if (j[maxIndex] < j[ii])
				maxIndex = ii;
		}
		addResult(sb, new FractionClassificationResult(truePositives[maxIndex], falsePositives[maxIndex], 0,
				n - truePositives[maxIndex]));

		sb.append(Utils.rounded(time / 1e6));

		// Calculate AUC (Average precision == Area Under Precision-Recall curve)
		final double auc = auc(p, r);
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
		final double auc2 = auc(maxp, r);

		sb.append("\t").append(Utils.rounded(auc));
		sb.append("\t").append(Utils.rounded(auc2));

		// Output the number of fit failures that must be processed to capture fractions of the true positives
		sb.append("\t").append(Utils.rounded(getFailures(cumul, 0.95)));
		sb.append("\t").append(Utils.rounded(getFailures(cumul, 0.99)));
		sb.append("\t").append(Utils.rounded(cumul[0][cumul[0].length - 1]));

		if (showPlot)
		{
			String title = TITLE + " Performance";
			Plot2 plot = new Plot2(title, "Spot Rank", "");
			plot.setLimits(0, rank.length, 0, 1.05);
			plot.setColor(Color.blue);
			plot.addPoints(rank, p, Plot2.LINE);
			//plot.addPoints(rank, maxp, Plot2.DOT);
			plot.setColor(Color.red);
			plot.addPoints(rank, r, Plot2.LINE);
			plot.setColor(Color.black);
			plot.addPoints(rank, j, Plot2.LINE);
			plot.setColor(Color.magenta);
			plot.drawLine(rank[fractionIndex], 0, rank[fractionIndex],
					Maths.max(p[fractionIndex], r[fractionIndex], j[fractionIndex]));
			plot.setColor(Color.pink);
			plot.drawLine(rank[maxIndex], 0, rank[maxIndex], Maths.max(p[maxIndex], r[maxIndex], j[maxIndex]));
			plot.setColor(Color.black);
			plot.addLabel(0, 0, "Precision=Blue, Recall=Red, Jaccard=Black");

			PlotWindow pw = Utils.display(title, plot);
			if (Utils.isNewWindow())
				idList[idCount++] = pw.getImagePlus().getID();

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
				idList[idCount++] = pw2.getImagePlus().getID();
		}

		summaryTable.append(sb.toString());
	}

	private double getFailures(double[][] cumul, double fraction)
	{
		int i = 0;
		final double[] sum = cumul[1];
		while (sum[i] < fraction && i < sum.length)
			i++;
		return i;
	}

	private void addResult(StringBuilder sb, FractionClassificationResult matchResult)
	{
		addCount(sb, matchResult.getTP());
		addCount(sb, matchResult.getFP());
		sb.append(Utils.rounded(matchResult.getRecall())).append("\t");
		sb.append(Utils.rounded(matchResult.getPrecision())).append("\t");
		sb.append(Utils.rounded(matchResult.getJaccard())).append("\t");
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
		sb.append("\t");
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
		sb.append("Type\tSearch\tBorder\tWidth\tFilter\tParam\tDescription\tA.Border\tMulti\tlower d\td\t");
		sb.append("TP\tFP\tRecall\tPrecision\tJaccard\t");
		sb.append("TP\tFP\tRecall\tPrecision\tJaccard\t");
		sb.append("TP\tFP\tRecall\tPrecision\tJaccard\t");
		sb.append("Time (ms)\t");
		sb.append("AUC\tAUC2\t");
		sb.append("Fail95\tFail99\tFail100");
		return sb.toString();
	}

	/**
	 * Calculates an estimate of the area under the precision-recall curve.
	 * <p>
	 * The estimate is computed using integration above the recall limit. Below the limit a simple linear interpolation
	 * is used from the given point to precision 1 at recall 0. This avoids noise in the lower recall section of the
	 * curve.
	 * <p>
	 * If no recall values are above the limit then the full integration is performed.
	 *
	 * @param precision
	 * @param recall
	 * @param recallLimit
	 *            Set to 0 to compute the full area.
	 * @return Area under the PR curve
	 */
	public static double auc(double[] precision, double[] recall, double recallLimit)
	{
		if (precision == null || recall == null)
			return 0;

		double area = 0.0;
		int k;

		if (recallLimit > 0)
		{
			// Move from high to low recall and find the first point below the limit
			k = recall.length - 1;
			while (k > 0 && recall[k] > recallLimit)
				k--;

			if (k > 0)
			{
				// Find the first point where precision was not 1
				int kk = 0;
				while (precision[kk + 1] == 1)
					kk++;

				// Full precision of 1 up to point kk
				area += (recall[kk] - recall[0]);

				// Interpolate from precision at kk to k
				area += (precision[k] + precision[kk]) * 0.5 * (recall[k] - recall[kk]);

				// Increment to start the remaining integral
				k++;
			}
		}
		else
		{
			// Complete integration from start
			k = 0;
		}

		// Integrate the rest
		double prevR = 0;
		double prevP = 1;
		if (recall[0] == 0)
			k++;

		for (; k < precision.length; k++)
		{
			final double delta = recall[k] - prevR;
			if (precision[k] == prevP)
				area += prevP * delta;
			else
				// Interpolate
				area += (precision[k] + prevP) * 0.5 * delta;
			prevR = recall[k];
			prevP = precision[k];
		}
		return area;
	}

	/**
	 * Calculates an estimate of the area under the precision-recall curve.
	 * <p>
	 * Assumes the first values in the two arrays are precision 1 at recall 0.
	 *
	 * @param precision
	 * @param recall
	 * @return Area under the PR curve
	 */
	private static double auc(double[] precision, double[] recall)
	{
		double area = 0.0;

		double prevR = 0;
		double prevP = 1;

		for (int k = 1; k < precision.length; k++)
		{
			final double delta = recall[k] - prevR;
			if (precision[k] == prevP)
				area += prevP * delta;
			else
				// Interpolate
				area += (precision[k] + prevP) * 0.5 * delta;
			prevR = recall[k];
			prevP = precision[k];
		}
		return area;
	}
}
