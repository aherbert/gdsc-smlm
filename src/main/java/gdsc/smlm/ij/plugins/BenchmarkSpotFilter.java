package gdsc.smlm.ij.plugins;

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
import gdsc.smlm.engine.FitEngine;
import gdsc.smlm.filters.MaximaSpotFilter;
import gdsc.smlm.filters.Spot;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.ImageConverter;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.match.BasePoint;
import gdsc.smlm.results.match.Coordinate;
import gdsc.smlm.results.match.MatchCalculator;
import gdsc.smlm.results.match.MatchResult;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.text.TextWindow;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

/**
 * Filters the benchmark spot image created by CreateData plugin to identify candidates and then assess the filter.
 */
public class BenchmarkSpotFilter implements PlugIn
{
	private static final String TITLE = "Benchmark Spot Filter";
	private static final String[] filterNames;
	private static final DataFilter[] filters;
	static
	{
		filters = DataFilter.values();
		filterNames = SettingsManager.getNames((Object[]) filters);
	}

	private static int search = 2;
	private static int border = 2;
	private static int dataFilter = 0;
	private static double smoothing = 1;
	private static boolean differenceFilter = false;
	private static int dataFilter2 = 0;
	private static double smoothing2 = 3;
	private static double distance = 1.5;
	private static boolean sDebug = false;
	private boolean extraOptions, debug = false;
	private long time = 0;

	private static TextWindow summaryTable = null;

	private ImagePlus imp;
	private MemoryPeakResults results;
	private CreateData.SimulationParameters simulationParameters;

	private class ScoredSpot implements Comparable<ScoredSpot>
	{
		final boolean match;
		final Spot spot;

		public ScoredSpot(boolean match, Spot spot)
		{
			this.match = match;
			this.spot = spot;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(ScoredSpot o)
		{
			if (spot.intensity > o.spot.intensity)
				return -1;
			if (spot.intensity < o.spot.intensity)
				return 1;
			return 0;
		}
	}

	private class FilterResult
	{
		final MatchResult result;
		final ScoredSpot[] spots;

		public FilterResult(MatchResult result, ScoredSpot[] spots)
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
		List<Coordinate> TP = new ArrayList<Coordinate>();
		List<Coordinate> FP = new ArrayList<Coordinate>();
		final HashMap<Integer, FilterResult> results;

		float[] data = null;
		long time = 0;

		public Worker(BlockingQueue<Integer> jobs, ImageStack stack, MaximaSpotFilter spotFilter,
				HashMap<Integer, ArrayList<Coordinate>> actualCoordinates, HashMap<Integer, FilterResult> results)
		{
			this.jobs = jobs;
			this.stack = stack;
			this.spotFilter = (MaximaSpotFilter) spotFilter.clone();
			this.actualCoordinates = actualCoordinates;
			this.results = results;
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

			ScoredSpot[] scoredSpots = new ScoredSpot[spots.length];
			MatchResult result;

			if (actual.length > 0)
			{
				Coordinate[] predicted = getCoordinates(spots);
				TP.clear();
				FP.clear();

				result = MatchCalculator.analyseResults2D(actual, predicted, distance, TP, FP, null, null);

				// Store the true and false positives
				int i = 0;
				for (Coordinate c : TP)
				{
					scoredSpots[i++] = new ScoredSpot(true, ((SpotCoordinate) c).spot);
				}
				for (Coordinate c : FP)
				{
					scoredSpots[i++] = new ScoredSpot(false, ((SpotCoordinate) c).spot);
				}
			}
			else
			{
				// All spots are false positives
				result = new MatchResult(0, spots.length, 0, 0);

				for (int i = 0; i < spots.length; i++)
				{
					scoredSpots[i] = new ScoredSpot(false, spots[i]);
				}
			}

			if (debug)
			{
				System.out.printf("Frame %d : N = %d, TP = %d, FP = %d, R = %.2f, P = %.2f\n", frame,
						result.getNumberActual(), result.getTruePositives(), result.getFalsePositives(),
						result.getRecall(), result.getPrecision());
			}

			results.put(frame, new FilterResult(result, scoredSpots));
		}

		private Coordinate[] getCoordinates(Spot[] spots)
		{
			Coordinate[] coords = new Coordinate[spots.length];
			for (int i = 0; i < spots.length; i++)
				coords[i] = new SpotCoordinate(spots[i]);
			return coords;
		}

		private class SpotCoordinate extends BasePoint
		{
			Spot spot;

			public SpotCoordinate(Spot spot)
			{
				super(spot.x, spot.y);
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
		extraOptions = Utils.isExtraOptions();

		simulationParameters = CreateData.simulationParameters;
		if (simulationParameters == null)
		{
			IJ.error(TITLE, "No benchmark spot parameters in memory");
			return;
		}
		imp = WindowManager.getImage(CreateData.CREATE_DATA_IMAGE_TITLE);
		if (imp == null)
		{
			IJ.error(TITLE, "No benchmark image");
			return;
		}
		results = MemoryPeakResults.getResults(CreateData.CREATE_DATA_IMAGE_TITLE + " (Create Data)");
		if (results == null)
		{
			IJ.error(TITLE, "No benchmark resuls in memory");
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

		gd.addMessage(String.format("Fits the benchmark image created by CreateData plugin.\nPSF width = %s",
				Utils.rounded(simulationParameters.s / simulationParameters.a)));

		gd.addSlider("Search", 0, 5, search);
		gd.addSlider("Border", 0, 5, border);
		gd.addChoice("Spot_filter", filterNames, filterNames[dataFilter]);
		gd.addSlider("Smoothing", 0, 4.5, smoothing);
		gd.addCheckbox("Difference_filter", differenceFilter);
		gd.addChoice("Spot_filter2", filterNames, filterNames[dataFilter2]);
		gd.addSlider("Smoothing2", 1.5, 6, smoothing2);
		gd.addSlider("Match_distance", 1, 3, distance);
		if (extraOptions)
			gd.addCheckbox("Debug", sDebug);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		search = Math.abs((int) gd.getNextNumber());
		border = Math.abs((int) gd.getNextNumber());
		dataFilter = gd.getNextChoiceIndex();
		smoothing = Math.abs(gd.getNextNumber());
		differenceFilter = gd.getNextBoolean();
		dataFilter2 = gd.getNextChoiceIndex();
		smoothing2 = Math.abs(gd.getNextNumber());
		distance = Math.abs(gd.getNextNumber());
		if (extraOptions)
			debug = sDebug = gd.getNextBoolean();

		if (gd.invalidNumber())
			return false;

		if (search < 1)
			search = 1;

		return true;
	}

	private void run()
	{
		MaximaSpotFilter spotFilter = FitEngine.createSpotFilter(search, border, filters[dataFilter], smoothing,
				differenceFilter, filters[dataFilter2], smoothing2);

		// Extract all the results in memory into a list per frame
		HashMap<Integer, ArrayList<Coordinate>> actualCoordinates = ResultsMatchCalculator.getCoordinates(results
				.getResults());

		HashMap<Integer, FilterResult> filterResults = new HashMap<Integer, BenchmarkSpotFilter.FilterResult>();

		final ImageStack stack = imp.getImageStack();

		// Create a pool of workers
		final int nThreads = Prefs.getThreads();
		BlockingQueue<Integer> jobs = new ArrayBlockingQueue<Integer>(nThreads * 2);
		List<Worker> workers = new LinkedList<Worker>();
		List<Thread> threads = new LinkedList<Thread>();
		for (int i = 0; i < nThreads; i++)
		{
			Worker worker = new Worker(jobs, stack, spotFilter, actualCoordinates, filterResults);
			Thread t = new Thread(worker);
			workers.add(worker);
			threads.add(t);
			t.start();
		}

		final int totalFrames = stack.getSize();

		// Fit the frames
		final int step = (totalFrames > 400) ? totalFrames / 200 : 2;
		for (int i = 1; i <= totalFrames; i++)
		{
			// TODO : Should we only process the frame if there were simulated spots?

			put(jobs, i);
			if (i % step == 0)
			{
				IJ.showProgress(i, totalFrames);
				IJ.showStatus("Frame: " + i + " / " + totalFrames);
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

		for (Worker w : workers)
			time += w.time;

		// Show a table of the results
		summariseResults(filterResults);

		IJ.showStatus("");
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

	private void summariseResults(HashMap<Integer, FilterResult> filterResults)
	{
		createTable();

		StringBuilder sb = new StringBuilder();

		double signal = (simulationParameters.minSignal + simulationParameters.maxSignal) * 0.5;

		// Create the benchmark settings and the fitting settings
		sb.append(imp.getStackSize()).append("\t");
		sb.append(imp.getWidth()).append("\t");
		sb.append(imp.getHeight()).append("\t");
		sb.append(simulationParameters.molecules).append("\t");
		double density = ((double) simulationParameters.molecules / imp.getStackSize()) /
				(imp.getWidth() * imp.getHeight()) / (simulationParameters.a * simulationParameters.a / 1e6);
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
		sb.append(Utils.rounded(signal / Math.sqrt(simulationParameters.b2))).append("\t");
		sb.append(search).append("\t");
		sb.append(border).append("\t");
		sb.append(filterNames[dataFilter]).append("\t");
		sb.append(Utils.rounded(smoothing)).append("\t");
		if (differenceFilter)
		{
			sb.append(filterNames[dataFilter2]).append("\t");
			sb.append(Utils.rounded(smoothing2)).append("\t");
		}
		else
			sb.append("-\t-\t");
		sb.append(Utils.rounded(distance)).append("\t");

		// Create the overall match score
		int tp = 0, fp = 0, fn = 0;
		ArrayList<ScoredSpot> allSpots = new ArrayList<BenchmarkSpotFilter.ScoredSpot>();
		for (FilterResult result : filterResults.values())
		{
			tp += result.result.getTruePositives();
			fp += result.result.getFalsePositives();
			fn += result.result.getFalseNegatives();
			allSpots.addAll(Arrays.asList(result.spots));
		}
		// Add the number of actual points
		sb.append(tp + fn).append("\t");
		addResult(sb, new MatchResult(tp, fp, fn, 0));

		// Rank the scored spots by intensity
		Collections.sort(allSpots);
		Collections.reverse(allSpots);

		// Output the match results when there are no more true positives
		int i = 0;
		while (i < allSpots.size() && !allSpots.get(i).match)
			i++;
		tp = fp = 0; // fn will be the same
		while (i < allSpots.size())
		{
			if (allSpots.get(i).match)
				tp++;
			else
				fp++;
			i++;
		}
		addResult(sb, new MatchResult(tp, fp, fn, 0));
		sb.append(Utils.rounded(time / 1e6));

		summaryTable.append(sb.toString());
	}

	private void addResult(StringBuilder sb, MatchResult matchResult)
	{
		sb.append(matchResult.getTruePositives()).append("\t");
		sb.append(matchResult.getFalsePositives()).append("\t");
		sb.append(Utils.rounded(matchResult.getRecall())).append("\t");
		sb.append(Utils.rounded(matchResult.getPrecision())).append("\t");
		sb.append(Utils.rounded(matchResult.getJaccard())).append("\t");
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
				"Frames\tW\tH\tMolecules\tDensity (um^-2)\tN\ts (nm)\ta (nm)\tDepth (nm)\tFixed\tGain\tReadNoise (ADUs)\tB (photons)\tb2 (photons)\tSNR\t");
		sb.append("Search\tBorder\tFilter\tParam\tFilter2\tParam2\td\t");
		sb.append("N\tTP\tFP\tRecall\tPrecision\tJaccard\tTP\tFP\tRecall\tPrecision\tJaccard\tTime (ms)");
		return sb.toString();
	}
}
