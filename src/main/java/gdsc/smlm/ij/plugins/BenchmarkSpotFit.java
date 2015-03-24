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

import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.fitting.FitResult;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.ij.plugins.BenchmarkSpotFilter.FilterResult;
import gdsc.smlm.ij.plugins.BenchmarkSpotFilter.ScoredSpot;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.utils.ImageConverter;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.match.Coordinate;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;
import ij.text.TextWindow;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

/**
 * Fits all the candidate spots identified by the benchmark spot filter plugin.
 */
public class BenchmarkSpotFit implements PlugIn
{
	private static final String TITLE = "Benchmark Spot Fit";

	private static FitConfiguration fitConfig;
	private static FitEngineConfiguration config;
	static
	{
		fitConfig = new FitConfiguration();
		config = new FitEngineConfiguration(fitConfig);
		// Set some default fit settings here
	}

	private static double fractionPositives = 100;
	private static double fractionNegativesAfterAllPositives = 50;
	private static int negativesAfterAllPositives = 10;
	private static double distance = 1.5;

	private boolean extraOptions = false;

	private static TextWindow summaryTable = null;

	private ImagePlus imp;
	private MemoryPeakResults results;
	private CreateData.SimulationParameters simulationParameters;

	private static HashMap<Integer, ArrayList<Coordinate>> actualCoordinates = null;
	private static HashMap<Integer, FilterCandidates> filterCandidates;

	private static int lastId = -1, lastFilterId = -1;
	private static double lastFractionPositives = -1;
	private static double lastFractionNegativesAfterAllPositives = -1;
	private static int lastNegativesAfterAllPositives = -1;

	private int idCount = 0;
	private int[] idList = new int[4];

	public class FilterCandidates implements Cloneable
	{
		final int p, n;
		final ScoredSpot[] spots;
		int tp, fp, tn, fn;
		FitResult[] fitResult;

		public FilterCandidates(int p, int n, ScoredSpot[] spots)
		{
			this.p = p;
			this.n = n;
			this.spots = spots;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Object#clone()
		 */
		public Object clone()
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

	/**
	 * Used to allow multi-threading of the fitting method
	 */
	private class Worker implements Runnable
	{
		volatile boolean finished = false;
		final BlockingQueue<Integer> jobs;
		final ImageStack stack;
		final HashMap<Integer, ArrayList<Coordinate>> actualCoordinates;
		final HashMap<Integer, FilterCandidates> filterCandidates;
		final HashMap<Integer, FilterCandidates> results;

		float[] data = null;

		public Worker(BlockingQueue<Integer> jobs, ImageStack stack,
				HashMap<Integer, ArrayList<Coordinate>> actualCoordinates,
				HashMap<Integer, FilterCandidates> filterCandidates)
		{
			this.jobs = jobs;
			this.stack = stack;
			this.actualCoordinates = actualCoordinates;
			this.filterCandidates = filterCandidates;
			this.results = new HashMap<Integer, FilterCandidates>();
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

			FilterCandidates candidates = filterCandidates.get(frame);
			int tp = 0, fp = 0, tn = 0, fn = 0;
			FitResult[] fitResult = new FitResult[candidates.spots.length];

			// Fit the candidates and store the results
			for (int i = 0; i < candidates.spots.length; i++)
			{
				ScoredSpot spot = candidates.spots[i];
				// TODO : Copy this from the BenchmarkFit plugin which uses minimal fit settings.

			}

			// TODO - Compute the matches of the fitted spots to the simulated positions
			// Copy this from the BenchmarkSpotFilter plugin

			// Mark the results 
			for (int i = 0; i < candidates.spots.length; i++)
			{
				ScoredSpot spot = candidates.spots[i];

				// Check if the fit succeeded
				boolean ok = true; //fitResult[i].getStatus() == FitStatus.OK;

				// Store if the candidate position was valid
				boolean match = spot.match;

				if (ok)
				{
					// Since the fit produced coordinates we can check if the final position matches 
					// an actual spot location
					// i.e. TP must start as valid candidates and finish as matches
					// TODO - use the match calculator results
					match = match && true;

					if (match)
						tp++;
					else
						fp++;
				}
				else
				{
					// Since the fit was not OK we use the original match label to set true or false negative					
					if (match)
						fn++;
					else
						tn++;
				}
			}

			// Store the results using a copy of the original
			candidates = (FilterCandidates) candidates.clone();
			candidates.tp = tp;
			candidates.fp = fp;
			candidates.tn = tn;
			candidates.fn = fn;
			candidates.fitResult = fitResult;
			results.put(frame, candidates);
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
			IJ.error(TITLE, "No benchmark results in memory");
			return;
		}
		if (BenchmarkSpotFilter.filterResults == null)
		{
			IJ.error(TITLE, "No benchmark spot candidates in memory");
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

		gd.addMessage(String
				.format("Fit candidate spots in the benchmark image created by CreateData plugin\nand identified by the Spot Filter plugin.\nPSF width = %s\n \nConfigure the fitting:",
						Utils.rounded(simulationParameters.s / simulationParameters.a)));

		gd.addSlider("Fraction_positives", 50, 100, fractionPositives);
		gd.addSlider("Fraction_negatives_after_positives", 0, 100, fractionNegativesAfterAllPositives);
		gd.addSlider("Min_negatives_after_positives", 0, 10, negativesAfterAllPositives);
		gd.addSlider("Match_distance", 1, 3, distance);

		// TODO: Collect options for fitting

		if (extraOptions)
		{
		}

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		fractionPositives = Math.abs(gd.getNextNumber());
		fractionNegativesAfterAllPositives = Math.abs(gd.getNextNumber());
		negativesAfterAllPositives = (int) Math.abs(gd.getNextNumber());
		distance = Math.abs(gd.getNextNumber());
		if (extraOptions)
		{
		}

		if (gd.invalidNumber())
			return false;

		GlobalSettings settings = new GlobalSettings();
		settings.setFitEngineConfiguration(config);
		if (!PeakFit.configureFitSolver(settings, null, extraOptions))
			return false;

		return true;
	}

	private void run()
	{
		// Extract all the results in memory into a list per frame. This can be cached
		if (lastId != simulationParameters.id)
		{
			actualCoordinates = ResultsMatchCalculator.getCoordinates(results.getResults(), true);
			lastId = simulationParameters.id;
		}

		// Extract all the candidates into a list per frame. This can be cached if the settings have not changed
		if (lastFilterId != BenchmarkSpotFilter.filterResultsId || lastFractionPositives != fractionPositives ||
				lastFractionNegativesAfterAllPositives != fractionNegativesAfterAllPositives ||
				lastNegativesAfterAllPositives != negativesAfterAllPositives)
		{
			filterCandidates = subsetFilterResults(BenchmarkSpotFilter.filterResults);

			lastFilterId = BenchmarkSpotFilter.filterResultsId;
			lastFractionPositives = fractionPositives;
			lastFractionNegativesAfterAllPositives = fractionNegativesAfterAllPositives;
			lastNegativesAfterAllPositives = negativesAfterAllPositives;
		}

		final ImageStack stack = imp.getImageStack();

		// Create a pool of workers
		final int nThreads = Prefs.getThreads();
		BlockingQueue<Integer> jobs = new ArrayBlockingQueue<Integer>(nThreads * 2);
		List<Worker> workers = new LinkedList<Worker>();
		List<Thread> threads = new LinkedList<Thread>();
		for (int i = 0; i < nThreads; i++)
		{
			Worker worker = new Worker(jobs, stack, actualCoordinates, filterCandidates);
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

		HashMap<Integer, FilterCandidates> results = new HashMap<Integer, FilterCandidates>();
		for (Worker w : workers)
		{
			results.putAll(w.results);
		}

		// TODO - Summarise the fitting results. N fits, N failures. 
		// Optimal match statistics if filtering is perfect (since fitting is not perfect).
		summariseResults(results);

		// Tile new windows
		if (idCount > 0)
		{
			idList = Arrays.copyOf(idList, idCount);
			new WindowOrganiser().tileWindows(idList);
		}

		IJ.showStatus("");
	}

	/**
	 * Extract all the filter candidates in order until the desired number of positives have been reached and the number
	 * of negatives matches the configured parameters.
	 * 
	 * @param filterResults
	 * @return The filter candidates
	 */
	private HashMap<Integer, FilterCandidates> subsetFilterResults(HashMap<Integer, FilterResult> filterResults)
	{
		// Convert fractions from percent 
		final double f1 = Math.min(1, fractionPositives / 100.0);
		final double f2 = fractionNegativesAfterAllPositives / 100.0;

		HashMap<Integer, FilterCandidates> subset = new HashMap<Integer, FilterCandidates>();
		for (Entry<Integer, FilterResult> result : filterResults.entrySet())
		{
			FilterResult r = result.getValue();

			// Determine the number of positives to find
			final int targetP = (int) Math.round(r.result.getTruePositives() * f1);
			// Count the number of positive & negatives
			int p = 0, n = 0;
			boolean reachedTarget = false;
			int nAfter = 0;

			int count = 0;
			for (ScoredSpot spot : r.spots)
			{
				if (spot.match)
				{
					p++;
					if (!reachedTarget)
					{
						reachedTarget = p >= targetP;
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

				count++;
			}

			// Debug
			System.out.printf("Frame %d : %d / (%d + %d). p=%d, n=%d, after=%d, f=%.1f\n", result.getKey().intValue(),
					r.result.getTruePositives(), r.result.getTruePositives(), r.result.getFalsePositives(), p, n,
					nAfter, (double) n / (n + p));

			subset.put(result.getKey(), new FilterCandidates(p, n, Arrays.copyOf(r.spots, count)));
		}
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

	private void summariseResults(HashMap<Integer, FilterCandidates> filterCandidates)
	{
		createTable();

		StringBuilder sb = new StringBuilder();
		summaryTable.append(sb.toString());
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
		StringBuilder sb = new StringBuilder("");
		return sb.toString();
	}
}
