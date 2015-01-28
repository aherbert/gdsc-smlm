package gdsc.smlm.engine;

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

import gdsc.smlm.filters.AverageDataProcessor;
import gdsc.smlm.filters.BlockAverageDataProcessor;
import gdsc.smlm.filters.DataProcessor;
import gdsc.smlm.filters.DoublePassSpotFilter;
import gdsc.smlm.filters.GaussianDataProcessor;
import gdsc.smlm.filters.MaximaSpotFilter;
import gdsc.smlm.filters.MedianDataProcessor;
import gdsc.smlm.filters.SinglePassSpotFilter;
import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.fitting.Gaussian2DFitter;
import gdsc.smlm.results.PeakResults;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import org.apache.commons.math3.util.FastMath;

/**
 * Fits local maxima using a 2D Gaussian.
 * <p>
 * Multi-threaded for speed. Uses a BlockingQueue to hold the ImageProcessor work which is then processed sequentially
 * by worker threads. The queue behaviour when the size is much greater than the number of worker threads can be
 * configured.
 */
public class FitEngine
{
	private final BlockingQueue<FitJob> jobs;
	private final List<FitWorker> workers;
	private final List<Thread> threads;
	private long time;
	private final FitQueue queueType;
	private final PeakResults results;
	private boolean isAlive = true;

	private FitJob sum = null;

	// Used by the FitWorkers 
	private double smooth;
	private double smooth2;
	private int border;
	private int search;
	private int fitting;
	private MaximaSpotFilter spotFilter;

	/**
	 * Return the actual smoothing window size calculated using the smoothing parameter and the configured peak widths.
	 * <p>
	 * Smoothing is done using a box filter of size (2n+1) if n >= 1, or a Gaussian approximation using a 3x3 kernel if
	 * n < 1. Disable smoothing using a smooth parameter of zero.
	 * 
	 * @return The size of the calculated smoothing window
	 */
	public double getSmooth()
	{
		return smooth;
	}

	/**
	 * Return the actual second smoothing window size calculated using the second smoothing parameter and the configured
	 * peak widths.
	 * <p>
	 * Smoothing is done using a box filter of size (2n+1) if n >= 1, or a Gaussian approximation using a 3x3 kernel if
	 * n < 1. Disable difference smoothing using a smooth2 parameter of zero.
	 * 
	 * @return The size of the calculated second smoothing window
	 */
	public double getSmooth2()
	{
		return smooth2;
	}

	/**
	 * Return the actual border size calculated using the smoothing parameter and the configured peak widths.
	 * <p>
	 * Alternatively the border can be specified in the constructor.
	 * 
	 * @return The size of the image border that is ignored
	 */
	public int getBorder()
	{
		return border;
	}

	/**
	 * Return the actual search window size (2n+1) calculated using the search parameter and the configured peak widths.
	 * 
	 * @return The size of the search window
	 */
	public int getSearch()
	{
		return search;
	}

	/**
	 * Return the actual fitting window size (2n+1) calculated using the fitting parameter and the configured peak
	 * widths.
	 * 
	 * @return The size of the fitting window
	 */
	public int getFitting()
	{
		return fitting;
	}

	/**
	 * @return The filter used for identifying candidate local maxima
	 */
	public MaximaSpotFilter getSpotFilter()
	{
		return (MaximaSpotFilter) spotFilter.clone();
	}

	/**
	 * Constructor
	 * 
	 * @param config
	 *            The fit configuration
	 * @param results
	 *            Output results
	 * @param threads
	 *            The number of threads to use (set to 1 if less than 1)
	 * @param queueType
	 *            Specify the queue behaviour
	 */
	public FitEngine(FitEngineConfiguration config, PeakResults results, int threads, FitQueue queueType)
	{
		this(config, results, threads, queueType, -1);
	}

	/**
	 * Constructor
	 * 
	 * @param config
	 *            The fit configuration
	 * @param results
	 *            Output results
	 * @param threads
	 *            The number of threads to use (set to 1 if less than 1)
	 * @param queueType
	 *            Specify the queue behaviour
	 * @param border
	 *            The border in pixels to ignore around the edge of the image
	 */
	public FitEngine(FitEngineConfiguration config, PeakResults results, int threads, FitQueue queueType, int border)
	{
		if (threads < 1)
			threads = 1;

		workers = new ArrayList<FitWorker>(threads);
		this.threads = new ArrayList<Thread>(threads);
		this.queueType = queueType;
		switch (queueType)
		{
			case BLOCKING:
			default:
				this.jobs = new ArrayBlockingQueue<FitJob>(threads * 3);
				break;
			case NON_BLOCKING:
			case IGNORE:
				this.jobs = new LinkedBlockingQueue<FitJob>();
				break;
		}
		this.results = results;

		initialiseWidths(config);
		if (border >= 0)
			this.border = border;

		createSpotFilter(config);

		// Create the workers
		for (int i = 0; i < threads; i++)
		{
			// Note - Clone the configuration and spot filter for each worker
			FitWorker worker = new FitWorker((FitEngineConfiguration) (config.clone()), results, jobs);
			worker.setSearchParameters(getSpotFilter(), fitting);
			Thread t = new Thread(worker);

			workers.add(worker);
			this.threads.add(t);

			t.start();
		}
	}

	/**
	 * Initialises the smoothing, border and search width for the fitting algorithm using the configured peak widths.
	 * 
	 * @param fitConfiguration
	 */
	private void initialiseWidths(FitEngineConfiguration config)
	{
		// TODO - This could be made a configuration option to be relative/absolute 

		// Use smallest width for smoothing and largest for the border and fitting region

		FitConfiguration fitConfiguration = config.getFitConfiguration();
		final double initialPeakStdDev0 = fitConfiguration.getInitialPeakStdDev0();
		final double initialPeakStdDev1 = fitConfiguration.getInitialPeakStdDev1();

		// Use 1 if zero to get at least a single pixel width
		double widthMin = (initialPeakStdDev0 > 0) ? initialPeakStdDev0 : 1;
		if (initialPeakStdDev1 > 0)
			widthMin = FastMath.min(initialPeakStdDev1, widthMin);
		double widthMax = (initialPeakStdDev0 > 0) ? initialPeakStdDev0 : 1;
		if (initialPeakStdDev1 > 0)
			widthMax = FastMath.max(initialPeakStdDev1, widthMax);

		// Get the half-width at half maximim
		final double hwhmMin = Gaussian2DFitter.sd2fwhm(widthMin) / 2;
		final double hwhmMax = Gaussian2DFitter.sd2fwhm(widthMax) / 2;

		// Width of smoothing box (can be zero)
		smooth = getSmoothingWindow(config.getSmooth(), hwhmMin);
		smooth2 = getSmoothingWindow(config.getSmooth2(), hwhmMin);

		// Removed this check now that there are different data filters
		//if (smooth2 <= smooth)
		//{
		//	// Cannot perform difference of smoothing
		//	smooth2 = 0;
		//}

		// Border where peaks are ignored
		border = (int) Math.floor(config.getSmooth() * hwhmMax);
		if (border < 1)
			border = 1;

		// Region for maxima finding
		search = (int) Math.ceil(config.getSearch() * hwhmMax);
		if (search < 1)
			search = 1;

		// Region for peak fitting
		fitting = (int) Math.ceil(config.getFitting() * widthMax);
		if (fitting < 2)
			fitting = 2;
	}

	private double getSmoothingWindow(double smoothingParameter, double hwhmMin)
	{
		//return BlockAverageDataProcessor.convert(smoothingParameter * hwhmMin);
		return smoothingParameter * hwhmMin;
	}

	private void createSpotFilter(FitEngineConfiguration config)
	{
		spotFilter = createSpotFilter(search, border, config.getDataFilter(), smooth, config.isDifferenceFilter(),
				config.getDataFilter2(), smooth2);
	}

	/**
	 * Create the spot filter that will be used in the fitting
	 * 
	 * @param search
	 * @param border
	 * @param dataFilter
	 * @param smooth
	 * @param differenceFilter
	 * @param dataFilter2
	 * @param smooth2
	 * @return the spot filter
	 */
	public static MaximaSpotFilter createSpotFilter(int search, int border, DataFilter dataFilter, double smooth,
			boolean differenceFilter, DataFilter dataFilter2, double smooth2)
	{
		DataProcessor processor1 = createDataProcessor(border, dataFilter, smooth);
		if (differenceFilter)
		{
			DataProcessor processor2 = createDataProcessor(border, dataFilter2, smooth2);
			return new DoublePassSpotFilter(search, border, processor1, processor2);
		}
		else
		{
			return new SinglePassSpotFilter(search, border, processor1);
		}
	}

	/**
	 * Create a data processor for the spot filter
	 * 
	 * @param border
	 * @param dataFilter
	 * @param parameter
	 * @return the data processor
	 */
	public static DataProcessor createDataProcessor(int border, DataFilter dataFilter, double parameter)
	{
		switch (dataFilter)
		{
			case MEAN:
				return new AverageDataProcessor(border, parameter);

			case BLOCK_MEAN:
				return new BlockAverageDataProcessor(border, parameter);

			case MEDIAN:
				return new MedianDataProcessor(border, parameter);

			case GAUSSIAN:
				return new GaussianDataProcessor(border, parameter);

			default:
				throw new RuntimeException("Not yet implemented: " + dataFilter.toString());
		}
	}

	/**
	 * Locate all the peaks in the given processor. Adds the work to the current queue.
	 * 
	 * @param job
	 *            The job
	 */
	public void run(FitJob job)
	{
		if (!isAlive || job == null || job.getData() == null)
			return;

		// Check the output is still OK. If no output then there is no point running any calculations.
		if (results.isActive())
		{
			// Allow the jobs to create a small backlog since some frames may process faster
			if (queueType == FitQueue.IGNORE && jobs.size() > threads.size() * 1.5)
				return;

			put(job);
		}
		else
		{
			isAlive = false;
		}
	}

	/**
	 * Adds the work to the current queue.
	 * 
	 * @param job
	 *            The job
	 */
	private void put(FitJob job)
	{
		try
		{
			jobs.put(job);
		}
		catch (InterruptedException e)
		{
			// TODO - Handle thread errors
			throw new RuntimeException("Unexpected interruption", e);
		}
	}

	/**
	 * Signal that no more fitting work will be added to the queue.
	 * <p>
	 * Ask all threads to end and wait. Returns when all threads have stopped running.
	 * 
	 * @param now
	 *            Stop the work immediately, otherwise finish all work in the queue
	 */
	public synchronized void end(boolean now)
	{
		if (threads.isEmpty())
			return;

		if (sum != null)
			put(sum); // Final frame

		time = 0;

		if (now)
		{
			// Request worker shutdown
			for (FitWorker worker : workers)
				worker.finish();

			// Workers may be waiting for a job. 
			// Add null jobs if the queue is not at capacity so they can be collected by alive workers.
			// If there are already jobs then the worker will stop due to the finish() signal.
			for (int i = 0; i < threads.size(); i++)
			{
				jobs.offer(new FitJob()); // non-blocking add to queue
			}
		}
		else
		{
			// Finish all the worker threads by passing in a null job
			for (int i = 0; i < threads.size(); i++)
			{
				put(new FitJob()); // blocking add to queue
			}
		}

		// Collect all the threads
		for (int i = 0; i < threads.size(); i++)
		{
			try
			{
				threads.get(i).join();
				time += workers.get(i).getTime();
			}
			catch (InterruptedException e)
			{
				// TODO - Handle thread errors
				e.printStackTrace();
			}
		}

		threads.clear();
	}

	/**
	 * @return the total fitting time
	 */
	public long getTime()
	{
		return time;
	}

	/**
	 * If false then the engine can be shutdown by using {@link #end(boolean)}
	 * 
	 * @return True if there are no worker threads
	 */
	public boolean isThreadsEmpty()
	{
		return threads.isEmpty();
	}

	/**
	 * @return True if there are no jobs queued
	 */
	public boolean isQueueEmpty()
	{
		return jobs.isEmpty();
	}
}
