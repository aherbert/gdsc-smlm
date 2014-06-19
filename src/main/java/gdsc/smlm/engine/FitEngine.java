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

import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.fitting.Gaussian2DFitter;
import gdsc.smlm.results.PeakResults;

import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

/**
 * Fits local maxima using a 2D Gaussian.
 * <p>
 * Multi-threaded for speed. Uses a BlockingQueue to hold the ImageProcessor work which is then processed sequentially
 * by worker threads. The queue behaviour when the size is much greater than the number of worker threads can be
 * configured.
 */
public class FitEngine
{
	private BlockingQueue<FitJob> jobs = null;
	private List<FitWorker> workers = new LinkedList<FitWorker>();
	private List<Thread> threads = new LinkedList<Thread>();
	private long time;
	private FitQueue queueType;
	private PeakResults results;
	private boolean isAlive = true;

	private FitJob sum = null;

	// Used by the FitWorkers 
	private double smooth;
	private double smooth2;
	private int border;
	private int search;

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
		init(config, results, threads, queueType, -1);
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
		init(config, results, threads, queueType, border);
	}

	private void init(FitEngineConfiguration config, PeakResults results, int threads, FitQueue queueType, int border)
	{
		if (threads < 1)
			threads = 1;

		createQueue(queueType, threads);
		this.results = results;

		initialiseWidths(config);
		if (border >= 0)
			this.border = border;

		// Create the workers
		for (int i = 0; i < threads; i++)
		{
			// Note - Clone the configuration for each worker
			FitWorker worker = new FitWorker((FitEngineConfiguration) (config.clone()), results, jobs);
			worker.smooth = smooth;
			worker.smooth2 = smooth2;
			worker.border = this.border;
			worker.search = search;
			Thread t = new Thread(worker);

			workers.add(worker);
			this.threads.add(t);

			t.start();
		}
	}

	private void createQueue(FitQueue queueType, int threads)
	{
		this.queueType = queueType;
		switch (queueType)
		{
			case BLOCKING:
				this.jobs = new ArrayBlockingQueue<FitJob>(threads * 3);
				break;
			case NON_BLOCKING:
			case IGNORE:
				this.jobs = new LinkedBlockingQueue<FitJob>();
				break;
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
		float initialPeakStdDev0 = fitConfiguration.getInitialPeakStdDev0();
		float initialPeakStdDev1 = fitConfiguration.getInitialPeakStdDev1();

		// Use 1 if zero to get at least a single pixel width
		float widthMin = (initialPeakStdDev0 > 0) ? initialPeakStdDev0 : 1;
		if (initialPeakStdDev1 > 0)
			widthMin = Math.min(initialPeakStdDev1, widthMin);
		float widthMax = (initialPeakStdDev0 > 0) ? initialPeakStdDev0 : 1;
		if (initialPeakStdDev1 > 0)
			widthMax = Math.max(initialPeakStdDev1, widthMax);

		// Get the half-width at half maximim
		float hwhmMin = Gaussian2DFitter.sd2fwhm(widthMin) / 2;
		float hwhmMax = Gaussian2DFitter.sd2fwhm(widthMax) / 2;

		// Width of smoothing box (can be zero)
		smooth = getSmoothingWindow(config.getSmooth(), hwhmMin);
		smooth2 = getSmoothingWindow(config.getSmooth2(), hwhmMin);

		if ((smooth > 0 && smooth < 1) && smooth2 < 1)
		{
			// Cannot perform difference of smoothing with both windows between 0 and 1.
			// (Only a 3x3 Gaussian is supported and so this would result in both filters being the same.)
			smooth2 = 0;
		}

		// Border where peaks are ignored
		border = (int) Math.floor(config.getSmooth() * hwhmMax);
		if (border < 1)
			border = 1;

		// TODO - What should this be to be robust?
		// 3 sigma should cover 99% of the Gaussian curve

		// Search region for peak fitting
		search = (int) Math.ceil(config.getSearch() * widthMax);
		if (search < 1)
			search = 1;
	}

	public double getSmoothingWindow(double smoothingParameter, float hwhmMin)
	{
		double smooth = smoothingParameter * hwhmMin;
		// Only allow non-integer amounts below 1 for Gaussian smoothing
		if (smooth > 1)
			smooth = (int) smooth;
		else if (smooth < 0)
			smooth = 0;
		return smooth;
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
