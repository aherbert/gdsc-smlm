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

import gdsc.smlm.filters.MaximaSpotFilter;
import gdsc.smlm.results.PeakResults;

import java.util.ArrayList;
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
	private final BlockingQueue<FitJob> jobs;
	private final List<FitWorker> workers;
	private final List<Thread> threads;
	private long time;
	private final FitQueue queueType;
	private final PeakResults results;
	private boolean isAlive = true;

	private FitJob sum = null;

	// Used by the FitWorkers 
	private int fitting;
	private MaximaSpotFilter spotFilter;

	/**
	 * Return the fitting window size calculated using the fitting parameter and the configured peak
	 * widths. The actual window is 2n+1 around the local maxima.
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

		fitting = config.getRelativeFitting();
		spotFilter = config.createSpotFilter(true);
		
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
