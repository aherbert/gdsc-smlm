package gdsc.smlm.engine;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import gdsc.core.logging.FileLogger;
import gdsc.core.logging.Logger;
import gdsc.smlm.filters.MaximaSpotFilter;
import gdsc.smlm.results.PeakResults;

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
	private Logger logger = null;
	private FileLogger fileLogger = null;
	private FitTypeCounter counter = null;

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
		return spotFilter.clone();
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
		this(config, results, threads, queueType, 3 * threads);
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
	 * @param queueSize
	 *            The size of the queue ({@link #queueType}
	 */
	public FitEngine(FitEngineConfiguration config, PeakResults results, int threads, FitQueue queueType, int queueSize)
	{
		if (threads < 1)
		{
			threads = 1;
			queueSize = 3;
		}

		workers = new ArrayList<FitWorker>(threads);
		this.threads = new ArrayList<Thread>(threads);
		this.queueType = queueType;
		switch (queueType)
		{
			case BLOCKING:
			default:
				this.jobs = new ArrayBlockingQueue<FitJob>(queueSize);
				break;
			case NON_BLOCKING:
			case IGNORE:
				this.jobs = new LinkedBlockingQueue<FitJob>();
				break;
		}
		this.results = results;

		fitting = config.getRelativeFitting();
		spotFilter = config.createSpotFilter(true);

		logger = config.getFitConfiguration().getLog();
		//// Allow debugging the fit process
		//try
		//{
		//	fileLogger = new FileLogger(
		//			String.format("/tmp/%s%s.log", config.getFitConfiguration().getFitSolver().getShortName(),
		//					(config.getFitConfiguration().isModelCamera()) ? "C" : ""));
		//}
		//catch (java.io.FileNotFoundException e)
		//{
		//}

		// Allow logging the type of fit
		if (logger != null)
			counter = new FitTypeCounter();

		// Create the workers
		for (int i = 0; i < threads; i++)
		{
			// Note - Clone the configuration and spot filter for each worker
			FitWorker worker = new FitWorker(config.clone(), results, jobs);
			worker.setSearchParameters(getSpotFilter(), fitting);
			worker.setLogger2(fileLogger);
			worker.setCounter(counter);
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

		if (fileLogger != null)
			fileLogger.close();

		// Output this to the log
		if (counter != null)
		{
			// Get the stats we want...
			//System.out.println(results.getName()); // Dataset name
			logger.info("Fitting paths...");
			final int total = counter.getTotal();
			final int single = counter.getUnset(FitType.NEIGHBOURS);
			report("Single", single, total);
			report("Neighbours", total - single, total);
			final int ok = counter.getSet(FitType.OK);
			report("OK", ok, total);
			final int multi = total - ok;
			report("Fail", multi, total);
			report("FailSingle", counter.getUnset(FitType.OK | FitType.NEIGHBOURS), single);
			report("FailMulti", counter.get(FitType.NEIGHBOURS, FitType.OK), multi);

			report("FitSingle", counter.get(FitType.OK, FitType.NEIGHBOURS), ok);
			report("FitSingleSingle", counter.get(FitType.OK, FitType.NEIGHBOURS | FitType.DOUBLET_OK), ok);
			report("FitSingleDoublet", counter.get(FitType.DOUBLET_OK, FitType.NEIGHBOURS), ok);
			report("FitMulti", counter.getSet(FitType.NEIGHBOURS_OK), ok);
			report("FailMultiFitSingle",
					counter.get(FitType.OK | FitType.NEIGHBOURS, FitType.NEIGHBOURS_OK | FitType.DOUBLET_OK), ok);
			report("FailMultiFitDoublet",
					counter.get(FitType.OK | FitType.NEIGHBOURS | FitType.DOUBLET_OK, FitType.NEIGHBOURS_OK), ok);
		}

		threads.clear();
	}

	private void report(String name, int count, int total)
	{
		logger.info("%s %d / %d = %.2f", name, count, total, (100.00 * count) / total);
		//System.out.printf("%s %d / %d = %.2f\n", name, count, total, (100.00 * count) / total);
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
