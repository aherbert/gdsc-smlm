/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */

package uk.ac.sussex.gdsc.smlm.engine;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang3.concurrent.ConcurrentRuntimeException;
import uk.ac.sussex.gdsc.core.logging.LoggerUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitEngineSettings;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.filters.MaximaSpotFilter;
import uk.ac.sussex.gdsc.smlm.results.PeakResults;

/**
 * Fits local maxima using a 2D Gaussian.
 *
 * <p>Multi-threaded for speed. Uses a BlockingQueue to hold the ImageProcessor work which is then
 * processed sequentially by worker threads. The queue behaviour when the size is much greater than
 * the number of worker threads can be configured.
 */
public class FitEngine {
  /** The empty job used as a shutdown signal. */
  private static final FitJob EMPTY_JOB = new FitJob();

  private final BlockingQueue<FitJob> jobs;
  private final List<FitWorker> workers;
  private List<Thread> threads;
  private long time;
  private final FitQueue queueType;
  private final PeakResults results;

  // Used by the FitWorkers
  private final int fitting;
  private final MaximaSpotFilter spotFilter;
  private final Logger logger;
  private FitTypeCounter counter;

  /**
   * Return the fitting window size calculated using the fitting parameter and the configured peak
   * widths. The actual window is 2n+1 around the local maxima.
   *
   * @return The size of the fitting window
   */
  public int getFitting() {
    return fitting;
  }

  /**
   * Gets the spot filter for identifying candidate local maxima.
   *
   * @return The filter used for identifying candidate local maxima.
   */
  public MaximaSpotFilter getSpotFilter() {
    return (MaximaSpotFilter) spotFilter.copy();
  }

  /**
   * Constructor.
   *
   * @param config The fit configuration
   * @param results Output results (must be thread safe if using multiple threads)
   * @param threads The number of threads to use (set to 1 if less than 1)
   * @param queueType Specify the queue behaviour
   * @param queueSize The size of the queue
   */
  private FitEngine(FitEngineConfiguration config, PeakResults results, int threads,
      FitQueue queueType, int queueSize) {
    workers = new ArrayList<>(threads);
    this.queueType = queueType;
    switch (queueType) {
      case NON_BLOCKING:
      case IGNORE:
        this.jobs = new LinkedBlockingQueue<>();
        break;
      case BLOCKING:
      default:
        this.jobs = new ArrayBlockingQueue<>(queueSize);
        break;
    }
    this.results = results;

    fitting = config.getFittingWidth();
    spotFilter = config.createSpotFilter();

    logger = config.getFitConfiguration().getLog();

    // Allow logging the type of fit
    if (logger != null) {
      counter = new FitTypeCounter();
    }

    // Create the workers

    // Note - Copy the configuration for each worker.
    // Set up for a direct copy of the settings.
    final FitEngineSettings fitEngineSettings = config.getFitEngineSettings();
    final FitConfiguration fitConfiguration = config.getFitConfiguration();
    final Calibration calibration = fitConfiguration.getCalibration();
    final PSF psf = fitConfiguration.getPsf();

    for (int i = 0; i < threads; i++) {
      final FitEngineConfiguration copy =
          new FitEngineConfiguration(fitEngineSettings, calibration, psf);
      // Copy anything else not in a proto object
      copy.getFitConfiguration().copySettings(fitConfiguration);

      final FitWorker worker = new FitWorker(copy, results, jobs);
      // Note - Copy the spot filter for each worker.
      worker.setSearchParameters(getSpotFilter(), fitting);
      worker.setCounter(counter);
      workers.add(worker);
    }
  }

  /**
   * Create a new FitEngine.
   *
   * @param config The fit configuration
   * @param results Output results (must be thread safe if using multiple threads)
   * @param threads The number of threads to use (set to 1 if less than 1)
   * @param queueType Specify the queue behaviour
   * @return the fit engine
   */
  public static FitEngine create(FitEngineConfiguration config, PeakResults results, int threads,
      FitQueue queueType) {
    return create(config, results, threads, queueType, 3 * threads);
  }

  /**
   * Create a new FitEngine.
   *
   * @param config The fit configuration
   * @param results Output results (must be thread safe if using multiple threads)
   * @param threads The number of threads to use (set to 1 if less than 1)
   * @param queueType Specify the queue behaviour
   * @param queueSize The size of the queue
   * @return the fit engine
   */
  public static FitEngine create(FitEngineConfiguration config, PeakResults results, int threads,
      FitQueue queueType, int queueSize) {
    if (threads < 1) {
      threads = 1;
      queueSize = 3;
    }
    final FitEngine fitEngine = new FitEngine(config, results, threads, queueType, queueSize);
    fitEngine.start();
    return fitEngine;
  }

  /**
   * Create the threads that run the workers.
   */
  private synchronized void start() {
    threads = new ArrayList<>(workers.size());
    for (final FitWorker worker : workers) {
      final Thread t = new Thread(worker);
      threads.add(t);
      t.start();
    }
  }

  /**
   * Adds the work to the current queue, waiting if necessary for space to become available.
   *
   * @param job The job
   * @throws ConcurrentRuntimeException if interrupted while waiting to add.
   */
  public void run(FitJob job) {
    if (job == null || job.getData() == null) {
      return;
    }

    // Check the output is still OK. If no output then there is no point running any calculations.
    if (results.isActive()) {
      // Allow the jobs to create a small backlog since some frames may process faster
      if (queueType == FitQueue.IGNORE && jobs.size() > workers.size() * 1.5) {
        return;
      }

      put(job);
    }
  }

  /**
   * Adds the work to the current queue.
   *
   * @param job The job
   */
  private void put(FitJob job) {
    try {
      jobs.put(job);
    } catch (final InterruptedException ex) {
      Thread.currentThread().interrupt();
      throw new ConcurrentRuntimeException("Unexpected interruption", ex);
    }
  }

  /**
   * Signal that no more fitting work will be added to the queue.
   *
   * <p>Ask all threads to end and wait. Returns when all threads have stopped running.
   *
   * @param now Stop the work immediately, otherwise finish all work in the queue
   * @throws ConcurrentRuntimeException if interrupted while waiting to add.
   */
  public synchronized void end(boolean now) {
    if (threads.isEmpty()) {
      return;
    }

    time = 0;

    if (now) {
      // Request worker shutdown
      for (final FitWorker worker : workers) {
        worker.finish();
      }

      // Workers may be waiting for a job.
      // Add null jobs if the queue is not at capacity so they can be collected by alive workers.
      // If there are already jobs then the worker will stop due to the finish() signal.
      for (int i = 0; i < threads.size(); i++) {
        // non-blocking add to queue
        if (!jobs.offer(EMPTY_JOB)) {
          // At capacity so stop adding more
          break;
        }
      }
    } else {
      // Finish all the worker threads by passing in a null job
      for (int i = 0; i < threads.size(); i++) {
        put(EMPTY_JOB); // blocking add to queue
      }
    }

    // Collect all the threads
    for (int i = 0; i < threads.size(); i++) {
      try {
        threads.get(i).join();
        time += workers.get(i).getTime();
      } catch (final InterruptedException ex) {
        Thread.currentThread().interrupt();
        Logger.getLogger(getClass().getName()).log(Level.SEVERE, "Unexpected interruption", ex);
        throw new ConcurrentRuntimeException(ex);
      }
    }

    // Output this to the log
    if (counter != null) {
      // Get the stats we want...

      // System.out.println(results.getName()); // Dataset name
      logger.info("Fitting paths...");
      final int total = counter.getTotal();
      final int single = counter.getUnset(FitType.MULTI);
      report("Single", single, total);
      report("Multi", total - single, total);
      final int ok = counter.getSet(FitType.OK);
      report("OK", ok, total);
      report("Fail", total - ok, total);
      final int multi = total - single;
      report("FailSingle", counter.getUnset(FitType.OK | FitType.MULTI), single);
      report("FailMulti", counter.get(FitType.MULTI, FitType.OK), multi);

      report("FitSingle", counter.get(FitType.OK, FitType.MULTI), ok);
      report("FitSingleSingle", counter.get(FitType.OK, FitType.MULTI | FitType.DOUBLET_OK), ok);
      report("FitSingleDoublet", counter.get(FitType.DOUBLET_OK, FitType.MULTI), ok);
      report("FitMulti", counter.getSet(FitType.OK | FitType.MULTI), ok);
      report("FitMultiSingle", counter.getSet(FitType.MULTI_OK), ok);
      report("FitMultiDoublet", counter.getSet(FitType.MULTI_DOUBLET_OK), ok);

      report("FailMultiFitSingle", counter.get(FitType.OK | FitType.MULTI,
          FitType.MULTI_OK | FitType.MULTI_DOUBLET_OK | FitType.DOUBLET_OK), ok);
      report("FailMultiFitDoublet", counter.get(FitType.OK | FitType.MULTI | FitType.DOUBLET_OK,
          FitType.MULTI_OK | FitType.MULTI_DOUBLET_OK), ok);
    }

    threads.clear();
  }

  private void report(String name, int count, int total) {
    LoggerUtils.log(logger, Level.INFO, "%s %d / %d = %.2f", name, count, total,
        (100.00 * count) / total);
  }

  /**
   * Gets the total fitting time.
   *
   * @return the total fitting time
   */
  public long getTime() {
    return time;
  }

  /**
   * If false then the engine can be shutdown by using {@link #end(boolean)}.
   *
   * @return True if there are no worker threads
   */
  public synchronized boolean isThreadsEmpty() {
    return threads.isEmpty();
  }

  /**
   * Checks if the jobs queue is empty.
   *
   * @return True if there are no jobs queued.
   */
  public boolean isQueueEmpty() {
    return jobs.isEmpty();
  }
}
