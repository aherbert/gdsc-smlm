package gdsc.smlm.results;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Class to allow fast comparison of peak results arrays without storing the array. Comparison is done using a staged
 * digest.
 */
public class PeakResultsSnapshot
{
	/**
	 * The Enum SnapshotState.
	 */
	public enum SnapshotState
	{
		/** No snapshot exists. */
		NONE,
		/** The snapshot is initialised */
		SNAPSHOT,
		/** A timeout occurred when computing the snapshot asynchronously. */
		TIMEOUT,
		/** An error occurred when computing the snapshot asynchronously. */
		ERROR;
	}

	/** The executor service. */
	private static ExecutorService executorService = null;

	/**
	 * Gets the executor service.
	 *
	 * @return the executor service
	 */
	private ExecutorService getExecutorService()
	{
		if (executorService == null)
			executorService = Executors.newFixedThreadPool(1);
		return executorService;
	}

	/**
	 * The default timeout in milliseconds. This should be enough to digest results sets of up to 10 million results.
	 */
	public long DEFAULT_TIMEOUT = 2000;

	/** The future. */
	private Future<?> future;

	/** The size. */
	private int size;

	/** The digest. */
	private PeakResultDigest digest = new PeakResultDigest();

	/** The values. */
	private String[] values;

	/**
	 * Instantiates a new peak results snapshot.
	 */
	public PeakResultsSnapshot()
	{
		size = -1;
	}

	/**
	 * Instantiates a new peak results snapshot with a snapshot of the results.
	 *
	 * @param peakResults
	 *            the peak results
	 */
	public PeakResultsSnapshot(PeakResult... peakResults)
	{
		snapshot(peakResults);
	}

	/**
	 * Take a snapshot of the results using the default timeout.
	 *
	 * @param peakResults
	 *            the peak results
	 * @return true, if successful
	 * @see #DEFAULT_TIMEOUT
	 */
	public boolean snapshot(PeakResult... peakResults)
	{
		size = -1;
		values = digest(peakResults);
		size = peakResults.length;
		return true;
	}

	/**
	 * Take a snapshot of the results.
	 *
	 * @param timeout
	 *            the timeout (in milliseconds)
	 * @param peakResults
	 *            the peak results
	 * @return true, if successful
	 */
	public boolean snapshot(long timeout, PeakResult... peakResults)
	{
		snapshotLater(peakResults);
		return waitForSnapshot(timeout);
	}

	/**
	 * Take a snapshot of the results asynchronously. Call {@link #waitForSnapshot(long)} to determine if a snapshot was
	 * taken.
	 *
	 * @param peakResults
	 *            the peak results
	 */
	public void snapshotLater(final PeakResult... peakResults)
	{
		snapshotLater(getExecutorService(), peakResults);
	}

	/**
	 * Take a snapshot of the results asynchronously. Call {@link #waitForSnapshot(long)} to determine if a snapshot was
	 * taken.
	 *
	 * @param executorService
	 *            the executor service used to execute the digest
	 * @param peakResults
	 *            the peak results
	 */
	public void snapshotLater(final ExecutorService executorService, final PeakResult... peakResults)
	{
		size = -1;
		future = executorService.submit(new Runnable()
		{
			public void run()
			{
				values = digest(peakResults);
				size = peakResults.length;
			}
		});
	}

	/**
	 * Do a staged digest of the results.
	 *
	 * @param peakResults
	 *            the peak results
	 * @return the string[] of digest results
	 */
	private String[] digest(PeakResult... peakResults)
	{
		// Snapshot digest of 3 values (at 2, 4 and all results). 
		// This should allow a fast compare of results sets with the same size.
		String[] values = new String[3];
		int j = 0;
		digest.reset();
		for (int i = 0, snap = 1; i < peakResults.length; i++)
		{
			digest.update(peakResults[i]);
			if (i == snap)
			{
				values[j++] = digest.digest();
				snap = (j == values.length - 1) ? peakResults.length : snap * 2;
			}
		}
		values[j++] = digest.digest();
		while (j < values.length)
			values[j++] = "";
		return values;
	}

	/**
	 * Wait for the result of {@link #snapshotLater(PeakResult...)} for the given time in milliseconds.
	 *
	 * @param timeout
	 *            the timeout (in milliseconds). Set to negative to wait indefinitely.
	 * @return true, if successful
	 */
	public boolean waitForSnapshot(long timeout)
	{
		if (future != null)
		{
			try
			{
				if (timeout > 0)
					future.get(timeout, TimeUnit.MILLISECONDS);
				else
					future.get();
			}
			catch (final InterruptedException ie)
			{
				size = -2;
			}
			catch (final Exception e)
			{
				size = -3;
				// Report this 
				e.printStackTrace();
			}
		}
		return hasSnapshot();
	}

	/**
	 * Checks for a snapshot.
	 *
	 * @return true, if successful
	 */
	public boolean hasSnapshot()
	{
		return size >= 0;
	}

	/**
	 * Gets the state of the snapshot.
	 *
	 * @return the state
	 */
	public SnapshotState getState()
	{
		if (hasSnapshot())
			return SnapshotState.SNAPSHOT;
		if (size == -1)
			return SnapshotState.NONE;
		if (size == -2)
			return SnapshotState.TIMEOUT;
		if (size == -3)
			return SnapshotState.ERROR;
		throw new IllegalStateException("Unknown state: " + size);
	}

	/**
	 * Check if the results match the snapshot.
	 *
	 * @param peakResults
	 *            the peak results
	 * @return true, if successful
	 */
	public boolean matches(final PeakResult... peakResults)
	{
		// Check the size. The size is -1 if we have no snapshot.
		if (size != peakResults.length)
			return false;

		// Same logic as digest
		int j = 0;
		digest.reset();
		for (int i = 0, snap = 1; i < peakResults.length; i++)
		{
			digest.update(peakResults[i]);
			if (i == snap)
			{
				if (!this.values[j++].equals(digest.digest()))
					return false;
				snap = (j == values.length - 1) ? peakResults.length : snap * 2;
			}
		}
		return (this.values[j].equals(digest.digest()));
	}
}
