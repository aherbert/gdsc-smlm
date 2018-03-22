package gdsc.smlm.results;

import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
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
public class PeakResultsDigest
{
	/**
	 * The default timeout in milliseconds. This should be enough to digest results sets of up to 10 million results.
	 */
	public static final long DEFAULT_TIMEOUT = 2000;

	/** The size. */
	private int size;

	/** The digest. */
	private PeakResultDigest digest = new PeakResultDigest();

	/** The digest of 1 result */
	private String value1;

	/** The digest of all results */
	private String value;

	/**
	 * Instantiates a new peak results digest.
	 */
	public PeakResultsDigest()
	{
		reset();
	}

	/**
	 * Instantiates a new peak results digest with a digest of the results.
	 *
	 * @param peakResults
	 *            the peak results
	 */
	public PeakResultsDigest(PeakResult... peakResults)
	{
		digest(peakResults);
	}

	/**
	 * Take a digest of the results using the default timeout.
	 *
	 * @param peakResults
	 *            the peak results
	 * @return true, if successful
	 */
	public boolean digest(PeakResult... peakResults)
	{
		reset();
		if (peakResults == null)
			return false;
		if (peakResults.length != 0)
		{
			digest.reset();
			digest.update(peakResults[0]);
			value1 = digest.snapshot();
			for (int i = 1; i < peakResults.length; i++)
			{
				digest.update(peakResults[i]);
			}
			value = digest.digest();
		}
		size = peakResults.length;
		return true;
	}

	/**
	 * Reset the state.
	 */
	private void reset()
	{
		size = -1;
		value1 = "";
		value = "";
	}

	/**
	 * Checks for a digest.
	 *
	 * @return true, if successful
	 */
	public boolean hasDigest()
	{
		return size >= 0;
	}

	/**
	 * Gets the digest.
	 * <p>
	 * Note that the digest is the empty string when the size is zero. This is different from using a MessageDigest on
	 * nothing which will still return a digest string, for example when using PeakResultDigest.getDigest().
	 *
	 * @return the digest
	 */
	public String getDigest()
	{
		if (size < 0)
			return null;
		if (size == 0)
			return "";
		return value;
	}

	/**
	 * Check if the digest matches the results.
	 *
	 * @param peakResults
	 *            the peak results
	 * @return true, if successful
	 */
	public boolean matches(final PeakResult... peakResults)
	{
		// Check the size. The size is -1 if we have no digest.
		if (peakResults == null)
			return size == -1;
		if (size != peakResults.length)
			return false;
		if (size == 0)
			// Nothing to digest
			return true;

		// Same logic as digest
		digest.reset();
		digest.update(peakResults[0]);
		if (value1.equals(digest.snapshot()))
			return true;
		for (int i = 1; i < peakResults.length; i++)
		{
			digest.update(peakResults[i]);
		}
		return (this.value.equals(digest.digest()));
	}

	/**
	 * Check if the digest matches the other digest.
	 * <p>
	 * Note: This will return true if there is no digest only if the state of the other
	 * digest is the same. It is left to the user to check that the state is valid.
	 *
	 * @param other
	 *            the other digest
	 * @return true, if successful
	 */
	public boolean matches(final PeakResultsDigest other)
	{
		// Check the size. We match the no digest state.
		if (size != other.size)
			return false;
		return value.equals(other.value);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode()
	{
		if (size <= 0)
			return -size; // Fast exit
		// Use the hashcode of the digest. String caches this.
		return value.hashCode();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj)
	{
		if (this == obj)
			return true;
		if (!(obj instanceof PeakResultsDigest))
			return false;
		return matches((PeakResultsDigest) obj);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		return String.format("size=%d;digest=%s", size, getDigest());
	}

	/**
	 * Take a digest of the results asynchronously. Call {@link #waitForDigest(long)} to determine if a digest was
	 * taken.
	 *
	 * @param executorService
	 *            the executor service used to execute the digest
	 * @param peakResults
	 *            the peak results
	 * @return the future
	 */
	public static Future<PeakResultsDigest> digestLater(final ExecutorService executorService,
			final PeakResult... peakResults)
	{
		return executorService.submit(new Callable<PeakResultsDigest>()
		{
			public PeakResultsDigest call() throws Exception
			{
				return new PeakResultsDigest(peakResults);
			}
		});
	}

	/**
	 * Wait for the result of {@link #digestLater(PeakResult...)} for the given time in milliseconds.
	 *
	 * @param future
	 *            the future
	 * @param timeout
	 *            the timeout (in milliseconds). Set to negative to wait indefinitely.
	 * @return the peak results digest (or null)
	 * @see #DEFAULT_TIMEOUT
	 */
	public static PeakResultsDigest waitForDigest(Future<PeakResultsDigest> future, long timeout)
	{
		try
		{
			if (timeout > 0)
				return future.get(timeout, TimeUnit.MILLISECONDS);
			else
				return future.get();
		}
		catch (final InterruptedException ie)
		{
		}
		catch (final Exception e)
		{
			// Report this 
			e.printStackTrace();
		}
		return null;
	}
}
