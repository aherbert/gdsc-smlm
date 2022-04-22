/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.results;

import java.util.concurrent.Callable;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import org.apache.commons.lang3.concurrent.ConcurrentRuntimeException;

/**
 * Class to allow fast comparison of peak results arrays without storing the array. Comparison is
 * done using a staged digest.
 */
public class PeakResultsDigest {
  /**
   * The default timeout in milliseconds. This should be enough to digest results sets of up to 10
   * million results.
   */
  public static final long DEFAULT_TIMEOUT = 2000;

  /** The size. */
  private int size;

  /** The digest. */
  private final PeakResultDigest digest = new PeakResultDigest();

  /** The digest of 1 result. */
  private String value1;

  /** The digest of all results. */
  private String value;

  /**
   * Creates a new peak results digest.
   */
  public PeakResultsDigest() {
    reset();
  }

  /**
   * Creates a new peak results digest with a digest of the results.
   *
   * @param peakResults the peak results
   * @return the peak results digest
   */
  public static PeakResultsDigest create(PeakResult... peakResults) {
    final PeakResultsDigest d = new PeakResultsDigest();
    d.digest(peakResults);
    return d;
  }

  /**
   * Take a digest of the results using the default timeout.
   *
   * @param peakResults the peak results
   * @return true, if successful
   */
  public boolean digest(PeakResult... peakResults) {
    reset();
    if (peakResults == null) {
      return false;
    }
    if (peakResults.length != 0) {
      digest.reset();
      digest.update(peakResults[0]);
      value1 = digest.snapshot();
      for (int i = 1; i < peakResults.length; i++) {
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
  private void reset() {
    size = -1;
    value1 = "";
    value = "";
  }

  /**
   * Checks for a digest.
   *
   * @return true, if successful
   */
  public boolean hasDigest() {
    return size >= 0;
  }

  /**
   * Gets the digest.
   *
   * <p>Note that the digest is the empty string when the size is zero. This is different from using
   * a MessageDigest on nothing which will still return a digest string, for example when using
   * PeakResultDigest.getDigest().
   *
   * @return the digest
   */
  public String getDigest() {
    if (size < 0) {
      return null;
    }
    if (size == 0) {
      return "";
    }
    return value;
  }

  /**
   * Check if the digest matches the results.
   *
   * @param peakResults the peak results
   * @return true, if successful
   */
  public boolean matches(final PeakResult... peakResults) {
    // Check the size. The size is -1 if we have no digest.
    if (peakResults == null) {
      return size == -1;
    }
    if (size != peakResults.length) {
      return false;
    }
    if (size == 0) {
      // Nothing to digest
      return true;
    }

    // Same logic as digest
    digest.reset();
    digest.update(peakResults[0]);
    if (value1.equals(digest.snapshot())) {
      return true;
    }
    for (int i = 1; i < peakResults.length; i++) {
      digest.update(peakResults[i]);
    }
    return (this.value.equals(digest.digest()));
  }

  /**
   * Check if the digest matches the other digest.
   *
   * <p>Note: This will return true if there is no digest only if the state of the other digest is
   * the same. It is left to the user to check that the state is valid.
   *
   * @param other the other digest
   * @return true, if successful
   */
  public boolean matches(final PeakResultsDigest other) {
    // Check the size. We match the no digest state.
    if (size != other.size) {
      return false;
    }
    return value.equals(other.value);
  }

  @Override
  public int hashCode() {
    if (size <= 0) {
      return -size; // Fast exit
    }
    // Use the hashcode of the digest. String caches this.
    return value.hashCode();
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) {
      return true;
    }
    if (!(obj instanceof PeakResultsDigest)) {
      return false;
    }
    return matches((PeakResultsDigest) obj);
  }

  @Override
  public String toString() {
    return String.format("size=%d;digest=%s", size, getDigest());
  }

  /**
   * Take a digest of the results asynchronously. Call {@link #waitForDigest(Future, long)} to
   * determine if a digest was taken.
   *
   * @param executorService the executor service used to execute the digest
   * @param peakResults the peak results
   * @return the future
   */
  public static Future<PeakResultsDigest> digestLater(final ExecutorService executorService,
      final PeakResult... peakResults) {
    return executorService.submit((Callable<PeakResultsDigest>) PeakResultsDigest::create);
  }

  /**
   * Wait for the result of {@link #digestLater(ExecutorService, PeakResult...)} for the given time
   * in milliseconds.
   *
   * <p>If a TimeoutException occurs this is ignored and the method returns null.
   *
   * @param future the future
   * @param timeout the timeout (in milliseconds). Set to negative to wait indefinitely.
   * @return the peak results digest (or null)
   * @throws ConcurrentRuntimeException if interrupted while waiting, the future was cancelled or an
   *         execution exception is thrown
   * @see #DEFAULT_TIMEOUT
   */
  public static PeakResultsDigest waitForDigest(Future<PeakResultsDigest> future, long timeout) {
    try {
      if (timeout > 0) {
        return future.get(timeout, TimeUnit.MILLISECONDS);
      }
      return future.get();
    } catch (final InterruptedException ex) {
      // Restore interrupted state...
      Thread.currentThread().interrupt();
      throw new ConcurrentRuntimeException(ex);
    } catch (final CancellationException | ExecutionException ex) {
      // Report this
      throw new ConcurrentRuntimeException(ex);
    } catch (final TimeoutException ignored) {
      // Ignore
    }
    return null;
  }
}
