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

package uk.ac.sussex.gdsc.smlm.function;

import org.apache.commons.math3.special.Gamma;

import java.util.Arrays;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

/**
 * Compute the log of {@code n!}.
 */
public class LogFactorial {
  /**
   * The master table containing {@code log(n!)}.
   *
   * <p>This should only ever be modified when holding the lock!
   */
  private static volatile double[] masterTable;

  /** All long-representable factorials. */
  static final long[] FACTORIALS = {1L, 1L, 2L, 6L, 24L, 120L, 720L, 5040L, 40320L, 362880L,
      3628800L, 39916800L, 479001600L, 6227020800L, 87178291200L, 1307674368000L, 20922789888000L,
      355687428096000L, 6402373705728000L, 121645100408832000L, 2432902008176640000L};

  static {
    masterTable = new double[FACTORIALS.length];
    // Since log(1) == 0 ignore the first two values
    for (int k = 2; k < FACTORIALS.length; k++) {
      masterTable[k] = Math.log(FACTORIALS[k]);
    }
  }

  /**
   * Main lock guarding all access to {@link #masterTable}.
   */
  private static final ReadWriteLock lock = new ReentrantReadWriteLock();

  /**
   * Gets the maximum N that is tabulated.
   *
   * <p>Note: This has synchronisation overhead.
   *
   * @return the max N
   */
  public static int getTableMaxN() {
    final ReadWriteLock rwl = LogFactorial.lock;
    rwl.readLock().lock();
    try {
      return masterTable.length - 1;
    } finally {
      // Unlock read
      rwl.readLock().unlock();
    }
  }

  /**
   * Increase the tabulated values up to a max n. Does nothing if already above the given n.
   *
   * <p>Note: This has synchronisation overhead.
   *
   * @param n the n
   */
  public static void increaseTableMaxN(int n) {
    n = getLowerLimitN(n);

    final ReadWriteLock rwl = LogFactorial.lock;
    rwl.readLock().lock();
    if (masterTable.length <= n) {
      // Must release read lock before acquiring write lock
      rwl.readLock().unlock();
      rwl.writeLock().lock();
      try {
        // Recheck state because another thread might have
        // acquired write lock and changed state before we did.
        if (masterTable.length <= n) {
          masterTable = changeSize(masterTable, n);
        }
        // Downgrade by acquiring read lock before releasing write lock
        rwl.readLock().lock();
      } finally {
        rwl.writeLock().unlock(); // Unlock write, still hold read
      }
    }
    // Unlock read
    rwl.readLock().unlock();
  }

  /**
   * Change the size of the table to tabulate up to a maximum of n.
   *
   * @param table the table
   * @param n the n
   * @return the new table (length = n+1)
   */
  private static double[] changeSize(double[] table, int n) {
    final double[] newTable = Arrays.copyOf(table, n + 1);

    // Using gamma for consistency with non-tabulated values
    int ki = table.length - 1;
    while (ki < n) {
      ki++;
      newTable[ki] = Gamma.logGamma(ki + 1.0);
    }

    return newTable;
  }

  /**
   * Reduces the tabulated values down to a max n. Does nothing if already below the given n.
   *
   * <p>Note: This has synchronisation overhead.
   *
   * @param n the new table max N
   */
  public static void reduceTableMaxN(int n) {
    final int n1 = getLowerLimitN(n) + 1;

    final ReadWriteLock rwl = LogFactorial.lock;
    rwl.readLock().lock();
    if (masterTable.length > n1) {
      // Must release read lock before acquiring write lock
      rwl.readLock().unlock();
      rwl.writeLock().lock();
      try {
        // Recheck state because another thread might have
        // acquired write lock and changed state before we did.
        if (masterTable.length > n1) {
          masterTable = changeSize(masterTable, n);
        }
        // Downgrade by acquiring read lock before releasing write lock
        rwl.readLock().lock();
      } finally {
        rwl.writeLock().unlock(); // Unlock write, still hold read
      }
    }

    // Unlock read
    rwl.readLock().unlock();
  }

  /**
   * Gets the lower limit on N in order to keep the representable factorials.
   *
   * @param n the n
   * @return the lower limit
   */
  private static int getLowerLimitN(int n) {
    return Math.max(n, FACTORIALS.length - 1);
  }

  /**
   * Compute the log of n!. Uses tabulated values or the gamma function if n is large.
   *
   * <p>Note: This has no synchronisation overhead.
   *
   * @param n the value n (must be positive)
   * @return log(n!)
   * @if n is negative
   */
  public static double logF(int n) {
    // This is not synchronized.
    // We read a copy of the table from main memory.
    // If it is too small then just compute using the Gamma function.
    final double[] logF = masterTable;
    if (n < logF.length) {
      return logF[n];
    }
    return Gamma.logGamma(n + 1.0);
  }

  /**
   * Compute the log of n!. Uses the gamma function.
   *
   * <p>Note: This has no synchronisation overhead.
   *
   * @param n the value n
   * @return log(n!)
   */
  public static double logF(double n) {
    if (n <= 1) {
      return 0;
    }
    return Gamma.logGamma(n + 1);
  }

  /////////////////////////////////////
  // Instance with pre-computed values
  /////////////////////////////////////

  private volatile double[] objectTable;

  /**
   * Main lock guarding all access to {@link #objectTable}.
   */
  private final ReadWriteLock objectLock = new ReentrantReadWriteLock();

  /**
   * Instantiates a new log factorial using the current static table.
   *
   * <p>The size of the static table can be obtained using {@link #getTableMaxN()}. Note that the
   * table may be changed by other threads. The size of the newly constructed object can be obtained
   * using {@link #getMaxN()}.
   */
  public LogFactorial() {
    // Copy the static values already present
    objectTable = LogFactorial.masterTable;
  }

  /**
   * Instantiates a new log factorial.
   *
   * @param objectTable the object table
   */
  private LogFactorial(double[] objectTable) {
    this.objectTable = objectTable;
  }

  /**
   * Instantiates a new log factorial using the given table size.
   *
   * @param n the n
   */
  public LogFactorial(int n) {
    // Copy the static values already present
    objectTable = changeSize(LogFactorial.masterTable, getLowerLimitN(n));
  }

  /**
   * Gets the maximum N that is tabulated.
   *
   * <p>Note: This has synchronisation overhead.
   *
   * @return the max N
   */
  public int getMaxN() {
    final ReadWriteLock rwl = this.objectLock;
    rwl.readLock().lock();
    try {
      return objectTable.length - 1;
    } finally {
      // Unlock read
      rwl.readLock().unlock();
    }
  }

  /**
   * Increase the tabulated values up to a max n. Does nothing if already above the given n.
   *
   * <p>Note: This has synchronisation overhead.
   *
   * @param n the n
   */
  public void increaseMaxN(int n) {
    n = getLowerLimitN(n);

    final ReadWriteLock rwl = LogFactorial.lock;
    rwl.readLock().lock();
    if (objectTable.length <= n) {
      // Must release read lock before acquiring write lock
      rwl.readLock().unlock();
      rwl.writeLock().lock();
      try {
        // Recheck state because another thread might have
        // acquired write lock and changed state before we did.
        if (objectTable.length <= n) {
          // Use the master table if it is bigger as that has all values pre-computed.
          final double[] mtable = LogFactorial.masterTable;
          if (mtable.length > objectTable.length) {
            objectTable = Arrays.copyOf(mtable, Math.min(n + 1, mtable.length));
          }

          if (objectTable.length <= n) {
            objectTable = changeSize(objectTable, n);
          }
        }
        // Downgrade by acquiring read lock before releasing write lock
        rwl.readLock().lock();
      } finally {
        rwl.writeLock().unlock(); // Unlock write, still hold read
      }
    }

    // Unlock read
    rwl.readLock().unlock();
  }

  /**
   * Ensure the table contains values for the specified range of N.
   *
   * <p>This can be called before using the object with a known range of n.
   *
   * <p>Note: This has synchronisation overhead.
   *
   * @param minN the min N
   * @param maxN the max N
   * @throws IllegalArgumentException If min is greater than max
   */
  public void ensureRange(int minN, int maxN) {
    // Validate the range
    if (minN < 0) {
      minN = 0;
    }
    maxN = getLowerLimitN(maxN);
    if (minN > maxN) {
      throw new IllegalArgumentException("Max must be greater than min");
    }

    final ReadWriteLock rwl = LogFactorial.lock;
    rwl.readLock().lock();
    if (objectTable.length <= maxN || !checkRange(minN, maxN)) {
      // Must release read lock before acquiring write lock
      rwl.readLock().unlock();
      rwl.writeLock().lock();
      try {
        // Recheck state because another thread might have
        // acquired write lock and changed state before we did.
        if (objectTable.length <= maxN) {
          // Use the master table if it is bigger as that has all values pre-computed.
          final double[] mtable = LogFactorial.masterTable;
          objectTable =
              Arrays.copyOf((mtable.length > objectTable.length) ? mtable : objectTable, maxN + 1);
        }

        // Check range has pre-computed values
        for (int n = Math.max(2, minN); n <= maxN; n++) {
          if (objectTable[n] == 0) {
            // Already hold the write lock so just write the values
            objectTable[n] = Gamma.logGamma(n + 1.0);
          }
        }

        // Downgrade by acquiring read lock before releasing write lock
        rwl.readLock().lock();
      } finally {
        rwl.writeLock().unlock(); // Unlock write, still hold read
      }
    }

    // Unlock read
    rwl.readLock().unlock();
  }

  /**
   * Check the table has computed values in the range minN to maxN inclusive.
   *
   * <p>Only call when holding the read lock.
   *
   * @param minN the min N
   * @param maxN the max N
   * @return true, if successful
   */
  private boolean checkRange(int minN, int maxN) {
    for (int n = Math.max(2, minN); n <= maxN; n++) {
      if (objectTable[n] == 0) {
        return false;
      }
    }
    return true;
  }

  /**
   * Get the log of n! using tabulated values.
   *
   * <p>Note: This has no synchronisation overhead.
   *
   * @param n the n (must be positive)
   * @return log(n!)
   * @throws ArrayIndexOutOfBoundsException if n is outside the table bounds
   * @see #getMaxN()
   */
  public double getLogF(int n) {
    final double value = objectTable[n];
    return (value == 0) ? LogFactorial.logF(n) : value;
  }

  /**
   * Get the log of n! using tabulated values.
   *
   * <p>Note: This has no synchronisation overhead.
   *
   * <p>No checks are made that the object table contains a pre-computed value, for instance in the
   * case where the table was created using {@link #ensureRange(int, int)} there may be values
   * outside the range that are zero.
   *
   * @param n the n (must be positive)
   * @return log(n!)
   * @throws ArrayIndexOutOfBoundsException if n is outside the table bounds
   * @see #getMaxN()
   */
  public double getLogFUnsafe(int n) {
    return objectTable[n];
  }

  /**
   * Copy the object.
   *
   * @return the log factorial
   */
  public LogFactorial copy() {
    return new LogFactorial(objectTable);
  }

  /**
   * Copy the object with the given table size.
   *
   * <p>If n is larger then the current size then the table will be expanded.
   *
   * @param n the n
   * @return the log factorial
   */
  public LogFactorial withN(int n) {
    // Copy the values already present
    final double[] mtable = LogFactorial.masterTable;
    double[] otable = this.objectTable;
    if (mtable.length > otable.length) {
      otable = mtable;
    }
    return new LogFactorial(changeSize(otable, getLowerLimitN(n)));
  }
}
