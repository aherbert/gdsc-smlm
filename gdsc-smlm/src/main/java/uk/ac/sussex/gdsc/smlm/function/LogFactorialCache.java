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

package uk.ac.sussex.gdsc.smlm.function;

import java.util.Arrays;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;
import org.apache.commons.math3.special.Gamma;

/**
 * Cache values of the natural logarithm of the factorial function {@code ln n!}.
 */
public class LogFactorialCache {
  /**
   * Table of log factorial values. Entries may be zero and are lazily computed when accessed.
   */
  private volatile double[] objectTable;

  /**
   * Main lock guarding all access to {@link #objectTable}.
   */
  private final ReadWriteLock objectLock = new ReentrantReadWriteLock();

  /**
   * Instantiates a new log factorial.
   *
   * <p>The size of the newly constructed object can be obtained using {@link #getMaxN()}.
   */
  public LogFactorialCache() {
    // ln(0!) = ln(1!) = 0
    objectTable = new double[2];
  }

  /**
   * Create a copy instance.
   *
   * @param objectTable the object table
   */
  private LogFactorialCache(double[] objectTable) {
    this.objectTable = objectTable;
  }

  /**
   * Instantiates a new log factorial using the given table size.
   *
   * @param n argument
   */
  public LogFactorialCache(int n) {
    objectTable = new double[getLowerLimitN(n) + 1];
  }

  /**
   * Gets the lower limit on N for the table size.
   *
   * @param n argument
   * @return the lower limit
   */
  private static int getLowerLimitN(int n) {
    return Math.max(n, 1);
  }

  /**
   * Compute the log of n!. Uses the gamma function.
   *
   * @param n the value n
   * @return log(n!)
   */
  private static double logF(double n) {
    if (n <= 1) {
      return 0;
    }
    return Gamma.logGamma(n + 1);
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
   * @param n argument
   */
  public void increaseMaxN(int n) {
    n = getLowerLimitN(n);

    final ReadWriteLock rwl = this.objectLock;
    rwl.readLock().lock();
    final double[] table = objectTable;
    if (table.length <= n) {
      // Must release read lock before acquiring write lock
      rwl.readLock().unlock();
      rwl.writeLock().lock();
      try {
        // Recheck state because another thread might have
        // acquired write lock and changed state before we did.
        if (table.length <= n) {
          objectTable = Arrays.copyOf(table, n + 1);
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

    final ReadWriteLock rwl = this.objectLock;
    rwl.readLock().lock();
    double[] table = objectTable;
    if (table.length <= maxN || !checkRange(minN, maxN, table)) {
      // Must release read lock before acquiring write lock
      rwl.readLock().unlock();
      rwl.writeLock().lock();
      try {
        // Recheck state because another thread might have
        // acquired write lock and changed state before we did.
        if (table.length <= maxN) {
          objectTable = table = Arrays.copyOf(table, maxN + 1);
        }
        computeRange(minN, maxN, table);

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
   * @param minN the min N
   * @param maxN the max N
   * @param table the table
   * @return true, if successful
   */
  private static boolean checkRange(int minN, int maxN, double[] table) {
    for (int n = Math.max(2, minN); n <= maxN; n++) {
      if (table[n] == 0) {
        return false;
      }
    }
    return true;
  }

  /**
   * Ensure the table has computed values in the range minN to maxN inclusive.
   *
   * @param minN the min N
   * @param maxN the max N
   * @param table the table
   */
  private static void computeRange(int minN, int maxN, double[] table) {
    for (int n = Math.max(2, minN); n <= maxN; n++) {
      if (table[n] == 0) {
        table[n] = logF(n);
      }
    }
  }

  /**
   * Get the log of n! using tabulated values.
   *
   * <p>Note: This has no synchronisation overhead.
   *
   * @param n argument (must be positive)
   * @return log(n!)
   * @throws ArrayIndexOutOfBoundsException if n is outside the table bounds
   * @see #getMaxN()
   */
  public double getLogFactorial(int n) {
    final double[] table = objectTable;
    double value = table[n];
    if (value == 0) {
      table[n] = value = logF(n);
    }
    return value;
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
   * @param n argument (must be positive)
   * @return log(n!)
   * @throws ArrayIndexOutOfBoundsException if n is outside the table bounds
   * @see #getMaxN()
   */
  public double getLogFactorialUnsafe(int n) {
    return objectTable[n];
  }

  /**
   * Copy the object.
   *
   * @return the log factorial cache
   */
  public LogFactorialCache copy() {
    return new LogFactorialCache(objectTable);
  }

  /**
   * Copy the object with the given table size.
   *
   * <p>If n is larger then the current size then the table will be expanded.
   *
   * @param n argument
   * @return the log factorial cache
   */
  public LogFactorialCache withN(int n) {
    // Copy the values already present
    final double[] table = this.objectTable;
    return new LogFactorialCache(Arrays.copyOf(table, getLowerLimitN(n) + 1));
  }
}
