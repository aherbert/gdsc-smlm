/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.util.Collection;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.Objects;
import java.util.Set;
import java.util.function.Predicate;
import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.data.VisibleForTesting;
import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.data.utils.Converter;
import uk.ac.sussex.gdsc.core.data.utils.IdentityTypeConverter;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationHelper;
import uk.ac.sussex.gdsc.smlm.data.config.ConfigurationException;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;
import uk.ac.sussex.gdsc.smlm.results.count.FrameCounter;
import uk.ac.sussex.gdsc.smlm.results.procedures.BResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.BirResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.BixyResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.BixyzResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.HResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.IResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.IxyResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.IxyrResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.IxyzResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.IxyzrResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.LsePrecisionBProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.LsePrecisionProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.MlePrecisionBProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.MlePrecisionProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedureX;
import uk.ac.sussex.gdsc.smlm.results.procedures.StoredPrecisionProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.TResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.TxyResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.WResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.WxWyResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyrResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyzResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyzrResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.ZResultProcedure;

/**
 * Stores peak results in memory.
 *
 * <p>The PeakResults interface add methods are not-thread safe. The results should be wrapped in a
 * SynchronizedPeakResults object if using on multiple threads.
 */
public class MemoryPeakResults extends AbstractPeakResults {
  private static final LinkedHashMap<String, MemoryPeakResults> resultsMap = new LinkedHashMap<>();

  /**
   * The memory size of a PeakResult object without deviations. Only the standard parameters are
   * included.
   */
  @VisibleForTesting
  static final int PEAK_RESULT_SIZE = 96;
  /**
   * The memory size of a PeakResult object with deviations. Only the standard parameters are
   * included.
   */
  @VisibleForTesting
  static final int PEAK_RESULT_SIZE_WITH_DEVIATIONS = 136;

  /** The preferred distance unit. */
  public static final DistanceUnit PREFERRED_DISTANCE_UNIT = DistanceUnit.PIXEL;
  /** The preferred intensity unit. */
  public static final IntensityUnit PREFERRED_INTENSITY_UNIT = IntensityUnit.PHOTON;
  /** The preferred angle unit. */
  public static final AngleUnit PREFERRED_ANGLE_UNIT = AngleUnit.RADIAN;

  private boolean sortAfterEnd;

  /////////////////////////////////////////////////////////////////
  // START OF RESULTS STORAGE METHODS
  /////////////////////////////////////////////////////////////////

  /**
   * The results. This is encapsulated to allow changing the data structure used to store the
   * results.
   */
  protected PeakResultStoreList results;

  /**
   * Gets the result.
   *
   * @param index the index
   * @return the peak result
   */
  public PeakResult get(int index) {
    if (index >= size()) {
      throw new IndexOutOfBoundsException("Index: " + index + ", Size: " + size());
    }
    return getfX(index);
  }

  /**
   * Gets the result. Note that this uses the get(int) method from the backing PeakResultStore which
   * may return stale data if index is outside of the current size.
   *
   * @param index the index
   * @return the peak result
   */
  PeakResult getf(int index) {
    return this.results.get(index);
  }

  /**
   * Gets the result for external use or modification. Note that this uses the get(int) method from
   * the backing PeakResultStore which may return stale data if index is outside of the current
   * size.
   *
   * @param index the index
   * @return the peak result
   */
  PeakResult getfX(int index) {
    return this.results.get(index);
  }

  @Override
  public int size() {
    return this.results.size();
  }

  /**
   * Add a result. Not synchronized.
   *
   * @param result the result
   */
  @Override
  public void add(PeakResult result) {
    this.results.add(result);
  }

  /**
   * Adds the results.
   *
   * @param results the results
   */
  public void add(MemoryPeakResults results) {
    this.results.addStore(results.results);
  }

  /**
   * Add a result.
   *
   * <p>Not synchronized. Use SynchronizedPeakResults to wrap this instance for use across threads.
   *
   * {@inheritDoc}
   */
  @Override
  public void add(int peak, int origX, int origY, float origValue, double error, float noise,
      float meanIntensity, float[] params, float[] paramsStdDev) {
    add(new PeakResult(peak, origX, origY, origValue, error, noise, meanIntensity, params,
        paramsStdDev));
  }

  /**
   * {@inheritDoc}
   *
   * <p>Not synchronized. Use SynchronizedPeakResults to wrap this instance for use across threads.
   */
  @Override
  public void addAll(Collection<PeakResult> results) {
    this.results.addCollection(results);
  }

  /**
   * Add all results.
   *
   * <p>Not synchronized. Use SynchronizedPeakResults to wrap this instance for use across threads.
   */
  @Override
  public void addAll(PeakResult[] results) {
    this.results.addArray(results);
  }

  /**
   * Add all results.
   *
   * <p>Not synchronized. Use SynchronizedPeakResults to wrap this instance for use across threads.
   */
  @Override
  public void addAll(PeakResultStore results) {
    this.results.addStore(results);
  }

  /**
   * Clear the results.
   */
  private void clear() {
    this.results.clear();
  }

  /**
   * Trims the capacity of this instance to be the current size. An application can use this
   * operation to minimize the storage of an instance.
   */
  public void trimToSize() {
    this.results.trimToSize();
  }

  /**
   * Sort the results. The sort order uses the frame in ascending order.
   */
  public void sort() {
    this.results.sort();
  }

  /**
   * Sort the results.
   *
   * @param comparator the comparator
   */
  public void sort(Comparator<PeakResult> comparator) {
    this.results.sort(comparator);
  }

  /**
   * Convert to an array. This is a new allocation of storage space.
   *
   * @return the peak result array
   */
  public PeakResult[] toArray() {
    return this.results.toArray();
  }

  /**
   * Removes the null results from the store.
   */
  public void removeNullResults() {
    this.results.removeIf(Objects::isNull);
  }

  /**
   * Removes the result if it matches the filter. If objects are removed then the order of elements
   * may change.
   *
   * @param filter the filter
   * @return true, if any were removed
   */
  public boolean removeIf(Predicate<PeakResult> filter) {
    return this.results.removeIf(filter);
  }

  /////////////////////////////////////////////////////////////////
  // END OF RESULTS STORAGE METHODS
  /////////////////////////////////////////////////////////////////

  /**
   * Instantiates a new memory peak results.
   *
   * @param store the backing storage implementation
   * @throws IllegalArgumentException If the store is null
   */
  public MemoryPeakResults(PeakResultStoreList store) {
    if (store == null) {
      throw new IllegalArgumentException("Store must not be null");
    }
    results = store;
  }

  /**
   * Instantiates a new memory peak results.
   */
  public MemoryPeakResults() {
    this(1000);
  }

  /**
   * Instantiates a new memory peak results.
   *
   * @param capacity the capacity
   */
  public MemoryPeakResults(int capacity) {
    // Use the fast and simple implementation for the store
    results = new ArrayPeakResultStore(capacity);
  }

  /**
   * Instantiates a new memory peak results.
   *
   * @param results the results
   */
  public MemoryPeakResults(Collection<PeakResult> results) {
    this();
    addAll(results);
  }

  /**
   * Instantiates a new memory peak results.
   *
   * @param psf the psf
   */
  public MemoryPeakResults(PSF psf) {
    this();
    setPsf(psf);
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   * @param copyResults Set to true to copy peak result objects
   */
  protected MemoryPeakResults(MemoryPeakResults source, boolean copyResults) {
    this.sortAfterEnd = source.sortAfterEnd;
    copySettings(source);
    // Deep copy the mutable bounds object
    setBounds(new Rectangle(source.getBounds()));
    // Copy the results
    results = (PeakResultStoreList) source.results.copy(copyResults);
  }

  /////////////////////////////////////////////////////////////////
  // START OF STATIC MEMORY STORAGE METHODS
  /////////////////////////////////////////////////////////////////

  /**
   * Gets the results.
   *
   * @param name The name of the results
   * @return Get the named results (or null if they do not exist)
   */
  public static MemoryPeakResults getResults(String name) {
    return resultsMap.get(name);
  }

  /**
   * Removes the results.
   *
   * @param name The name of the results
   * @return The removed results (or null if they do not exist)
   */
  public static MemoryPeakResults removeResults(String name) {
    return resultsMap.remove(name);
  }

  /**
   * Add the results to memory. The name is taken from the results.
   *
   * @param results the results
   */
  public static void addResults(MemoryPeakResults results) {
    if (results == null) {
      throw new NullPointerException("Results must not be null");
    }
    results.trimToSize();
    resultsMap.put(results.getName(), results);
  }

  /**
   * Gets the result names.
   *
   * @return A set of the available named results held in memory
   */
  public static Set<String> getResultNames() {
    return resultsMap.keySet();
  }

  /**
   * Gets the all results.
   *
   * @return A collection of the results held in memory
   */
  public static Collection<MemoryPeakResults> getAllResults() {
    return resultsMap.values();
  }

  /**
   * Count the number of result sets in memory.
   *
   * @return the results memory size
   */
  public static int getResultsMemorySize() {
    return resultsMap.size();
  }

  /**
   * Return true if there are no non-empty results in memory.
   *
   * @return true, if is memory empty
   */
  public static boolean isMemoryEmpty() {
    if (resultsMap.isEmpty()) {
      return true;
    }
    for (final MemoryPeakResults r : resultsMap.values()) {
      if (!r.isEmpty()) {
        return false;
      }
    }
    return true;
  }

  /**
   * Performs the specified test on each set of results in memory; any results must pass to return
   * {@code true}. Returns {@code false} if there are no results in memory.
   *
   * @param test the test
   * @return true, if the test passes on any result
   */
  public static boolean isAnyInMemory(Predicate<MemoryPeakResults> test) {
    for (final MemoryPeakResults r : resultsMap.values()) {
      if (test.test(r)) {
        return true;
      }
    }
    return false;
  }

  /**
   * Performs the specified test on each set of results in memory; all results must pass to return
   * {@code true}. Returns {@code false} if there are no results in memory.
   *
   * @param test the test
   * @return true, if the test passes on all results
   */
  public static boolean isAllInMemory(Predicate<MemoryPeakResults> test) {
    if (resultsMap.isEmpty()) {
      return false;
    }
    for (final MemoryPeakResults r : resultsMap.values()) {
      if (!test.test(r)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Count the total number of results in memory.
   *
   * @return the int
   */
  public static int countMemorySize() {
    int size = 0;
    for (final MemoryPeakResults r : resultsMap.values()) {
      size += r.size();
    }
    return size;
  }

  /**
   * Clear the results from memory.
   */
  public static void clearMemory() {
    resultsMap.clear();
  }

  /**
   * Estimate the total size of results in memory.
   *
   * @return the long
   */
  public static long estimateMemorySize() {
    long memorySize = 0;
    for (final MemoryPeakResults r : resultsMap.values()) {
      memorySize += estimateMemorySize(r);
    }
    return memorySize;
  }

  /**
   * Return an estimate of the memory size taken by PeakResult objects.
   *
   * <p>Note: This is just a guess based on measured sizes for the objects in memory.
   *
   * @param results the results
   * @return The memory size
   */
  public static long estimateMemorySize(MemoryPeakResults results) {
    long memorySize = 0;
    if (results != null && results.size() > 0) {
      final boolean includeDeviations = results.getf(0).hasParameterDeviations();
      memorySize = estimateMemorySize(results.size(), includeDeviations);
    }
    return memorySize;
  }

  /**
   * Return an estimate of the memory size taken by {@link PeakResult} objects.
   *
   * <p>Note: This is just a guess based on measured sizes for the objects in memory.
   *
   * @param size the number of results
   * @param includeDeviations Set to true if the results have deviations
   * @return The memory size
   */
  public static long estimateMemorySize(int size, boolean includeDeviations) {
    return (long) size
        * ((includeDeviations) ? PEAK_RESULT_SIZE : PEAK_RESULT_SIZE_WITH_DEVIATIONS);
  }

  /////////////////////////////////////////////////////////////////
  // END OF STATIC MEMORY STORAGE METHODS
  /////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////
  // START OF PeakResults interface METHODS
  // Note: Most of the methods are in the section for
  // storage of peak results
  /////////////////////////////////////////////////////////////////

  /**
   * {@inheritDoc}
   *
   * <p>This clears the current results but does not reduce storage allocation. This can be done
   * with {@link #trimToSize()}.
   *
   * @see #trimToSize()
   */
  @Override
  public void begin() {
    clear();
  }

  @Override
  public void end() {
    if (isSortAfterEnd()) {
      sort();
    }
  }

  /////////////////////////////////////////////////////////////////
  // END OF PeakResults interface METHODS
  /////////////////////////////////////////////////////////////////

  /**
   * Gets the bounds. These are returned in Pixel units as the bounds are defined in the PeakResults
   * interface as the bounds used to create the results.
   *
   * @param calculate Set to true to calculate the bounds if they are null or zero width/height
   * @return the bounds of the result coordinates
   * @throws DataException if conversion to pixel units is not possible
   */
  public Rectangle getBounds(boolean calculate) {
    Rectangle bounds = getBounds();
    if ((bounds == null || bounds.width == 0 || bounds.height == 0) && calculate) {
      bounds = new Rectangle();
      // Note: The bounds should be in pixels
      final Rectangle2D.Float b = getDataBounds(DistanceUnit.PIXEL);

      // Round to integer
      bounds.x = (int) Math.floor(b.x);
      bounds.y = (int) Math.floor(b.y);

      final int maxX = (int) Math.ceil(b.x + b.width);
      final int maxY = (int) Math.ceil(b.y + b.height);

      bounds.width = maxX - bounds.x;
      bounds.height = maxY - bounds.y;

      setBounds(bounds);
    }
    return bounds;
  }

  /**
   * Gets the data bounds.
   *
   * @param distanceUnit the distance unit (if null then the data bounds will be in native units)
   * @return the bounds of the result coordinates
   * @throws DataException if conversion to the required units is not possible
   */
  public Rectangle2D.Float getDataBounds(DistanceUnit distanceUnit) {
    if (isEmpty()) {
      return new Rectangle2D.Float();
    }

    // Create this first to throw an exception if invalid
    final TypeConverter<DistanceUnit> c;
    if (distanceUnit == null) {
      c = new IdentityTypeConverter<>(null);
    } else {
      c = CalibrationHelper.getDistanceConverter(getCalibration(), distanceUnit);
    }

    // Get the native bounds
    float minX = getf(0).getXPosition();
    float maxX = minX;
    float minY = getf(0).getYPosition();
    float maxY = minY;
    for (int i = 1, size = size(); i < size; i++) {
      final PeakResult p = getf(i);
      final float x = p.getXPosition();
      final float y = p.getYPosition();
      if (minX > x) {
        minX = x;
      } else if (maxX < x) {
        maxX = x;
      }
      if (minY > y) {
        minY = y;
      } else if (maxY < y) {
        maxY = y;
      }
    }

    // Convert the results
    //@formatter:off
    return new Rectangle2D.Float(
        c.convert(minX),
        c.convert(minY),
        c.convert(maxX - minX),
        c.convert(maxY - minY));
    //@formatter:on
  }

  @Override
  public boolean isActive() {
    return true;
  }

  /**
   * Sets the sort after end.
   *
   * @param sortAfterEnd True if the results should be sorted after the {@link #end()} method
   */
  public void setSortAfterEnd(boolean sortAfterEnd) {
    this.sortAfterEnd = sortAfterEnd;
  }

  /**
   * Checks if is sort after end.
   *
   * @return True if the results should be sorted after the {@link #end()} method
   */
  public boolean isSortAfterEnd() {
    return sortAfterEnd;
  }

  /**
   * Copy the results. Create new objects for the properties (avoiding a shallow copy) but does not
   * deep copy all of the peak results. Allows results to be resorted but not modified.
   *
   * @return the memory peak results
   */
  public MemoryPeakResults copy() {
    return copy(false);
  }

  /**
   * Copy the results. Create new objects for the properties (avoiding a shallow copy) and
   * optionally a deep copy all of the peak results. Copying the peak results allows modification of
   * their properties.
   *
   * @param copyResults Set to true to copy peak result objects
   * @return the memory peak results
   */
  public MemoryPeakResults copy(boolean copyResults) {
    return new MemoryPeakResults(this, copyResults);
  }

  /**
   * Checks if the results are {@code null} or empty.
   *
   * @param results the results
   * @return true if {@code null} or empty
   */
  public static boolean isEmpty(MemoryPeakResults results) {
    return results == null || results.isEmpty();
  }

  /**
   * Checks if is empty.
   *
   * @return True if empty
   */
  public boolean isEmpty() {
    return size() == 0;
  }

  /**
   * Checks if is not empty.
   *
   * @return True if not empty
   */
  public boolean isNotEmpty() {
    return size() != 0;
  }

  /**
   * Checks for background. At least one result must have a non-zero background.
   *
   * @return true, if successful
   */
  public boolean hasBackground() {
    for (int i = 0, size = size(); i < size; i++) {
      if (getf(i).getBackground() != 0) {
        return true;
      }
    }
    return false;
  }

  /**
   * Checks for noise. At least one result must have a positive noise.
   *
   * @return true, if successful
   */
  public boolean hasNoise() {
    for (int i = 0, size = size(); i < size; i++) {
      if (getf(i).hasNoise()) {
        return true;
      }
    }
    return false;
  }

  /**
   * Checks for noise. At least one result must have a positive mean intensity.
   *
   * @return true, if successful
   */
  public boolean hasMeanIntensity() {
    for (int i = 0, size = size(); i < size; i++) {
      if (getf(i).hasMeanIntensity()) {
        return true;
      }
    }
    return false;
  }

  /**
   * Checks for intensity. At least one result must have a positive intensity.
   *
   * @return true, if successful
   */
  public boolean hasIntensity() {
    for (int i = 0, size = size(); i < size; i++) {
      if (getf(i).getIntensity() > 0) {
        return true;
      }
    }
    return false;
  }

  /**
   * Checks for deviations. All results must have deviations.
   *
   * @return true, if successful
   */
  public boolean hasDeviations() {
    for (int i = 0, size = size(); i < size; i++) {
      if (!getf(i).hasParameterDeviations()) {
        return false;
      }
    }
    return !isEmpty();
  }

  /**
   * Checks for end frame. At least one result must have an end frame.
   *
   * @return true, if successful
   */
  public boolean hasEndFrame() {
    for (int i = 0, size = size(); i < size; i++) {
      if (getf(i).getFrame() != getf(i).getEndFrame()) {
        return true;
      }
    }
    return false;
  }

  /**
   * Checks for id. At least one result must have a non-zero ID.
   *
   * @return true, if successful
   */
  public boolean hasId() {
    for (int i = 0, size = size(); i < size; i++) {
      if (getf(i).getId() != 0) {
        return true;
      }
    }
    return false;
  }

  /**
   * Checks for id. At least two results must have different non-zero Ids. If the number of results
   * is 1 then the id must not be zero.
   *
   * @return true, if successful
   */
  public boolean hasMultipleId() {
    if (isEmpty()) {
      return false;
    }
    if (size() == 1) {
      return getf(0).getId() != 0;
    }
    for (int i = 0, size = size(); i < size; i++) {
      final int id = getf(i).getId();
      if (id != 0) {
        while (++i < size) {
          final int id2 = getf(i).getId();
          if (id2 != 0 && id2 != id) {
            return true;
          }
        }
        break;
      }
    }
    return false;
  }

  /**
   * Checks for precision. All results must have a stored precision value.
   *
   * @return true, if successful
   */
  public boolean hasPrecision() {
    for (int i = 0, size = size(); i < size; i++) {
      if (!getf(i).hasPrecision()) {
        return false;
      }
    }
    return true;
  }

  /**
   * Checks for 3D results. At least 1 result must have a non-zero z-coordinate.
   *
   * @return true, if successful
   */
  public boolean is3D() {
    for (int i = 0, size = size(); i < size; i++) {
      if (getf(i).getZPosition() != 0) {
        return true;
      }
    }
    return false;
  }

  /**
   * Gets the first result.
   *
   * @return the first result
   * @throws IllegalStateException If the size is zero
   */
  public PeakResult getFirst() {
    checkNotEmpty();
    return getfX(0);
  }

  private void checkNotEmpty() {
    if (isEmpty()) {
      throw new IllegalStateException("Empty");
    }
  }

  /**
   * Gets the last result.
   *
   * @return the last result
   * @throws IllegalStateException If the size is zero
   */
  public PeakResult getLast() {
    checkNotEmpty();
    return getfX(size() - 1);
  }

  /**
   * Gets the first frame.
   *
   * <p>This may be different from {@link #getMinFrame()} if the results are not sorted by frame.
   *
   * @return the first frame
   * @throws IllegalStateException If the size is zero
   */
  public int getFirstFrame() {
    checkNotEmpty();
    return getf(0).getFrame();
  }

  /**
   * Gets the last frame.
   *
   * <p>This may be different from {@link #getMaxFrame()} if the results are not sorted by frame.
   *
   * @return the last frame
   * @throws IllegalStateException If the size is zero
   */
  public int getLastFrame() {
    checkNotEmpty();
    return getf(size() - 1).getEndFrame();
  }

  /**
   * Gets the minimum frame.
   *
   * @return the min frame
   * @throws IllegalStateException If the size is zero
   */
  public int getMinFrame() {
    checkNotEmpty();
    int min = getf(0).getFrame();
    for (int i = 1, size = size(); i < size; i++) {
      if (min > getf(i).getFrame()) {
        min = getf(i).getFrame();
      }
    }
    return min;
  }

  /**
   * Gets the maximum frame.
   *
   * @return the max frame
   * @throws IllegalStateException If the size is zero
   */
  public int getMaxFrame() {
    checkNotEmpty();
    int max = getf(0).getEndFrame();
    for (int i = 1, size = size(); i < size; i++) {
      if (max < getf(i).getEndFrame()) {
        max = getf(i).getEndFrame();
      }
    }
    return max;
  }

  /**
   * Checks for null results in the store.
   *
   * @return true, if null PeakResult object(s) exist
   */
  public boolean hasNullResults() {
    for (int i = 0; i < size(); i++) {
      if (getf(i) == null) {
        return true;
      }
    }
    return false;
  }

  /**
   * Checks if is distance in preferred units.
   *
   * @return true, if is distance in preferred units
   */
  public boolean isDistanceInPreferredUnits() {
    return getDistanceUnit() == PREFERRED_DISTANCE_UNIT;
  }

  /**
   * Checks if is intensity in preferred units.
   *
   * @return true, if is intensity in preferred units
   */
  public boolean isIntensityInPreferredUnits() {
    return (getIntensityUnit() == PREFERRED_INTENSITY_UNIT);
  }

  /**
   * Checks if is angle in preferred units.
   *
   * @return true, if is angle in preferred units
   */
  public boolean isAngleInPreferredUnits() {
    return getAngleUnit() == PREFERRED_ANGLE_UNIT;
  }

  /**
   * Convert to preferred units.
   *
   * @return true, if the data is now stored in the preferred units.
   */
  public boolean convertToPreferredUnits() {
    return convertToUnits(PREFERRED_DISTANCE_UNIT, PREFERRED_INTENSITY_UNIT, PREFERRED_ANGLE_UNIT);
  }

  /**
   * Convert to the specified units. If the units are null they will remain unchanged.
   *
   * @param distanceUnit the distance unit
   * @param intensityUnit the intensity unit
   * @param angleUnit the angle unit
   * @return true, if the data is now stored in the preferred units.
   */
  public boolean convertToUnits(DistanceUnit distanceUnit, IntensityUnit intensityUnit,
      AngleUnit angleUnit) {
    if (!hasCalibration()) {
      return false;
    }

    final PeakResultConversionHelper helper =
        new PeakResultConversionHelper(getCalibration(), getPsf());
    helper.setIntensityUnit(intensityUnit);
    helper.setDistanceUnit(distanceUnit);
    if (PsfHelper.hasAngleParameters(getPsf())) {
      helper.setAngleUnit(angleUnit);
    }
    final Converter[] converters = helper.getConverters();

    if (!helper.isCalibrationChanged()) {
      // Check if already in the specified units
      return helper.isValidConversion();
    }

    // Update the calibration
    setCalibration(helper.getCalibration());

    // We must convert the noise and mean intensity
    final Converter intensityConverter = converters[PeakResult.INTENSITY];

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult p = getfX(i);
      p.setNoise(intensityConverter.convert(p.getNoise()));
      p.setMeanIntensity(intensityConverter.convert(p.getMeanIntensity()));
      if (p.hasParameterDeviations()) {
        for (int j = 0; j < converters.length; j++) {
          p.setParameter(j, converters[j].convert(p.getParameter(j)));
          p.setParameterDeviation(j, converters[j].convert(p.getParameterDeviation(j)));
        }
      } else {
        for (int j = 0; j < converters.length; j++) {
          p.setParameter(j, converters[j].convert(p.getParameter(j)));
        }
      }
    }

    return true;
  }

  /////////////////////////////////////////////////////////////////
  // START OF PROCEDURE METHODS
  // Note the converters are always created (and not cached) to
  // support thread safety, i.e. accessing the results in different
  // units across threads.
  /////////////////////////////////////////////////////////////////

  /**
   * For the first result execute the procedure.
   *
   * @param procedure the procedure
   */
  public void forFirst(PeakResultProcedure procedure) {
    if (isEmpty()) {
      return;
    }
    procedure.execute(getfX(0));
  }

  /**
   * For the first result execute the procedure.
   *
   * <p>Warning: Results with be in their native units since no unit conversion is performed.
   *
   * @param procedure the procedure
   */
  public void forFirstNative(BixyzResultProcedure procedure) {
    if (isEmpty()) {
      return;
    }
    final PeakResult r = getf(0);
    procedure.executeBixyz(r.getBackground(), r.getIntensity(), r.getXPosition(), r.getYPosition(),
        r.getZPosition());
  }

  /**
   * For the first result execute the procedure.
   *
   * <p>Warning: Results with be in their native units since no unit conversion is performed.
   *
   * @param procedure the procedure
   */
  public void forFirstNative(BResultProcedure procedure) {
    if (isEmpty()) {
      return;
    }
    final PeakResult r = getf(0);
    procedure.executeB(r.getBackground());
  }

  /**
   * For the first result execute the procedure.
   *
   * <p>Warning: Results with be in their native units since no unit conversion is performed.
   *
   * @param procedure the procedure
   */
  public void forFirstNative(IResultProcedure procedure) {
    if (isEmpty()) {
      return;
    }
    final PeakResult r = getf(0);
    procedure.executeI(r.getIntensity());
  }

  /**
   * For the first result execute the procedure.
   *
   * <p>Warning: Results with be in their native units since no unit conversion is performed.
   *
   * @param procedure the procedure
   */
  public void forFirstNative(XyzResultProcedure procedure) {
    if (isEmpty()) {
      return;
    }
    final PeakResult r = getf(0);
    procedure.executeXyz(r.getXPosition(), r.getYPosition(), r.getZPosition());
  }

  /**
   * For each result execute the procedure.
   *
   * <p>Warning: Results with be in their native units since no unit conversion is performed.
   *
   * @param procedure the procedure
   */
  public void forEachNative(BixyzResultProcedure procedure) {
    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);
      procedure.executeBixyz(r.getBackground(), r.getIntensity(), r.getXPosition(),
          r.getYPosition(), r.getZPosition());
    }
  }

  /**
   * For each result execute the procedure.
   *
   * <p>Warning: Results with be in their native units since no unit conversion is performed.
   *
   * @param procedure the procedure
   */
  public void forEachNative(BResultProcedure procedure) {
    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);
      procedure.executeB(r.getBackground());
    }
  }

  /**
   * For each result execute the procedure.
   *
   * <p>Warning: Results with be in their native units since no unit conversion is performed.
   *
   * @param procedure the procedure
   */
  public void forEachNative(IResultProcedure procedure) {
    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);
      procedure.executeI(r.getIntensity());
    }
  }

  /**
   * For each result execute the procedure.
   *
   * <p>Warning: Results with be in their native units since no unit conversion is performed.
   *
   * @param procedure the procedure
   */
  public void forEachNative(XyzResultProcedure procedure) {
    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);
      procedure.executeXyz(r.getXPosition(), r.getYPosition(), r.getZPosition());
    }
  }

  /**
   * For each result execute the procedure.
   *
   * <p>Warning: Results with be in their native units since no unit conversion is performed.
   *
   * @param procedure the procedure
   */
  public void forEach(PeakResultProcedure procedure) {
    for (int i = 0, size = size(); i < size; i++) {
      procedure.execute(getfX(i));
    }
  }

  /**
   * For each result execute the fast-exit procedure.
   *
   * @param procedure the procedure
   * @return true, if a fast exit occurred
   */
  public boolean forEach(PeakResultProcedureX procedure) {
    for (int i = 0, size = size(); i < size; i++) {
      if (procedure.execute(getfX(i))) {
        return true;
      }
    }
    return false;
  }

  /**
   * For each result execute the procedure using the specified units.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param intensityUnit the intensity unit
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(IntensityUnit intensityUnit, BResultProcedure procedure) {
    final TypeConverter<IntensityUnit> ic = getIntensityConverter(intensityUnit);

    for (int i = 0, size = size(); i < size; i++) {
      procedure.executeB(ic.convert(getf(i).getBackground()));
    }
  }

  /**
   * For each result execute the procedure using the specified units.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param intensityUnit the intensity unit
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(IntensityUnit intensityUnit, BirResultProcedure procedure) {
    final TypeConverter<IntensityUnit> ic = getIntensityConverter(intensityUnit);

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);
      //@formatter:off
      procedure.executeBir(
          ic.convert(r.getBackground()),
          ic.convert(r.getIntensity()),
          r);
      //@formatter:on
    }
  }

  /**
   * For each result execute the procedure using the specified units.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param intensityUnit the intensity unit
   * @param distanceUnit the distance unit
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(IntensityUnit intensityUnit, DistanceUnit distanceUnit,
      BixyResultProcedure procedure) {
    final TypeConverter<IntensityUnit> ic = getIntensityConverter(intensityUnit);
    final TypeConverter<DistanceUnit> dc =
        getCalibrationReader().getDistanceConverter(distanceUnit);

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);
      //@formatter:off
      procedure.executeBixy(
          ic.convert(r.getBackground()),
          ic.convert(r.getIntensity()),
          dc.convert(r.getXPosition()),
          dc.convert(r.getYPosition()));
      //@formatter:on
    }
  }

  /**
   * For each result execute the procedure using the specified units.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param intensityUnit the intensity unit
   * @param distanceUnit the distance unit
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(IntensityUnit intensityUnit, DistanceUnit distanceUnit,
      BixyzResultProcedure procedure) {
    final TypeConverter<IntensityUnit> ic = getIntensityConverter(intensityUnit);
    final TypeConverter<DistanceUnit> dc =
        getCalibrationReader().getDistanceConverter(distanceUnit);

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);
      //@formatter:off
      procedure.executeBixyz(
          ic.convert(r.getBackground()),
          ic.convert(r.getIntensity()),
          dc.convert(r.getXPosition()),
          dc.convert(r.getYPosition()),
          dc.convert(r.getZPosition()));
      //@formatter:on
    }
  }

  /**
   * For each result execute the procedure using the specified units.
   *
   * <p>This procedure is exclusive to data fit with a Gaussian2D PSF as it computes the height of
   * the Gaussian from the integral (Intensity) and the widths.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param intensityUnit the intensity unit
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(IntensityUnit intensityUnit, HResultProcedure procedure) {
    checkCalibration();

    final int[] indices = PsfHelper.getGaussian2DWxWyIndices(getPsf());

    final int isx = indices[0];
    final int isy = indices[1];
    final double twoPi = 2 * Math.PI;

    final TypeConverter<IntensityUnit> ic =
        getCalibrationReader().getIntensityConverter(intensityUnit);
    final TypeConverter<DistanceUnit> dc =
        getCalibrationReader().getDistanceConverter(DistanceUnit.PIXEL);

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);

      // Convert the widths to pixels
      final float sx = dc.convert(r.getParameter(isx));
      final float sy = dc.convert(r.getParameter(isy));

      //@formatter:off
      procedure.executeH(
          (float)(ic.convert((double)r.getIntensity()) / (twoPi * sx * sy)));
      //@formatter:on
    }
  }

  /**
   * For each result execute the procedure using the specified units.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param intensityUnit the intensity unit
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(IntensityUnit intensityUnit, IResultProcedure procedure) {
    final TypeConverter<IntensityUnit> ic = getIntensityConverter(intensityUnit);

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);
      //@formatter:off
      procedure.executeI(
          ic.convert(r.getIntensity()));
      //@formatter:on
    }
  }

  /**
   * For each result execute the procedure using the specified units.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param intensityUnit the intensity unit
   * @param distanceUnit the distance unit
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(IntensityUnit intensityUnit, DistanceUnit distanceUnit,
      IxyResultProcedure procedure) {
    final TypeConverter<IntensityUnit> ic = getIntensityConverter(intensityUnit);
    final TypeConverter<DistanceUnit> dc =
        getCalibrationReader().getDistanceConverter(distanceUnit);

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);
      //@formatter:off
      procedure.executeIxy(
          ic.convert(r.getIntensity()),
          dc.convert(r.getXPosition()),
          dc.convert(r.getYPosition()));
      //@formatter:on
    }
  }

  /**
   * For each result execute the procedure using the specified units.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param intensityUnit the intensity unit
   * @param distanceUnit the distance unit
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(IntensityUnit intensityUnit, DistanceUnit distanceUnit,
      IxyrResultProcedure procedure) {
    final TypeConverter<IntensityUnit> ic = getIntensityConverter(intensityUnit);
    final TypeConverter<DistanceUnit> dc =
        getCalibrationReader().getDistanceConverter(distanceUnit);

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getfX(i);
      //@formatter:off
      procedure.executeIxyr(
          ic.convert(r.getIntensity()),
          dc.convert(r.getXPosition()),
          dc.convert(r.getYPosition()),
          r);
      //@formatter:on
    }
  }

  /**
   * For each result execute the procedure using the specified units.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param intensityUnit the intensity unit
   * @param distanceUnit the distance unit
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(IntensityUnit intensityUnit, DistanceUnit distanceUnit,
      IxyzResultProcedure procedure) {
    final TypeConverter<IntensityUnit> ic = getIntensityConverter(intensityUnit);
    final TypeConverter<DistanceUnit> dc =
        getCalibrationReader().getDistanceConverter(distanceUnit);

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);
      //@formatter:off
      procedure.executeIxyz(
          ic.convert(r.getIntensity()),
          dc.convert(r.getXPosition()),
          dc.convert(r.getYPosition()),
          dc.convert(r.getZPosition()));
      //@formatter:on
    }
  }

  /**
   * For each result execute the procedure using the specified units.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param intensityUnit the intensity unit
   * @param distanceUnit the distance unit
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(IntensityUnit intensityUnit, DistanceUnit distanceUnit,
      IxyzrResultProcedure procedure) {
    final TypeConverter<IntensityUnit> ic = getIntensityConverter(intensityUnit);
    final TypeConverter<DistanceUnit> dc =
        getCalibrationReader().getDistanceConverter(distanceUnit);

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getfX(i);
      //@formatter:off
      procedure.executeIxyzr(
          ic.convert(r.getIntensity()),
          dc.convert(r.getXPosition()),
          dc.convert(r.getYPosition()),
          dc.convert(r.getZPosition()),
          r);
      //@formatter:on
    }
  }

  /**
   * For each result execute the procedure using the specified units.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(TResultProcedure procedure) {
    checkCalibration();

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);
      //@formatter:off
      procedure.executeT(
          r.getFrame());
      //@formatter:on
    }
  }

  /**
   * For each result execute the procedure using the specified units.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param distanceUnit the distance unit
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(DistanceUnit distanceUnit, TxyResultProcedure procedure) {
    final TypeConverter<DistanceUnit> dc = getDistanceConverter(distanceUnit);

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);
      //@formatter:off
      procedure.executeTxy(
          r.getFrame(),
          dc.convert(r.getXPosition()),
          dc.convert(r.getYPosition()));
      //@formatter:on
    }
  }

  /**
   * For each result execute the procedure using the specified units.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param distanceUnit the distance unit
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(DistanceUnit distanceUnit, WResultProcedure procedure) {
    checkCalibration();

    // Note that in the future we may support more than just Gaussian2D PSF
    // so this may have to change

    final int[] indices = PsfHelper.getGaussian2DWxWyIndices(getPsf());

    final int isx = indices[0];
    final int isy = indices[1];

    final TypeConverter<DistanceUnit> dc =
        getCalibrationReader().getDistanceConverter(distanceUnit);

    if (isx == isy) {
      for (int i = 0, size = size(); i < size; i++) {
        final PeakResult r = getf(i);
        procedure.executeW(dc.convert(r.getParameter(isx)));
      }
    } else {
      for (int i = 0, size = size(); i < size; i++) {
        final PeakResult r = getf(i);
        // Convert the separate widths into a single width
        final double s = Gaussian2DPeakResultHelper.getStandardDeviation(r.getParameter(isx),
            r.getParameter(isy));
        procedure.executeW((float) dc.convert(s));
      }
    }
  }

  /**
   * For each result execute the procedure using the specified units.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param distanceUnit the distance unit
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(DistanceUnit distanceUnit, WxWyResultProcedure procedure) {
    checkCalibration();

    // Note that in the future we may support more than just Gaussian2D PSF
    // so this may have to change

    final int[] indices = PsfHelper.getGaussian2DWxWyIndices(getPsf());

    final int isx = indices[0];
    final int isy = indices[1];

    final TypeConverter<DistanceUnit> dc =
        getCalibrationReader().getDistanceConverter(distanceUnit);

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);
      //@formatter:off
      procedure.executeWxWy(
          dc.convert(r.getParameter(isx)),
          dc.convert(r.getParameter(isy)));
      //@formatter:on
    }
  }

  /**
   * For each result execute the procedure using the specified units.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param distanceUnit the distance unit
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(DistanceUnit distanceUnit, XyResultProcedure procedure) {
    final TypeConverter<DistanceUnit> dc = getDistanceConverter(distanceUnit);

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);
      //@formatter:off
      procedure.executeXy(
          dc.convert(r.getXPosition()),
          dc.convert(r.getYPosition()));
      //@formatter:on
    }
  }

  /**
   * For each result execute the procedure using the specified units.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * <p>Warning: The peak result with be in native units.
   *
   * @param distanceUnit the distance unit
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(DistanceUnit distanceUnit, XyrResultProcedure procedure) {
    final TypeConverter<DistanceUnit> dc = getDistanceConverter(distanceUnit);

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getfX(i);
      //@formatter:off
      procedure.executeXyr(
          dc.convert(r.getXPosition()),
          dc.convert(r.getYPosition()),
          r);
      //@formatter:on
    }
  }

  /**
   * For each result execute the procedure using the specified units.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param distanceUnit the distance unit
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(DistanceUnit distanceUnit, XyzResultProcedure procedure) {
    final TypeConverter<DistanceUnit> dc = getDistanceConverter(distanceUnit);

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);
      //@formatter:off
      procedure.executeXyz(
          dc.convert(r.getXPosition()),
          dc.convert(r.getYPosition()),
          dc.convert(r.getZPosition()));
      //@formatter:on
    }
  }

  /**
   * For each result execute the procedure using the specified units.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param distanceUnit the distance unit
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(DistanceUnit distanceUnit, XyzrResultProcedure procedure) {
    final TypeConverter<DistanceUnit> dc = getDistanceConverter(distanceUnit);

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getfX(i);
      //@formatter:off
      procedure.executeXyzr(
          dc.convert(r.getXPosition()),
          dc.convert(r.getYPosition()),
          dc.convert(r.getZPosition()),
          r);
      //@formatter:on
    }
  }

  /**
   * For each result execute the procedure using the specified units.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param distanceUnit the distance unit
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(DistanceUnit distanceUnit, ZResultProcedure procedure) {
    final TypeConverter<DistanceUnit> dc = getDistanceConverter(distanceUnit);

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);
      procedure.executeZ(dc.convert(r.getZPosition()));
    }
  }

  /**
   * For each result execute the procedure.
   *
   * <p>Note the precision may not be stored in the results. The default precision for a result is
   * NaN.
   *
   * @param procedure the procedure
   */
  public void forEach(StoredPrecisionProcedure procedure) {
    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);
      procedure.executeStoredPrecision(r.getPrecision());
    }
  }

  /**
   * For each result execute the procedure.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(LsePrecisionProcedure procedure) {
    checkCalibration();

    final Gaussian2DPeakResultCalculator calculator = Gaussian2DPeakResultHelper.create(getPsf(),
        getCalibration(), Gaussian2DPeakResultHelper.LSE_PRECISION);

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);
      procedure.executeLsePrecision(calculator.getLsePrecision(r.getParameters(), r.getNoise()));
    }
  }

  /**
   * For each result execute the procedure.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(LsePrecisionBProcedure procedure) {
    checkCalibration();

    final Gaussian2DPeakResultCalculator calculator = Gaussian2DPeakResultHelper.create(getPsf(),
        getCalibration(), Gaussian2DPeakResultHelper.LSE_PRECISION_X);

    for (int i = 0, size = size(); i < size; i++) {
      procedure.executeLsePrecisionB(calculator.getLsePrecision(getf(i).getParameters()));
    }
  }

  /**
   * For each result execute the procedure.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(MlePrecisionProcedure procedure) {
    checkCalibration();

    final Gaussian2DPeakResultCalculator calculator = Gaussian2DPeakResultHelper.create(getPsf(),
        getCalibration(), Gaussian2DPeakResultHelper.MLE_PRECISION);

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);
      procedure.executeMlePrecision(calculator.getMlePrecision(r.getParameters(), r.getNoise()));
    }
  }

  /**
   * For each result execute the procedure.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param procedure the procedure
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void forEach(MlePrecisionBProcedure procedure) {
    checkCalibration();

    final Gaussian2DPeakResultCalculator calculator = Gaussian2DPeakResultHelper.create(getPsf(),
        getCalibration(), Gaussian2DPeakResultHelper.MLE_PRECISION_X);

    for (int i = 0, size = size(); i < size; i++) {
      procedure.executeMlePrecisionB(calculator.getMlePrecision(getf(i).getParameters()));
    }
  }

  /**
   * Apply a pixel translation to the results. The original coordinates are updated using a raw
   * pixel shift. The current XY coordinates are updated by converting the pixel shift to the
   * current distance units. If there are current bounds then the origin will be updated with a raw
   * pixel shift.
   *
   * @param x the x pixel shift
   * @param y the y pixel shift
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void translate(int x, int y) {
    if (x == 0 && y == 0) {
      return;
    }

    checkCalibration();

    final Rectangle bounds = getBounds();
    if (bounds != null) {
      bounds.x += x;
      bounds.y += y;
      setBounds(bounds);
    }

    // Convert the pixel shift to the units of the results
    final TypeConverter<DistanceUnit> dc =
        getCalibrationReader().getDistanceConverter(getDistanceUnit());
    final float xx = dc.convert(x);
    final float yy = dc.convert(y);

    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getf(i);
      r.setOrigX(r.getOrigX() + x);
      r.setOrigY(r.getOrigY() + y);
      r.setXPosition(r.getXPosition() + xx);
      r.setYPosition(r.getYPosition() + yy);
    }
  }

  /////////////////////////////////////////////////////////////////
  // END OF PROCEDURE METHODS
  /////////////////////////////////////////////////////////////////

  /**
   * Gets the distance unit.
   *
   * @return the distance unit
   */
  public DistanceUnit getDistanceUnit() {
    if (hasCalibration() && getCalibrationReader().hasDistanceUnit()) {
      return getCalibrationReader().getDistanceUnit();
    }
    return null;
  }

  /**
   * Gets the intensity unit.
   *
   * @return the intensity unit
   */
  public IntensityUnit getIntensityUnit() {
    if (hasCalibration() && getCalibrationReader().hasIntensityUnit()) {
      return getCalibrationReader().getIntensityUnit();
    }
    return null;
  }

  /**
   * Gets the angle unit.
   *
   * @return the angle unit
   */
  public AngleUnit getAngleUnit() {
    if (hasCalibration() && getCalibrationReader().hasAngleUnit()) {
      return getCalibrationReader().getAngleUnit();
    }
    return null;
  }

  /**
   * Gets the distance converter.
   *
   * @param distanceUnit the distance unit
   * @return the distance converter
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public TypeConverter<DistanceUnit> getDistanceConverter(DistanceUnit distanceUnit) {
    checkCalibration();
    return getCalibrationReader().getDistanceConverter(distanceUnit);
  }

  /**
   * Gets the intensity converter.
   *
   * @param intensityUnit the intensity unit
   * @return the intensity converter
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public TypeConverter<IntensityUnit> getIntensityConverter(IntensityUnit intensityUnit) {
    checkCalibration();
    return getCalibrationReader().getIntensityConverter(intensityUnit);
  }

  /**
   * Gets the time converter.
   *
   * @param timeUnit the time unit
   * @return the time converter
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public TypeConverter<TimeUnit> getTimeConverter(TimeUnit timeUnit) {
    checkCalibration();
    return getCalibrationReader().getTimeConverter(timeUnit);
  }

  /**
   * Check calibration exits.
   *
   * @throws ConfigurationException if the configuration is invalid
   */
  private void checkCalibration() {
    if (!hasCalibration()) {
      throw new ConfigurationException("No calibration");
    }
  }

  /**
   * Fix zero background to the given background.
   *
   * <p>This will fail if the calibration is missing information to convert the units.
   *
   * @param intensityUnit the intensity unit
   * @param newBackground the new background
   * @throws ConversionException if the conversion is not possible
   * @throws ConfigurationException if the configuration is invalid
   */
  public void setZeroBackground(IntensityUnit intensityUnit, float newBackground) {
    final TypeConverter<IntensityUnit> ic = getIntensityConverter(intensityUnit);
    newBackground = ic.convertBack(newBackground);
    for (int i = 0, size = size(); i < size; i++) {
      final PeakResult r = getfX(i);
      if (r.getBackground() == 0) {
        r.setBackground(newBackground);
      }
    }
  }

  /**
   * Creates the frame counter. It will be initialised with the value of {@link #getFirstFrame()} -
   * 1. This ensures that the frame from the first result will be recognised as a new frame.
   *
   * @return the frame counter
   */
  public FrameCounter newFrameCounter() {
    return new FrameCounter((isEmpty()) ? 0 : getFirstFrame());
  }

  /**
   * Gets a view of the results.
   *
   * @return the view
   */
  public PeakResultView getPeakResultView() {
    return new DynamicPeakResultView(results);
  }

  /**
   * Gets a snapshot view of the results. The view can cache the results so changes to the results
   * may not be reflected in the view. Use in a read-only context.
   *
   * @return the view
   */
  public PeakResultView getSnapshotPeakResultView() {
    return new CachedPeakResultView(results);
  }

  /**
   * Get a subset of the results that match the filter.
   *
   * @param filter the filter
   * @return the results
   */
  public PeakResult[] getSubset(Predicate<PeakResult> filter) {
    return results.subset(filter);
  }

  /**
   * Find the index of the given result. More formally, returns the lowest index <tt>i</tt> such
   * that <tt>(result==null&nbsp;?&nbsp;get(i)==null&nbsp;:&nbsp;result.equals(get(i)))</tt>, or -1
   * if there is no such index.
   *
   * @param result the result
   * @return the index (or -1)
   */
  public int indexOf(PeakResult result) {
    return results.indexOf(result);
  }
}
