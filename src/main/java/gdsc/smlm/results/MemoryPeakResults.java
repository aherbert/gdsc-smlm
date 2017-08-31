package gdsc.smlm.results;

import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.util.Collection;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.Set;

import gdsc.core.data.DataException;
import gdsc.core.data.utils.ConversionException;
import gdsc.core.data.utils.Converter;
import gdsc.core.data.utils.IdentityTypeConverter;
import gdsc.core.data.utils.TypeConverter;
import gdsc.smlm.data.config.CalibrationHelper;
import gdsc.smlm.data.config.ConfigurationException;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.UnitProtos.AngleUnit;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import gdsc.smlm.results.procedures.BIXYResultProcedure;
import gdsc.smlm.results.procedures.BIXYZResultProcedure;
import gdsc.smlm.results.procedures.BResultProcedure;
import gdsc.smlm.results.procedures.HResultProcedure;
import gdsc.smlm.results.procedures.IResultProcedure;
import gdsc.smlm.results.procedures.IXYRResultProcedure;
import gdsc.smlm.results.procedures.IXYResultProcedure;
import gdsc.smlm.results.procedures.IXYZResultProcedure;
import gdsc.smlm.results.procedures.LSEPrecisionBProcedure;
import gdsc.smlm.results.procedures.LSEPrecisionProcedure;
import gdsc.smlm.results.procedures.MLEPrecisionBProcedure;
import gdsc.smlm.results.procedures.MLEPrecisionProcedure;
import gdsc.smlm.results.procedures.PeakResultProcedure;
import gdsc.smlm.results.procedures.PeakResultProcedureX;
import gdsc.smlm.results.procedures.TResultProcedure;
import gdsc.smlm.results.procedures.TXYResultProcedure;
import gdsc.smlm.results.procedures.WResultProcedure;
import gdsc.smlm.results.procedures.WxWyResultProcedure;
import gdsc.smlm.results.procedures.XYRResultProcedure;
import gdsc.smlm.results.procedures.XYResultProcedure;
import gdsc.smlm.results.procedures.XYZResultProcedure;

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

/**
 * Stores peak results in memory.
 * <p>
 * The PeakResults interface add methods are not-thread safe. The results should be wrapped in a SynchronizedPeakResults
 * object if using on multiple threads.
 */
public class MemoryPeakResults extends AbstractPeakResults implements Cloneable
{
	private static LinkedHashMap<String, MemoryPeakResults> resultsMap = new LinkedHashMap<String, MemoryPeakResults>();
	private static final Runtime s_runtime = Runtime.getRuntime();
	private static int byteSize = 0;
	private static int byteSizeWithDeviations = 0;

	private static final int DEFAULT_SIZE = 96;
	private static final int DEFAULT_SIZE_WITH_DEVIATIONS = 144;

	private boolean sortAfterEnd;

	/////////////////////////////////////////////////////////////////
	// START OF RESULTS STORAGE METHODS 
	/////////////////////////////////////////////////////////////////

	/**
	 * The results.
	 * This is encapsulated to allow changing the data structure used to store the results.
	 */
	private PeakResultStore results;

	/**
	 * Gets the result.
	 *
	 * @param index
	 *            the index
	 * @return the peak result
	 */
	public PeakResult get(int index)
	{
		if (index >= size())
			throw new IndexOutOfBoundsException("Index: " + index + ", Size: " + size());
		return getf(index);
	}

	/**
	 * Gets the result. Note that this uses the get(int) method from the backing PeakResultStore which may return stale
	 * data if index is outside of the current size.
	 *
	 * @param index
	 *            the index
	 * @return the peak result
	 */
	PeakResult getf(int index)
	{
		return this.results.get(index);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#size()
	 */
	public int size()
	{
		return this.results.size();
	}

	/**
	 * Add a result. Not synchronized.
	 *
	 * @param result
	 *            the result
	 */
	public void add(PeakResult result)
	{
		this.results.add(result);
	}

	/**
	 * Add all results.
	 * <p>
	 * Not synchronized. Use SynchronizedPeakResults to wrap this instance for use across threads.
	 * 
	 * {@inheritDoc}
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#addCollection(java.util.Collection)
	 */
	public void addAll(Collection<PeakResult> results)
	{
		this.results.addCollection(results);
	}

	/**
	 * Add all results.
	 * <p>
	 * Not synchronized. Use SynchronizedPeakResults to wrap this instance for use across threads.
	 * 
	 * @see gdsc.smlm.results.PeakResults#addAll(gdsc.smlm.results.PeakResult[])
	 */
	public void addAll(PeakResult[] results)
	{
		this.results.addArray(results);
	}

	/**
	 * Adds the results.
	 *
	 * @param results
	 *            the results
	 */
	public void add(MemoryPeakResults results)
	{
		this.results.addStore(results.results);
	}

	/**
	 * Clear the results.
	 */
	private void clear()
	{
		this.results.clear();
	}

	/**
	 * Trims the capacity of this instance to be the current size. An application can use this operation to minimize
	 * the storage of an instance.
	 */
	public void trimToSize()
	{
		this.results.trimToSize();
	}

	/**
	 * Sort the results.
	 */
	public void sort()
	{
		this.results.sort();
	}

	/**
	 * Sort the results.
	 *
	 * @param comparator
	 *            the comparator
	 */
	public void sort(Comparator<PeakResult> comparator)
	{
		this.results.sort(comparator);
	}

	/**
	 * Convert to an array. This is a new allocation of storage space.
	 *
	 * @return the peak result array
	 */
	public PeakResult[] toArray()
	{
		return this.results.toArray();
	}

	/**
	 * Removes the null results from the store.
	 */
	public void removeNullResults()
	{
		this.results.removeIf(new PeakResultPredicate()
		{
			public boolean test(PeakResult t)
			{
				return t == null;
			}
		});
	}

	/**
	 * Removes the result if it matches the filter. If objects are removed then the order of elements may change.
	 *
	 * @param filter
	 *            the filter
	 * @return true, if any were removed
	 */
	public boolean removeIf(PeakResultPredicate filter)
	{
		return this.results.removeIf(filter);
	}

	/////////////////////////////////////////////////////////////////
	// END OF RESULTS STORAGE METHODS 
	/////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////
	// START OF STATIC MEMORY STORAGE METHODS 
	/////////////////////////////////////////////////////////////////

	/**
	 * Instantiates a new memory peak results.
	 *
	 * @param store
	 *            the backing storage implementation
	 * @throws IllegalArgumentException
	 *             If the store is null
	 */
	public MemoryPeakResults(PeakResultStore store) throws IllegalArgumentException
	{
		if (store == null)
			throw new IllegalArgumentException("Store must not be null");
		results = store;
	}

	/**
	 * Instantiates a new memory peak results.
	 */
	public MemoryPeakResults()
	{
		this(1000);
	}

	/**
	 * Instantiates a new memory peak results.
	 *
	 * @param capacity
	 *            the capacity
	 */
	public MemoryPeakResults(int capacity)
	{
		// Use the fast and simple implementation for the store
		results = new ArrayPeakResultStore(capacity);
	}

	/**
	 * Instantiates a new memory peak results.
	 *
	 * @param results
	 *            the results
	 */
	public MemoryPeakResults(Collection<PeakResult> results)
	{
		this();
		addAll(results);
	}

	/**
	 * Instantiates a new memory peak results.
	 *
	 * @param psf
	 *            the psf
	 */
	public MemoryPeakResults(PSF psf)
	{
		this();
		setPSF(psf);
	}

	/**
	 * Gets the results.
	 *
	 * @param name
	 *            The name of the results
	 * @return Get the named results (or null if they do not exist)
	 */
	public static MemoryPeakResults getResults(String name)
	{
		return resultsMap.get(name);
	}

	/**
	 * Removes the results.
	 *
	 * @param name
	 *            The name of the results
	 * @return The removed results (or null if they do not exist)
	 */
	public static MemoryPeakResults removeResults(String name)
	{
		return resultsMap.remove(name);
	}

	/**
	 * Add the results to memory. The name is taken from the results.
	 *
	 * @param results
	 *            the results
	 */
	public static void addResults(MemoryPeakResults results)
	{
		if (results == null)
			throw new NullPointerException("Results must not be null");
		results.trimToSize();
		resultsMap.put(results.getName(), results);
	}

	/**
	 * Gets the result names.
	 *
	 * @return A set of the available named results held in memory
	 */
	public static Set<String> getResultNames()
	{
		return resultsMap.keySet();
	}

	/**
	 * Gets the all results.
	 *
	 * @return A collection of the results held in memory
	 */
	public static Collection<MemoryPeakResults> getAllResults()
	{
		return resultsMap.values();
	}

	/**
	 * Count the number of result sets in memory.
	 *
	 * @return the results memory size
	 */
	public static int getResultsMemorySize()
	{
		return resultsMap.size();
	}

	/**
	 * Return true if there are no non-empty results in memory.
	 *
	 * @return true, if is memory empty
	 */
	public static boolean isMemoryEmpty()
	{
		if (resultsMap.isEmpty())
			return true;
		for (MemoryPeakResults r : resultsMap.values())
			if (!r.isEmpty())
				return false;
		return true;
	}

	/**
	 * Count the total number of results in memory.
	 *
	 * @return the int
	 */
	public static int countMemorySize()
	{
		int size = 0;
		for (MemoryPeakResults r : resultsMap.values())
		{
			size += r.size();
		}
		return size;
	}

	/**
	 * Clear the results from memory.
	 */
	public static void clearMemory()
	{
		resultsMap.clear();
	}

	/**
	 * Estimate the total size of results in memory.
	 *
	 * @return the long
	 */
	public static long estimateMemorySize()
	{
		long memorySize = 0;
		for (MemoryPeakResults r : resultsMap.values())
		{
			memorySize += estimateMemorySize(r);
		}
		return memorySize;
	}

	/**
	 * Convert the size in bytes into a string.
	 *
	 * @param memorySize
	 *            the memory size
	 * @return The memory size string
	 */
	public static String memorySizeString(long memorySize)
	{
		return memorySize < 10000 * 1024 ? memorySize / 1024L + "K" : memorySize / 1048576L + "MB";
	}

	/**
	 * Return an estimate of the memory size taken by PeakResult objects.
	 * <p>
	 * Note: This is just a guess based on measured sizes for the objects in memory.
	 *
	 * @param r
	 *            the r
	 * @return The memory size
	 */
	public static long estimateMemorySize(MemoryPeakResults r)
	{
		long memorySize = 0;
		if (r != null && r.size() > 0)
		{
			boolean includeDeviations = r.getf(0).paramStdDevs != null;
			memorySize = MemoryPeakResults.estimateMemorySize(r.size(), includeDeviations);
		}
		return memorySize;
	}

	/**
	 * Return an estimate of the memory size taken by PeakResult objects.
	 * <p>
	 * Note: This is just a guess based on measured sizes for the objects in memory.
	 *
	 * @param size
	 *            the size
	 * @param includeDeviations
	 *            the include deviations
	 * @return The memory size
	 */
	public static long estimateMemorySize(int size, boolean includeDeviations)
	{
		if (byteSize == 0)
		{
			// Comment out to speed up the code
			//byteSize = (int) (measureSize(10000, false) / 10000);
			//byteSizeWithDeviations = (int) (measureSize(10000, true) / 10000);
			//System.out.printf("Size = %d,  Size with deviations = %d", byteSize, byteSizeWithDeviations);

			// Check just in case the estimate is bad
			if (byteSize <= 0)
				byteSize = DEFAULT_SIZE;
			if (byteSizeWithDeviations <= 0)
				byteSizeWithDeviations = DEFAULT_SIZE_WITH_DEVIATIONS;
		}
		return size * ((includeDeviations) ? byteSize : byteSizeWithDeviations);
	}

	// The following code can be used to determine the memory size of an object.
	// Taken from: http://www.javaworld.com/javaworld/javatips/jw-javatip130.html?page=1

	/**
	 * Measure size.
	 *
	 * @param size
	 *            the size
	 * @param includeDeviations
	 *            the include deviations
	 * @return the long
	 */
	public static long measureSize(int size, boolean includeDeviations)
	{
		// Warm up all classes/methods we will use
		runGC();
		usedMemory();
		// Array to keep strong references to allocated objects
		final int count = 1000;
		Object[] objects = new Object[count];

		long heap1 = 0;
		// Allocate count+1 objects, discard the first one
		for (int i = -1; i < count; ++i)
		{
			Object object = null;

			// Instantiate your data here and assign it to object

			object = new PeakResult(0, 1, 2, 3.0f, 4.0, 5.0f, new float[7], (includeDeviations) ? new float[7] : null);

			if (i >= 0)
				objects[i] = object;
			else
			{
				object = null; // Discard the warm up object
				runGC();
				heap1 = usedMemory(); // Take a before heap snapshot
			}
		}
		runGC();
		long heap2 = usedMemory(); // Take an after heap snapshot:

		long memorySize = Math.round(((double) (heap2 - heap1)) / count);
		//System.out.println("'before' heap: " + heap1 + ", 'after' heap: " + heap2);
		//System.out.println("heap delta: " + (heap2 - heap1) + ", {" + objects[0].getClass() + "} size = " + memorySize +
		//		" bytes");
		for (int i = 0; i < count; ++i)
			objects[i] = null;
		objects = null;
		runGC();

		return memorySize * size;
	}

	/**
	 * Run the garbage collector multiple times to free memory.
	 */
	public static void runGC()
	{
		// It helps to call Runtime.gc()
		// using several method calls:
		for (int r = 0; r < 4; ++r)
			_runGC();
	}

	/**
	 * Run GC.
	 */
	private static void _runGC()
	{
		long usedMem1 = usedMemory(), usedMem2 = Long.MAX_VALUE;
		for (int i = 0; (usedMem1 < usedMem2) && (i < 500); ++i)
		{
			runGCOnce();
			Thread.currentThread();
			Thread.yield();

			usedMem2 = usedMem1;
			usedMem1 = usedMemory();
		}
	}

	/**
	 * Run GC once.
	 */
	public static void runGCOnce()
	{
		s_runtime.runFinalization();
		s_runtime.gc();
	}

	/**
	 * Used memory.
	 *
	 * @return the long
	 */
	public static long usedMemory()
	{
		return s_runtime.totalMemory() - s_runtime.freeMemory();
	}

	/**
	 * Total memory.
	 *
	 * @return the long
	 */
	public static long totalMemory()
	{
		return s_runtime.totalMemory();
	}

	/**
	 * Free memory.
	 *
	 * @return the long
	 */
	public static long freeMemory()
	{
		return s_runtime.freeMemory();
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
	 * <p>
	 * This clears the current results but does not reduce storage allocation. This can be done with
	 * {@link #trimToSize()}.
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#begin()
	 */
	public void begin()
	{
		clear();
	}

	/**
	 * Add a result.
	 * <p>
	 * Not synchronized. Use SynchronizedPeakResults to wrap this instance for use across threads.
	 * 
	 * {@inheritDoc}
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#add(int, int, int, float, double, float, float[], float[])
	 */
	public void add(int peak, int origX, int origY, float origValue, double chiSquared, float noise, float[] params,
			float[] paramsStdDev)
	{
		add(new PeakResult(peak, origX, origY, origValue, chiSquared, noise, params, paramsStdDev));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#end()
	 */
	public void end()
	{
		if (isSortAfterEnd())
			sort();
	}

	/////////////////////////////////////////////////////////////////
	// END OF PeakResults interface METHODS 
	/////////////////////////////////////////////////////////////////

	/**
	 * Gets the bounds. These are returned in Pixel units as the bounds are defined in the PeakResults interface as the
	 * bounds used to create the results.
	 *
	 * @param calculate
	 *            Set to true to calculate the bounds if they are null or zero width/height
	 * @return the bounds of the result coordinates
	 * @throws DataException
	 *             if conversion to pixel units is not possible
	 */
	public Rectangle getBounds(boolean calculate) throws DataException
	{
		Rectangle bounds = getBounds();
		if ((bounds == null || bounds.width == 0 || bounds.height == 0) && calculate)
		{
			bounds = new Rectangle();
			// Note: The bounds should be in pixels
			Rectangle2D.Float b = getDataBounds(DistanceUnit.PIXEL);

			// Round to integer
			bounds.x = (int) Math.floor(b.x);
			bounds.y = (int) Math.floor(b.y);

			int maxX = (int) Math.ceil(b.x + b.width);
			int maxY = (int) Math.ceil(b.y + b.height);

			// For compatibility with drawing images add one to the limits if they are integers
			// Q. Is this still necessary since drawing images has been re-written to handle edge cases?
			//if (maxX == b.x + b.width)
			//	maxX += 1;
			//if (maxY == b.y + b.height)
			//	maxY += 1;

			bounds.width = maxX - bounds.x;
			bounds.height = maxY - bounds.y;

			setBounds(bounds);
		}
		return bounds;
	}

	/**
	 * Gets the data bounds.
	 *
	 * @param distanceUnit
	 *            the distance unit (if null then the data bounds will be in native units)
	 * @return the bounds of the result coordinates
	 * @throws DataException
	 *             if conversion to the required units is not possible
	 */
	public Rectangle2D.Float getDataBounds(DistanceUnit distanceUnit) throws DataException
	{
		if (isEmpty())
			return new Rectangle2D.Float();

		// Create this first to throw an exception if invalid
		final TypeConverter<DistanceUnit> c;
		if (distanceUnit == null)
		{
			c = new IdentityTypeConverter<DistanceUnit>(null);
		}
		else
		{
			c = CalibrationHelper.getDistanceConverter(getCalibration(), distanceUnit);
		}

		// Get the native bounds
		float minX = getf(0).getXPosition(), maxX = minX;
		float minY = getf(0).getYPosition(), maxY = minY;
		for (int i = 1, size = size(); i < size; i++)
		{
			PeakResult p = getf(i);
			float x = p.getXPosition();
			float y = p.getYPosition();
			if (minX > x)
				minX = x;
			else if (maxX < x)
				maxX = x;
			if (minY > y)
				minY = y;
			else if (maxY < y)
				maxY = y;
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#isActive()
	 */
	public boolean isActive()
	{
		return true;
	}

	/**
	 * Sets the sort after end.
	 *
	 * @param sortAfterEnd
	 *            True if the results should be sorted after the {@link #end()} method
	 */
	public void setSortAfterEnd(boolean sortAfterEnd)
	{
		this.sortAfterEnd = sortAfterEnd;
	}

	/**
	 * Checks if is sort after end.
	 *
	 * @return True if the results should be sorted after the {@link #end()} method
	 */
	public boolean isSortAfterEnd()
	{
		return sortAfterEnd;
	}

	/**
	 * Shallow copy this set of results. To create new object references use {@link #copy()}.
	 *
	 * @return the memory peak results
	 * @see java.lang.Object#clone()
	 */
	@Override
	public MemoryPeakResults clone()
	{
		try
		{
			return (MemoryPeakResults) super.clone();
		}
		catch (CloneNotSupportedException e)
		{
			// This should not happen so ignore
		}
		return null;
	}

	/**
	 * Copy the results. Create new objects for the properties (avoiding a shallow copy) but does not
	 * deep copy all of the peak results. Allows results to be resorted but not modified.
	 *
	 * @return the memory peak results
	 */
	public MemoryPeakResults copy()
	{
		MemoryPeakResults copy = clone();
		if (copy != null)
		{
			// Deep copy the objects that are not immutable
			if (getBounds() != null)
				copy.setBounds(new Rectangle(getBounds()));
			copy.results = results.copy();
		}
		return copy;
	}

	/**
	 * Checks if is empty.
	 *
	 * @return True if empty
	 */
	public boolean isEmpty()
	{
		return size() == 0;
	}

	/**
	 * Checks if is not empty.
	 *
	 * @return True if not empty
	 */
	public boolean isNotEmpty()
	{
		return size() != 0;
	}

	/**
	 * Checks for background. At least one result must have a non-zero background.
	 *
	 * @return true, if successful
	 */
	public boolean hasBackground()
	{
		for (int i = 0, size = size(); i < size; i++)
		{
			if (getf(i).getBackground() != 0)
				return true;
		}
		return false;
	}

	/**
	 * Checks for noise. At least one result must have a positive noise.
	 *
	 * @return true, if successful
	 */
	public boolean hasNoise()
	{
		for (int i = 0, size = size(); i < size; i++)
		{
			if (getf(i).noise > 0)
				return true;
		}
		return false;
	}

	/**
	 * Checks for intensity. At least one result must have a positive intensity.
	 *
	 * @return true, if successful
	 */
	public boolean hasIntensity()
	{
		for (int i = 0, size = size(); i < size; i++)
		{
			if (getf(i).getSignal() > 0)
				return true;
		}
		return false;
	}

	/**
	 * Checks for deviations. All results must have deviations.
	 *
	 * @return true, if successful
	 */
	public boolean hasDeviations()
	{
		for (int i = 0, size = size(); i < size; i++)
		{
			if (getf(i).paramStdDevs == null)
				return false;
		}
		return !isEmpty();
	}

	/**
	 * Checks for end frame. At least one result must have an end frame.
	 *
	 * @return true, if successful
	 */
	public boolean hasEndFrame()
	{
		for (int i = 0, size = size(); i < size; i++)
		{
			if (getf(i).getFrame() != getf(i).getEndFrame())
				return true;
		}
		return false;
	}

	/**
	 * Checks for id. At least one result must have a non-zero ID.
	 *
	 * @return true, if successful
	 */
	public boolean hasId()
	{
		for (int i = 0, size = size(); i < size; i++)
		{
			if (getf(i).getId() != 0)
				return true;
		}
		return false;
	}

	/**
	 * Checks for id. At least two results must have different non-zero Ids. If the number of results is 1 then the id
	 * must not be zero.
	 *
	 * @return true, if successful
	 */
	public boolean hasMultipleId()
	{
		if (isEmpty())
			return false;
		if (size() == 1)
			return getf(0).getId() != 0;
		for (int i = 0, size = size(); i < size; i++)
		{
			int id = getf(i).getId();
			if (id != 0)
			{
				while (++i < size)
				{
					int id2 = getf(i).getId();
					if (id2 != 0 && id2 != id)
						return true;
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
	public boolean hasPrecision()
	{
		for (int i = 0, size = size(); i < size; i++)
			if (!getf(i).hasPrecision())
				return false;
		return true;
	}

	/**
	 * Checks for 3D results. At least 1 result must have a non-zero z-coordinate.
	 *
	 * @return true, if successful
	 */
	public boolean is3D()
	{
		for (int i = 0, size = size(); i < size; i++)
			if (getf(i).getZPosition() != 0)
				return true;
		return false;
	}

	/**
	 * Gets the first frame.
	 * <p>
	 * This may be different from {@link #getMinFrame()} if the results are not sorted by frame.
	 *
	 * @return the first frame
	 * @throws IllegalStateException
	 *             If the size is zero
	 */
	public int getFirstFrame()
	{
		if (isEmpty())
			throw new IllegalStateException("Empty");
		return getf(0).getFrame();
	}

	/**
	 * Gets the last frame.
	 * <p>
	 * This may be different from {@link #getMaxFrame()} if the results are not sorted by frame.
	 *
	 * @return the last frame
	 * @throws IllegalStateException
	 *             If the size is zero
	 */
	public int getLastFrame()
	{
		if (isEmpty())
			throw new IllegalStateException("Empty");
		return getf(size() - 1).getEndFrame();
	}

	/**
	 * Gets the minimum frame.
	 *
	 * @return the min frame
	 * @throws IllegalStateException
	 *             If the size is zero
	 */
	public int getMinFrame()
	{
		if (isEmpty())
			throw new IllegalStateException("Empty");
		int min = getf(0).getFrame();
		for (int i = 1, size = size(); i < size; i++)
			if (min > getf(i).getFrame())
				min = getf(i).getFrame();
		return min;
	}

	/**
	 * Gets the maximum frame.
	 *
	 * @return the max frame
	 * @throws IllegalStateException
	 *             If the size is zero
	 */
	public int getMaxFrame()
	{
		if (isEmpty())
			throw new IllegalStateException("Empty");
		int max = getf(0).getEndFrame();
		for (int i = 1, size = size(); i < size; i++)
			if (max < getf(i).getEndFrame())
				max = getf(i).getEndFrame();
		return max;
	}

	/**
	 * Checks if all results have a stored precision value.
	 *
	 * @return true, if all results have a stored precision value
	 */
	public boolean hasStoredPrecision()
	{
		for (int i = 0; i < size(); i++)
		{
			if (!getf(i).hasPrecision())
				return false;
		}
		return true;
	}

	/**
	 * Checks for null results in the store.
	 *
	 * @return true, if null PeakResult object(s) exist
	 */
	public boolean hasNullResults()
	{
		for (int i = 0; i < size(); i++)
		{
			if (getf(i) == null)
				return true;
		}
		return false;
	}

	/** The preferred distance unit */
	public static final DistanceUnit PREFERRED_DISTANCE_UNIT = DistanceUnit.PIXEL;
	/** The preferred intensity unit */
	public static final IntensityUnit PREFERRED_INTENSITY_UNIT = IntensityUnit.PHOTON;
	/** The preferred angle unit */
	public static final AngleUnit PREFERRED_ANGLE_UNIT = AngleUnit.RADIAN;

	/**
	 * Checks if is distance in preferred units.
	 *
	 * @return true, if is distance in preferred units
	 */
	public boolean isDistanceInPreferredUnits()
	{
		return getDistanceUnit() == PREFERRED_DISTANCE_UNIT;
	}

	/**
	 * Checks if is intensity in preferred units.
	 *
	 * @return true, if is intensity in preferred units
	 */
	public boolean isIntensityInPreferredUnits()
	{
		return (getIntensityUnit() == PREFERRED_INTENSITY_UNIT);
	}

	/**
	 * Checks if is angle in preferred units.
	 *
	 * @return true, if is angle in preferred units
	 */
	public boolean isAngleInPreferredUnits()
	{
		return getAngleUnit() == PREFERRED_ANGLE_UNIT;
	}

	/**
	 * Convert to preferred units.
	 *
	 * @return true, if the data is now stored in the preferred units.
	 */
	public boolean convertToPreferredUnits()
	{
		return convertToUnits(PREFERRED_DISTANCE_UNIT, PREFERRED_INTENSITY_UNIT, PREFERRED_ANGLE_UNIT);
	}

	/**
	 * Convert to the specified units. If the units are null they will remain unchanged.
	 *
	 * @param distanceUnit
	 *            the distance unit
	 * @param intensityUnit
	 *            the intensity unit
	 * @param angleUnit
	 *            the angle unit
	 * @return true, if the data is now stored in the preferred units.
	 */
	public boolean convertToUnits(DistanceUnit distanceUnit, IntensityUnit intensityUnit, AngleUnit angleUnit)
	{
		if (!hasCalibration())
			return false;

		PeakResultsHelper helper = new PeakResultsHelper(getCalibration(), getPSF());
		helper.setIntensityUnit(intensityUnit);
		helper.setDistanceUnit(distanceUnit);
		helper.setAngleUnit(angleUnit);
		final Converter[] converters = helper.getConverters();

		if (!helper.isCalibrationChanged())
		{
			// Check if already in the specified units
			return helper.isValidConversion();
		}

		// Update the calibration
		setCalibration(helper.getCalibration());

		// We must convert the noise
		Converter noiseConverter = converters[PeakResult.INTENSITY];

		for (int i = 0, size = size(); i < size; i++)
		{
			PeakResult p = getf(i);
			p.noise = noiseConverter.convert(p.noise);
			final float[] params = p.params;
			final float[] paramsStdDev = p.paramStdDevs;
			if (paramsStdDev == null)
			{
				for (int j = 0; j < converters.length; j++)
					params[j] = converters[j].convert(params[j]);
			}
			else
			{
				for (int j = 0; j < converters.length; j++)
				{
					params[j] = converters[j].convert(params[j]);
					paramsStdDev[j] = converters[j].convert(paramsStdDev[j]);
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
	 * For each result execute the procedure.
	 * <p>
	 * Warning: Results with be in their native units since no unit conversion is performed.
	 *
	 * @param procedure
	 *            the procedure
	 */
	public void forEach(PeakResultProcedure procedure)
	{
		for (int i = 0, size = size(); i < size; i++)
		{
			procedure.execute(getf(i));
		}
	}

	/**
	 * For the first result execute the procedure.
	 *
	 * @param procedure
	 *            the procedure
	 */
	public void forFirst(PeakResultProcedure procedure)
	{
		if (isEmpty())
			return;
		procedure.execute(getf(0));
	}

	/**
	 * For each result execute the fast-exit procedure.
	 *
	 * @param procedure
	 *            the procedure
	 * @return true, if a fast exit occurred
	 */
	public boolean forEach(PeakResultProcedureX procedure)
	{
		for (int i = 0, size = size(); i < size; i++)
		{
			if (procedure.execute(getf(i)))
				return true;
		}
		return false;
	}

	/**
	 * For each result execute the procedure.
	 * <p>
	 * Warning: Results with be in their native units since no unit conversion is performed.
	 *
	 * @param procedure
	 *            the procedure
	 */
	public void forEachNative(BIXYZResultProcedure procedure)
	{
		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = getf(i);
			procedure.executeBIXYZ(r.getBackground(), r.getSignal(), r.getXPosition(), r.getYPosition(),
					r.getZPosition());
		}
	}

	/**
	 * For the first result execute the procedure.
	 * <p>
	 * Warning: Results with be in their native units since no unit conversion is performed.
	 *
	 * @param procedure
	 *            the procedure
	 */
	public void forFirstNative(BIXYZResultProcedure procedure)
	{
		if (isEmpty())
			return;
		final PeakResult r = getf(0);
		procedure.executeBIXYZ(r.getBackground(), r.getSignal(), r.getXPosition(), r.getYPosition(), r.getZPosition());
	}

	/**
	 * For each result execute the procedure using the specified units.
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 *
	 * @param intensityUnit
	 *            the intensity unit
	 * @param procedure
	 *            the procedure
	 * @throws ConversionException
	 *             if the conversion is not possible
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public void forEach(IntensityUnit intensityUnit, BResultProcedure procedure)
			throws ConversionException, ConfigurationException
	{
		if (!hasCalibration())
			throw new ConfigurationException("No calibration");

		TypeConverter<IntensityUnit> ic = getCalibrationReader().getIntensityConverter(intensityUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			procedure.executeB(ic.convert(getf(i).getBackground()));
		}
	}

	/**
	 * For each result execute the procedure using the specified units.
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 *
	 * @param intensityUnit
	 *            the intensity unit
	 * @param distanceUnit
	 *            the distance unit
	 * @param procedure
	 *            the procedure
	 * @throws ConversionException
	 *             if the conversion is not possible
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public void forEach(IntensityUnit intensityUnit, DistanceUnit distanceUnit, BIXYResultProcedure procedure)
			throws ConversionException, ConfigurationException
	{
		if (!hasCalibration())
			throw new ConfigurationException("No calibration");

		TypeConverter<IntensityUnit> ic = getCalibrationReader().getIntensityConverter(intensityUnit);
		TypeConverter<DistanceUnit> dc = getCalibrationReader().getDistanceConverter(distanceUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = getf(i);
			//@formatter:off
			procedure.executeBIXY(
					ic.convert(r.getBackground()), 
					ic.convert(r.getSignal()), 
					dc.convert(r.getXPosition()),
					dc.convert(r.getYPosition()));
			//@formatter:on
		}
	}

	/**
	 * For each result execute the procedure using the specified units.
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 *
	 * @param intensityUnit
	 *            the intensity unit
	 * @param distanceUnit
	 *            the distance unit
	 * @param procedure
	 *            the procedure
	 * @throws ConversionException
	 *             if the conversion is not possible
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public void forEach(IntensityUnit intensityUnit, DistanceUnit distanceUnit, BIXYZResultProcedure procedure)
			throws ConversionException, ConfigurationException
	{
		if (!hasCalibration())
			throw new ConfigurationException("No calibration");

		TypeConverter<IntensityUnit> ic = getCalibrationReader().getIntensityConverter(intensityUnit);
		TypeConverter<DistanceUnit> dc = getCalibrationReader().getDistanceConverter(distanceUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = getf(i);
			//@formatter:off
			procedure.executeBIXYZ(
					ic.convert(r.getBackground()), 
					ic.convert(r.getSignal()), 
					dc.convert(r.getXPosition()),
					dc.convert(r.getYPosition()), 
					dc.convert(r.getZPosition()));
			//@formatter:on
		}
	}

	/**
	 * For each result execute the procedure using the specified units.
	 * <p>
	 * This procedure is exclusive to data fit with a Gaussian2D PSF as it computes the height of the Gaussian from the
	 * integral (Intensity) and the widths.
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 *
	 * @param intensityUnit
	 *            the intensity unit
	 * @param procedure
	 *            the procedure
	 * @throws ConversionException
	 *             if the conversion is not possible
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public void forEach(IntensityUnit intensityUnit, HResultProcedure procedure)
			throws ConversionException, ConfigurationException
	{
		if (!hasCalibration())
			throw new ConfigurationException("No calibration");

		int[] indices = PSFHelper.getGaussian2DWxWyIndices(getPSF());

		final int isx = indices[0];
		final int isy = indices[1];
		final double twoPi = 2 * Math.PI;

		TypeConverter<IntensityUnit> ic = getCalibrationReader().getIntensityConverter(intensityUnit);
		TypeConverter<DistanceUnit> dc = getCalibrationReader().getDistanceConverter(DistanceUnit.PIXEL);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = getf(i);

			// Convert the widths to pixels
			float sx = dc.convert(r.getParameter(isx));
			float sy = dc.convert(r.getParameter(isy));

			//@formatter:off
			procedure.executeH(
					(float)(ic.convert((double)r.getSignal()) / (twoPi * sx * sy)));
			//@formatter:on
		}
	}

	/**
	 * For each result execute the procedure using the specified units.
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 *
	 * @param intensityUnit
	 *            the intensity unit
	 * @param procedure
	 *            the procedure
	 * @throws ConversionException
	 *             if the conversion is not possible
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public void forEach(IntensityUnit intensityUnit, IResultProcedure procedure)
			throws ConversionException, ConfigurationException
	{
		if (!hasCalibration())
			throw new ConfigurationException("No calibration");

		TypeConverter<IntensityUnit> ic = getCalibrationReader().getIntensityConverter(intensityUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = getf(i);
			//@formatter:off
			procedure.executeI(
					ic.convert(r.getSignal()));
			//@formatter:on
		}
	}

	/**
	 * For each result execute the procedure using the specified units.
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 *
	 * @param intensityUnit
	 *            the intensity unit
	 * @param distanceUnit
	 *            the distance unit
	 * @param procedure
	 *            the procedure
	 * @throws ConversionException
	 *             if the conversion is not possible
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public void forEach(IntensityUnit intensityUnit, DistanceUnit distanceUnit, IXYResultProcedure procedure)
			throws ConversionException, ConfigurationException
	{
		if (!hasCalibration())
			throw new ConfigurationException("No calibration");

		TypeConverter<IntensityUnit> ic = getCalibrationReader().getIntensityConverter(intensityUnit);
		TypeConverter<DistanceUnit> dc = getCalibrationReader().getDistanceConverter(distanceUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = getf(i);
			//@formatter:off
			procedure.executeIXY(
					ic.convert(r.getSignal()), 
					dc.convert(r.getXPosition()),
					dc.convert(r.getYPosition()));
			//@formatter:on
		}
	}

	/**
	 * For each result execute the procedure using the specified units.
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 *
	 * @param intensityUnit
	 *            the intensity unit
	 * @param distanceUnit
	 *            the distance unit
	 * @param procedure
	 *            the procedure
	 * @throws ConversionException
	 *             if the conversion is not possible
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public void forEach(IntensityUnit intensityUnit, DistanceUnit distanceUnit, IXYRResultProcedure procedure)
			throws ConversionException, ConfigurationException
	{
		if (!hasCalibration())
			throw new ConfigurationException("No calibration");

		TypeConverter<IntensityUnit> ic = getCalibrationReader().getIntensityConverter(intensityUnit);
		TypeConverter<DistanceUnit> dc = getCalibrationReader().getDistanceConverter(distanceUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = getf(i);
			//@formatter:off
			procedure.executeIXYR(
					ic.convert(r.getSignal()), 
					dc.convert(r.getXPosition()),
					dc.convert(r.getYPosition()),
					r);
			//@formatter:on
		}
	}

	/**
	 * For each result execute the procedure using the specified units.
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 *
	 * @param intensityUnit
	 *            the intensity unit
	 * @param distanceUnit
	 *            the distance unit
	 * @param procedure
	 *            the procedure
	 * @throws ConversionException
	 *             if the conversion is not possible
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public void forEach(IntensityUnit intensityUnit, DistanceUnit distanceUnit, IXYZResultProcedure procedure)
			throws ConversionException, ConfigurationException
	{
		if (!hasCalibration())
			throw new ConfigurationException("No calibration");

		TypeConverter<IntensityUnit> ic = getCalibrationReader().getIntensityConverter(intensityUnit);
		TypeConverter<DistanceUnit> dc = getCalibrationReader().getDistanceConverter(distanceUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = getf(i);
			//@formatter:off
			procedure.executeIXYZ(
					ic.convert(r.getSignal()), 
					dc.convert(r.getXPosition()),
					dc.convert(r.getYPosition()),
					dc.convert(r.getZPosition()));
			//@formatter:on
		}
	}

	/**
	 * For each result execute the procedure using the specified units.
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 *
	 * @param procedure
	 *            the procedure
	 * @throws ConversionException
	 *             if the conversion is not possible
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public void forEach(TResultProcedure procedure) throws ConfigurationException
	{
		if (!hasCalibration())
			throw new ConfigurationException("No calibration");

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = getf(i);
			//@formatter:off
			procedure.executeT(
					r.getFrame());
			//@formatter:on
		}
	}

	/**
	 * For each result execute the procedure using the specified units.
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 *
	 * @param distanceUnit
	 *            the distance unit
	 * @param procedure
	 *            the procedure
	 * @throws ConversionException
	 *             if the conversion is not possible
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public void forEach(DistanceUnit distanceUnit, TXYResultProcedure procedure)
			throws ConversionException, ConfigurationException
	{
		if (!hasCalibration())
			throw new ConfigurationException("No calibration");

		TypeConverter<DistanceUnit> dc = getCalibrationReader().getDistanceConverter(distanceUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = getf(i);
			//@formatter:off
			procedure.executeTXY(
					r.getFrame(),
					dc.convert(r.getXPosition()),
					dc.convert(r.getYPosition()));
			//@formatter:on
		}
	}

	/**
	 * For each result execute the procedure using the specified units.
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 *
	 * @param distanceUnit
	 *            the distance unit
	 * @param procedure
	 *            the procedure
	 * @throws ConversionException
	 *             if the conversion is not possible
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public void forEach(DistanceUnit distanceUnit, WResultProcedure procedure)
			throws ConversionException, ConfigurationException
	{
		if (!hasCalibration())
			throw new ConfigurationException("No calibration");

		// Note that in the future we may support more than just Gaussian2D PSF
		// so this may have to change

		int[] indices = PSFHelper.getGaussian2DWxWyIndices(getPSF());

		final int isx = indices[0];
		final int isy = indices[1];

		TypeConverter<DistanceUnit> dc = getCalibrationReader().getDistanceConverter(distanceUnit);

		if (isx == isy)
		{
			for (int i = 0, size = size(); i < size; i++)
			{
				final PeakResult r = getf(i);
				//@formatter:off
    			procedure.executeW(
    					dc.convert(r.getParameter(isx)));
    			//@formatter:on
			}
		}
		else
		{
			for (int i = 0, size = size(); i < size; i++)
			{
				final PeakResult r = getf(i);
				// Convert the separate widths into a single width
				double s = Gaussian2DPeakResultHelper.getStandardDeviation(r.getParameter(isx), r.getParameter(isy));
				//@formatter:off
    			procedure.executeW(
    					(float)dc.convert(s));
    			//@formatter:on
			}
		}
	}

	/**
	 * For each result execute the procedure using the specified units.
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 *
	 * @param distanceUnit
	 *            the distance unit
	 * @param procedure
	 *            the procedure
	 * @throws ConversionException
	 *             if the conversion is not possible
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public void forEach(DistanceUnit distanceUnit, WxWyResultProcedure procedure)
			throws ConversionException, ConfigurationException
	{
		if (!hasCalibration())
			throw new ConfigurationException("No calibration");

		// Note that in the future we may support more than just Gaussian2D PSF
		// so this may have to change

		int[] indices = PSFHelper.getGaussian2DWxWyIndices(getPSF());

		final int isx = indices[0];
		final int isy = indices[1];

		TypeConverter<DistanceUnit> dc = getCalibrationReader().getDistanceConverter(distanceUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
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
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 *
	 * @param distanceUnit
	 *            the distance unit
	 * @param procedure
	 *            the procedure
	 * @throws ConversionException
	 *             if the conversion is not possible
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public void forEach(DistanceUnit distanceUnit, XYResultProcedure procedure)
			throws ConversionException, ConfigurationException
	{
		if (!hasCalibration())
			throw new ConfigurationException("No calibration");

		TypeConverter<DistanceUnit> dc = getCalibrationReader().getDistanceConverter(distanceUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = getf(i);
			//@formatter:off
			procedure.executeXY(
					dc.convert(r.getXPosition()),
					dc.convert(r.getYPosition()));
			//@formatter:on
		}
	}

	/**
	 * For each result execute the procedure using the specified units.
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 * <p>
	 * Warning: The peak result with be in native units.
	 *
	 * @param distanceUnit
	 *            the distance unit
	 * @param procedure
	 *            the procedure
	 * @throws ConversionException
	 *             if the conversion is not possible
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public void forEach(DistanceUnit distanceUnit, XYRResultProcedure procedure)
			throws ConversionException, ConfigurationException
	{
		if (!hasCalibration())
			throw new ConfigurationException("No calibration");

		TypeConverter<DistanceUnit> dc = getCalibrationReader().getDistanceConverter(distanceUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = getf(i);
			//@formatter:off
			procedure.executeXYR(
					dc.convert(r.getXPosition()),
					dc.convert(r.getYPosition()),
					r);
			//@formatter:on
		}
	}

	/**
	 * For each result execute the procedure using the specified units.
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 *
	 * @param distanceUnit
	 *            the distance unit
	 * @param procedure
	 *            the procedure
	 * @throws ConversionException
	 *             if the conversion is not possible
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public void forEach(DistanceUnit distanceUnit, XYZResultProcedure procedure)
			throws ConversionException, ConfigurationException
	{
		if (!hasCalibration())
			throw new ConfigurationException("No calibration");

		TypeConverter<DistanceUnit> dc = getCalibrationReader().getDistanceConverter(distanceUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = getf(i);
			//@formatter:off
			procedure.executeXYZ(
					dc.convert(r.getXPosition()),
					dc.convert(r.getYPosition()),
					dc.convert(r.getZPosition()));
			//@formatter:on
		}
	}

	/**
	 * For each result execute the procedure
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 *
	 * @param procedure
	 *            the procedure
	 * @throws ConversionException
	 *             if the conversion is not possible
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public void forEach(LSEPrecisionProcedure procedure) throws ConversionException, ConfigurationException
	{
		if (!hasCalibration())
			throw new ConfigurationException("No calibration");

		Gaussian2DPeakResultCalculator calculator = Gaussian2DPeakResultHelper.create(getPSF(), getCalibration(),
				Gaussian2DPeakResultHelper.LSE_PRECISION);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = getf(i);
			procedure.executeLSEPrecision(calculator.getLSEPrecision(r.params, r.noise));
		}
	}

	/**
	 * For each result execute the procedure
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 *
	 * @param procedure
	 *            the procedure
	 * @throws ConversionException
	 *             if the conversion is not possible
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public void forEach(LSEPrecisionBProcedure procedure) throws ConversionException, ConfigurationException
	{
		if (!hasCalibration())
			throw new ConfigurationException("No calibration");

		Gaussian2DPeakResultCalculator calculator = Gaussian2DPeakResultHelper.create(getPSF(), getCalibration(),
				Gaussian2DPeakResultHelper.LSE_PRECISION_X);

		for (int i = 0, size = size(); i < size; i++)
		{
			procedure.executeLSEPrecisionB(calculator.getLSEPrecision(getf(i).params));
		}
	}

	/**
	 * For each result execute the procedure
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 *
	 * @param procedure
	 *            the procedure
	 * @throws ConversionException
	 *             if the conversion is not possible
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public void forEach(MLEPrecisionProcedure procedure) throws ConversionException, ConfigurationException
	{
		if (!hasCalibration())
			throw new ConfigurationException("No calibration");

		Gaussian2DPeakResultCalculator calculator = Gaussian2DPeakResultHelper.create(getPSF(), getCalibration(),
				Gaussian2DPeakResultHelper.MLE_PRECISION);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = getf(i);
			procedure.executeMLEPrecision(calculator.getMLEPrecision(r.params, r.noise));
		}
	}

	/**
	 * For each result execute the procedure
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 *
	 * @param procedure
	 *            the procedure
	 * @throws ConversionException
	 *             if the conversion is not possible
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public void forEach(MLEPrecisionBProcedure procedure) throws ConversionException, ConfigurationException
	{
		if (!hasCalibration())
			throw new ConfigurationException("No calibration");

		Gaussian2DPeakResultCalculator calculator = Gaussian2DPeakResultHelper.create(getPSF(), getCalibration(),
				Gaussian2DPeakResultHelper.MLE_PRECISION_X);

		for (int i = 0, size = size(); i < size; i++)
		{
			procedure.executeMLEPrecisionB(calculator.getMLEPrecision(getf(i).params));
		}
	}

	/**
	 * Apply a pixel translation to the results. The original coordinates are updated using a raw pixel shift. The
	 * current XY coordinates are updated by converting the pixel shift to the current distance units. If there are
	 * current bounds then the origin will be updated with a raw pixel shift.
	 *
	 * @param x
	 *            the x pixel shift
	 * @param y
	 *            the y pixel shift
	 * @throws ConversionException
	 *             if the conversion is not possible
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public void translate(int x, int y) throws ConversionException, ConfigurationException
	{
		if (x == 0 && y == 0)
			return;

		if (!hasCalibration())
			throw new ConfigurationException("No calibration");

		Rectangle bounds = getBounds();
		if (bounds != null)
		{
			bounds.x += x;
			bounds.y += y;
			setBounds(bounds);
		}

		// Convert the pixel shift to the units of the results
		TypeConverter<DistanceUnit> dc = getCalibrationReader().getDistanceConverter(getDistanceUnit());
		float xx = dc.convert(x);
		float yy = dc.convert(y);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = getf(i);
			r.origX += x;
			r.origY += y;
			r.params[PeakResult.X] += xx;
			r.params[PeakResult.Y] += yy;
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
	public DistanceUnit getDistanceUnit()
	{
		if (hasCalibration())
		{
			if (getCalibrationReader().hasDistanceUnit())
				return getCalibrationReader().getDistanceUnit();
		}
		return null;
	}

	/**
	 * Gets the intensity unit.
	 *
	 * @return the intensity unit
	 */
	public IntensityUnit getIntensityUnit()
	{
		if (hasCalibration())
		{
			if (getCalibrationReader().hasIntensityUnit())
				return getCalibrationReader().getIntensityUnit();
		}
		return null;
	}

	/**
	 * Gets the angle unit.
	 *
	 * @return the angle unit
	 */
	public AngleUnit getAngleUnit()
	{
		if (hasCalibration())
		{
			if (getCalibrationReader().hasAngleUnit())
				return getCalibrationReader().getAngleUnit();
		}
		return null;
	}

	/**
	 * Fix zero background to the given background.
	 *
	 * @param newBackground
	 *            the new background
	 */
	public void setZeroBackground(float newBackground)
	{
		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = getf(i);
			if (r.params[0] == 0)
				r.params[0] = newBackground;
		}
	}

	/**
	 * Creates the frame counter. It will be initialised with the value of {@link #getFirstFrame()} - 1. This ensures
	 * that the frame from the first result will be recognised as a new frame.
	 *
	 * @return the frame counter
	 */
	public FrameCounter newFrameCounter()
	{
		return new FrameCounter((isEmpty()) ? 0 : getFirstFrame());
	}
}
