package gdsc.smlm.results;

import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Set;

import gdsc.smlm.data.config.CalibrationHelper;
import gdsc.smlm.data.config.ConfigurationException;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.SMLMSettings.AngleUnit;
import gdsc.smlm.data.config.SMLMSettings.DistanceUnit;
import gdsc.smlm.data.config.SMLMSettings.IntensityUnit;
import gdsc.core.data.utils.ConversionException;
import gdsc.core.data.utils.TypeConverter;
import gdsc.smlm.data.config.UnitConverterFactory;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.results.procedures.IResultProcedure;
import gdsc.smlm.results.procedures.IXYResultProcedure;
import gdsc.smlm.results.procedures.IXYZResultProcedure;
import gdsc.smlm.results.procedures.LSEPrecisionProcedure;
import gdsc.smlm.results.procedures.PeakResultProcedure;
import gdsc.smlm.results.procedures.PeakResultProcedureX;
import gdsc.smlm.results.procedures.StandardResultProcedure;
import gdsc.smlm.results.procedures.TXYResultProcedure;
import gdsc.smlm.results.procedures.WResultProcedure;
import gdsc.smlm.results.procedures.WxWyResultProcedure;
import gdsc.smlm.results.procedures.XYRResultProcedure;
import gdsc.smlm.results.procedures.BIXYResultProcedure;
import gdsc.smlm.results.procedures.BIXYZResultProcedure;
import gdsc.smlm.results.procedures.HResultProcedure;
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
	private ArrayList<PeakResult> results;

	/**
	 * Gets the result.
	 *
	 * @param index
	 *            the index
	 * @return the peak result
	 */
	private PeakResult get(int index)
	{
		return results.get(index);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#size()
	 */
	public int size()
	{
		return results.size();
	}

	/**
	 * Add a result. Not synchronized.
	 *
	 * @param result
	 *            the result
	 */
	public void add(PeakResult result)
	{
		add(result);
	}

	/**
	 * Add all results.
	 * <p>
	 * Not synchronized. Use SynchronizedPeakResults to wrap this instance for use across threads.
	 * 
	 * {@inheritDoc}
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#addAll(java.util.Collection)
	 */
	public void addAll(Collection<PeakResult> results)
	{
		this.results.addAll(results);
	}

	/**
	 * Clear the results.
	 */
	private void clear()
	{
		results.clear();
	}

	/**
	 * Trims the capacity of this instance to be the current size. An application can use this operation to minimize
	 * the storage of an instance.
	 */
	public void trimToSize()
	{
		results.trimToSize();
	}

	/**
	 * Sort the results.
	 */
	public void sort()
	{
		Collections.sort(results);
	}

	/**
	 * Convert to an array.
	 *
	 * @return the peak result array
	 */
	public PeakResult[] toArray()
	{
		return results.toArray(new PeakResult[size()]);
	}

	/**
	 * Convert to an array reusing the space if provided.
	 *
	 * @param array
	 *            the array (can be null)
	 * @return the peak result array
	 */
	public PeakResult[] toArray(PeakResult[] array)
	{
		if (array == null || array.length < size())
			return toArray();
		return results.toArray(array);
	}

	/**
	 * Copy results the the copy instance.
	 *
	 * @param copy
	 *            the copy
	 */
	private void copyResults(MemoryPeakResults copy)
	{
		copy.results = new ArrayList<PeakResult>(results);
	}

	/**
	 * Removes the null results from the store.
	 */
	public void removeNullResults()
	{
		ArrayList<PeakResult> list = new ArrayList<PeakResult>(size());
		for (int i = 0, size = size(); i < size; i++)
		{
			PeakResult p = get(i);
			if (p != null)
				list.add(p);
		}
		this.results = list;
	}

	/////////////////////////////////////////////////////////////////
	// END OF RESULTS STORAGE METHODS 
	/////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////
	// START OF STATIC MEMORY STORAGE METHODS 
	/////////////////////////////////////////////////////////////////

	/**
	 * Instantiates a new memory peak results.
	 */
	public MemoryPeakResults()
	{
		results = new ArrayList<PeakResult>(1000);
	}

	/**
	 * Instantiates a new memory peak results.
	 *
	 * @param capacity
	 *            the capacity
	 */
	public MemoryPeakResults(int capacity)
	{
		results = new ArrayList<PeakResult>(capacity);
	}

	/**
	 * Instantiates a new memory peak results.
	 *
	 * @param results
	 *            the results
	 */
	public MemoryPeakResults(Collection<PeakResult> results)
	{
		addAll(results);
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
			boolean includeDeviations = r.get(0).paramsStdDev != null;
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
	 * bounds used to create the resuls.
	 *
	 * @param calculate
	 *            Set to true to calculate the bounds if they are null or zero width/height
	 * @return the bounds of the result coordinates
	 * @throws DataException
	 *             if conversion to pixel units is not possible
	 */
	public Rectangle getBounds(boolean calculate)
	{
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
		}
		return bounds;
	}

	/**
	 * Gets the data bounds.
	 *
	 * @param distanceUnit
	 *            the distance unit
	 * @return the bounds of the result coordinates
	 * @throws DataException
	 *             if conversion to the required units is not possible
	 */
	public Rectangle2D.Float getDataBounds(DistanceUnit distanceUnit)
	{
		if (isEmpty())
			return new Rectangle2D.Float();

		StandardResultProcedure p = new StandardResultProcedure(this, distanceUnit);
		p.getXY();

		float minX = p.x[0], maxX = minX;
		float minY = p.y[0], maxY = minY;
		for (int i = 1, size = size(); i < size; i++)
		{
			float x = p.x[i];
			float y = p.y[i];
			if (minX > x)
				minX = x;
			else if (maxX < x)
				maxX = x;
			if (minY > y)
				minY = y;
			else if (maxY < y)
				maxY = y;
		}
		return new Rectangle2D.Float(minX, minY, maxX - minX, maxY - minY);
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
			// Deep copy the objects
			if (bounds != null)
				copy.bounds = new Rectangle(bounds);
			if (calibration != null)
				copy.calibration = calibration.clone();
			copyResults(copy);
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
	 * Checks for background. At least one result must have a non-zero background.
	 *
	 * @return true, if successful
	 */
	public boolean hasBackground()
	{
		for (int i = 0, size = size(); i < size; i++)
		{
			if (get(i).getBackground() != 0)
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
			if (get(i).noise > 0)
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
			if (get(i).getSignal() > 0)
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
			if (get(i).paramsStdDev == null)
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
			if (get(i).getFrame() != get(i).getEndFrame())
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
			if (get(i).getId() != 0)
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
			return get(0).getId() != 0;
		for (int i = 0, size = size(); i < size; i++)
		{
			int id = get(i).getId();
			if (id != 0)
			{
				while (++i < size)
				{
					int id2 = get(i).getId();
					if (id2 != 0 && id2 != id)
						return true;
				}
				break;
			}
		}
		return false;
	}

	/**
	 * Gets the first frame.
	 *
	 * @return the first frame
	 * @throws IllegalStateException
	 *             If the size is zero
	 */
	public int getFirstFrame()
	{
		if (isEmpty())
			throw new IllegalStateException("Empty");
		return get(0).getFrame();
	}

	/**
	 * Gets the last frame.
	 *
	 * @return the last frame
	 * @throws IllegalStateException
	 *             If the size is zero
	 */
	public int getLastFrame()
	{
		if (isEmpty())
			throw new IllegalStateException("Empty");
		return get(size() - 1).getEndFrame();
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
			if (!get(i).hasPrecision())
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
			if (get(i) == null)
				return true;
		}
		return false;
	}

	/**
	 * Convert to preferred units.
	 *
	 * @return true, if the data is now stored in the preferred units.
	 */
	public boolean convertToPreferredUnits()
	{
		if (calibration == null)
			return false;
		// TODO - use the helper
		CalibrationHelper helper = null; //new CalibrationHelper(calibration);
		boolean success = convertDistanceToPixelUnits(helper);
		success &= convertIntensityToPhotonUnits(helper);
		success &= convertAngleToRadianUnits(helper);
		//setCalibration(helper.getCalibration());
		return success;
	}

	/**
	 * Checks if distance is in pixel units.
	 *
	 * @return true, if distance is in pixels
	 */
	public boolean isDistanceInPixelUnits()
	{
		return getDistanceUnit() == DistanceUnit.PIXEL;
	}

	/**
	 * Convert the distance units to pixels. Requires the calibration to have distance units and nm/pixel.
	 *
	 * @param helper
	 *            the helper
	 * @return true, if the distance units are now in pixels
	 */
	private boolean convertDistanceToPixelUnits(CalibrationHelper helper)
	{
		if (isDistanceInPixelUnits())
			return true;

		if (calibration.hasNmPerPixel())
		{
			try
			{
				TypeConverter<DistanceUnit> c = UnitConverterFactory.createConverter(getDistanceUnit(),
						DistanceUnit.PIXEL, calibration.getNmPerPixel());
				// Convert data
				for (int i = 0, size = size(); i < size; i++)
				{
					PeakResult p = get(i);
					// Leave the original positions
					//p.origX
					//p.origY
					convertDistance(p.params, c);
					if (p.paramsStdDev != null)
						convertDistance(p.paramsStdDev, c);
				}
				calibration.setDistanceUnit(DistanceUnit.PIXEL);
				return true;
			}
			catch (ConversionException e)
			{
				// Gracefully fail so ignore this
			}
		}
		return false;
	}

	/** The Constant offsetYSD. */
	private final static int offsetY, offsetXSD, offsetYSD;
	static
	{
		offsetY = Gaussian2DFunction.Y_POSITION - Gaussian2DFunction.X_POSITION;
		offsetXSD = Gaussian2DFunction.X_SD - Gaussian2DFunction.X_POSITION;
		offsetYSD = Gaussian2DFunction.Y_SD - Gaussian2DFunction.X_POSITION;
	}

	/**
	 * Convert distance.
	 *
	 * @param params
	 *            the params
	 * @param c
	 *            the c
	 */
	private void convertDistance(float[] params, TypeConverter<DistanceUnit> c)
	{
		for (int i = Gaussian2DFunction.X_POSITION; i < params.length; i += 6)
		{
			params[i] = c.convert(params[i]);
			params[i + offsetY] = c.convert(params[i + offsetY]);
			params[i + offsetXSD] = c.convert(params[i + offsetXSD]);
			params[i + offsetYSD] = c.convert(params[i + offsetYSD]);
		}
	}

	/**
	 * Checks if intensity is in photons units.
	 *
	 * @return true, if intensity is in photons
	 */
	public boolean isIntensityInPhotonUnits()
	{
		return (getIntensityUnit() == IntensityUnit.PHOTON);
	}

	/**
	 * Convert the intensity units to photons. Requires the calibration to have intensity units, gain and bias.
	 *
	 * @param helper
	 *            the helper
	 * @return true, if the intensity units are now in photons
	 */
	private boolean convertIntensityToPhotonUnits(CalibrationHelper helper)
	{
		if (isIntensityInPhotonUnits())
			return true;

		if (calibration.hasGain() && calibration.hasBias())
		{
			try
			{
				TypeConverter<IntensityUnit> bc = UnitConverterFactory.createConverter(getIntensityUnit(),
						IntensityUnit.PHOTON, calibration.getBias(), calibration.getGain());
				TypeConverter<IntensityUnit> c = UnitConverterFactory.createConverter(getIntensityUnit(),
						IntensityUnit.PHOTON, calibration.getGain());
				// Convert data
				for (int i = 0, size = size(); i < size; i++)
				{
					PeakResult p = get(i);
					// Leave the original value
					//p.origValue
					p.noise = (float) c.convert(p.noise);
					// Background must account for the bias
					p.params[Gaussian2DFunction.BACKGROUND] = (float) bc
							.convert(p.params[Gaussian2DFunction.BACKGROUND]);
					convertIntensity(p.params, c);
					if (p.paramsStdDev != null)
					{
						// Standard deviations so do not subtract the bias from the background
						p.paramsStdDev[Gaussian2DFunction.BACKGROUND] = (float) c
								.convert(p.paramsStdDev[Gaussian2DFunction.BACKGROUND]);
						convertIntensity(p.paramsStdDev, c);
					}
				}
				calibration.setIntensityUnit(IntensityUnit.PHOTON);
				return true;
			}
			catch (ConversionException e)
			{
				// Gracefully fail so ignore this
			}
		}
		return false;
	}

	/**
	 * Convert intensity.
	 *
	 * @param params
	 *            the params
	 * @param c
	 *            the c
	 */
	private void convertIntensity(float[] params, TypeConverter<IntensityUnit> c)
	{
		for (int i = Gaussian2DFunction.SIGNAL; i < params.length; i += 6)
		{
			params[Gaussian2DFunction.SIGNAL] = c.convert(params[Gaussian2DFunction.SIGNAL]);
		}
	}

	/**
	 * Checks if angle is in radian units.
	 *
	 * @return true, if angle is in radians
	 */
	public boolean isAngleInRadianUnits()
	{
		return getAngleUnit() == AngleUnit.RADIAN;
	}

	/**
	 * Convert the angle units to radians. Requires the calibration to have angle units.
	 *
	 * @param helper
	 *            the helper
	 * @return true, if the angle units are now in radians
	 */
	private boolean convertAngleToRadianUnits(CalibrationHelper helper)
	{
		if (isAngleInRadianUnits())
			return true;

		try
		{
			TypeConverter<AngleUnit> c = UnitConverterFactory.createConverter(getAngleUnit(), AngleUnit.RADIAN);
			// Convert data
			for (int i = 0, size = size(); i < size; i++)
			{
				PeakResult p = get(i);
				convertAngle(p.params, c);
				if (p.paramsStdDev != null)
					convertAngle(p.paramsStdDev, c);
			}
			calibration.setAngleUnit(AngleUnit.RADIAN);
			return true;
		}
		catch (ConversionException e)
		{
			// Gracefully fail so ignore this
		}
		return false;
	}

	/**
	 * Convert angle.
	 *
	 * @param params
	 *            the params
	 * @param c
	 *            the c
	 */
	private void convertAngle(float[] params, TypeConverter<AngleUnit> c)
	{
		for (int i = Gaussian2DFunction.SHAPE; i < params.length; i += 6)
		{
			params[i] = c.convert(params[i]);
		}
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
			procedure.execute(get(i));
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
		procedure.execute(get(0));
	}

	/**
	 * For each result execute the procedure.
	 *
	 * @param procedure
	 *            the procedure
	 */
	public void forEach(PeakResultProcedureX procedure)
	{
		for (int i = 0, size = size(); i < size; i++)
		{
			if (procedure.execute(get(i)))
				return;
		}
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
			final PeakResult r = get(i);
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
		final PeakResult r = get(0);
		procedure.executeBIXYZ(r.getBackground(), r.getSignal(), r.getXPosition(), r.getYPosition(), r.getZPosition());
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
		if (calibration == null)
			throw new ConfigurationException("No calibration");

		ArrayList<TypeConverter<IntensityUnit>> list = calibration.getDualIntensityConverter(intensityUnit);
		TypeConverter<IntensityUnit> ic = list.get(0);
		TypeConverter<IntensityUnit> bic = list.get(1);
		TypeConverter<DistanceUnit> dc = calibration.getDistanceConverter(distanceUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = get(i);
			//@formatter:off
			procedure.executeBIXY(
					bic.convert(r.getBackground()), 
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
		if (calibration == null)
			throw new ConfigurationException("No calibration");

		ArrayList<TypeConverter<IntensityUnit>> list = calibration.getDualIntensityConverter(intensityUnit);
		TypeConverter<IntensityUnit> ic = list.get(0);
		TypeConverter<IntensityUnit> bic = list.get(1);
		TypeConverter<DistanceUnit> dc = calibration.getDistanceConverter(distanceUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = get(i);
			//@formatter:off
			procedure.executeBIXYZ(
					bic.convert(r.getBackground()), 
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
		if (calibration == null)
			throw new ConfigurationException("No calibration");

		int[] indices = PSFHelper.getGaussian2DWxWyIndices(psf);

		final int ix = indices[0] + PeakResult.STANDARD_PARAMETERS;
		final int iy = indices[1] + PeakResult.STANDARD_PARAMETERS;
		final double twoPi = 2 * Math.PI;

		TypeConverter<IntensityUnit> ic = calibration.getIntensityConverter(intensityUnit);
		TypeConverter<DistanceUnit> dc = calibration.getDistanceConverter(DistanceUnit.PIXEL);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = get(i);

			// Convert the widths to pixels
			float sx = dc.convert(r.getParameter(ix));
			float sy = dc.convert(r.getParameter(iy));

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
		if (calibration == null)
			throw new ConfigurationException("No calibration");

		TypeConverter<IntensityUnit> ic = calibration.getIntensityConverter(intensityUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = get(i);
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
		if (calibration == null)
			throw new ConfigurationException("No calibration");

		TypeConverter<IntensityUnit> ic = calibration.getIntensityConverter(intensityUnit);
		TypeConverter<DistanceUnit> dc = calibration.getDistanceConverter(distanceUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = get(i);
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
	public void forEach(IntensityUnit intensityUnit, DistanceUnit distanceUnit, IXYZResultProcedure procedure)
			throws ConversionException, ConfigurationException
	{
		if (calibration == null)
			throw new ConfigurationException("No calibration");

		TypeConverter<IntensityUnit> ic = calibration.getIntensityConverter(intensityUnit);
		TypeConverter<DistanceUnit> dc = calibration.getDistanceConverter(distanceUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = get(i);
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
		if (calibration == null)
			throw new ConfigurationException("No calibration");

		TypeConverter<DistanceUnit> dc = calibration.getDistanceConverter(distanceUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = get(i);
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
		if (calibration == null)
			throw new ConfigurationException("No calibration");

		// Note that in the future we may support more than just Gaussian2D PSF
		// so this may have to change
		
		int[] indices = PSFHelper.getGaussian2DWxWyIndices(psf);

		final int ix = indices[0] + PeakResult.STANDARD_PARAMETERS;
		final int iy = indices[1] + PeakResult.STANDARD_PARAMETERS;

		TypeConverter<DistanceUnit> dc = calibration.getDistanceConverter(distanceUnit);

		if (ix == iy)
		{
			for (int i = 0, size = size(); i < size; i++)
			{
				final PeakResult r = get(i);
				//@formatter:off
    			procedure.executeW(
    					dc.convert(r.getParameter(ix)));
    			//@formatter:on
			}
		}
		else
		{
			for (int i = 0, size = size(); i < size; i++)
			{
				final PeakResult r = get(i);
				// Convert the separate widths into a single width
				double s = PeakResultHelper.getGaussian2DStandardDeviation(r.getParameter(ix), r.getParameter(iy));
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
		if (calibration == null)
			throw new ConfigurationException("No calibration");

		// Note that in the future we may support more than just Gaussian2D PSF
		// so this may have to change
		
		int[] indices = PSFHelper.getGaussian2DWxWyIndices(psf);

		final int ix = indices[0] + PeakResult.STANDARD_PARAMETERS;
		final int iy = indices[1] + PeakResult.STANDARD_PARAMETERS;

		TypeConverter<DistanceUnit> dc = calibration.getDistanceConverter(distanceUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = get(i);
			//@formatter:off
			procedure.executeWxWy(
					dc.convert(r.getParameter(ix)),
					dc.convert(r.getParameter(iy)));
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
		if (calibration == null)
			throw new ConfigurationException("No calibration");

		TypeConverter<DistanceUnit> dc = calibration.getDistanceConverter(distanceUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = get(i);
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
		if (calibration == null)
			throw new ConfigurationException("No calibration");

		TypeConverter<DistanceUnit> dc = calibration.getDistanceConverter(distanceUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = get(i);
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
		if (calibration == null)
			throw new ConfigurationException("No calibration");

		TypeConverter<DistanceUnit> dc = calibration.getDistanceConverter(distanceUnit);

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = get(i);
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
		if (calibration == null)
			throw new ConfigurationException("No calibration");

		// TODO - do the conversion using the calibration helper.

		// TODO - Check if this is a Gaussian2DFunction and throw an error if not
		// Otherwise determine the PSF fields to obtain the distance
		if (!isCCDCamera())
			throw new ConfigurationException("Not a CCD camera");
		final boolean emCCD = isEMCCD();

		ArrayList<TypeConverter<IntensityUnit>> list = calibration.getDualIntensityConverter(IntensityUnit.PHOTON);
		TypeConverter<IntensityUnit> ic = list.get(0);
		TypeConverter<IntensityUnit> bic = list.get(1);
		TypeConverter<DistanceUnit> dc = calibration.getDistanceConverter(DistanceUnit.NM);

		// This will be fine if the intensity converter was created
		final double nmPerPixel = getNmPerPixel();

		for (int i = 0, size = size(); i < size; i++)
		{
			final PeakResult r = get(i);
			float s = r.getSD();
			procedure.executeLSEPrecision(PeakResult.getPrecision(nmPerPixel, dc.convert(s), ic.convert(r.getSignal()),
					bic.convert(r.getBackground()), emCCD));
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
		if (calibration != null)
		{
			if (calibration.hasDistanceUnit())
				return calibration.getDistanceUnit();
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
		if (calibration != null)
		{
			if (calibration.hasIntensityUnit())
				return calibration.getIntensityUnit();
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
		if (calibration != null)
		{
			if (calibration.hasAngleUnit())
				return calibration.getAngleUnit();
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
			final PeakResult r = get(i);
			if (r.params[0] == 0)
				r.params[0] = newBackground;
		}
	}
}
