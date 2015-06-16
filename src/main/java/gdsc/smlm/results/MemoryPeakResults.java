package gdsc.smlm.results;

import gdsc.smlm.function.gaussian.Gaussian2DFunction;

import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.ListIterator;
import java.util.Set;

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
 * Stores peak results in memory. The PeakResults interface add methods are thread safe.
 */
public class MemoryPeakResults extends AbstractPeakResults implements Iterable<PeakResult>
{
	private static LinkedHashMap<String, MemoryPeakResults> resultsMap = new LinkedHashMap<String, MemoryPeakResults>();
	private static final Runtime s_runtime = Runtime.getRuntime();
	private static int byteSize = 0;
	private static int byteSizeWithDeviations = 0;

	private static final int DEFAULT_SIZE = 96;
	private static final int DEFAULT_SIZE_WITH_DEVIATIONS = 144;

	public static final double DEFAULT_NM_PER_PIXEL = 100;
	public static final double DEFAULT_GAIN = 1;
	public static final boolean DEFAULT_EMCCD = true;

	private ArrayList<PeakResult> results;
	private boolean sortAfterEnd;

	public MemoryPeakResults()
	{
		results = new ArrayList<PeakResult>(1000);
	}

	public MemoryPeakResults(int capacity)
	{
		results = new ArrayList<PeakResult>(capacity);
	}

	/**
	 * @param name
	 *            The name of the results
	 * @return Get the named results (or null if they do not exist)
	 */
	public static MemoryPeakResults getResults(String name)
	{
		return resultsMap.get(name);
	}

	/**
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
	 */
	public static void addResults(MemoryPeakResults results)
	{
		if (results == null)
			throw new NullPointerException("Results must not be null");
		results.trimToSize();
		resultsMap.put(results.getName(), results);
	}

	/**
	 * @return A set of the available named results held in memory
	 */
	public static Set<String> getResultNames()
	{
		return resultsMap.keySet();
	}

	/**
	 * @return A collection of the results held in memory
	 */
	public static Collection<MemoryPeakResults> getAllResults()
	{
		return resultsMap.values();
	}

	/**
	 * Count the total number of results in memory
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
	 * Estimate the total size of results in memory
	 */
	public static long estimateMemorySize()
	{
		long memorySize = 0;
		for (MemoryPeakResults r : resultsMap.values())
		{
			memorySize += estimateMemorySize(r.getResults());
		}
		return memorySize;
	}

	/**
	 * Clear the results from memory
	 */
	public static void clearMemory()
	{
		resultsMap.clear();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#begin()
	 */
	public void begin()
	{
		// Q. Should a new array be allocated if the previous one was very large?
		if (results.size() > 10000)
			results = new ArrayList<PeakResult>(1000);
		else
			results.clear();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#add(int, int, int, float, double, float, float[], float[])
	 */
	public synchronized void add(int peak, int origX, int origY, float origValue, double chiSquared, float noise,
			float[] params, float[] paramsStdDev)
	{
		results.add(new PeakResult(peak, origX, origY, origValue, chiSquared, noise, params, paramsStdDev));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#addAll(java.util.Collection)
	 */
	public synchronized void addAll(Collection<PeakResult> results)
	{
		this.results.addAll(results);
	}

	/**
	 * Add a result. Not synchronized.
	 * 
	 * @param result
	 */
	public void add(PeakResult result)
	{
		results.add(result);
	}

	/**
	 * Add a result. Synchronized.
	 * 
	 * @param result
	 */
	public synchronized void addSync(PeakResult result)
	{
		results.add(result);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#size()
	 */
	public int size()
	{
		results.trimToSize();
		return results.size();
	}

	/**
	 * @see {@link java.util.ArrayList#trimToSize() }
	 */
	public void trimToSize()
	{
		results.trimToSize();
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

	/**
	 * Sort the results
	 */
	public void sort()
	{
		Collections.sort(results);
	}

	/**
	 * @return The peak results
	 */
	public List<PeakResult> getResults()
	{
		return results;
	}

	/**
	 * @param calculate
	 *            Set to true to calculate the bounds if they are null or zero width/height
	 * @return the bounds of the result coordinates
	 */
	public Rectangle getBounds(boolean calculate)
	{
		if ((bounds == null || bounds.width == 0 || bounds.height == 0) && calculate)
		{
			bounds = new Rectangle();
			float minX = Float.POSITIVE_INFINITY, minY = Float.POSITIVE_INFINITY;
			float maxX = Float.NEGATIVE_INFINITY, maxY = Float.NEGATIVE_INFINITY;
			for (PeakResult result : results)
			{
				if (minX > result.params[Gaussian2DFunction.X_POSITION])
					minX = result.params[Gaussian2DFunction.X_POSITION];
				if (maxX < result.params[Gaussian2DFunction.X_POSITION])
					maxX = result.params[Gaussian2DFunction.X_POSITION];
				if (minY > result.params[Gaussian2DFunction.Y_POSITION])
					minY = result.params[Gaussian2DFunction.Y_POSITION];
				if (maxY < result.params[Gaussian2DFunction.Y_POSITION])
					maxY = result.params[Gaussian2DFunction.Y_POSITION];
			}
			// Round to integer
			bounds.x = (int) Math.floor(minX);
			bounds.y = (int) Math.floor(minY);

			// For compatibility with drawing images add one to the limits if they are integers
			if (maxX == (int) maxX)
				maxX += 1;
			if (maxY == (int) maxY)
				maxY += 1;

			bounds.width = (int) Math.ceil(maxX) - bounds.x;
			bounds.height = (int) Math.ceil(maxY) - bounds.y;
		}
		return bounds;
	}

	/**
	 * Returns a list iterator over the elements in this list (in proper
	 * sequence).
	 * 
	 * <p>
	 * The returned list iterator is <a href="#fail-fast"><i>fail-fast</i></a>.
	 * 
	 * @see java.util.List#listIterator()
	 */
	public ListIterator<PeakResult> listIterator()
	{
		return results.listIterator();
	}

	/**
	 * Returns an iterator over the elements in this list in proper sequence.
	 * 
	 * <p>
	 * The returned iterator is <a href="#fail-fast"><i>fail-fast</i></a>.
	 * 
	 * @see java.util.List#iterator()
	 * @return an iterator over the elements in this list in proper sequence
	 */
	public Iterator<PeakResult> iterator()
	{
		return results.iterator();
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
	 * @param sortAfterEnd
	 *            True if the results should be sorted after the {@link #end()} method
	 */
	public void setSortAfterEnd(boolean sortAfterEnd)
	{
		this.sortAfterEnd = sortAfterEnd;
	}

	/**
	 * @return True if the results should be sorted after the {@link #end()} method
	 */
	public boolean isSortAfterEnd()
	{
		return sortAfterEnd;
	}

	/**
	 * Convert the size in bytes into a string
	 * 
	 * @param memorySize
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
	 * @param size
	 * @param includeDeviations
	 * @return The memory size
	 */
	public static long estimateMemorySize(List<PeakResult> results)
	{
		long memorySize = 0;
		if (results != null && results.size() > 0)
		{
			boolean includeDeviations = results.get(0).paramsStdDev != null;
			memorySize = MemoryPeakResults.estimateMemorySize(results.size(), includeDeviations);
		}
		return memorySize;
	}

	/**
	 * Return an estimate of the memory size taken by PeakResult objects.
	 * <p>
	 * Note: This is just a guess based on measured sizes for the objects in memory.
	 * 
	 * @param size
	 * @param includeDeviations
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

	private static void runGC()
	{
		// It helps to call Runtime.gc()
		// using several method calls:
		for (int r = 0; r < 4; ++r)
			_runGC();
	}

	private static void _runGC()
	{
		long usedMem1 = usedMemory(), usedMem2 = Long.MAX_VALUE;
		for (int i = 0; (usedMem1 < usedMem2) && (i < 500); ++i)
		{
			s_runtime.runFinalization();
			s_runtime.gc();
			Thread.currentThread();
			Thread.yield();

			usedMem2 = usedMem1;
			usedMem1 = usedMemory();
		}
	}

	private static long usedMemory()
	{
		return s_runtime.totalMemory() - s_runtime.freeMemory();
	}

	/**
	 * Get the nm-per-pixel from the calibration, or if not available, return the {@link #DEFAULT_NM_PER_PIXEL}
	 * 
	 * @return the nmPerPixel
	 */
	public double getNmPerPixel()
	{
		return (calibration != null) ? calibration.nmPerPixel : DEFAULT_NM_PER_PIXEL;
	}

	/**
	 * Get the gain from the calibration, or if not available, return the {@link #DEFAULT_GAIN}
	 * 
	 * @return the gain
	 */
	public double getGain()
	{
		return (calibration != null) ? calibration.gain : DEFAULT_GAIN;
	}

	/**
	 * Get the EMCCD flag from the calibration, or if not available, return the {@link #DEFAULT_EMCCD}
	 * 
	 * @return the gain
	 */
	public boolean isEMCCD()
	{
		return (calibration != null) ? calibration.emCCD : DEFAULT_EMCCD;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.AbstractPeakResults#setCalibration(gdsc.smlm.results.Calibration)
	 */
	@Override
	public void setCalibration(Calibration calibration)
	{
		super.setCalibration(calibration);
	}

	/**
	 * Copy the results. Create new objects for the properties (avoiding a shallow copy) but does not
	 * deep copy all of the peak results. Allows results to be resorted but not modified.
	 */
	public MemoryPeakResults copy()
	{
		MemoryPeakResults copy;
		try
		{
			copy = (MemoryPeakResults) super.clone();

			// Deep copy the objects
			if (bounds != null)
				copy.bounds = new Rectangle(bounds);
			if (calibration != null)
				copy.calibration = calibration.clone();
			if (results != null)
				copy.results = new ArrayList<PeakResult>(results);

			return copy;
		}
		catch (CloneNotSupportedException e)
		{
			// This should not happen so ignore
		}
		return null;
	}
}
