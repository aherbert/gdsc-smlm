package gdsc.smlm.results;

import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.ListIterator;
import java.util.Set;

import gdsc.smlm.data.config.SMLMSettings.AngleUnit;
import gdsc.smlm.data.config.SMLMSettings.DistanceUnit;
import gdsc.smlm.data.config.SMLMSettings.IntensityUnit;
import gdsc.core.data.utils.ConversionException;
import gdsc.core.data.utils.TypeConverter;
import gdsc.smlm.data.config.UnitConverterFactory;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;

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
 * Stores peak results in memory. The PeakResults interface add methods are thread safe as they are synchronized. There
 * are equivalent non-.synchronized methods.
 */
public class MemoryPeakResults extends AbstractPeakResults implements Cloneable, Iterable<PeakResult>
{
	private static LinkedHashMap<String, MemoryPeakResults> resultsMap = new LinkedHashMap<String, MemoryPeakResults>();
	private static final Runtime s_runtime = Runtime.getRuntime();
	private static int byteSize = 0;
	private static int byteSizeWithDeviations = 0;

	private static final int DEFAULT_SIZE = 96;
	private static final int DEFAULT_SIZE_WITH_DEVIATIONS = 144;

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
	 * Count the number of result sets in memory
	 */
	public static int getResultsMemorySize()
	{
		return resultsMap.size();
	}

	/**
	 * Return true if there are no non-empty results in memory
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

	/**
	 * Add a result. Synchronized.
	 * 
	 * {@inheritDoc}
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#add(int, int, int, float, double, float, float[], float[])
	 */
	public synchronized void add(int peak, int origX, int origY, float origValue, double chiSquared, float noise,
			float[] params, float[] paramsStdDev)
	{
		results.add(new PeakResult(peak, origX, origY, origValue, chiSquared, noise, params, paramsStdDev));
	}

	/**
	 * Add a result. Not synchronized.
	 *
	 * @param peak
	 *            the peak (e.g. frame number)
	 * @param origX
	 *            the original X position
	 * @param origY
	 *            the original Y position
	 * @param origValue
	 *            the original value
	 * @param error
	 *            the error
	 * @param noise
	 *            the noise
	 * @param params
	 *            the parameters
	 * @param paramsStdDev
	 *            the parameters standard deviation (or null)
	 */
	public void addf(int peak, int origX, int origY, float origValue, double error, float noise, float[] params,
			float[] paramsStdDev)
	{
		results.add(new PeakResult(peak, origX, origY, origValue, error, noise, params, paramsStdDev));
	}

	/**
	 * Add all results. Synchronized.
	 * 
	 * {@inheritDoc}
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#addAll(java.util.Collection)
	 */
	public synchronized void addAll(Collection<PeakResult> results)
	{
		this.results.addAll(results);
	}

	/**
	 * Add all results. Not synchronized.
	 * 
	 */
	public void addAllf(Collection<PeakResult> results)
	{
		this.results.addAll(results);
	}

	/**
	 * Add a result. Not synchronized.
	 *
	 * @param result
	 *            the result
	 */
	@Override
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
	 * Gets the bounds.
	 *
	 * @param calculate
	 *            Set to true to calculate the bounds if they are null or zero width/height
	 * @return the bounds of the result coordinates
	 */
	public Rectangle getBounds(boolean calculate)
	{
		if ((bounds == null || bounds.width == 0 || bounds.height == 0) && calculate)
		{
			bounds = new Rectangle();
			Rectangle2D.Float b = getDataBounds();

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
	 * @return the bounds of the result coordinates
	 */
	public Rectangle2D.Float getDataBounds()
	{
		if (isEmpty())
			return new Rectangle2D.Float();

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
		return new Rectangle2D.Float(minX, minY, maxX - minX, maxY - minY);
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

	/**
	 * Run the garbage collector multiple times to free memory
	 */
	public static void runGC()
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
			runGCOnce();
			Thread.currentThread();
			Thread.yield();

			usedMem2 = usedMem1;
			usedMem1 = usedMemory();
		}
	}

	public static void runGCOnce()
	{
		s_runtime.runFinalization();
		s_runtime.gc();
	}

	public static long usedMemory()
	{
		return s_runtime.totalMemory() - s_runtime.freeMemory();
	}

	public static long totalMemory()
	{
		return s_runtime.totalMemory();
	}

	public static long freeMemory()
	{
		return s_runtime.freeMemory();
	}

	/**
	 * Shallow copy this set of results. To create new object references use {@link #copy()}.
	 * 
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
			if (results != null)
				copy.results = new ArrayList<PeakResult>(results);
		}
		return copy;
	}

	/**
	 * @return True if empty
	 */
	public boolean isEmpty()
	{
		return results.isEmpty();
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
	 * Gets the head position in the set of results.
	 *
	 * @return the head
	 */
	public PeakResult getHead()
	{
		if (isEmpty())
			return null;
		return results.get(0);
	}

	/**
	 * Gets the tail position in the set of results.
	 *
	 * @return the tail
	 */
	public PeakResult getTail()
	{
		if (isEmpty())
			return null;
		return results.get(size() - 1);
	}

	/**
	 * Checks if all results have a stored precision value.
	 *
	 * @return true, if all results have a stored precision value
	 */
	public boolean hasStoredPrecision()
	{
		for (int i = 0; i < results.size(); i++)
		{
			if (!results.get(i).hasPrecision())
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
		for (int i = 0; i < results.size(); i++)
		{
			if (results.get(i) == null)
				return true;
		}
		return false;
	}

	/**
	 * Removes the null results from the store.
	 */
	public void removeNullResults()
	{
		ArrayList<PeakResult> list = new ArrayList<PeakResult>(size());
		// Use an iterator to take advantage of the fail-fast thread safety
		for (PeakResult e : results)
		{
			if (e != null)
				list.add(e);
		}
		this.results = list;
	}

	/**
	 * Checks if distance is in pixel units.
	 *
	 * @return true, if distance is in pixels
	 */
	public boolean isDistanceInPixelUnits()
	{
		if (calibration != null && calibration.hasDistanceUnit())
		{
			if (calibration.getDistanceUnit() == DistanceUnit.PIXEL)
				return true;
		}
		return false;
	}

	/**
	 * Convert the distance units to pixels. Requires the calibration to have distance units and nm/pixel.
	 *
	 * @return true, if the distance units are now in pixels
	 */
	public boolean convertDistanceToPixelUnits()
	{
		if (calibration != null && calibration.hasDistanceUnit())
		{
			if (calibration.getDistanceUnit() == DistanceUnit.PIXEL)
				return true;

			if (calibration.hasNmPerPixel())
			{
				try
				{
					TypeConverter<DistanceUnit> c = UnitConverterFactory.createConverter(calibration.getDistanceUnit(),
							DistanceUnit.PIXEL, calibration.getNmPerPixel());
					// Convert data
					for (PeakResult p : results)
					{
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
		}
		return false;
	}

	private final static int offsetY, offsetXSD, offsetYSD;
	static
	{
		offsetY = Gaussian2DFunction.Y_POSITION - Gaussian2DFunction.X_POSITION;
		offsetXSD = Gaussian2DFunction.X_SD - Gaussian2DFunction.X_POSITION;
		offsetYSD = Gaussian2DFunction.Y_SD - Gaussian2DFunction.X_POSITION;
	}

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
		if (calibration != null && calibration.hasIntensityUnit())
		{
			if (calibration.getIntensityUnit() == IntensityUnit.PHOTON)
				return true;
		}
		return false;
	}

	/**
	 * Convert the intensity units to photons. Requires the calibration to have intensity units, gain and bias.
	 *
	 * @return true, if the intensity units are now in photons
	 */
	public boolean convertIntensityToPhotonUnits()
	{
		if (calibration != null && calibration.hasIntensityUnit())
		{
			if (calibration.getIntensityUnit() == IntensityUnit.PHOTON)
				return true;

			if (calibration.hasGain() && calibration.hasBias())
			{
				try
				{
					TypeConverter<IntensityUnit> bc = UnitConverterFactory.createConverter(
							calibration.getIntensityUnit(), IntensityUnit.PHOTON, calibration.getBias(),
							calibration.getGain());
					TypeConverter<IntensityUnit> c = UnitConverterFactory.createConverter(
							calibration.getIntensityUnit(), IntensityUnit.PHOTON, calibration.getGain());
					// Convert data
					for (PeakResult p : results)
					{
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
		}
		return false;
	}

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
		if (calibration != null && calibration.hasAngleUnit())
		{
			if (calibration.getAngleUnit() == AngleUnit.RADIAN)
				return true;
		}
		return false;
	}

	/**
	 * Convert the angle units to radians. Requires the calibration to have angle units.
	 *
	 * @return true, if the angle units are now in radians
	 */
	public boolean convertAngleToRadianUnits()
	{
		if (calibration != null && calibration.hasAngleUnit())
		{
			if (calibration.getAngleUnit() == AngleUnit.RADIAN)
				return true;

			try
			{
				TypeConverter<AngleUnit> c = UnitConverterFactory.createConverter(calibration.getAngleUnit(),
						AngleUnit.RADIAN);
				// Convert data
				for (PeakResult p : results)
				{
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
		}
		return false;
	}

	private void convertAngle(float[] params, TypeConverter<AngleUnit> c)
	{
		for (int i = Gaussian2DFunction.SHAPE; i < params.length; i += 6)
		{
			params[i] = c.convert(params[i]);
		}
	}

	/**
	 * Convert to preference units.
	 *
	 * @return true, if successful
	 */
	public boolean convertToPreferenceUnits()
	{
		// TODO - do the conversion using the calibration helper.
		// Remove the conversion methods from this class.
		// CalibrationHelper helper = new CalibrationHelper(null);

		boolean success = convertDistanceToPixelUnits();
		success &= convertIntensityToPhotonUnits();
		success &= convertAngleToRadianUnits();
		return success;
	}

	/**
	 * For each result execute the procedure.
	 *
	 * @param procedure
	 *            the procedure
	 */
	public void forEach(ResultProcedure procedure)
	{
		for (int i = 0, size = results.size(); i < size; i++)
		{
			final PeakResult r = results.get(i);
			procedure.execute(r.getBackground(), r.getSignal(), r.getXPosition(), r.getYPosition(),
					r.getXPosition());
		}
	}

	/**
	 * For each result execute the procedure using the specified units.
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 *
	 * @param procedure
	 *            the procedure
	 * @param intensityUnit
	 *            the intensity unit
	 * @param distanceUnit
	 *            the distance unit
	 * @throws ConversionException
	 *             if the conversion is not possible
	 */
	public void forEach(ResultProcedure procedure, IntensityUnit intensityUnit, DistanceUnit distanceUnit)
	{
		if (calibration == null)
			throw new ConversionException("No calibration");

		ArrayList<TypeConverter<IntensityUnit>> list = calibration.getIntensityConverter(intensityUnit);
		TypeConverter<IntensityUnit> ic = list.get(0);
		TypeConverter<IntensityUnit> bic = list.get(1);
		TypeConverter<DistanceUnit> dc = calibration.getDistanceConverter(distanceUnit);

		for (int i = 0, size = results.size(); i < size; i++)
		{
			final PeakResult r = results.get(i);
			//@formatter:off
			procedure.execute(
					bic.convert(r.getBackground()), 
					ic.convert(r.getSignal()), 
					dc.convert(r.getXPosition()),
					dc.convert(r.getYPosition()), 
					dc.convert(r.getZPosition()));
			//@formatter:on
		}
	}

	/**
	 * For each result execute the procedure.
	 * <p>
	 * This will fail if the calibration is missing information to convert the units.
	 *
	 * @param procedure
	 *            the procedure
	 * @throws ConversionException
	 *             if the conversion is not possible
	 */
	public void forEach(LSEPrecisionProcedure procedure)
	{
		if (calibration == null)
			throw new ConversionException("No calibration");

		// TODO - do the conversion using the calibration helper.
		
		// TODO - Check if this is a Gaussian2DFunction and throw an error if not
		// Otherwise determine the PSF fields to obtain the distance
		if (!isCCDCamera())
			throw new ConversionException("Not a CCD camera");
		final boolean emCCD = isEMCCD();

		ArrayList<TypeConverter<IntensityUnit>> list = calibration.getIntensityConverter(IntensityUnit.PHOTON);
		TypeConverter<IntensityUnit> ic = list.get(0);
		TypeConverter<IntensityUnit> bic = list.get(1);
		TypeConverter<DistanceUnit> dc = calibration.getDistanceConverter(DistanceUnit.NM);

		// This will be fine if the intensity converter was created
		final double nmPerPixel = getNmPerPixel();
		
		for (int i = 0, size = results.size(); i < size; i++)
		{
			final PeakResult r = results.get(i);
			float s = r.getSD();
			procedure.execute(PeakResult.getPrecision(nmPerPixel, dc.convert(s), ic.convert(r.getSignal()), bic.convert(r.getBackground()), emCCD));
		}
	}
}
