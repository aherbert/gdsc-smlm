package gdsc.smlm.results;

import java.util.LinkedList;

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
 * Define a cluster of localisations from different frames that represent a single molecule trace
 */
public class Trace extends Cluster
{
	private int nBlinks = -1;
	private int[] onTimes, offTimes;

	public Trace()
	{
		super();
	}
	
	public Trace(PeakResult peakResult)
	{
		super(peakResult);
	}

	@Override
	public void add(PeakResult result)
	{
		super.add(result);
		nBlinks = -1; // Invalidate the analysis
	}

	private void analyse()
	{
		if (nBlinks == -1)
		{
			if (results.isEmpty())
			{
				nBlinks = 0;
				onTimes = offTimes = null;
				return;
			}
			if (results.size() == 1)
			{
				nBlinks = 1;
				onTimes = new int[] { 1 };
				offTimes = null;
				return;
			}

			// Ensure in the correct time-order
			sort();
			LinkedList<Integer> on = new LinkedList<Integer>();
			LinkedList<Integer> off = new LinkedList<Integer>();

			nBlinks = 1;
			int t1 = results.get(0).peak;
			int onStart = t1;
			for (int i = 0; i < results.size() - 1; i++)
			{
				int t2 = results.get(i + 1).peak;
				int diff = t2 - t1;
				if (diff > 1)
				{
					off.add(diff - 1);
					on.add(t1 - onStart + 1);
					nBlinks++;
					onStart = t2;
				}
				t1 = t2;
			}
			on.add(t1 - onStart + 1);
			
			onTimes = toArray(on);
			offTimes = toArray(off);
		}
	}

	private int[] toArray(LinkedList<Integer> data)
	{
		if (data.isEmpty())
			return null;
		int[] array = new int[data.size()];
		int i=0;
		for (int value : data)
		{
			array[i++] = value;
		}
		return array;
	}

	/**
	 * @return The number of times the molecule blinked
	 */
	public int getNBlinks()
	{
		analyse();
		return nBlinks;
	}

	/**
	 * @return The average on time for the molecule
	 */
	public double getOnTime()
	{
		analyse();
		return getAverage(onTimes);
	}

	/**
	 * @return The average off time for the molecule
	 */
	public double getOffTime()
	{
		analyse();
		return getAverage(offTimes);
	}

	private double getAverage(int[] times)
	{
		if (times != null)
		{
			double av = 0;
			for (int t : times)
				av += t;
			return av / times.length;
		}
		return 0;
	}

	/**
	 * @return the on-times
	 */
	public int[] getOnTimes()
	{
		analyse();
		return onTimes;
	}

	/**
	 * @return the off-times
	 */
	public int[] getOffTimes()
	{
		analyse();
		return offTimes;
	}
}
