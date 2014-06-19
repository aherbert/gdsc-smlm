package gdsc.smlm.model;

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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Contains a model for a blinking fluorophore.
 * <p>
 * Based on the work of Coltharp et al (2012) Accurate Construction of photoactivated localization microscopy images for
 * quantitative measurements. PLOS One 7, Issue 12, pp 1-15
 */
public abstract class FluorophoreSequenceModel extends MoleculeModel implements Comparable<FluorophoreSequenceModel>
{
	public FluorophoreSequenceModel(int id, double x, double y, double z)
	{
		super(id, x, y, z);
	}

	public FluorophoreSequenceModel(int id, double[] xyz)
	{
		super(id, xyz);
	}

	/**
	 * The number of times the molecule went into the dark state
	 */
	private int blinks = 0;
	/**
	 * A sequence of fluorescent bursts in pairs of {on,off} times. The burst sequence will be length = 2 * (blinks+1)
	 */
	private double[] burstSequence = new double[] { 0, 0 };

	protected void setBurstSequence(double[] sequence)
	{
		if (sequence != null && sequence.length > 1)
		{
			blinks = (sequence.length / 2) - 1;

			// Ensure the sequence array is an even number in length 
			final int length = 2 * (blinks + 1);
			if (sequence.length == length)
				burstSequence = sequence;
			else
				burstSequence = Arrays.copyOf(sequence, length);
		}
	}

	/**
	 * @return The number of times the fluorophore blinked
	 */
	public int getNumberOfBlinks()
	{
		return blinks;
	}

	/**
	 * Get the start time, i.e. when the molecule activated.
	 * <p>
	 * Note that a molecule will always have a start time even if it has no blinks. This models a molecule that turns on
	 * and then bleaches immediately.
	 * 
	 * @return The start time
	 */
	public double getStartTime()
	{
		return burstSequence[0];
	}

	/**
	 * Get the end time, i.e. when the molecule bleached.
	 * 
	 * @return The end time
	 */
	public double getEndTime()
	{
		return burstSequence[burstSequence.length - 1];
	}

	/**
	 * @return Fluorescent bursts arranged as list of on/off times: {onT,offT}
	 */
	public List<double[]> getBurstSequence()
	{
		ArrayList<double[]> data = new ArrayList<double[]>(blinks + 1);
		for (int i = 0; i <= blinks; i++)
		{
			data.add(new double[] { burstSequence[i * 2], burstSequence[i * 2 + 1] });
		}
		return data;
	}

	/**
	 * @return Fluorescent bursts arranged as list of on/off times in integer sampling intervals: {onT,offT}
	 */
	public List<int[]> getSampledBurstSequence()
	{
		ArrayList<int[]> data = new ArrayList<int[]>(blinks + 1);
		for (int i = 0; i <= blinks; i++)
		{
			data.add(new int[] { (int) (burstSequence[i * 2]), (int) (burstSequence[i * 2 + 1]) });
		}
		return data;
	}

	/**
	 * Order by time ascending
	 * 
	 * @param o
	 *            The other fluorophore
	 * @return -1,0,1
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(FluorophoreSequenceModel o)
	{
		final double s1 = getStartTime();
		final double s2 = o.getStartTime();
		return (s1 < s2) ? -1 : (s1 == s2) ? 0 : 1;
	}

	/**
	 * @return The duration of the on times
	 */
	public double[] getOnTimes()
	{
		double[] onTimes = new double[blinks + 1];
		for (int i = 0; i <= blinks; i++)
		{
			onTimes[i] = burstSequence[i * 2 + 1] - burstSequence[i * 2];
		}
		return onTimes;
	}

	/**
	 * @return The duration of the off times
	 */
	public double[] getOffTimes()
	{
		if (blinks < 1)
			return new double[0];

		double[] offTimes = new double[blinks];
		for (int i = 1; i <= blinks; i++)
		{
			offTimes[i - 1] = burstSequence[i * 2] - burstSequence[i * 2 - 1];
		}
		return offTimes;
	}

	/**
	 * @return The duration of the on times if sampled at integer time intervals
	 */
	public int[] getSampledOnTimes()
	{
		if (blinks == 0)
			return new int[] { end(burstSequence[1]) - start(burstSequence[0]) };

		// Process all blinks. Join together blinks with an off-time that would not be noticed, 
		// i.e. where the molecule was on in consecutive frames.
		int[] onTimes = new int[blinks + 1];
		int n = 0;
		int tStart = (int) burstSequence[0];
		for (int i = 0; i < blinks; i++)
		{
			int end1 = end(burstSequence[i * 2 + 1]);
			int start2 = start(burstSequence[(i + 1) * 2]);

			if (start2 - end1 > 0)
			{
				onTimes[n++] = end1 - tStart;
				tStart = start2;
			}
		}
		onTimes[n++] = end(getEndTime()) - tStart;

		return Arrays.copyOf(onTimes, n);
	}

	private int start(double t)
	{
		return (int) t;
	}

	private int end(double t)
	{
		return (int) (Math.ceil(t));
	}

	/**
	 * @return The duration of the off times if sampled at integer time intervals
	 */
	public int[] getSampledOffTimes()
	{
		if (blinks == 0)
			return new int[0];

		// Process all blinks. Join together blinks with an off-time that would not be noticed, 
		// i.e. where the molecule was on in consecutive frames.
		int[] offTimes = new int[blinks];
		int n = 0;
		for (int i = 0; i < blinks; i++)
		{
			int end1 = end(burstSequence[i * 2 + 1]);
			int start2 = start(burstSequence[(i + 1) * 2]);

			if (start2 - end1 > 0)
			{
				offTimes[n++] = start2 - end1;
			}
		}

		return Arrays.copyOf(offTimes, n);
	}

	/**
	 * @return An array of frames when the molecule was on
	 */
	public int[] getOnFrames()
	{
		int sequenceStartT = (int) getStartTime();
		int sequenceEndT = (int) getEndTime();

		int n = 0;
		int[] onFrames = new int[sequenceEndT - sequenceStartT + 1];
		for (int i = 0; i <= blinks; i++)
		{
			int on = (int) (burstSequence[i * 2]);
			int off = (int) (burstSequence[i * 2 + 1]);

			for (int t = on; t <= off; t++)
			{
				onFrames[n++] = t;
			}
		}

		return Arrays.copyOf(onFrames, n);
	}

	/**
	 * Scale the times using the specified factor. Allows adjusting the relative time of the sequence.
	 * 
	 * @param scale
	 */
	public void adjustTime(double scale)
	{
		if (scale < 0)
			throw new IllegalArgumentException("Scale factor must be above zero");
		for (int i = 0; i < burstSequence.length; i++)
			burstSequence[i] *= scale;
	}
}
