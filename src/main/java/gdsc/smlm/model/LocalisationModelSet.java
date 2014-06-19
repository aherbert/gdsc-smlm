package gdsc.smlm.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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
 * Contains a collection of localisations associated with a discrete-time. This can be used to model a moving
 * localisation
 * simulated on a more refined time-scale to a single localisation.
 */
public class LocalisationModelSet implements Comparable<LocalisationModelSet>
{
	private int id, time;
	private List<LocalisationModel> localisations = new ArrayList<LocalisationModel>();
	private double[] data = null;
	private LocalisationModelSet previous, next;

	/**
	 * Create a new localisation
	 * 
	 * @param id
	 * @param time
	 */
	public LocalisationModelSet(int id, int time)
	{
		this.id = id;
		this.time = time;
	}

	public void add(LocalisationModel l)
	{
		localisations.add(l);
	}

	public int size()
	{
		return localisations.size();
	}

	public List<LocalisationModel> getLocalisations()
	{
		return localisations;
	}

	/**
	 * @return The Id
	 */
	public int getId()
	{
		return id;
	}

	/**
	 * Allow the package to set the id
	 * 
	 * @param id
	 *            The Id
	 */
	void setId(int id)
	{
		this.id = id;
	}

	/**
	 * @return The time
	 */
	public int getTime()
	{
		return time;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(LocalisationModelSet o)
	{
		return time - o.time;
	}

	/**
	 * @return True if this localisation is on for the entire duration of the time interval
	 */
	public boolean isContinuous()
	{
		if (localisations.isEmpty())
			return false;

		// All localisations must be continuous and consecutive in time
		int[] t = new int[localisations.size()];
		int c = 0;
		for (LocalisationModel l : localisations)
		{
			t[c++] = l.getTime();
			if (!l.isContinuous())
				return false;
		}

		// Check consecutive in time
		Arrays.sort(t);

		if ((t[t.length - 1] - t[0] + 1) > c)
			return false;

		//System.out.printf("continuous size = %d\n",  localisations.size());
		return true;
	}

	/**
	 * @return the data
	 */
	public double[] getData()
	{
		return data;
	}

	/**
	 * @param data
	 *            the data to set
	 */
	public void setData(double[] data)
	{
		this.data = data;
	}

	/**
	 * Convert the set of localisations to a single localisation with the combined signal and the centroid location
	 * (centre-of-mass weighted by intensity).
	 * 
	 * @return
	 */
	public LocalisationModel toLocalisation()
	{
		double intensity = 0;
		double[] xyz = new double[3];
		for (LocalisationModel l : localisations)
		{
			final double s = l.getIntensity();
			intensity += s;
			final double[] xyz2 = l.getCoordinates();
			for (int i = 0; i < 3; i++)
				xyz[i] += xyz2[i] * s;
		}
		if (!localisations.isEmpty())
			for (int i = 0; i < 3; i++)
				xyz[i] /= intensity;

		LocalisationModel l = new LocalisationModel(id, time, xyz, intensity,
				isContinuous() ? LocalisationModel.CONTINUOUS : LocalisationModel.SINGLE);
		l.setData(data);
		return l;
	}

	/**
	 * @return The total intensity
	 */
	public double getIntensity()
	{
		double intensity = 0;
		for (LocalisationModel l : localisations)
			intensity += l.getIntensity();
		return intensity;
	}

	/**
	 * @return the previous
	 */
	public LocalisationModelSet getPrevious()
	{
		return previous;
	}

	/**
	 * @param previous
	 *            the previous to set
	 */
	public void setPrevious(LocalisationModelSet previous)
	{
		this.previous = previous;
		if (previous != null)
			previous.next = this;
	}

	/**
	 * @return the next
	 */
	public LocalisationModelSet getNext()
	{
		return next;
	}

	/**
	 * @param next
	 *            the next to set
	 */
	public void setNext(LocalisationModelSet next)
	{
		this.next = next;
		if (next != null)
			next.previous = this;
	}

	/**
	 * @return True if either of the previous/next pointers are not null
	 */
	public boolean hasNeighbour()
	{
		return next != null || previous != null;
	}
}
