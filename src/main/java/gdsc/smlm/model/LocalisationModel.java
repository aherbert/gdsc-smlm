package gdsc.smlm.model;

import java.util.Arrays;

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
 * Contains a discrete-time model for the position and intensity of a localisation.
 */
public class LocalisationModel implements Comparable<LocalisationModel>
{
	private double[] xyz;
	private int id;
	private int time;
	private double intensity;
	private int state;
	private LocalisationModel previous, next;
	private double[] data = null;
	private int label;

	public static final int SINGLE = 0;
	public static final int PREVIOUS = 1;
	public static final int NEXT = 2;
	public static final int NEIGHBOUR = PREVIOUS | NEXT;
	public static final int CONTINUOUS = 4 | PREVIOUS | NEXT;

	/**
	 * Create a new localisation
	 * 
	 * @param id
	 * @param time
	 * @param x
	 * @param y
	 * @param z
	 * @param intensity
	 * @param state
	 */
	public LocalisationModel(int id, int time, double x, double y, double z, double intensity, int state)
	{
		init(id, time, new double[] { x, y, z }, intensity, state, null, null);
	}

	/**
	 * Create a new localisation
	 * 
	 * @param id
	 * @param time
	 * @param xyz
	 *            Coordinates (will be deep copied)
	 * @param intensity
	 * @param state
	 */
	public LocalisationModel(int id, int time, double[] xyz, double intensity, int state)
	{
		init(id, time, Arrays.copyOf(xyz, xyz.length), intensity, state, null, null);
	}

	/**
	 * Create a new localisation
	 * 
	 * @param id
	 * @param time
	 * @param x
	 * @param y
	 * @param z
	 * @param intensity
	 * @param state
	 * @param previous
	 *            The previous localisation in this pulse (should be continuous time points)
	 * @param next
	 *            The next localisation in this pulse (should be continuous time points)
	 */
	public LocalisationModel(int id, int time, double x, double y, double z, double intensity, int state,
			LocalisationModel previous, LocalisationModel next)
	{
		init(id, time, new double[] { x, y, z }, intensity, state, previous, next);
	}

	/**
	 * Create a new localisation
	 * 
	 * @param id
	 * @param time
	 * @param xyz
	 *            Coordinates (will be deep copied)
	 * @param intensity
	 * @param state
	 * @param previous
	 *            The previous localisation in this pulse (should be continuous time points)
	 * @param next
	 *            The next localisation in this pulse (should be continuous time points)
	 */
	public LocalisationModel(int id, int time, double[] xyz, double intensity, int state, LocalisationModel previous,
			LocalisationModel next)
	{
		init(id, time, Arrays.copyOf(xyz, xyz.length), intensity, state, previous, next);
	}

	private void init(int id, int time, double[] xyz, double intensity, int state, LocalisationModel previous,
			LocalisationModel next)
	{
		this.xyz = xyz;
		this.id = id;
		this.time = time;
		this.intensity = intensity;
		this.state = state;
		this.previous = previous;
		this.next = next;
	}

	public double getX()
	{
		return xyz[0];
	}

	public double getY()
	{
		return xyz[1];
	}

	public double getZ()
	{
		return xyz[2];
	}

	/**
	 * @return The coordinates (x,y,z)
	 */
	public double[] getCoordinates()
	{
		return xyz;
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

	/**
	 * @return the intensity
	 */
	public double getIntensity()
	{
		return intensity;
	}

	/**
	 * @param intensity
	 *            The new intensity
	 */
	public void setIntensity(double intensity)
	{
		this.intensity = intensity;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(LocalisationModel o)
	{
		if (time == o.time)
		{
			if (intensity > o.intensity)
				return -1;
			if (intensity < o.intensity)
				return 1;
			return 0;
		}
		return time - o.time;
	}

	/**
	 * @return True if this localisation is on in the previous time interval
	 */
	public boolean hasPrevious()
	{
		return (state & PREVIOUS) == PREVIOUS;
	}

	/**
	 * @return True if this localisation is on in the next time interval
	 */
	public boolean hasNext()
	{
		return (state & NEXT) == NEXT;
	}

	/**
	 * @return True if this localisation is on in the previous or next time interval
	 */
	public boolean hasNeighbour()
	{
		return (state & (NEIGHBOUR)) != 0;
	}

	/**
	 * @return True if this localisation is on for the entire duration of the time interval
	 */
	public boolean isContinuous()
	{
		return (state & CONTINUOUS) == CONTINUOUS;
	}

	/**
	 * @return the previous
	 */
	public LocalisationModel getPrevious()
	{
		return previous;
	}

	/**
	 * @param previous
	 *            the previous to set
	 */
	public void setPrevious(LocalisationModel previous)
	{
		this.previous = previous;
		if (previous != null)
			previous.next = this;
	}

	/**
	 * @return the next
	 */
	public LocalisationModel getNext()
	{
		return next;
	}

	/**
	 * @param next
	 *            the next to set
	 */
	public void setNext(LocalisationModel next)
	{
		this.next = next;
		if (next != null)
			next.previous = this;
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
	 * Gets the label.
	 *
	 * @return the label
	 */
	public int getLabel()
	{
		return label;
	}

	/**
	 * Sets the label. This can be used to identify subsets of molecules.
	 *
	 * @param label
	 *            the new label
	 */
	public void setLabel(int label)
	{
		this.label = label;
	}
}
