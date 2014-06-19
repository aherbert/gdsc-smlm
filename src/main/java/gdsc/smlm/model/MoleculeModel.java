package gdsc.smlm.model;

import org.apache.commons.math3.random.RandomGenerator;

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
 * Contains a model for a moving molecule.
 */
public class MoleculeModel 
{
	private int id;
	protected double[] xyz;
	protected double mass;

	/**
	 * Create a new molecule
	 * 
	 * @param id
	 * @param xyz
	 *            [x,y,z]
	 */
	public MoleculeModel(int id, double[] xyz)
	{
		this.id = id;
		this.xyz = xyz;
	}

	/**
	 * Create a new molecule
	 * 
	 * @param id
	 * @param x
	 * @param y
	 * @param z
	 */
	public MoleculeModel(int id, double x, double y, double z)
	{
		this.id = id;
		xyz = new double[] { x, y, z };
	}

	/**
	 * Create a new molecule
	 * 
	 * @param mass
	 * @param xyz
	 *            [x,y,z]
	 */
	public MoleculeModel(double mass, double[] xyz)
	{
		this.mass = mass;
		this.xyz = xyz;
	}

	/**
	 * Create a new molecule
	 * 
	 * @param mass
	 * @param x
	 * @param y
	 * @param z
	 */
	public MoleculeModel(double mass, double x, double y, double z)
	{
		this.mass = mass;
		xyz = new double[] { x, y, z };
	}

	/**
	 * Create a new molecule
	 * 
	 * @param id
	 * @param mass
	 * @param xyz
	 *            [x,y,z]
	 */
	public MoleculeModel(int id, double mass, double[] xyz)
	{
		this.id = id;
		this.mass = mass;
		this.xyz = xyz;
	}

	/**
	 * Create a new molecule
	 * 
	 * @param id
	 * @param mass
	 * @param x
	 * @param y
	 * @param z
	 */
	public MoleculeModel(int id, double mass, double x, double y, double z)
	{
		this.id = id;
		this.mass = mass;
		xyz = new double[] { x, y, z };
	}

	/**
	 * @return the x
	 */
	public double getX()
	{
		return xyz[0];
	}

	/**
	 * @return the y
	 */
	public double getY()
	{
		return xyz[1];
	}

	/**
	 * @return the z
	 */
	public double getZ()
	{
		return xyz[2];
	}

	/**
	 * @return the id
	 */
	public int getId()
	{
		return id;
	}

	/**
	 * Package level set method to allow renumbering
	 * @param id
	 */
	void setId(int id)
	{
		this.id = id;
	}

	/**
	 * @return the mass
	 */
	public double getMass()
	{
		return mass;
	}
	
	/**
	 * @return The coordinates (x,y,z)
	 */
	public double[] getCoordinates()
	{
		return xyz;
	}
	
	/**
	 * Move the molecule using a random Gaussian shift with standard deviation of the given diffusion rate.
	 * <p>
	 * Note: The array provided by {@link #getCoordinates()} is updated and returned.
	 * 
	 * @param diffusionRate
	 * @param random
	 * @return The new coordinates
	 */
	public double[] move(double diffusionRate, RandomGenerator random)
	{
		double[] xyz = getCoordinates();
		if (diffusionRate > 0)
		{
			for (int i = 0; i < 3; i++)
			{
				final double shift = random.nextGaussian() * diffusionRate;
				// Clip the movement
				//if (shift > 5*diffusionRate)
				//	xyz[i] += 5*diffusionRate;
				//else
				xyz[i] += shift;
			}
		}
		return xyz;
	}

	/**
	 * Move the molecule using a random walk with the given step size
	 * <p>
	 * Note: The array provided by {@link #getCoordinates()} is updated and returned.
	 * 
	 * @param stepSize
	 * @param random
	 * @return The new coordinates
	 */
	public double[] walk(double stepSize, RandomGenerator random)
	{
		double[] xyz = getCoordinates();
		if (stepSize > 0)
		{
			for (int i = 0; i < 3; i++)
			{
				if (random.nextDouble() < 0.5)
					xyz[i] += stepSize;
				else
					xyz[i] -= stepSize;
			}
		}
		return xyz;
	}
}
