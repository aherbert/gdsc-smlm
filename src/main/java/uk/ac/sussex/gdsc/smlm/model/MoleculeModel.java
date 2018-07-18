/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
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
package uk.ac.sussex.gdsc.smlm.model;

import org.apache.commons.math3.random.RandomGenerator;

/**
 * Contains a model for a moving molecule.
 */
public class MoleculeModel
{
	/** The id. */
	private int id;

	/** The xyz. */
	protected double[] xyz;

	/** The mass. */
	protected double mass;

	/** The label. */
	private int label;

	/**
	 * Create a new molecule.
	 *
	 * @param id
	 *            the id
	 * @param xyz
	 *            [x,y,z]
	 */
	public MoleculeModel(int id, double[] xyz)
	{
		this.id = id;
		this.xyz = xyz;
	}

	/**
	 * Create a new molecule.
	 *
	 * @param id
	 *            the id
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param z
	 *            the z
	 */
	public MoleculeModel(int id, double x, double y, double z)
	{
		this.id = id;
		xyz = new double[] { x, y, z };
	}

	/**
	 * Create a new molecule.
	 *
	 * @param mass
	 *            the mass
	 * @param xyz
	 *            [x,y,z]
	 */
	public MoleculeModel(double mass, double[] xyz)
	{
		this.mass = mass;
		this.xyz = xyz;
	}

	/**
	 * Create a new molecule.
	 *
	 * @param mass
	 *            the mass
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param z
	 *            the z
	 */
	public MoleculeModel(double mass, double x, double y, double z)
	{
		this.mass = mass;
		xyz = new double[] { x, y, z };
	}

	/**
	 * Create a new molecule.
	 *
	 * @param id
	 *            the id
	 * @param mass
	 *            the mass
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
	 * Create a new molecule.
	 *
	 * @param id
	 *            the id
	 * @param mass
	 *            the mass
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param z
	 *            the z
	 */
	public MoleculeModel(int id, double mass, double x, double y, double z)
	{
		this.id = id;
		this.mass = mass;
		xyz = new double[] { x, y, z };
	}

	/**
	 * Gets the x.
	 *
	 * @return the x
	 */
	public double getX()
	{
		return xyz[0];
	}

	/**
	 * Gets the y.
	 *
	 * @return the y
	 */
	public double getY()
	{
		return xyz[1];
	}

	/**
	 * Gets the z.
	 *
	 * @return the z
	 */
	public double getZ()
	{
		return xyz[2];
	}

	/**
	 * Gets the id.
	 *
	 * @return the id
	 */
	public int getId()
	{
		return id;
	}

	/**
	 * Package level set method to allow renumbering.
	 *
	 * @param id
	 *            the new id
	 */
	void setId(int id)
	{
		this.id = id;
	}

	/**
	 * Gets the mass.
	 *
	 * @return the mass
	 */
	public double getMass()
	{
		return mass;
	}

	/**
	 * Gets the coordinates.
	 *
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
	 *            Diffusion rate for each dimension
	 * @param random
	 *            Random generator
	 * @return The new coordinates
	 */
	public double[] move(double diffusionRate, RandomGenerator random)
	{
		final double[] xyz = getCoordinates();
		if (diffusionRate > 0)
			for (int i = 0; i < 3; i++)
			{
				final double shift = random.nextGaussian() * diffusionRate;
				// Clip the movement
				//if (shift > 5*diffusionRate)
				//	xyz[i] += 5*diffusionRate;
				//else if (shift < -5*diffusionRate)
				//	xyz[i] -= 5*diffusionRate;
				//else
				xyz[i] += shift;
			}
		return xyz;
	}

	/**
	 * Move the molecule using a random Gaussian shift with standard deviation of the given diffusion rate.
	 * <p>
	 * Note: The array provided by {@link #getCoordinates()} is updated and returned.
	 *
	 * @param diffusionRate
	 *            Diffusion rate for each dimension
	 * @param random
	 *            Random generator (one per dimension)
	 * @return The new coordinates
	 */
	public double[] move(double diffusionRate, RandomGenerator[] random)
	{
		final double[] xyz = getCoordinates();
		if (diffusionRate > 0)
			for (int i = 0; i < 3; i++)
			{
				final double shift = random[i].nextGaussian() * diffusionRate;
				// Clip the movement
				//if (shift > 5*diffusionRate)
				//	xyz[i] += 5*diffusionRate;
				//else if (shift < -5*diffusionRate)
				//	xyz[i] -= 5*diffusionRate;
				//else
				xyz[i] += shift;
			}
		return xyz;
	}

	/**
	 * Move the molecule using a random walk with the given step size
	 * <p>
	 * Note: The array provided by {@link #getCoordinates()} is updated and returned.
	 *
	 * @param stepSize
	 *            Step size for each dimension
	 * @param random
	 *            Random generator
	 * @return The new coordinates
	 */
	public double[] walk(double stepSize, RandomGenerator random)
	{
		final double[] xyz = getCoordinates();
		if (stepSize > 0)
			for (int i = 0; i < 3; i++)
				if (random.nextDouble() < 0.5)
					xyz[i] += stepSize;
				else
					xyz[i] -= stepSize;
		return xyz;
	}

	/**
	 * Move the molecule using a random walk with the given step size
	 * <p>
	 * Note: The array provided by {@link #getCoordinates()} is updated and returned.
	 *
	 * @param stepSize
	 *            Step size for each dimension
	 * @param random
	 *            Random generator (one per dimension)
	 * @return The new coordinates
	 */
	public double[] walk(double stepSize, RandomGenerator[] random)
	{
		final double[] xyz = getCoordinates();
		if (stepSize > 0)
			for (int i = 0; i < 3; i++)
				if (random[i].nextDouble() < 0.5)
					xyz[i] += stepSize;
				else
					xyz[i] -= stepSize;
		return xyz;
	}

	/**
	 * Slide the molecule along a unit vector using a random Gaussian shift with standard deviation of the given
	 * diffusion rate.
	 * <p>
	 * Note: The array provided by {@link #getCoordinates()} is updated and returned.
	 *
	 * @param diffusionRate
	 *            Diffusion rate for 3D diffusion
	 * @param axis
	 *            The linear axis to move along (must be a unit vector)
	 * @param random
	 *            Random number generator
	 * @return The new coordinates
	 */
	public double[] slide(double diffusionRate, double[] axis, RandomGenerator random)
	{
		final double[] xyz = getCoordinates();
		if (diffusionRate > 0)
		{
			final double shift;
			// Sample from a Gaussian - This may only be relevant for 1D diffusion
			shift = random.nextGaussian() * diffusionRate;

			// Sample from the cumulative probability distribution for the MSD.
			// Then get a square root to find the shift and assign a direction
			//RandomDataGenerator r = new RandomDataGenerator(random);
			//shift = ((random.nextDouble() < 0.5) ? 1 : -1) * Math.sqrt(r.nextExponential(diffusionRate*diffusionRate));

			// Clip the movement
			//if (shift > 5*diffusionRate)
			//	shift = 5*diffusionRate;
			//else if (shift < -5*diffusionRate)
			//	shift = -5*diffusionRate;
			for (int i = 0; i < 3; i++)
				xyz[i] += shift * axis[i];
		}
		return xyz;
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
