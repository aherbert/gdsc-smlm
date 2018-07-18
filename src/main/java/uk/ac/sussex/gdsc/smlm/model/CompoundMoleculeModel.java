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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.random.RandomGenerator;

/**
 * Contains a model for a compound moving molecule (contains multiple molecules held in a fixed configuration).
 * <p>
 * The coordinates of the object represent the current centre of mass. The coordinates of each molecule can be retrieved
 * using the {@link #getCoordinates(int)} method;
 */
public class CompoundMoleculeModel extends MoleculeModel
{
	/** The constant to convert degrees to radians: <code>180.0 / Math.PI</code>. */
	static final double DEGREES_TO_RADIANS = 180.0 / Math.PI;

	/** Define the Z-axis */
	public static final double[] Z_AXIS = new double[] { 0, 0, 1 };

	private int label;

	/**
	 * Identity matrix for no rotation
	 */
	private static final double[] I = new double[] { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

	private List<? extends MoleculeModel> molecules;

	/**
	 * The fraction of a population that this molecule represents.
	 * <p>
	 * This is used when creating a particle distribution of different types of molecules. This is included so that the
	 * population description can be serialised.
	 */
	private double fraction = 1;

	/**
	 * The diffusion rate for the molecule
	 */
	private double diffusionRate = 0;
	/**
	 * The diffusion type for the molecule
	 */
	private DiffusionType diffusionType = DiffusionType.RANDOM_WALK;

	/**
	 * Create a new molecule
	 * <p>
	 * Note: molecules mass may be updated, see {@link #getMass()}.
	 *
	 * @param id
	 *            the id
	 * @param xyz
	 *            [x,y,z]
	 * @param molecules
	 *            the molecules
	 * @throws IllegalArgumentException
	 *             If the input list contains nulls
	 */
	public CompoundMoleculeModel(int id, double[] xyz, List<? extends MoleculeModel> molecules)
			throws IllegalArgumentException
	{
		this(id, xyz, molecules, true);
	}

	/**
	 * Create a new molecule
	 * <p>
	 * Note: molecules mass may be updated, see {@link #getMass()}.
	 *
	 * @param id
	 *            the id
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param z
	 *            the z
	 * @param molecules
	 *            the molecules
	 * @throws IllegalArgumentException
	 *             If the input list contains nulls
	 */
	public CompoundMoleculeModel(int id, double x, double y, double z, List<? extends MoleculeModel> molecules)
			throws IllegalArgumentException
	{
		this(id, x, y, z, molecules, true);
	}

	/**
	 * Create a new molecule
	 * <p>
	 * Note: molecules mass may be updated, see {@link #getMass()}.
	 *
	 * @param id
	 *            the id
	 * @param xyz
	 *            [x,y,z]
	 * @param molecules
	 *            the molecules
	 * @param centre
	 *            Shift molecules to have a centre-of-mass at 0,0,0.
	 * @throws IllegalArgumentException
	 *             If the input list contains nulls
	 */
	public CompoundMoleculeModel(int id, double[] xyz, List<? extends MoleculeModel> molecules, boolean centre)
			throws IllegalArgumentException
	{
		super(id, xyz);
		setMolecules(molecules, centre);
	}

	/**
	 * Create a new molecule
	 * <p>
	 * Note: molecules mass may be updated, see {@link #getMass()}.
	 *
	 * @param id
	 *            the id
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param z
	 *            the z
	 * @param molecules
	 *            the molecules
	 * @param centre
	 *            Shift molecules to have a centre-of-mass at 0,0,0.
	 * @throws IllegalArgumentException
	 *             If the input list contains nulls
	 */
	public CompoundMoleculeModel(int id, double x, double y, double z, List<? extends MoleculeModel> molecules,
			boolean centre) throws IllegalArgumentException
	{
		super(id, x, y, z);
		setMolecules(molecules, centre);
	}

	private void setMolecules(List<? extends MoleculeModel> molecules, boolean centre) throws IllegalArgumentException
	{
		if (molecules == null)
			molecules = new ArrayList<>(0);
		this.molecules = molecules;
		checkMass();

		if (centre)
			centre();
	}

	/**
	 * Check that all the molecules have a valid mass. If any are below zero then set to zero. If all are below zero
	 * then set to 1 so that the centre-of-mass is valid.
	 */
	public void checkMass()
	{
		int invalidMass = 0;
		for (final MoleculeModel m : molecules)
		{
			if (m == null)
				throw new IllegalArgumentException("Input list contains null molecules");
			if (m.mass <= 0 || Double.isNaN(m.mass))
				invalidMass++;
				//throw new IllegalArgumentException("Input list contains molecules with invalid mass: " + m.mass);
		}

		this.mass = 0;
		if (invalidMass > 0)
		{
			final double resetMass = (invalidMass == molecules.size()) ? 1 : 0;
			for (final MoleculeModel m : molecules)
			{
				if (m.mass <= 0 || Double.isNaN(m.mass))
					m.mass = resetMass;
				this.mass += m.mass;
			}
		}
		else
			for (final MoleculeModel m : molecules)
				this.mass += m.mass;
	}

	/**
	 * Shift molecules to have a centre-of-mass at 0,0,0.
	 * Coordinates are therefore relative to the origin.
	 */
	public void centre()
	{
		if (molecules.isEmpty())
			return;
		if (this.mass == 0)
			checkMass();

		final double[] com = new double[3];
		for (final MoleculeModel m : molecules)
		 for (int i = 0; i < 3; i++)
				com[i] += m.xyz[i] * m.mass;

		for (int i = 0; i < 3; i++)
			com[i] /= this.mass;

		for (final MoleculeModel m : molecules)
			for (int i = 0; i < 3; i++)
				m.xyz[i] -= com[i];
	}

	/**
	 * Rotate the molecule using a random axis and a random rotation angle. The rotation is around the centre-of-mass.
	 *
	 * @param maxAngle
	 *            The maximum angle to rotate (in either direction) in degrees
	 * @param random
	 *            the random
	 */
	public void rotateRandom(double maxAngle, RandomGenerator random)
	{
		if (maxAngle == 0 || getSize() < 2)
			return;
		if (this.mass == 0)
			checkMass();

		final double[] axis = new double[3];
		for (int i = 0; i < 3; i++)
			axis[i] = random.nextGaussian();

		final double angle = (-maxAngle + random.nextDouble() * 2.0 * maxAngle);
		if (angle == 0)
			return;

		rotateMolecules(axis, angle);
	}

	/**
	 * Rotate the molecule using a specified axis and a random rotation angle. The rotation is around the
	 * centre-of-mass.
	 *
	 * @param axis
	 *            The axis to rotate around
	 * @param maxAngle
	 *            The maximum angle to rotate (in either direction) in degrees
	 * @param random
	 *            the random
	 */
	public void rotateRandomAngle(double[] axis, double maxAngle, RandomGenerator random)
	{
		if (maxAngle == 0 || getSize() < 2)
			return;
		if (this.mass == 0)
			checkMass();

		final double angle = (-maxAngle + random.nextDouble() * 2.0 * maxAngle);
		if (angle == 0)
			return;

		rotateMolecules(axis, angle);
	}

	/**
	 * Rotate the molecule using a random axis and a specified rotation angle. The rotation is around the
	 * centre-of-mass.
	 *
	 * @param angle
	 *            The angle to rotate (in either direction) in degrees
	 * @param random
	 *            the random
	 */
	public void rotateRandomAxis(double angle, RandomGenerator random)
	{
		if (angle == 0 || getSize() < 2)
			return;
		if (this.mass == 0)
			checkMass();

		final double[] axis = new double[3];
		for (int i = 0; i < 3; i++)
			axis[i] = random.nextGaussian();

		rotateMolecules(axis, angle);
	}

	/**
	 * Rotate the molecule using a specified axis and rotation angle. The rotation is around the centre-of-mass.
	 *
	 * @param angle
	 *            The angle to rotate (in degrees)
	 * @param axis
	 *            The axis to rotate around
	 */
	public void rotate(double[] axis, double angle)
	{
		if (angle == 0 || getSize() < 2 || axis == null || axis.length < 3)
			return;
		if (this.mass == 0)
			checkMass();

		rotateMolecules(axis, angle);
	}

	/**
	 * Rotate the molecule using a specified axis and rotation angle. The rotation is around the centre-of-mass.
	 *
	 * @param angle
	 *            The angle to rotate (in degrees)
	 * @param axis
	 *            The axis to rotate around
	 */
	private void rotateMolecules(double[] axis, double angle)
	{
		final double[] r = getRotationMatrix(axis, angle);
		if (r == I)
			return;

		// Use the actual centre of mass of the molecules to avoid drift over time
		// (i.e. do not assume the centre is 0,0,0)
		final double[] com = new double[3];
		for (final MoleculeModel m : molecules)
			for (int i = 0; i < 3; i++)
				com[i] += m.xyz[i] * m.mass;
		for (int i = 0; i < 3; i++)
			com[i] /= this.mass;

		for (final MoleculeModel m : molecules)
		{
			final double xtmp = m.xyz[0] - com[0];
			final double ytmp = m.xyz[1] - com[1];
			final double ztmp = m.xyz[2] - com[2];
			// Do not add back COM since the molecules should be centred around the origin
			m.xyz[0] = xtmp * r[0] + ytmp * r[1] + ztmp * r[2];
			m.xyz[1] = xtmp * r[3] + ytmp * r[4] + ztmp * r[5];
			m.xyz[2] = xtmp * r[6] + ytmp * r[7] + ztmp * r[8];
		}
	}

	/**
	 * Get the rotation matrix for a rotation around an axis
	 *
	 * @param axis
	 *            axis
	 * @param angle
	 *            angle (in degrees)
	 * @return The rotation matrix
	 */
	private static double[] getRotationMatrix(double[] axis, double angle)
	{
		/* Set to unit length */
		double length = 0;
		for (int i = 0; i < 3; i++)
			length += axis[i] * axis[i];
		if (length == 0 || Double.isNaN(length))
			return I;
		length = Math.sqrt(length);
		for (int i = 0; i < 3; i++)
		 axis[i] /= length;

		/* Store the components and their squares */
		final double u = axis[0];
		final double u2 = u * u;
		final double v = axis[1];
		final double v2 = v * v;
		final double w = axis[2];
		final double w2 = w * w;
		final double uv = u * v;
		final double uw = u * w;
		final double vw = v * w;

		angle /= DEGREES_TO_RADIANS;
		final double cost = Math.cos(angle);
		final double sint = Math.sin(angle);

		final double[] r = new double[9];
		/* Set the rotation matrix */
		r[0] = u2 + ((v2 + w2) * cost);
		r[1] = uv - uv * cost - w * sint;
		r[2] = uw - uw * cost + v * sint;
		r[3] = uv - uv * cost + w * sint;
		r[4] = v2 + ((u2 + w2) * cost);
		r[5] = vw - vw * cost - u * sint;
		r[6] = -1.0 * (uw * (-1.0 + cost)) - v * sint;
		r[7] = -1.0 * (vw * (-1.0 + cost)) + u * sint;
		r[8] = w2 + ((u2 + v2) * cost);
		/*
		 * r[0] = u2 + ((v2 + w2)*cost);
		 * r[3] = uv - uv*cost - w*sint;
		 * r[6] = uw - uw*cost + v*sint;
		 * r[1] = uv - uv*cost + w*sint;
		 * r[4] = v2 + ((u2 + w2)*cost);
		 * r[7] = vw - vw*cost - u*sint;
		 * r[2] = -1.0*(uw*(-1.0 + cost)) - v*sint;
		 * r[5] = -1.0*(vw*(-1.0 + cost)) + u*sint;
		 * r[8] = w2 + ((u2 + v2)*cost);
		 */
		return r;
	}

	/**
	 * @return The number of molecules in the compound molecule
	 */
	public int getSize()
	{
		return molecules.size();
	}

	/**
	 * Get the current coordinates of the nth molecule in the compound molecule
	 *
	 * @param n
	 *            The requested molecule (0 <= n < {@link #getSize()})
	 * @return The xyz coordinates
	 * @see #getSize()
	 * @throws IndexOutOfBoundsException
	 *             If the requested molecule does not exist
	 */
	public double[] getCoordinates(int n) throws IndexOutOfBoundsException
	{
		final double[] xyz = Arrays.copyOf(this.xyz, 3);
		final MoleculeModel m = molecules.get(n);
		for (int i = 0; i < 3; i++)
			xyz[i] += m.xyz[i];
		return xyz;
	}

	/**
	 * Get the current coordinates of the nth molecule in the compound molecule relative to the centre-of-mass
	 *
	 * @param n
	 *            The requested molecule (0 <= n < {@link #getSize()})
	 * @return The xyz coordinates
	 * @see #getSize()
	 * @throws IndexOutOfBoundsException
	 *             If the requested molecule does not exist
	 */
	public double[] getRelativeCoordinates(int n) throws IndexOutOfBoundsException
	{
		final MoleculeModel m = molecules.get(n);
		return Arrays.copyOf(m.xyz, 3);
	}

	/**
	 * Get the nth molecule in the compound molecule
	 * <p>
	 * Note that the molecule coordinates are relative the centre-of-mass of the compound
	 *
	 * @param n
	 *            The requested molecule (0 <= n < {@link #getSize()})
	 * @return The molecule
	 * @see #getSize()
	 * @throws IndexOutOfBoundsException
	 *             If the requested molecule does not exist
	 */
	public MoleculeModel getMolecule(int n) throws IndexOutOfBoundsException
	{
		return molecules.get(n);
	}

	/**
	 * @return The fraction of a population that this molecule represents
	 */
	public double getFraction()
	{
		return fraction;
	}

	/**
	 * Set the fraction of a population that this molecule represents
	 *
	 * @param fraction
	 *            the fraction to set
	 */
	public void setFraction(double fraction)
	{
		this.fraction = fraction;
	}

	/**
	 * Return the of all the molecules.
	 * <p>
	 * The mass is calculated once during initialisation. If any molecule has a mass less than zero it will be set to
	 * zero. If all molecules have a mass of zero then the mass for each will be reset to one. This allows the
	 * centre-of-mass calculation to function correctly. If a molecule is part of the compound and has no mass it will
	 * be rotated and moved but will not contribute to the COM.
	 *
	 * @return The mass of all the molecules
	 */
	@Override
	public double getMass()
	{
		return mass;
	}

	/**
	 * Scale the molecules relative coordinates by the given factor.
	 *
	 * @param factor
	 *            the factor
	 */
	public void scale(double factor)
	{
		for (final MoleculeModel m : molecules)
			for (int i = 0; i < 3; i++)
				m.xyz[i] *= factor;
	}

	/**
	 * @return the diffusionRate
	 */
	public double getDiffusionRate()
	{
		return diffusionRate;
	}

	/**
	 * @param diffusionRate
	 *            the diffusionRate to set
	 */
	public void setDiffusionRate(double diffusionRate)
	{
		this.diffusionRate = diffusionRate;
	}

	/**
	 * Get the diffusion type
	 *
	 * @return The diffusion type
	 */
	public DiffusionType getDiffusionType()
	{
		return diffusionType;
	}

	/**
	 * Set the diffusion type
	 *
	 * @param diffusionType
	 *            The diffusion type
	 */
	public void setDiffusionType(DiffusionType diffusionType)
	{
		this.diffusionType = diffusionType;
	}

	/**
	 * Gets the label.
	 *
	 * @return the label
	 */
	@Override
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
	@Override
	public void setLabel(int label)
	{
		this.label = label;
	}
}
