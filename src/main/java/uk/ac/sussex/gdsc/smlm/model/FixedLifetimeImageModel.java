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

/**
 * Contains a model for an image of fixed lifetime fluorophores. All fluorphores will have the same on time. The
 * activation time will be the incremented by the on-time plus the dark time between fluorophores.
 */
public class FixedLifetimeImageModel extends ImageModel
{
	private double next = 0;

	/**
	 * Construct a new image model
	 *
	 * @param tOn
	 *            Fixed on-state time
	 * @param tOff
	 *            Dark time between successive fluorophores
	 */
	public FixedLifetimeImageModel(double tOn, double tOff)
	{
		super(tOn, tOff, 0, 0, 0);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.ImageModel#createActivationTime(double[])
	 */
	@Override
	protected double createActivationTime(double[] xyz)
	{
		final double tAct = next + getRandom().getRandomGenerator().nextDouble();
		// Ensure at least tOff full dark frames between lifetimes:
		// Frames:    |      |      |      |
		//         ------|
		//              end             |--------
		//                              start
		//                     tOff
		final int unit = (int) Math.ceil(tOff);
		final int endT = (int) Math.ceil((tAct + tOn) / unit);
		next = (endT + 1) * unit;
		return tAct;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.ImageModel#createFluorophore(int, double[], double)
	 */
	@Override
	protected FluorophoreSequenceModel createFluorophore(int id, double[] xyz, double tAct)
	{
		return new SimpleFluorophoreSequenceModel(id, xyz, tAct, tOn);
	}
}
