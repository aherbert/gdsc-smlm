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
 * Contains a model for an image of blinking fluorophores under pulsed activation illumination. Activation energy is
 * sampled from an exponential distribution. Fluorophores are created when the activation energy has been achieved under
 * the given illumination.
 * <p>
 * Based on the work of Coltharp et al (2012) Accurate Construction of photoactivated localization microscopy images for
 * quantitative measurements. PLOS One 7, Issue 12, pp 1-15
 */
public class ActivationEnergyImageModel extends ImageModel
{
	private double eAct;
	private SpatialIllumination illumination;

	/**
	 * Construct a new image model
	 *
	 * @param eAct
	 *            Average energy for activation
	 * @param illumination
	 *            The illumination model
	 * @param tOn
	 *            Average on-state time
	 * @param tOff
	 *            Average off-state time for the first dark state
	 * @param tOff2
	 *            Average off-state time for the second dark state
	 * @param nBlinks
	 *            Average number of blinks int the first dark state (used for each burst between second dark states)
	 * @param nBlinks2
	 *            Average number of blinks into the second dark state
	 */
	public ActivationEnergyImageModel(double eAct, SpatialIllumination illumination, double tOn, double tOff,
			double tOff2, double nBlinks, double nBlinks2)
	{
		super(tOn, tOff, tOff2, nBlinks, nBlinks2);
		init(eAct, illumination);
	}

	private void init(double eAct, SpatialIllumination illumination)
	{
		checkParameter("eAct", eAct);
		if (illumination == null)
			throw new IllegalArgumentException("SpatialIllumination is null");
		this.eAct = eAct;
		this.illumination = illumination;
	}

	/**
	 * @return the average energy for activation
	 */
	public double getActivationEnergy()
	{
		return eAct;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.ImageModel#createActivationTime(double[])
	 */
	@Override
	protected double createActivationTime(double[] xyz)
	{
		return getActivationTime(xyz, frameLimit);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.ImageModel#createFluorophore(int, double[], double)
	 */
	@Override
	protected FluorophoreSequenceModel createFluorophore(int id, double[] xyz, double tAct)
	{
		return new StandardFluorophoreSequenceModel(id, xyz, tAct, tOn, tOff, tOff2, nBlinks, nBlinks2,
				isUseGeometricDistribution(), getRandom());
	}

	private double getActivationTime(double[] xyz, int frames)
	{
		final double activation = getRandom().nextExponential(eAct);
		double e = 0;
		for (int t = 0; t < frames; t++)
		{
			// Q. Should the molecule be moving during the activation phase?
			final double[] photons = illumination.getPulsedPhotons(xyz, t + 1);

			e += photons[0]; // pulse energy
			if (e > activation)
				return t;

			e += photons[1]; // during energy
			if (e > activation)
				// Interpolate
				return t + 1 - (e - activation) / photons[1];
		}
		return frames; // default to the number of frames.
	}
}
