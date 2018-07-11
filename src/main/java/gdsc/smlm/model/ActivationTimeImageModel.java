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
package gdsc.smlm.model;

/**
 * Contains a model for an image of blinking fluorophores under constant activation illumination.
 * <p>
 * Based on the work of Coltharp et al (2012) Accurate Construction of photoactivated localization microscopy images for
 * quantitative measurements. PLOS One 7, Issue 12, pp 1-15
 */
public class ActivationTimeImageModel extends ImageModel
{
	private double tAct;

	/**
	 * Construct a new image model
	 *
	 * @param tAct
	 *            Average time for activation
	 * @param tOn
	 *            Average on-state time
	 * @param tOn
	 *            Average on-state time
	 * @param tOff
	 *            Average off-state time for the first dark state
	 * @param tOff
	 *            Average off-state time for the second dark state
	 * @param nBlinks
	 *            Average number of blinks int the first dark state (used for each burst between second dark states)
	 * @param nBlinks2
	 *            Average number of blinks into the second dark state
	 */
	public ActivationTimeImageModel(double tAct, double tOn, double tOff, double tOff2, double nBlinks, double nBlinks2)
	{
		super(tOn, tOff, tOff2, nBlinks, nBlinks2);
		init(tAct);
	}

	private void init(double tAct)
	{
		checkParameter("tAct", tAct);
		this.tAct = tAct;
	}

	/**
	 * @return the tAct
	 */
	public double gettAct()
	{
		return tAct;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.ImageModel#createActivationTime(double[])
	 */
	@Override
	protected double createActivationTime(double[] xyz)
	{
		return getRandom().nextExponential(tAct);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.ImageModel#createFluorophore(int, double[], double)
	 */
	@Override
	protected FluorophoreSequenceModel createFluorophore(int id, double[] xyz, double tAct)
	{
		return new StandardFluorophoreSequenceModel(id, xyz, tAct, tOn, tOff, tOff2, nBlinks, nBlinks2,
				isUseGeometricDistribution(), getRandom());
	}
}
