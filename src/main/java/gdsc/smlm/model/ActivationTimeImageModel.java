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
	 * @see gdsc.smlm.model.ImageModel#createFluorophore(int, double[], int)
	 */
	@Override
	protected FluorophoreSequenceModel createFluorophore(int id, double[] xyz, int frames)
	{
		double t = getRandom().nextExponential(tAct);
		if (t >= frames)
			return null;
		return new StandardFluorophoreSequenceModel(id, xyz, tAct, tOn, tOff, tOff2, nBlinks,
				nBlinks2, isUseGeometricDistribution(), getRandom());
	}
}
