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
 * Contains a model for an image of fixed lifetime fluorophores. Al fluorphores will have the same on time. The
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
	 * @see gdsc.smlm.model.ImageModel#createFluorophore(int, double[], int)
	 */
	@Override
	protected FluorophoreSequenceModel createFluorophore(int id, double[] xyz, int frames)
	{
		// Activate randomly within the start frame
		final double tAct = next + getRandom().getRandomGenerator().nextDouble();
		if (tAct >= frames)
			return null;
		// Round up to next frame
		next = Math.ceil(tAct + tOn + tOff);
		return new SimpleFluorophoreSequenceModel(id, xyz, tAct, tOn);
	}
}
