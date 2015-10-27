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
 * Contains a single-time model for a non-blinking fluorophore
 */
public class SimpleFluorophoreSequenceModel extends FluorophoreSequenceModel
{
	/**
	 * Construct a new flourophore
	 * 
	 * @param id
	 *            The identifier
	 * @param xyz
	 *            The [x,y,z] coordinates
	 * @param tAct
	 *            The time the fluorophore turned on
	 * @param tOn
	 *            The time the fluorophore was on
	 */
	public SimpleFluorophoreSequenceModel(int id, double[] xyz, double tAct, double tOn)
	{
		super(id, xyz);
		if (tOn < 0)
			tOn = 0;
		setBurstSequence(new double[] { tAct, tAct + tOn });
	}
}
