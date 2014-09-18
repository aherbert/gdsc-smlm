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
	 * @param time
	 *            The time the fluorophore is on
	 */
	public SimpleFluorophoreSequenceModel(int id, double[] xyz, int time)
	{
		super(id, xyz);
		setBurstSequence(new double[] { time, time });
	}
}
