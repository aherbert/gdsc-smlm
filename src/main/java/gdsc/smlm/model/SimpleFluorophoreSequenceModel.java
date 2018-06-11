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
