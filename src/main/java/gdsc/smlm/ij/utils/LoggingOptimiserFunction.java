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
package gdsc.smlm.ij.utils;

import gdsc.smlm.function.OptimiserFunction;
import ij.IJ;

/**
 * Allow progress tracking of the Apache Commons Math 3 Optimiser in ImageJ
 */
public abstract class LoggingOptimiserFunction extends OptimiserFunction
{
	private boolean logging = false;
	private int evalCount = 0;
	protected String name = "Optimiser";

	public LoggingOptimiserFunction(String name)
	{
		this.name = name;
	}

	/**
	 * Log the count of evaluations to the ImageJ status bar
	 *
	 * @param b
	 */
	public void setLogging(boolean b)
	{
		logging = b;
		if (b)
			evalCount = 0;
	}

	public void increment()
	{
		evalCount++;
		if (logging)
			IJ.showStatus(name + " Evaluation " + evalCount);
	}

	/**
	 * @return The function name
	 */
	public String getName()
	{
		return name;
	}
}
