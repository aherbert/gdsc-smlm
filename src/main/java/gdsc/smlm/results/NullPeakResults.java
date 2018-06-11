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
package gdsc.smlm.results;


/**
 * Does nothing for any of the PeakResults methods
 */
public class NullPeakResults extends AbstractPeakResults implements ThreadSafePeakResults
{
	public void begin()
	{
	}

	public void add(int peak, int origX, int origY, float origValue, double error, float noise, float meanSignal,
			float[] params, float[] paramsStdDev)
	{
	}

	public void add(PeakResult result)
	{
	}

	public void addAll(PeakResult[] results)
	{
	}

	public int size()
	{
		return 0;
	}

	public void end()
	{
	}

	public boolean isActive()
	{
		return true;
	}
}
