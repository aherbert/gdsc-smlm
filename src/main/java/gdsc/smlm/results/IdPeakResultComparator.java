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

import java.util.Comparator;

/**
 * Compares the results by Id
 */
public class IdPeakResultComparator implements Comparator<PeakResult>
{
	/** An instance of the comparator */
	public static final IdPeakResultComparator INSTANCE = new IdPeakResultComparator();

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
	 */
	@Override
	public int compare(PeakResult o1, PeakResult o2)
	{
		int id1 = o1.getId();
		int id2 = o2.getId();
		if (id1 < id2)
			return -1;
		if (id1 == id2)
			// Sort by frame
			return Integer.compare(o1.getFrame(), o2.getFrame());
		return 1;
	}
}
