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
package gdsc.smlm.results.sort;

import java.util.Comparator;

import gdsc.smlm.results.PeakResult;

/**
 * Compares the results by frame, end frame and then intensity descending
 */
public class FrameEndFrameIntensityPeakResultComparator implements Comparator<PeakResult>
{
	/** An instance of the comparator */
	public static final FrameEndFrameIntensityPeakResultComparator INSTANCE = new FrameEndFrameIntensityPeakResultComparator();

	/*
	 * (non-Javadoc)
	 *
	 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
	 */
	@Override
	public int compare(PeakResult o1, PeakResult o2)
	{
		int f1 = o1.getFrame();
		int f2 = o2.getFrame();
		if (f1 < f2)
			return -1;
		if (f1 > f2)
			return 1;
		f1 = o1.getEndFrame();
		f2 = o2.getEndFrame();
		if (f1 < f2)
			return -1;
		if (f1 > f2)
			return 1;
		return Float.compare(o2.getIntensity(), o1.getIntensity());
	}
}
