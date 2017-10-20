package gdsc.smlm.results;

import java.util.Comparator;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Compares the results by Id
 */
public class IdPeakResultComparator implements Comparator<PeakResult>
{
	/*
	 * (non-Javadoc)
	 * 
	 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
	 */
	public int compare(PeakResult o1, PeakResult o2)
	{
		return Integer.compare(o1.getId(), o2.getId());
	}
}
