package gdsc.smlm.results.filter;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Contains a set of components of the multi filter.
 */
public class MultiFilterComponentSetFactory
{
	public static MultiFilterComponentSet create(MultiFilterComponent[] components, int size)
	{
		switch (size)
		{
			//@formatter:off
			case 1: return new MultiFilterComponentSet1(components); 
			case 2: return new MultiFilterComponentSet2(components); 
			case 3: return new MultiFilterComponentSet3(components); 
			case 4: return new MultiFilterComponentSet4(components); 
			case 5: return new MultiFilterComponentSet5(components); 
			case 6: return new MultiFilterComponentSet6(components); 
			//@formatter:on
		}
		return new MultiFilterComponentSet0(components);
	}
}