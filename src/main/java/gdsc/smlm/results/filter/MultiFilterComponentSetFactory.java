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
package gdsc.smlm.results.filter;

import java.util.Arrays;


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
			case 0: return new MultiFilterComponentSet0(components); 
			case 1: return new MultiFilterComponentSet1(components); 
			case 2: return new MultiFilterComponentSet2(components); 
			case 3: return new MultiFilterComponentSet3(components); 
			case 4: return new MultiFilterComponentSet4(components); 
			case 5: return new MultiFilterComponentSet5(components); 
			case 6: return new MultiFilterComponentSet6(components); 
			case 7: return new MultiFilterComponentSet7(components); 
			//@formatter:on
		}
		return new MultiFilterComponentSetDefault(Arrays.copyOf(components, size));
	}
}
