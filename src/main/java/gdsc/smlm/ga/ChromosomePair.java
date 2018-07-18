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
package gdsc.smlm.ga;

/**
 * Define a pair of chromosomes.
 *
 * @param <T>
 *            the generic type
 */
public class ChromosomePair<T extends Comparable<T>>
{
	/** The first chromosome. */
	final Chromosome<T> c1;
	/** The second chromosome. */
	final Chromosome<T> c2;

	/**
	 * Create a pair.
	 *
	 * @param c1
	 *            the first chromosome
	 * @param c2
	 *            the second chromosome
	 * @throws IllegalArgumentException
	 *             if either chromosome is null
	 */
	public ChromosomePair(Chromosome<T> c1, Chromosome<T> c2)
	{
		if (c1 == null || c2 == null)
			throw new IllegalArgumentException("Chromosomes must not be null");
		this.c1 = c1;
		this.c2 = c2;
	}
}
