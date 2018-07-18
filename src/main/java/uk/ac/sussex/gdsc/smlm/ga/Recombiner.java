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
package uk.ac.sussex.gdsc.smlm.ga;

/**
 * Defines recombination crossover of a chromosome pair.
 *
 * @param <T>
 *            the generic type
 */
public interface Recombiner<T extends Comparable<T>>
{
	/**
	 * Crossover the provided chromosomes to produce one or more new sequences.
	 *
	 * @param chromosome1
	 *            the chromosome 1
	 * @param chromosome2
	 *            the chromosome 2
	 * @return one or more new sequences
	 */
	public Chromosome<T>[] cross(Chromosome<T> chromosome1, Chromosome<T> chromosome2);
}
