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

import java.util.Collection;

import gdsc.core.data.DataException;
import gdsc.smlm.results.predicates.PeakResultPredicate;
import gdsc.smlm.results.procedures.PeakResultProcedure;

/**
 * Stores peak results and prevents modification.
 */
public class ImmutablePeakResultStore implements PeakResultStore
{
	private final PeakResultStore store;

	/**
	 * Instantiates a new immutable peak result store.
	 *
	 * @param store
	 *            the store
	 */
	public ImmutablePeakResultStore(PeakResultStore store)
	{
		if (store == null)
			throw new NullPointerException("Store must not be null");
		this.store = store;
	}

	@Override
	public int size()
	{
		return store.size();
	}

	@Override
	public boolean add(PeakResult result)
	{
		throw new DataException("This result store is immutable");
	}

	@Override
	public boolean addCollection(Collection<PeakResult> results)
	{
		throw new DataException("This result store is immutable");
	}

	@Override
	public boolean addArray(PeakResult[] results)
	{
		throw new DataException("This result store is immutable");
	}

	@Override
	public boolean addStore(PeakResultStore results)
	{
		throw new DataException("This result store is immutable");
	}

	@Override
	public boolean remove(PeakResult result)
	{
		throw new DataException("This result store is immutable");
	}

	@Override
	public boolean removeCollection(Collection<PeakResult> results)
	{
		throw new DataException("This result store is immutable");
	}

	@Override
	public boolean removeArray(PeakResult[] results)
	{
		throw new DataException("This result store is immutable");
	}

	@Override
	public boolean removeStore(PeakResultStore results)
	{
		throw new DataException("This result store is immutable");
	}

	@Override
	public boolean retainCollection(Collection<PeakResult> results)
	{
		throw new DataException("This result store is immutable");
	}

	@Override
	public boolean retainArray(PeakResult[] results)
	{
		throw new DataException("This result store is immutable");
	}

	@Override
	public boolean retainStore(PeakResultStore results)
	{
		throw new DataException("This result store is immutable");
	}

	@Override
	public void clear()
	{
		throw new DataException("This result store is immutable");
	}

	@Override
	public void trimToSize()
	{
		store.trimToSize();
	}

	@Override
	public PeakResult[] toArray()
	{
		return makeImmutable(store.toArray());
	}

	/**
	 * Make the array a collection of immutable peak result objects.
	 *
	 * @param array
	 *            the array
	 * @return the array
	 */
	private static PeakResult[] makeImmutable(PeakResult[] array)
	{
		for (int i = 0; i < array.length; i++)
			array[i] = new ImmutablePeakResult(array[i]);
		return array;
	}

	@Override
	public PeakResultStore copy()
	{
		return new ImmutablePeakResultStore(store.copy());
	}

	@Override
	public PeakResultStore copy(boolean deepCopy)
	{
		return new ImmutablePeakResultStore(store.copy(deepCopy));
	}

	@Override
	public boolean removeIf(PeakResultPredicate filter)
	{
		throw new DataException("This result store is immutable");
	}

	/**
	 * Used to wrap the results to make them immutable
	 */
	private static class ImmutablePeakResultProcedure implements PeakResultProcedure
	{
		PeakResultProcedure procedure;

		public ImmutablePeakResultProcedure(PeakResultProcedure procedure)
		{
			this.procedure = procedure;
		}

		@Override
		public void execute(PeakResult peakResult)
		{
			procedure.execute(new ImmutablePeakResult(peakResult));
		}
	}

	@Override
	public void forEach(PeakResultProcedure procedure)
	{
		store.forEach(new ImmutablePeakResultProcedure(procedure));
	}

	@Override
	public PeakResult[] subset(PeakResultPredicate filter)
	{
		return makeImmutable(store.subset(filter));
	}

	@Override
	public boolean contains(PeakResult result)
	{
		return store.contains(result);
	}
}
