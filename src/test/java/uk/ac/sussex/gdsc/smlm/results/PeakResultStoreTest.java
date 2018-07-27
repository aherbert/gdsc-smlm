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
package uk.ac.sussex.gdsc.smlm.results;

import java.util.Arrays;
import java.util.Comparator;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;

import uk.ac.sussex.gdsc.core.utils.Random;
import uk.ac.sussex.gdsc.smlm.results.predicates.PeakResultPredicate;
import uk.ac.sussex.gdsc.smlm.results.sort.FrameIdPeakResultComparator;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;

@SuppressWarnings({ "javadoc" })
public class PeakResultStoreTest
{
	int capacity = 1;

	@SeededTest
	public void canStoreResultsUsingArrayList(RandomSeed seed)
	{
		canStoreResults(seed, new ArrayListPeakResultStore(capacity));
	}

	@SeededTest
	public void canStoreResultsUsingTurboList(RandomSeed seed)
	{
		canStoreResults(seed, new TurboListPeakResultStore(capacity));
	}

	@SeededTest
	public void canStoreResultsUsingArray(RandomSeed seed)
	{
		canStoreResults(seed, new ArrayPeakResultStore(capacity));
	}

	@SeededTest
	public void canStoreResultsUsingHashSet(RandomSeed seed)
	{
		canStoreResults(seed, new SetPeakResultStore(capacity));
	}

	@SuppressWarnings("null")
	private static void canStoreResults(RandomSeed seed, PeakResultStore store)
	{
		final boolean isList = store instanceof PeakResultStoreList;
		final PeakResultStoreList storeList = (isList) ? (PeakResultStoreList) store : null;
		PeakResult result;
		final UniformRandomProvider r = TestSettings.getRandomGenerator(seed.getSeed());

		PeakResult[] list = new PeakResult[20];
		int size = 0;

		Assertions.assertEquals(size, store.size(), "Not empty on construction");
		Assertions.assertEquals(size, store.toArray().length, "Not empty list");

		// Can store data in order
		for (int i = 0; i < 10; i++)
		{
			result = create(r);
			list[size++] = result;
			store.add(result);
		}

		assertEquals(list, size, store);

		// Can sort
		if (isList)
		{
			Arrays.sort(list, 0, size, FrameIdPeakResultComparator.INSTANCE);
			storeList.sort();

			for (int i = 0; i < size; i++)
				Assertions.assertTrue(list[i] == storeList.get(i), "List entry not same reference");

			final Comparator<PeakResult> c = new Comparator<PeakResult>()
			{
				@Override
				public int compare(PeakResult o1, PeakResult o2)
				{
					return o1.getOrigX() - o2.getOrigX();
				}
			};
			Arrays.sort(list, 0, size, c);
			storeList.sort(c);

			for (int i = 0; i < size; i++)
				Assertions.assertTrue(list[i] == storeList.get(i), "List entry not same reference after sort");
		}

		// Can trim to size
		store.trimToSize();

		assertEquals(list, size, store);

		// Can remove null results
		store.add(null);
		result = create(r);
		list[size++] = result;
		store.add(result);
		store.add(null);
		for (int i = 0; i < 3; i++)
		{
			result = create(r);
			list[size++] = result;
			store.add(result);
		}
		Assertions.assertTrue(size != store.size(), "Same size after adding null results");
		store.removeIf(new PeakResultPredicate()
		{
			@Override
			public boolean test(PeakResult t)
			{
				return t == null;
			}
		});
		assertEquals(list, size, store);

		// Can clear
		size = 0;
		store.clear();
		assertEquals(list, size, store);

		// Can add collection
		for (int i = 0; i < 10; i++)
		{
			result = create(r);
			list[size++] = result;
		}
		store.addCollection(Arrays.asList(Arrays.copyOf(list, size)));
		assertEquals(list, size, store);

		// Can add array
		size = 0;
		store.clear();
		for (int i = 0; i < 10; i++)
		{
			result = create(r);
			list[size++] = result;
		}
		store.addArray(Arrays.copyOf(list, size));
		assertEquals(list, size, store);

		// Can add store
		size = 0;
		store.clear();
		for (int i = 0; i < 10; i++)
		{
			result = create(r);
			list[size++] = result;
		}
		final PeakResultStore store2 = store.copy();
		store2.addArray(Arrays.copyOf(list, size));
		// Ensure the store is not at full capacity
		store2.remove(list[--size]);
		store.addStore(store2);
		assertEquals(list, size, store);

		// Can copy
		final PeakResultStore copy = store.copy();
		//Assertions.assertNotEquals(copy, store);
		Assertions.assertTrue(copy != store, "Copy is the same reference");
		assertEquals(list, size, store);

		// Can remove single
		for (int j = 0; j < 5; j++)
		{
			size = 0;
			store.clear();
			PeakResult toRemove = null;
			for (int i = 0; i < 5; i++)
			{
				result = create(r);
				if (j != i)
					list[size++] = result;
				else
					toRemove = result;
				store.add(result);
			}
			Assertions.assertTrue(size != store.size(), "Same size after adding extra single result");
			store.remove(toRemove);
			assertEquals(list, size, store);
		}

		// Can remove collection
		for (int j = 0; j < 5; j++)
		{
			size = 0;
			store.clear();
			for (int i = 0; i < 20; i++)
			{
				result = create(r);
				list[size++] = result;
				store.add(result);
			}

			final int[] indices = Random.sample(3, size, r);
			Arrays.sort(indices);
			final PeakResult[] list1 = new PeakResult[size - indices.length];
			final PeakResult[] toRemove = new PeakResult[indices.length];
			for (int i = 0; i < indices.length; i++)
				toRemove[i] = list[indices[i]];
			for (int i = 0, k = 0; i < size; i++)
				if (Arrays.binarySearch(indices, i) < 0)
					list1[k++] = list[i];
			store.removeCollection(Arrays.asList(toRemove));
			assertEquals(list1, list1.length, store);
		}

		// Can remove range
		if (!isList)
			return;

		size = 0;
		store.clear();
		for (int i = 0; i < 20; i++)
		{
			result = create(r);
			list[size++] = result;
			store.add(result);
		}

		// From the end
		storeList.remove(size - 1, size - 1);
		assertEquals(list, --size, store);

		storeList.remove(size - 2, size - 1);
		size -= 2;
		assertEquals(list, size, store);

		// From the start
		storeList.remove(0, 0);
		list = Arrays.copyOfRange(list, 1, size);
		size = list.length;
		assertEquals(list, size, store);

		storeList.remove(0, 1);
		list = Arrays.copyOfRange(list, 2, size);
		size = list.length;
		assertEquals(list, size, store);

		// From the middle
		storeList.remove(3, 3);
		System.arraycopy(list, 4, list, 3, size - 4);
		size--;
		assertEquals(list, size, store);

		storeList.remove(3, 4);
		System.arraycopy(list, 5, list, 3, size - 5);
		size -= 2;
		assertEquals(list, size, store);
	}

	private static PeakResult create(UniformRandomProvider r)
	{
		return new PeakResult(r.nextInt(), r.nextInt(), r.nextInt(), r.nextFloat(), r.nextDouble(), r.nextFloat(),
				r.nextFloat(),
				PeakResult.createParams(r.nextFloat(), r.nextFloat(), r.nextFloat(), r.nextFloat(), r.nextFloat()),
				null);
	}

	private static void assertEquals(PeakResult[] list, int size, PeakResultStore store)
	{
		Assertions.assertEquals(size, store.size(), "Not the same size");

		final boolean isList = store instanceof PeakResultStoreList;
		if (isList)
		{
			final PeakResultStoreList storeList = (PeakResultStoreList) store;
			for (int i = 0; i < size; i++)
			{
				Assertions.assertTrue(list[i] == storeList.get(i), "Not the same list index reference");
				Assertions.assertEquals(i, storeList.indexOf(storeList.get(i)), "indexOf finds wrong item");
			}
		}
		// Set equals
		final PeakResult[] list2 = store.toArray();
		Assertions.assertEquals(size, list2.length, "toArray() creates wrong size");
		for (int i = 0; i < size; i++)
			Assertions.assertTrue(contains(list, size, list2[i]), "Cannot find item in the array");
	}

	private static boolean contains(PeakResult[] list, int size, PeakResult peakResult)
	{
		for (int i = 0; i < size; i++)
			if (list[i] == peakResult)
				return true;
		return false;
	}
}
