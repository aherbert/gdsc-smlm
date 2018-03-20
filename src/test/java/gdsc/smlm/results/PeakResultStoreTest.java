package gdsc.smlm.results;

import java.util.Arrays;
import java.util.Comparator;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.Random;
import gdsc.smlm.results.predicates.PeakResultPredicate;

public class PeakResultStoreTest
{
	int capacity = 1;

	@Test
	public void canStoreResultsUsingArrayList()
	{
		canStoreResults(new ArrayListPeakResultStore(capacity));
	}

	@Test
	public void canStoreResultsUsingTurboList()
	{
		canStoreResults(new TurboListPeakResultStore(capacity));
	}

	@Test
	public void canStoreResultsUsingArray()
	{
		canStoreResults(new ArrayPeakResultStore(capacity));
	}

	@Test
	public void canStoreResultsUsingHashSet()
	{
		canStoreResults(new SetPeakResultStore(capacity));
	}

	private void canStoreResults(PeakResultStore store)
	{
		final boolean isList = store instanceof PeakResultStoreList;
		PeakResultStoreList storeList = (isList) ? (PeakResultStoreList) store : null;
		PeakResult result;
		RandomGenerator r = new Well19937c(30051977);

		PeakResult[] list = new PeakResult[20];
		int size = 0;

		Assert.assertEquals(size, store.size());
		Assert.assertEquals(size, store.toArray().length);

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
			Arrays.sort(list, 0, size);
			storeList.sort();

			for (int i = 0; i < size; i++)
				Assert.assertTrue(list[i] == storeList.get(i));

			Comparator<PeakResult> c = new Comparator<PeakResult>()
			{
				public int compare(PeakResult o1, PeakResult o2)
				{
					return o1.getOrigX() - o2.getOrigX();
				}
			};
			Arrays.sort(list, 0, size, c);
			storeList.sort(c);

			for (int i = 0; i < size; i++)
				Assert.assertTrue(list[i] == storeList.get(i));
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
		Assert.assertNotEquals(size, store.size());
		store.removeIf(new PeakResultPredicate()
		{
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

		// Can copy
		PeakResultStore copy = store.copy();
		Assert.assertNotEquals(copy, store);
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
			Assert.assertNotEquals(size, store.size());
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

			int[] indices = Random.sample(3, size, r);
			Arrays.sort(indices);
			PeakResult[] list1 = new PeakResult[size - indices.length];
			PeakResult[] toRemove = new PeakResult[indices.length];
			for (int i = 0; i < indices.length; i++)
			{
				toRemove[i] = list[indices[i]];
			}
			for (int i = 0, k = 0; i < size; i++)
			{
				if (Arrays.binarySearch(indices, i) < 0)
					list1[k++] = list[i];
			}
			store.removeCollection(Arrays.asList(toRemove));
			assertEquals(list1, list1.length, store);
		}
	}

	private PeakResult create(RandomGenerator r)
	{
		return new PeakResult(r.nextInt(), r.nextInt(), r.nextInt(), r.nextFloat(), r.nextDouble(), r.nextFloat(),
				PeakResult.createParams(r.nextFloat(), r.nextFloat(), r.nextFloat(), r.nextFloat(), r.nextFloat()),
				null);
	}

	private void assertEquals(PeakResult[] list, int size, PeakResultStore store)
	{
		Assert.assertEquals(size, store.size());

		final boolean isList = store instanceof PeakResultStoreList;
		if (isList)
		{
			PeakResultStoreList storeList = (PeakResultStoreList) store;
			for (int i = 0; i < size; i++)
			{
				Assert.assertTrue(list[i] == storeList.get(i));
				Assert.assertEquals(i, storeList.indexOf(storeList.get(i)));
			}
		}
		// Set equals
		PeakResult[] list2 = store.toArray();
		Assert.assertEquals(size, list2.length);
		for (int i = 0; i < size; i++)
			Assert.assertTrue(contains(list, size, list2[i]));
	}

	private boolean contains(PeakResult[] list, int size, PeakResult peakResult)
	{
		for (int i = 0; i < size; i++)
			if (list[i] == peakResult)
				return true;
		return false;
	}
}
