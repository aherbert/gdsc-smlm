package gdsc.smlm.results;

import java.util.Arrays;
import java.util.Comparator;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

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

	private void canStoreResults(PeakResultStore store)
	{
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

			Assert.assertEquals(size, store.size());
			Assert.assertTrue(result == store.get(i));

			// Can convert to an array
			PeakResult[] array = store.toArray();
			Assert.assertEquals(size, array.length);
			for (int j = 0; j < size; j++)
				Assert.assertTrue(list[i] == array[i]);
		}

		// Can sort
		Arrays.sort(list, 0, size);
		store.sort();

		for (int i = 0; i < size; i++)
			Assert.assertTrue(list[i] == store.get(i));

		Comparator<PeakResult> c = new Comparator<PeakResult>()
		{
			public int compare(PeakResult o1, PeakResult o2)
			{
				return o1.origX - o2.origX;
			}
		};
		Arrays.sort(list, 0, size, c);
		store.sort(c);

		for (int i = 0; i < size; i++)
			Assert.assertTrue(list[i] == store.get(i));

		// Can trim to size
		store.trimToSize();
		Assert.assertEquals(size, store.size());

		for (int i = 0; i < size; i++)
			Assert.assertTrue(list[i] == store.get(i));

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
		Assert.assertEquals(size, store.size());

		for (int i = 0; i < size; i++)
			Assert.assertTrue(list[i] == store.get(i));

		// Can clear
		size = 0;
		store.clear();
		Assert.assertEquals(size, store.size());

		// Can add collection
		for (int i = 0; i < 10; i++)
		{
			result = create(r);
			list[size++] = result;
		}
		store.addCollection(Arrays.asList(Arrays.copyOf(list, size)));
		Assert.assertEquals(size, store.size());
		for (int i = 0; i < size; i++)
			Assert.assertTrue(list[i] == store.get(i));

		// Can add array
		size = 0;
		store.clear();
		for (int i = 0; i < 10; i++)
		{
			result = create(r);
			list[size++] = result;
		}
		store.addArray(Arrays.copyOf(list, size));
		Assert.assertEquals(size, store.size());
		for (int i = 0; i < size; i++)
			Assert.assertTrue(list[i] == store.get(i));

		// Can copy
		PeakResultStore copy = store.copy();
		Assert.assertNotEquals(copy, store);
		Assert.assertEquals(size, copy.size());
		for (int i = 0; i < size; i++)
			Assert.assertTrue(list[i] == copy.get(i));
	}

	private PeakResult create(RandomGenerator r)
	{
		return new PeakResult(r.nextInt(), r.nextInt(), r.nextInt(), r.nextFloat(), r.nextDouble(), r.nextFloat(),
				PeakResult.createParams(r.nextFloat(), r.nextFloat(), r.nextFloat(), r.nextFloat(), r.nextFloat()),
				null);
	}
}
