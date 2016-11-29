package gdsc.smlm.results.filter;

import org.junit.Test;

import org.junit.Assert;

public class CoordinateStoreTest
{
	@Test
	public void canCreateStore()
	{
		CoordinateStore s;
		s = CoordinateStoreFactory.create(10, 10, 0);
		Assert.assertTrue(s instanceof NullCoordinateStore);
		s = CoordinateStoreFactory.create(10, 10, 0.2);
		Assert.assertTrue(s instanceof GridCoordinateStore);
		s = CoordinateStoreFactory.create(10, 10, 0.5);
		Assert.assertTrue(s instanceof GridCoordinateStore1);
		s = CoordinateStoreFactory.create(10, 10, 1);
		Assert.assertTrue(s instanceof GridCoordinateStore1);
		s = CoordinateStoreFactory.create(10, 10, 1.5);
		Assert.assertTrue(s instanceof GridCoordinateStore);
		s = CoordinateStoreFactory.create(10, 10, 2);
		Assert.assertTrue(s instanceof GridCoordinateStore);
	}
	
	@Test
	public void canDetectDuplicates()
	{
		double[] datax = { 0.1, 4.1 };
		double[] datay = { 3.1, 7.1 };
		double[] resolution = { 0.3, 0.5, 1.5 };

		for (int i = 0; i < resolution.length; i++)
		{
			CoordinateStore s = CoordinateStoreFactory.create(10, 10, resolution[i]);
			for (int j = 0; j < datax.length; j++)
				s.add(datax[j], datay[j]);

			for (int j = 0; j < datax.length; j++)
			{
				String msg = resolution[i] + " [" + j + "]"; 
				Assert.assertTrue(msg, s.contains(datax[j], datay[j]));
				Assert.assertTrue(msg, s.contains(datax[j] + resolution[i] * 0.99, datay[j]));
				Assert.assertTrue(msg, s.contains(datax[j], datay[j] + resolution[i] * 0.99));
				Assert.assertFalse(msg, s.contains(datax[j] + resolution[i], datay[j] + resolution[i]));
			}
		}
	}

	@Test
	public void canQueueToGrid()
	{
		double[] datax = { 0.1, 4.1 };
		double[] datay = { 3.1, 7.1 };
		double[] resolution = { 0.3, 0.5, 1.5 };

		for (int i = 0; i < resolution.length; i++)
		{
			CoordinateStore s = CoordinateStoreFactory.create(10, 10, resolution[i]);
			for (int j = 0; j < datax.length; j++)
			{
				s.addToQueue(datax[j], datay[j]);
				Assert.assertFalse(s.contains(datax[j], datay[j]));
				for (int k = 0; k < j; k++)
					Assert.assertTrue(s.contains(datax[k], datay[k]));
				s.flush();
				for (int k = 0; k <= j; k++)
					Assert.assertTrue(s.contains(datax[k], datay[k]));
			}
		}
	}

	@Test
	public void canClearGrid()
	{
		double[] datax = { 0.1, 4.1 };
		double[] datay = { 3.1, 7.1 };
		double[] resolution = { 0.3, 0.5, 1.5 };

		for (int i = 0; i < resolution.length; i++)
		{
			CoordinateStore s = CoordinateStoreFactory.create(10, 10, resolution[i]);

			// Add then clear
			for (int j = 0; j < datax.length; j++)
			{
				s.add(datax[j], datay[j]);
				Assert.assertTrue(s.contains(datax[j], datay[j]));
			}

			s.clear();

			for (int j = 0; j < datax.length; j++)
				Assert.assertFalse(s.contains(datax[j], datay[j]));

			// Queue then flush then clear
			for (int j = 0; j < datax.length; j++)
			{
				s.addToQueue(datax[j], datay[j]);
				Assert.assertFalse(s.contains(datax[j], datay[j]));
			}
			s.flush();
			for (int j = 0; j < datax.length; j++)
				Assert.assertTrue(s.contains(datax[j], datay[j]));

			s.clear();

			for (int j = 0; j < datax.length; j++)
				Assert.assertFalse(s.contains(datax[j], datay[j]));
		}
	}

	@Test(expected = IndexOutOfBoundsException.class)
	public void cannotAddOutsideGrid1XLow()
	{
		CoordinateStore s = CoordinateStoreFactory.create(10, 10, 1);
		s.add(-1, 0);
	}

	@Test(expected = IndexOutOfBoundsException.class)
	public void cannotAddOutsideGridXHigh()
	{
		CoordinateStore s = CoordinateStoreFactory.create(10, 10, 1);
		s.add(11, 0);
	}

	@Test(expected = IndexOutOfBoundsException.class)
	public void cannotAddOutsideGrid1YLow()
	{
		CoordinateStore s = CoordinateStoreFactory.create(10, 10, 1);
		s.add(0, -1);
	}

	@Test(expected = IndexOutOfBoundsException.class)
	public void cannotAddOutsideGridYHigh()
	{
		CoordinateStore s = CoordinateStoreFactory.create(10, 10, 1);
		s.add(0, 11);
	}
	
	@Test
	public void canSafeAddOutsideGrid()
	{
		GridCoordinateStore s = new GridCoordinateStore(10, 10, 1);
		s.safeAdd(-1, 0);
		s.safeAdd(11, 0);
		s.safeAdd(-1, 0);
		s.safeAdd(0, 11);
	}
	
	@Test
	public void containsOutsideGridIsFalse()
	{
		GridCoordinateStore s = new GridCoordinateStore(10, 10, 1);
		s.safeAdd(-1, 0);
		s.safeAdd(11, 0);
		s.safeAdd(-1, 0);
		s.safeAdd(0, 11);
		Assert.assertFalse(s.contains(-1, 0));
		Assert.assertFalse(s.contains(11, 0));
		Assert.assertFalse(s.contains(-1, 0));
		Assert.assertFalse(s.contains(0, 11));
	}
	
	@Test
	public void findOutsideGridIsNull()
	{
		GridCoordinateStore s = new GridCoordinateStore(10, 10, 1);
		s.safeAdd(-1, 0);
		s.safeAdd(11, 0);
		s.safeAdd(-1, 0);
		s.safeAdd(0, 11);
		Assert.assertNull(s.find(-1, 0));
		Assert.assertNull(s.find(11, 0));
		Assert.assertNull(s.find(-1, 0));
		Assert.assertNull(s.find(0, 11));
	}

	@Test
	public void canChangeResolution()
	{
		double[] datax = { 0.1, 4.1 };
		double[] datay = { 3.1, 7.1 };
		double[] resolution = { 0.3, 0.5, 1.5 };

		GridCoordinateStore s = new GridCoordinateStore(10, 10, 0);
		for (int i = 0; i < resolution.length; i++)
		{
			s.changeResolution(resolution[i]);

			for (int j = 0; j < datax.length; j++)
				s.add(datax[j], datay[j]);

			for (int j = 0; j < datax.length; j++)
			{
				String msg = resolution[i] + " [" + j + "]"; 
				Assert.assertTrue(msg, s.contains(datax[j], datay[j]));
				Assert.assertTrue(msg, s.contains(datax[j] + resolution[i] * 0.99, datay[j]));
				Assert.assertTrue(msg, s.contains(datax[j], datay[j] + resolution[i] * 0.99));
				Assert.assertFalse(msg, s.contains(datax[j] + resolution[i], datay[j] + resolution[i]));
			}
		}
	}
}
