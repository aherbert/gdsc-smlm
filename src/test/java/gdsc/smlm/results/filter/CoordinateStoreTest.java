package gdsc.smlm.results.filter;

import org.junit.Test;

import org.junit.Assert;

public class CoordinateStoreTest
{
	@Test
	public void canCreateStore()
	{
		CoordinateStore s;
		s = CoordinateStoreFactory.create(10, 10, -1, -1);
		Assert.assertTrue(s instanceof NullCoordinateStore);
		s = CoordinateStoreFactory.create(10, 10, 0, -1);
		Assert.assertTrue(s instanceof GridCoordinateStore);
		s = CoordinateStoreFactory.create(10, 10, 0.2, -1);
		Assert.assertTrue(s instanceof GridCoordinateStore);
		s = CoordinateStoreFactory.create(10, 10, 0.5, -1);
		Assert.assertTrue(s instanceof GridCoordinateStore1);
		s = CoordinateStoreFactory.create(10, 10, 1, -1);
		Assert.assertTrue(s instanceof GridCoordinateStore1);
		s = CoordinateStoreFactory.create(10, 10, 1.5, -1);
		Assert.assertTrue(s instanceof GridCoordinateStore);
		s = CoordinateStoreFactory.create(10, 10, 2, -1);
		Assert.assertTrue(s instanceof GridCoordinateStore);
	}

	@Test
	public void canDetectXYDuplicates()
	{
		double[] datax = { 0.1, 4.1 };
		double[] datay = { 3.1, 7.1 };
		double[] dataz = { 0, 0.1 };
		double[] resolution = { 0, 0.3, 0.5, 1.5 };

		for (int i = 0; i < resolution.length; i++)
		{
			CoordinateStore s = CoordinateStoreFactory.create(10, 10, resolution[i], -1);
			for (int j = 0; j < datax.length; j++)
				s.add(datax[j], datay[j], dataz[j]);

			for (int j = 0; j < datax.length; j++)
			{
				String msg = resolution[i] + " [" + j + "]";
				//@formatter:off
				Assert.assertTrue(msg, s.contains(datax[j], datay[j], dataz[j]));
				Assert.assertTrue(msg, s.contains(datax[j] + resolution[i] * 0.99, datay[j], dataz[j]));
				Assert.assertTrue(msg, s.contains(datax[j], datay[j] + resolution[i] * 0.99, dataz[j]));
				Assert.assertFalse(msg,	s.contains(datax[j] + increase(resolution[i], 1.01), datay[j] + resolution[i], dataz[j]));
				Assert.assertFalse(msg,	s.contains(datax[j] + resolution[i], datay[j] + increase(resolution[i], 1.01), dataz[j]));
				//@formatter:on
			}
		}
	}

	private double increase(double value, double delta)
	{
		return (value == 0) ? delta : value * delta;
	}

	@Test
	public void canDetectZDuplicates()
	{
		double x = 3.1, y = 4.3, z = 1.1;
		double[] resolution = { 0.3, 0.5, 1.5 };

		double zResolution;
		CoordinateStore s;
		String msg;

		for (int i = 0; i < resolution.length; i++)
		{
			// No 3D
			zResolution = -1;
			s = CoordinateStoreFactory.create(10, 10, resolution[i], zResolution);
			s.add(x, y, z);

			msg = resolution[i] + "," + zResolution;
			Assert.assertTrue(msg, s.contains(x, y, z));

			// 3D exact match
			zResolution = 0;
			s = CoordinateStoreFactory.create(10, 10, resolution[i], zResolution);
			s.add(x, y, z);

			msg = resolution[i] + "," + zResolution;
			Assert.assertTrue(msg, s.contains(x, y, z));
			Assert.assertFalse(msg, s.contains(x, y, z + 0.01));
			Assert.assertFalse(msg, s.contains(x, y, z - 0.01));

			// 3D match within z-resolution
			zResolution = 1;
			s = CoordinateStoreFactory.create(10, 10, resolution[i], zResolution);
			s.add(x, y, z);

			msg = resolution[i] + "," + zResolution;
			Assert.assertTrue(msg, s.contains(x, y, z));
			Assert.assertTrue(msg, s.contains(x, y, z + zResolution));
			Assert.assertTrue(msg, s.contains(x, y, z - zResolution));
			Assert.assertFalse(msg, s.contains(x, y, z + zResolution * 1.01));
			Assert.assertFalse(msg, s.contains(x, y, z - zResolution * 1.01));
		}
	}

	@Test
	public void canQueueToGrid()
	{
		double[] datax = { 0.1, 4.1 };
		double[] datay = { 3.1, 7.1 };
		double[] dataz = { 0, 0.1 };
		double[] resolution = { 0.3, 0.5, 1.5 };

		for (int i = 0; i < resolution.length; i++)
		{
			CoordinateStore s = CoordinateStoreFactory.create(10, 10, resolution[i], -1);
			for (int j = 0; j < datax.length; j++)
			{
				s.addToQueue(datax[j], datay[j], dataz[j]);
				Assert.assertFalse(s.contains(datax[j], datay[j], dataz[j]));
				for (int k = 0; k < j; k++)
					Assert.assertTrue(s.contains(datax[k], datay[k], dataz[j]));
				s.flush();
				for (int k = 0; k <= j; k++)
					Assert.assertTrue(s.contains(datax[k], datay[k], dataz[j]));
			}
		}
	}

	@Test
	public void canClearGrid()
	{
		double[] datax = { 0.1, 4.1 };
		double[] datay = { 3.1, 7.1 };
		double[] dataz = { 0, 0.1 };
		double[] resolution = { 0.3, 0.5, 1.5 };

		for (int i = 0; i < resolution.length; i++)
		{
			CoordinateStore s = CoordinateStoreFactory.create(10, 10, resolution[i], -1);

			// Add then clear
			for (int j = 0; j < datax.length; j++)
			{
				s.add(datax[j], datay[j], dataz[j]);
				Assert.assertTrue(s.contains(datax[j], datay[j], dataz[j]));
			}

			s.clear();

			for (int j = 0; j < datax.length; j++)
				Assert.assertFalse(s.contains(datax[j], datay[j], dataz[j]));

			// Queue then flush then clear
			for (int j = 0; j < datax.length; j++)
			{
				s.addToQueue(datax[j], datay[j], dataz[j]);
				Assert.assertFalse(s.contains(datax[j], datay[j], dataz[j]));
			}
			s.flush();
			for (int j = 0; j < datax.length; j++)
				Assert.assertTrue(s.contains(datax[j], datay[j], dataz[j]));

			s.clear();

			for (int j = 0; j < datax.length; j++)
				Assert.assertFalse(s.contains(datax[j], datay[j], dataz[j]));
		}
	}

	@Test(expected = IndexOutOfBoundsException.class)
	public void cannotAddOutsideGrid1XLow()
	{
		CoordinateStore s = CoordinateStoreFactory.create(10, 10, 1, 0);
		s.add(-1, 0, 0);
	}

	@Test(expected = IndexOutOfBoundsException.class)
	public void cannotAddOutsideGridXHigh()
	{
		CoordinateStore s = CoordinateStoreFactory.create(10, 10, 1, 0);
		s.add(11, 0, 0);
	}

	@Test(expected = IndexOutOfBoundsException.class)
	public void cannotAddOutsideGrid1YLow()
	{
		CoordinateStore s = CoordinateStoreFactory.create(10, 10, 1, 0);
		s.add(0, -1, 0);
	}

	@Test(expected = IndexOutOfBoundsException.class)
	public void cannotAddOutsideGridYHigh()
	{
		CoordinateStore s = CoordinateStoreFactory.create(10, 10, 1, 0);
		s.add(0, 11, 0);
	}

	@Test
	public void canSafeAddOutsideGrid()
	{
		GridCoordinateStore s = new GridCoordinateStore(10, 10, 1, 0.0);
		s.safeAdd(-1, 0, 0);
		s.safeAdd(11, 0, 0);
		s.safeAdd(-1, 0, 0);
		s.safeAdd(0, 11, 0);
	}

	@Test
	public void containsOutsideGridIsFalse()
	{
		GridCoordinateStore s = new GridCoordinateStore(10, 10, 1, 0.0);
		s.safeAdd(-1, 0, 0);
		s.safeAdd(11, 0, 0);
		s.safeAdd(-1, 0, 0);
		s.safeAdd(0, 11, 0);
		Assert.assertFalse(s.contains(-1, 0, 0));
		Assert.assertFalse(s.contains(11, 0, 0));
		Assert.assertFalse(s.contains(-1, 0, 0));
		Assert.assertFalse(s.contains(0, 11, 0));
	}

	@Test
	public void findOutsideGridIsNull()
	{
		GridCoordinateStore s = new GridCoordinateStore(10, 10, 1, 0.0);
		s.safeAdd(-1, 0, 0);
		s.safeAdd(11, 0, 0);
		s.safeAdd(-1, 0, 0);
		s.safeAdd(0, 11, 0);
		Assert.assertNull(s.find(-1, 0, 0));
		Assert.assertNull(s.find(11, 0, 0));
		Assert.assertNull(s.find(-1, 0, 0));
		Assert.assertNull(s.find(0, 11, 0));
	}

	@Test
	public void canChangeXYResolution()
	{
		double[] datax = { 0.1, 4.1 };
		double[] datay = { 3.1, 7.1 };
		double[] dataz = { 0, 0.1 };
		double[] resolution = { 0, 0.3, 0.5, 1.5 };

		GridCoordinateStore s = new GridCoordinateStore(10, 10, 0, 0.0);
		for (int i = 0; i < resolution.length; i++)
		{
			s.changeXYResolution(resolution[i]);

			for (int j = 0; j < datax.length; j++)
				s.add(datax[j], datay[j], dataz[j]);

			for (int j = 0; j < datax.length; j++)
			{
				String msg = resolution[i] + " [" + j + "]";
				//@formatter:off
				Assert.assertTrue(msg, s.contains(datax[j], datay[j], dataz[j]));
				Assert.assertTrue(msg, s.contains(datax[j] + resolution[i] * 0.99, datay[j], dataz[j]));
				Assert.assertTrue(msg, s.contains(datax[j], datay[j] + resolution[i] * 0.99, dataz[j]));
				Assert.assertFalse(msg,	s.contains(datax[j] + increase(resolution[i], 1.01), datay[j] + resolution[i], dataz[j]));
				Assert.assertFalse(msg,	s.contains(datax[j] + resolution[i], datay[j] + increase(resolution[i], 1.01), dataz[j]));
				//@formatter:on
			}
		}
	}

	@Test
	public void canChangeXYResolutionOnFixedStore()
	{
		double[] datax = { 0.1, 4.1 };
		double[] datay = { 3.1, 7.1 };
		double[] dataz = { 0, 0.1 };
		double[] resolution = { 0, 0.3, 0.5 };

		GridCoordinateStore1 s = new GridCoordinateStore1(10, 10, 0, 0.0);
		for (int i = 0; i < resolution.length; i++)
		{
			s.changeXYResolution(resolution[i]);

			for (int j = 0; j < datax.length; j++)
				s.add(datax[j], datay[j], dataz[j]);

			for (int j = 0; j < datax.length; j++)
			{
				String msg = resolution[i] + " [" + j + "]";
				//@formatter:off
				Assert.assertTrue(msg, s.contains(datax[j], datay[j], dataz[j]));
				Assert.assertTrue(msg, s.contains(datax[j] + resolution[i] * 0.99, datay[j], dataz[j]));
				Assert.assertTrue(msg, s.contains(datax[j], datay[j] + resolution[i] * 0.99, dataz[j]));
				Assert.assertFalse(msg,	s.contains(datax[j] + increase(resolution[i], 1.01), datay[j] + resolution[i], dataz[j]));
				Assert.assertFalse(msg,	s.contains(datax[j] + resolution[i], datay[j] + increase(resolution[i], 1.01), dataz[j]));
				//@formatter:on
			}
		}
	}

	@Test(expected = IllegalArgumentException.class)
	public void cannotChangeToBadXYResolutionOnFixedStore()
	{
		GridCoordinateStore1 s = new GridCoordinateStore1(10, 10, 0, 0.0);
		s.changeXYResolution(1.1);
	}
}
