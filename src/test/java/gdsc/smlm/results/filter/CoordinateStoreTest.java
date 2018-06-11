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

import org.junit.Assert;
import org.junit.Test;

public class CoordinateStoreTest
{
	// TODO - test for a crop store ...

	@Test
	public void canCreateStore()
	{
		CoordinateStore s;
		s = CoordinateStoreFactory.create(1, 2, 10, 11, -1, -1);
		Assert.assertTrue(s instanceof NullCoordinateStore);
		s = CoordinateStoreFactory.create(1, 2, 10, 11, 0, -1);
		Assert.assertTrue(s instanceof GridCoordinateStore);
		s = CoordinateStoreFactory.create(1, 2, 10, 11, 0.2, -1);
		Assert.assertTrue(s instanceof GridCoordinateStore);
		s = CoordinateStoreFactory.create(1, 2, 10, 11, 0.5, -1);
		Assert.assertTrue(s instanceof GridCoordinateStore1);
		s = CoordinateStoreFactory.create(1, 2, 10, 11, 1, -1);
		Assert.assertTrue(s instanceof GridCoordinateStore1);
		s = CoordinateStoreFactory.create(1, 2, 10, 11, 1.5, -1);
		Assert.assertTrue(s instanceof GridCoordinateStore);
		s = CoordinateStoreFactory.create(1, 2, 10, 11, 2, -1);
		Assert.assertTrue(s instanceof GridCoordinateStore);
	}

	@Test
	public void canDetectXYDuplicates()
	{
		double[] datax = { 1.1, 4.1 };
		double[] datay = { 3.1, 7.1 };
		double[] dataz = { 0, 0.1 };
		double[] resolution = { 0, 0.3, 0.5, 1.5 };

		for (int i = 0; i < resolution.length; i++)
		{
			CoordinateStore s = CoordinateStoreFactory.create(1, 2, 10, 11, resolution[i], -1);
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
			s = CoordinateStoreFactory.create(1, 2, 10, 11, resolution[i], zResolution);
			s.add(x, y, z);

			msg = resolution[i] + "," + zResolution;
			Assert.assertTrue(msg, s.contains(x, y, z));

			// 3D exact match
			zResolution = 0;
			s = CoordinateStoreFactory.create(1, 2, 10, 11, resolution[i], zResolution);
			s.add(x, y, z);

			msg = resolution[i] + "," + zResolution;
			Assert.assertTrue(msg, s.contains(x, y, z));
			Assert.assertFalse(msg, s.contains(x, y, z + 0.01));
			Assert.assertFalse(msg, s.contains(x, y, z - 0.01));

			// 3D match within z-resolution
			zResolution = 1;
			s = CoordinateStoreFactory.create(1, 2, 10, 11, resolution[i], zResolution);
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
		double[] datax = { 1.1, 4.1 };
		double[] datay = { 3.1, 7.1 };
		double[] dataz = { 0, 0.1 };
		double[] resolution = { 0.3, 0.5, 1.5 };

		for (int i = 0; i < resolution.length; i++)
		{
			CoordinateStore s = CoordinateStoreFactory.create(1, 2, 10, 11, resolution[i], -1);
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
		double[] datax = { 1.1, 4.1 };
		double[] datay = { 3.1, 7.1 };
		double[] dataz = { 0, 0.1 };
		double[] resolution = { 0.3, 0.5, 1.5 };

		for (int i = 0; i < resolution.length; i++)
		{
			CoordinateStore s = CoordinateStoreFactory.create(1, 2, 10, 11, resolution[i], -1);

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
		CoordinateStore s = CoordinateStoreFactory.create(1, 2, 10, 11, 1, 0);
		s.add(s.getMinX() - 1, s.getMinY(), 0);
	}

	@Test(expected = IndexOutOfBoundsException.class)
	public void cannotAddOutsideGridXHigh()
	{
		CoordinateStore s = CoordinateStoreFactory.create(1, 2, 10, 11, 1, 0);
		s.add(s.getMinX() + s.getWidth() + 1, s.getMinY(), 0);
	}

	@Test(expected = IndexOutOfBoundsException.class)
	public void cannotAddOutsideGrid1YLow()
	{
		CoordinateStore s = CoordinateStoreFactory.create(1, 2, 10, 11, 1, 0);
		s.add(s.getMinX(), s.getMinY() - 1, 0);
	}

	@Test(expected = IndexOutOfBoundsException.class)
	public void cannotAddOutsideGridYHigh()
	{
		CoordinateStore s = CoordinateStoreFactory.create(1, 2, 10, 11, 1, 0);
		s.add(s.getMinX(), s.getMinY() + s.getHeight() + 1, 0);
	}

	@Test
	public void canSafeAddOutsideGrid()
	{
		GridCoordinateStore s = new GridCoordinateStore(1, 2, 10, 11, 1, 0.0);
		s.safeAdd(s.getMinX() - 1, s.getMinY(), 0);
		s.safeAdd(s.getMinX() + s.getWidth() + 1, s.getMinY(), 0);
		s.safeAdd(s.getMinX(), s.getMinY() - 1, 0);
		s.safeAdd(s.getMinX(), s.getMinY() + s.getHeight() + 1, 0);
	}

	@Test
	public void containsOutsideGridIsFalse()
	{
		GridCoordinateStore s = new GridCoordinateStore(1, 2, 10, 11, 1, 0.0);
		Assert.assertFalse(addAndFind(s, s.getMinX() - 1, s.getMinY(), 0));
		Assert.assertFalse(addAndFind(s, s.getMinX() + s.getWidth() + 1, s.getMinY(), 0));
		Assert.assertFalse(addAndFind(s, s.getMinX(), s.getMinY() - 1, 0));
		Assert.assertFalse(addAndFind(s, s.getMinX(), s.getMinY() + s.getHeight() + 1, 0));
	}

	private boolean addAndFind(GridCoordinateStore s, double x, double y, double z)
	{
		s.safeAdd(x, y, z);
		return s.contains(x, y, z);
	}

	@Test
	public void findOutsideGridIsNull()
	{
		GridCoordinateStore s = new GridCoordinateStore(1, 2, 10, 11, 1, 0.0);
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
		double[] datax = { 1.1, 4.1 };
		double[] datay = { 3.1, 7.1 };
		double[] dataz = { 0, 0.1 };
		double[] resolution = { 0, 0.3, 0.5, 1.5 };

		GridCoordinateStore s = new GridCoordinateStore(1, 2, 10, 11, 0, 0.0);
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
		double[] datax = { 1.1, 4.1 };
		double[] datay = { 3.1, 7.1 };
		double[] dataz = { 0, 0.1 };
		double[] resolution = { 0, 0.3, 0.5 };

		GridCoordinateStore1 s = new GridCoordinateStore1(1, 2, 10, 11, 0, 0.0);
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
		GridCoordinateStore1 s = new GridCoordinateStore1(1, 2, 10, 11, 0, 0.0);
		s.changeXYResolution(1.1);
	}
}
