package gdsc.smlm.results.filter;

import org.junit.Test;

import org.junit.Assert;

public class GridCoordinateStoreTest
{
	@Test
	public void canDetectDuplicates()
	{
		double[] datax = { 0.1, 3.1 };
		double[] datay = { 3.1, 5.1 };
		double[] resolution = { 0.5, 1 };

		for (int i = 0; i < resolution.length; i++)
		{
			GridCoordinateStore s = new GridCoordinateStore(10, 10, resolution[i]);
			for (int j = 0; j < datax.length; j++)
				s.add(datax[j], datay[j]);

			for (int j = 0; j < datax.length; j++)
			{
				Assert.assertTrue(s.contains(datax[j], datay[j]));
				Assert.assertTrue(s.contains(datax[j] + resolution[i] * 0.99, datay[j]));
				Assert.assertTrue(s.contains(datax[j], datay[j] + resolution[i] * 0.99));
				Assert.assertFalse(s.contains(datax[j] + resolution[i], datay[j] + resolution[i]));
			}
		}
	}

	@Test
	public void canQueueToGrid()
	{
		double[] datax = { 0.1, 3.1 };
		double[] datay = { 3.1, 5.1 };
		double[] resolution = { 0.5, 1 };

		for (int i = 0; i < resolution.length; i++)
		{
			GridCoordinateStore s = new GridCoordinateStore(10, 10, resolution[i]);
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
		double[] datax = { 0.1, 3.1 };
		double[] datay = { 3.1, 5.1 };
		double[] resolution = { 0.5, 1 };

		for (int i = 0; i < resolution.length; i++)
		{
			GridCoordinateStore s = new GridCoordinateStore(10, 10, resolution[i]);
			
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
}
