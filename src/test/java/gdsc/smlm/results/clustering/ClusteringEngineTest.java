package gdsc.smlm.results.clustering;

import gdsc.smlm.ij.utils.Utils;

import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

public class ClusteringEngineTest
{
	//private RandomGenerator rand = new Well19937c(System.currentTimeMillis() + System.identityHashCode(this));
	private RandomGenerator rand = new Well19937c(30051977);

	// Store the closest pair of clusters
	int ii, jj;

	@Test
	public void canClusterClusterPointsAtDifferentDensitiesUsingClosest()
	{
		for (double radius : new double[] { 5, 10, 20 })
		{
			for (int size : new int[] { 2000, 1000, 500, 400, 300, 200, 100 })
			{
				testClusting(ClusteringAlgorithm.Closest, radius, 100, size);
			}
		}
	}

	@Test
	public void canClusterClusterPointsAtDifferentDensitiesUsingPairwiseWithoutNeighbours()
	{
		for (double radius : new double[] { 5, 10, 20 })
		{
			for (int size : new int[] { 2000, 1000, 500, 400, 300, 200, 100 })
			{
				testClusting(ClusteringAlgorithm.PairwiseWithoutNeighbours, radius, 100, size);
			}
		}
	}

	@Test
	public void pairwiseWithoutNeighboursIsFasterAtLowDensities()
	{
		int Repeats = 10;
		double radius = 50;
		Object[] points = new Object[Repeats];
		for (int i = 0; i < Repeats; i++)
			points[i] = createClusters(50, 1000, 2, radius / 2);

		long t1 = runSpeedTest(points, ClusteringAlgorithm.Closest, radius);
		long t2 = runSpeedTest(points, ClusteringAlgorithm.PairwiseWithoutNeighbours, radius);

		System.out.printf("Closest %d, PairwiseWithoutNeighbours %d = %fx faster\n", t1, t2, (double) t1 / t2);
		Assert.assertTrue(t2 < t1);
	}

	@Test
	public void pairwiseWithoutNeighboursIsSlowerAtHighDensities()
	{
		int Repeats = 10;
		double radius = 50;
		Object[] points = new Object[Repeats];
		for (int i = 0; i < Repeats; i++)
			points[i] = createClusters(500, 1000, 2, radius / 2);

		long t1 = runSpeedTest(points, ClusteringAlgorithm.Closest, radius);
		long t2 = runSpeedTest(points, ClusteringAlgorithm.PairwiseWithoutNeighbours, radius);

		System.out.printf("Closest %d, PairwiseWithoutNeighbours %d = %fx faster\n", t1, t2, (double) t1 / t2);
		Assert.assertTrue(t1 < t2);
	}

	@Test
	public void pairwiseIsFaster()
	{
		int Repeats = 20;
		Object[] points = new Object[Repeats];
		for (int i = 0; i < Repeats; i++)
			points[i] = createPoints(500, 1000);
		double radius = 50;

		long t1 = runSpeedTest(points, ClusteringAlgorithm.Closest, radius);
		long t2 = runSpeedTest(points, ClusteringAlgorithm.Pairwise, radius);

		System.out.printf("Closest %d, Pairwise %d = %fx faster\n", t1, t2, (double) t1 / t2);
		Assert.assertTrue(t2 < t1);
	}

	private long runSpeedTest(Object[] points, ClusteringAlgorithm algorithm, double radius)
	{
		ClusteringEngine engine = new ClusteringEngine();
		engine.setClusteringAlgorithm(algorithm);

		// Initialise
		engine.findClusters((ArrayList<ClusterPoint>) points[0], radius);

		long start = System.nanoTime();
		for (int i = 0; i < points.length; i++)
			engine.findClusters((ArrayList<ClusterPoint>) points[i], radius);
		return System.nanoTime() - start;
	}

	private void testClusting(ClusteringAlgorithm algorithm, double radius, int n, int size)
	{
		ClusteringEngine engine = new ClusteringEngine();
		engine.setClusteringAlgorithm(algorithm);
		ArrayList<ClusterPoint> points = createPoints(n, size);

		// Report density of the clustering we are testing. Size/radius are in nm
		System.out.printf("Testing n=%d, Size=%d, Density=%s um^-2, Radius=%s nm\n", n, size,
				Utils.rounded(n * 1e6 / (size * size)), Utils.rounded(radius));

		ArrayList<Cluster> exp = findClusters(points, radius);
		ArrayList<Cluster> obs = engine.findClusters(points, radius);
		Collections.sort(exp);
		Collections.sort(obs);

		try
		{
			Assert.assertEquals("# clusters is different", exp.size(), obs.size());
			for (int i = 0; i < exp.size(); i++)
			{
				assertEqual(i, exp.get(i), obs.get(i));
			}
		}
		catch (AssertionError e)
		{
			print("Expected", exp);
			print("Observed", obs);
			throw e;
		}
	}

	private void print(String name, ArrayList<Cluster> clusters)
	{
		System.out.printf(name + " : size=%d\n", clusters.size());
		for (int i = 0; i < clusters.size(); i++)
		{
			Cluster c = clusters.get(i);
			System.out.printf("[%d] : head=%d, n=%d, cx=%g, cy=%g\n", i, c.head.id, c.n, c.x, c.y);
		}
	}

	private void assertEqual(int i, Cluster cluster, Cluster cluster2)
	{
		Assert.assertEquals(i + " cluster: Size is different", cluster.n, cluster2.n);
		Assert.assertEquals(i + " cluster: X is different", cluster.x, cluster2.x, 1e-4);
		Assert.assertEquals(i + " cluster: Y is different", cluster.y, cluster2.y, 1e-4);
		// Q. Should we check each cluster member is the same ?
	}

	/**
	 * Perform centroid-linkage clustering up to the given radius
	 * 
	 * @param points
	 * @param radius
	 * @return The clusters
	 */
	private ArrayList<Cluster> findClusters(ArrayList<ClusterPoint> points, double radius)
	{
		// Initialise all clusters with one molecule
		ArrayList<Cluster> clusters = new ArrayList<Cluster>(points.size());
		for (int i = 0; i < points.size(); i++)
		{
			final ClusterPoint m = points.get(i);
			clusters.add(new Cluster(new ClusterPoint(i, m.x, m.y)));
		}

		// Iteratively find the closest pair
		while (findClosest(clusters, radius))
		{
			clusters.get(ii).add(clusters.get(jj));
			clusters.remove(jj);
		}

		return clusters;
	}

	/**
	 * Implement and all-vs-all search for the closest pair of clusters within the given radius. Set the class level
	 * variables ii and jj to the indices of the closest pair.
	 * 
	 * @param clusters
	 * @param radius
	 * @return True if a pair was found
	 */
	private boolean findClosest(ArrayList<Cluster> clusters, double radius)
	{
		double minD = radius * radius;
		ii = -1;
		for (int i = 0; i < clusters.size(); i++)
		{
			Cluster c1 = clusters.get(i);
			for (int j = i + 1; j < clusters.size(); j++)
			{
				final double d2 = c1.distance2(clusters.get(j));
				if (d2 < minD)
				{
					ii = i;
					jj = j;
					minD = d2;
				}
			}
		}

		return ii > -1;
	}

	/**
	 * Create n points in a 2D distribution of size * size.
	 * 
	 * @param n
	 * @param size
	 * @return The points
	 */
	private ArrayList<ClusterPoint> createPoints(int n, int size)
	{
		ArrayList<ClusterPoint> points = new ArrayList<ClusterPoint>(n);
		while (n-- > 0)
			points.add(new ClusterPoint(n, rand.nextDouble() * size, rand.nextDouble() * size));
		return points;
	}

	/**
	 * Create n clusters of m points in a 2D distribution of size * size. Clusters will be spread in a radius*radius
	 * square.
	 * 
	 * @param n
	 * @param size
	 * @param m
	 * @param radius
	 * @return The points
	 */
	private ArrayList<ClusterPoint> createClusters(int n, int size, int m, double radius)
	{
		ArrayList<ClusterPoint> points = new ArrayList<ClusterPoint>(n);
		int id = 0;
		while (n-- > 0)
		{
			double x = rand.nextDouble() * size;
			double y = rand.nextDouble() * size;
			for (int i = m; i-- > 0;)
			{
				points.add(new ClusterPoint(id++, x + rand.nextDouble() * radius, y + rand.nextDouble() * radius));
			}
		}
		return points;
	}
}
