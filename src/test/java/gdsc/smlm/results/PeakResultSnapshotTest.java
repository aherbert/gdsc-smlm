package gdsc.smlm.results;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Assume;
import org.junit.Test;

public class PeakResultSnapshotTest
{
	@Test
	public void sameResultsAreEqual()
	{
		final RandomGenerator r = new Well19937c();
		PeakResult[] r1 = createResults(r, 10, 5, false, false, false, false);
		PeakResultsSnapshot snap = newPeakResultsSnapshot(r1);
		Assert.assertTrue(snap.matches(r1));
		Assert.assertTrue(snap.matches(snap));
	}

	@Test
	public void sameEmptyResultsAreEqual()
	{
		final RandomGenerator r = new Well19937c();
		PeakResult[] r1 = createResults(r, 0, 5, false, false, false, false);
		PeakResultsSnapshot snap = newPeakResultsSnapshot(r1);
		Assert.assertTrue(snap.matches(r1));
		Assert.assertTrue(snap.matches(snap));
	}

	@Test
	public void sameResultsAreEqualWithDeviation()
	{
		final RandomGenerator r = new Well19937c();
		PeakResult[] r1 = createResults(r, 10, 5, true, false, false, false);
		PeakResultsSnapshot snap = newPeakResultsSnapshot(r1);
		Assert.assertTrue(snap.matches(r1));
		Assert.assertTrue(snap.matches(snap));
	}

	@Test
	public void sameResultsAreEqualWithId()
	{
		final RandomGenerator r = new Well19937c();
		PeakResult[] r1 = createResults(r, 10, 5, false, true, false, false);
		PeakResultsSnapshot snap = newPeakResultsSnapshot(r1);
		Assert.assertTrue(snap.matches(r1));
		Assert.assertTrue(snap.matches(snap));
	}

	@Test
	public void sameResultsAreEqualWithEndFrame()
	{
		final RandomGenerator r = new Well19937c();
		PeakResult[] r1 = createResults(r, 10, 5, false, false, true, false);
		PeakResultsSnapshot snap = newPeakResultsSnapshot(r1);
		Assert.assertTrue(snap.matches(r1));
		Assert.assertTrue(snap.matches(snap));
	}

	@Test
	public void sameResultsAreEqualWithPrecision()
	{
		final RandomGenerator r = new Well19937c();
		PeakResult[] r1 = createResults(r, 10, 5, false, false, false, true);
		PeakResultsSnapshot snap = newPeakResultsSnapshot(r1);
		Assert.assertTrue(snap.matches(r1));
		Assert.assertTrue(snap.matches(snap));
	}

	@Test
	public void differentResultsAreNotEqual()
	{
		final RandomGenerator r = new Well19937c();
		PeakResult[] r1 = createResults(r, 10, 5, false, false, false, false);
		PeakResultsSnapshot snap = newPeakResultsSnapshot(r1);
		for (int size : new int[] { 10, 1, 0 })
		{
			PeakResult[] r2 = createResults(r, size, 5, false, false, false, false);
			Assert.assertFalse(snap.matches(r2));
			Assert.assertFalse(snap.matches(newPeakResultsSnapshot(r2)));
		}
	}

	@Test
	public void timeDigest()
	{
		Assume.assumeTrue(false);

		final RandomGenerator r = new Well19937c();
		PeakResultsSnapshot snap = new PeakResultsSnapshot();
		int N = 5;
		for (int size = 1000; size < 2000000; size *= 2)
		{
			PeakResult[] r1 = createResults(r, size, 5, false, false, false, false);
			long time = System.nanoTime();
			for (int i = N; i-- > 0;)
				snap.snapshot(-1, r1);
			time = System.nanoTime() - time;
			System.out.printf("size = %d, time = %g ms\n", size, (1e-6 * time) / N);
		}
	}

	private PeakResult[] createResults(RandomGenerator r, int size, int n, boolean withDeviations, boolean withId,
			boolean withEndFrame, boolean withPrecision)
	{
		final ArrayPeakResultStore store = new ArrayPeakResultStore(10);
		while (size-- > 0)
		{
			float[] params = createParams(n, r);
			float[] paramsDev = (withDeviations) ? createParams(n, r) : null;
			AttributePeakResult p = new AttributePeakResult(r.nextInt(), r.nextInt(), r.nextInt(), r.nextFloat(),
					r.nextDouble(), r.nextFloat(), params, paramsDev);
			if (withId)
				p.setId(r.nextInt());
			if (withEndFrame)
				//p.setEndFrame(p.getFrame() +  1 + r.nextInt(5));
				p.setEndFrame(r.nextInt());
			if (withPrecision)
				p.setPrecision(r.nextDouble());
			store.add(p);
		}
		return store.toArray();
	}

	private float[] createParams(int n, RandomGenerator r)
	{
		float[] p = new float[n];
		while (n-- > 0)
			p[n] = r.nextFloat();
		return p;
	}

	private PeakResultsSnapshot newPeakResultsSnapshot(PeakResult[] peakResults)
	{
		// For debug testing do not use a timeout
		PeakResultsSnapshot snap = new PeakResultsSnapshot();
		snap.snapshot(-1, peakResults);
		return snap;
	}
}
