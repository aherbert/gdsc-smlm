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

import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;

import java.util.logging.Logger;

@SuppressWarnings({ "javadoc" })
public class PeakResultDigestTest
{
    @SeededTest
    public void sameResultsAreEqual(RandomSeed seed)
    {
        final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
        final PeakResult[] r1 = createResults(r, 10, 5, false, false, false, false);
        final PeakResultsDigest digest = new PeakResultsDigest(r1);
        Assertions.assertTrue(digest.matches(r1));
        Assertions.assertTrue(digest.matches(digest));
    }

    @SeededTest
    public void sameSize1ResultsAreEqual(RandomSeed seed)
    {
        final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
        final PeakResult[] r1 = createResults(r, 1, 5, false, false, false, false);
        final PeakResultsDigest digest = new PeakResultsDigest(r1);
        Assertions.assertTrue(digest.matches(r1));
        Assertions.assertTrue(digest.matches(digest));
    }

    @SeededTest
    public void sameEmptyResultsAreEqual(RandomSeed seed)
    {
        final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
        final PeakResult[] r1 = createResults(r, 0, 5, false, false, false, false);
        final PeakResultsDigest digest = new PeakResultsDigest(r1);
        Assertions.assertTrue(digest.matches(r1));
        Assertions.assertTrue(digest.matches(digest));
    }

    @SeededTest
    public void sameResultsAreEqualWithDeviation(RandomSeed seed)
    {
        final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
        final PeakResult[] r1 = createResults(r, 10, 5, true, false, false, false);
        final PeakResultsDigest digest = new PeakResultsDigest(r1);
        Assertions.assertTrue(digest.matches(r1));
        Assertions.assertTrue(digest.matches(digest));
    }

    @SeededTest
    public void sameResultsAreEqualWithId(RandomSeed seed)
    {
        final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
        final PeakResult[] r1 = createResults(r, 10, 5, false, true, false, false);
        final PeakResultsDigest digest = new PeakResultsDigest(r1);
        Assertions.assertTrue(digest.matches(r1));
        Assertions.assertTrue(digest.matches(digest));
    }

    @SeededTest
    public void sameResultsAreEqualWithEndFrame(RandomSeed seed)
    {
        final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
        final PeakResult[] r1 = createResults(r, 10, 5, false, false, true, false);
        final PeakResultsDigest digest = new PeakResultsDigest(r1);
        Assertions.assertTrue(digest.matches(r1));
        Assertions.assertTrue(digest.matches(digest));
    }

    @SeededTest
    public void sameResultsAreEqualWithPrecision(RandomSeed seed)
    {
        final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
        final PeakResult[] r1 = createResults(r, 10, 5, false, false, false, true);
        final PeakResultsDigest digest = new PeakResultsDigest(r1);
        Assertions.assertTrue(digest.matches(r1));
        Assertions.assertTrue(digest.matches(digest));
    }

    @SeededTest
    public void differentResultsAreNotEqual(RandomSeed seed)
    {
        final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
        final PeakResult[] r1 = createResults(r, 10, 5, false, false, false, false);
        final PeakResultsDigest digest = new PeakResultsDigest(r1);
        for (final int size : new int[] { 10, 1, 0 })
        {
            final PeakResult[] r2 = createResults(r, size, 5, false, false, false, false);
            Assertions.assertFalse(digest.matches(r2));
            Assertions.assertFalse(digest.matches(new PeakResultsDigest(r2)));
        }
    }

    @SeededTest
    public void digestMatchesPeakResultDigest(RandomSeed seed)
    {
        final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
        for (int size = 1; size < 5; size++)
        {
            final PeakResult[] r1 = createResults(r, size, 5, false, false, false, false);
            final PeakResultsDigest digest = new PeakResultsDigest(r1);

            final PeakResultDigest d = new PeakResultDigest();
            for (final PeakResult rr : r1)
                d.update(rr);
            Assertions.assertEquals(d.digest(), digest.getDigest());
        }
    }

    @SeededTest
    public void digestIsEmptyStringWhenSizeIsZero()
    {
        Assertions.assertEquals("", new PeakResultsDigest(new PeakResult[0]).getDigest());
    }

    @SeededTest
    public void digestHandlesNull()
    {
        final PeakResult[] r1 = null;
        final PeakResult[] r0 = new PeakResult[0];
        final PeakResultsDigest digest = new PeakResultsDigest(r1);
        Assertions.assertTrue(digest.matches(r1));
        Assertions.assertFalse(digest.matches(r0));
    }

    @SeededTest
    public void digestHandlesEmptyArray()
    {
        final PeakResult[] r1 = null;
        final PeakResult[] r0 = new PeakResult[0];
        final PeakResultsDigest digest = new PeakResultsDigest(r0);
        Assertions.assertTrue(digest.matches(r0));
        Assertions.assertFalse(digest.matches(r1));
    }

    @SeededTest
    public void timeDigest(RandomSeed seed)
    {
        Assumptions.assumeTrue(false);
        final Logger logger = Logger.getLogger(PeakResultDigestTest.class.getName());

        final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
        final PeakResultsDigest digest = new PeakResultsDigest();
        final int N = 5;
        for (int size = 1000; size < 2000000; size *= 2)
        {
            final PeakResult[] r1 = createResults(r, size, 5, false, false, false, false);
            long time = System.nanoTime();
            for (int i = N; i-- > 0;)
                digest.digest(r1);
            time = System.nanoTime() - time;
            logger.info(FunctionUtils.getSupplier("size = %d, time = %g ms", size, (1e-6 * time) / N));
        }
    }

    private static PeakResult[] createResults(UniformRandomProvider r, int size, int n, boolean withDeviations,
            boolean withId, boolean withEndFrame, boolean withPrecision)
    {
        final ArrayPeakResultStore store = new ArrayPeakResultStore(10);
        while (size-- > 0)
        {
            final float[] params = createParams(n, r);
            final float[] paramsDev = (withDeviations) ? createParams(n, r) : null;
            final AttributePeakResult p = new AttributePeakResult(r.nextInt(), r.nextInt(), r.nextInt(), r.nextFloat(),
                    r.nextDouble(), r.nextFloat(), r.nextFloat(), params, paramsDev);
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

    private static float[] createParams(int n, UniformRandomProvider r)
    {
        final float[] p = new float[n];
        while (n-- > 0)
            p[n] = r.nextFloat();
        return p;
    }
}
