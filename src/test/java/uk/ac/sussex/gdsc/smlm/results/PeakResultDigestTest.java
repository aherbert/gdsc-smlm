/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

import java.util.logging.Logger;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;

@SuppressWarnings({"javadoc"})
class PeakResultDigestTest {
  @SeededTest
  void sameResultsAreEqual(RandomSeed seed) {
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    final PeakResult[] r1 = createResults(r, 10, 5, false, false, false, false);
    final PeakResultsDigest digest = new PeakResultsDigest(r1);
    Assertions.assertTrue(digest.matches(r1));
    Assertions.assertTrue(digest.matches(digest));
  }

  @SeededTest
  void sameSize1ResultsAreEqual(RandomSeed seed) {
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    final PeakResult[] r1 = createResults(r, 1, 5, false, false, false, false);
    final PeakResultsDigest digest = new PeakResultsDigest(r1);
    Assertions.assertTrue(digest.matches(r1));
    Assertions.assertTrue(digest.matches(digest));
  }

  @SeededTest
  void sameEmptyResultsAreEqual(RandomSeed seed) {
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    final PeakResult[] r1 = createResults(r, 0, 5, false, false, false, false);
    final PeakResultsDigest digest = new PeakResultsDigest(r1);
    Assertions.assertTrue(digest.matches(r1));
    Assertions.assertTrue(digest.matches(digest));
  }

  @SeededTest
  void sameResultsAreEqualWithDeviation(RandomSeed seed) {
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    final PeakResult[] r1 = createResults(r, 10, 5, true, false, false, false);
    final PeakResultsDigest digest = new PeakResultsDigest(r1);
    Assertions.assertTrue(digest.matches(r1));
    Assertions.assertTrue(digest.matches(digest));
  }

  @SeededTest
  void sameResultsAreEqualWithId(RandomSeed seed) {
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    final PeakResult[] r1 = createResults(r, 10, 5, false, true, false, false);
    final PeakResultsDigest digest = new PeakResultsDigest(r1);
    Assertions.assertTrue(digest.matches(r1));
    Assertions.assertTrue(digest.matches(digest));
  }

  @SeededTest
  void sameResultsAreEqualWithEndFrame(RandomSeed seed) {
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    final PeakResult[] r1 = createResults(r, 10, 5, false, false, true, false);
    final PeakResultsDigest digest = new PeakResultsDigest(r1);
    Assertions.assertTrue(digest.matches(r1));
    Assertions.assertTrue(digest.matches(digest));
  }

  @SeededTest
  void sameResultsAreEqualWithPrecision(RandomSeed seed) {
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    final PeakResult[] r1 = createResults(r, 10, 5, false, false, false, true);
    final PeakResultsDigest digest = new PeakResultsDigest(r1);
    Assertions.assertTrue(digest.matches(r1));
    Assertions.assertTrue(digest.matches(digest));
  }

  @SeededTest
  void differentResultsAreNotEqual(RandomSeed seed) {
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    final PeakResult[] r1 = createResults(r, 10, 5, false, false, false, false);
    final PeakResultsDigest digest = new PeakResultsDigest(r1);
    for (final int size : new int[] {10, 1, 0}) {
      final PeakResult[] r2 = createResults(r, size, 5, false, false, false, false);
      Assertions.assertFalse(digest.matches(r2));
      Assertions.assertFalse(digest.matches(new PeakResultsDigest(r2)));
    }
  }

  @SeededTest
  void digestMatchesPeakResultDigest(RandomSeed seed) {
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    for (int size = 1; size < 5; size++) {
      final PeakResult[] r1 = createResults(r, size, 5, false, false, false, false);
      final PeakResultsDigest digest = new PeakResultsDigest(r1);

      final PeakResultDigest d = new PeakResultDigest();
      for (final PeakResult rr : r1) {
        d.update(rr);
      }
      Assertions.assertEquals(d.digest(), digest.getDigest());
    }
  }

  @SeededTest
  void digestIsEmptyStringWhenSizeIsZero() {
    Assertions.assertEquals("", new PeakResultsDigest(new PeakResult[0]).getDigest());
  }

  @SeededTest
  void digestHandlesNull() {
    final PeakResult[] r1 = null;
    final PeakResult[] r0 = new PeakResult[0];
    final PeakResultsDigest digest = new PeakResultsDigest(r1);
    Assertions.assertTrue(digest.matches(r1));
    Assertions.assertFalse(digest.matches(r0));
  }

  @SeededTest
  void digestHandlesEmptyArray() {
    final PeakResult[] r1 = null;
    final PeakResult[] r0 = new PeakResult[0];
    final PeakResultsDigest digest = new PeakResultsDigest(r0);
    Assertions.assertTrue(digest.matches(r0));
    Assertions.assertFalse(digest.matches(r1));
  }

  @SeededTest
  void timeDigest(RandomSeed seed) {
    Assumptions.assumeTrue(false);
    final Logger logger = Logger.getLogger(PeakResultDigestTest.class.getName());

    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    final PeakResultsDigest digest = new PeakResultsDigest();
    final int N = 5;
    for (int size = 1000; size < 2000000; size *= 2) {
      final PeakResult[] r1 = createResults(r, size, 5, false, false, false, false);
      long time = System.nanoTime();
      for (int i = N; i-- > 0;) {
        digest.digest(r1);
      }
      time = System.nanoTime() - time;
      logger.info(FunctionUtils.getSupplier("size = %d, time = %g ms", size, (1e-6 * time) / N));
    }
  }

  private static PeakResult[] createResults(UniformRandomProvider rng, int size, int np,
      boolean withDeviations, boolean withId, boolean withEndFrame, boolean withPrecision) {
    final ArrayPeakResultStore store = new ArrayPeakResultStore(10);
    while (size-- > 0) {
      final float[] params = createParams(np, rng);
      final float[] paramsDev = (withDeviations) ? createParams(np, rng) : null;
      final AttributePeakResult p =
          new AttributePeakResult(rng.nextInt(), rng.nextInt(), rng.nextInt(), rng.nextFloat(),
              rng.nextDouble(), rng.nextFloat(), rng.nextFloat(), params, paramsDev);
      if (withId) {
        p.setId(rng.nextInt());
      }
      if (withEndFrame) {
        // p.setEndFrame(p.getFrame() + 1 + rng.nextInt(5));
        p.setEndFrame(rng.nextInt());
      }
      if (withPrecision) {
        p.setPrecision(rng.nextDouble());
      }
      store.add(p);
    }
    return store.toArray();
  }

  private static float[] createParams(int np, UniformRandomProvider rng) {
    final float[] params = new float[np];
    while (np-- > 0) {
      params[np] = rng.nextFloat();
    }
    return params;
  }
}
