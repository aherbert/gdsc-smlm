/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;

@SuppressWarnings({"javadoc"})
public class PeakResultTest {
  @SeededTest
  public void sameResultIsEqual(RandomSeed seed) {
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    final PeakResult[] r1 = createResults(r, 1, 5, false, false, false, false);
    final PeakResult p = r1[0];
    Assertions.assertTrue(PeakResult.equals(p, p), "Same object");
    Assertions.assertTrue(PeakResult.equals(null, null), "Null object");
  }

  @SeededTest
  public void sameResultsAreEqual(RandomSeed seed) {
    UniformRandomProvider rng;
    final int size = 10;
    final int n = 5;
    final boolean[] both = {true, false};
    for (final boolean withDeviations : both) {
      for (final boolean withId : both) {
        for (final boolean withEndFrame : both) {
          for (final boolean withPrecision : both) {
            rng = RngUtils.create(seed.getSeed());
            final PeakResult[] r1 =
                createResults(rng, size, n, withDeviations, withId, withEndFrame, withPrecision);
            rng = RngUtils.create(seed.getSeed());
            final PeakResult[] r2 =
                createResults(rng, size, n, withDeviations, withId, withEndFrame, withPrecision);
            for (int i = 0; i < r1.length; i++) {
              Assertions.assertTrue(PeakResult.equals(r1[i], r2[i]));
            }
          }
        }
      }
    }
  }

  @SeededTest
  public void differentResultsAreNotEqual(RandomSeed seed) {
    UniformRandomProvider rng;
    final int size = 1;
    final int n = 5;
    final boolean[] both = {true, false};
    for (final boolean withDeviations : both) {
      for (final boolean withId : both) {
        for (final boolean withEndFrame : both) {
          for (final boolean withPrecision : both) {
            rng = RngUtils.create(seed.getSeed());
            final PeakResult[] r1 =
                createResults(rng, size, n, withDeviations, withId, withEndFrame, withPrecision);
            final PeakResult[] r2 =
                createResults(rng, size, n, withDeviations, withId, withEndFrame, withPrecision);
            for (int i = 0; i < r1.length; i++) {
              Assertions.assertFalse(PeakResult.equals(r1[i], r2[i]));
            }
          }
        }
      }
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
    final float[] p = new float[np];
    while (np-- > 0) {
      p[np] = rng.nextFloat();
    }
    return p;
  }
}
