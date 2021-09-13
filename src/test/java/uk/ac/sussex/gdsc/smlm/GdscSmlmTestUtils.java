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

package uk.ac.sussex.gdsc.smlm;

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.DiscreteSampler;
import org.apache.commons.rng.sampling.distribution.PoissonSamplerCache;

/**
 * Provide testing utilities.
 */
public class GdscSmlmTestUtils {
  private static final PoissonSamplerCache POISSON_CACHE = new PoissonSamplerCache(0, 1000);

  /**
   * Creates the poisson sampler.
   *
   * @param rng the rng
   * @param mean the mean
   * @return the continuous sampler
   */
  public static DiscreteSampler createPoissonSampler(UniformRandomProvider rng, double mean) {
    return POISSON_CACHE.createSharedStateSampler(rng, mean);
  }

  /**
   * Generate an array of double values in the interval {@code [0, 1)}.
   *
   * @param n the number of generate
   * @param rng the random generator
   * @return the values
   */
  public static double[] generateDoubles(int n, UniformRandomProvider rng) {
    final double[] data = new double[n];
    for (int i = 0; i < n; i++) {
      data[i] = rng.nextDouble();
    }
    return data;
  }
}
