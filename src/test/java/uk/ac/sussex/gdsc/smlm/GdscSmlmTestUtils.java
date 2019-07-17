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
    return POISSON_CACHE.createPoissonSampler(rng, mean);
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
