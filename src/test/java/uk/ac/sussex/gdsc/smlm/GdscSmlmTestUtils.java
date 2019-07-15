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
}
