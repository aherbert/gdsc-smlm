/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.model;

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.ContinuousSampler;
import uk.ac.sussex.gdsc.core.utils.rng.PoissonSamplers;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;

/**
 * Contains a continuous-time model for a blinking fluorophore. Assumes a constant activation laser
 * and a simple exponential model for the average activation time.
 *
 * <p>Based on the work of Coltharp et al (2012) Accurate Construction of photoactivated
 * localization microscopy images for quantitative measurements. PLOS One 7, Issue 12, pp 1-15.
 */
public class StandardFluorophoreSequenceModel extends FluorophoreSequenceModel {
  /**
   * Construct a new flourophore.
   *
   * @param averageActivationTime Average time for activation
   * @param id The identifier
   * @param xyz The [x,y,z] coordinates
   * @param onTime Average on-state time
   * @param offTime Average off-state time for the first dark state
   * @param offTime2 Average off-state time for the second dark state
   * @param blinks1 Average number of blinks int the first dark state (used for each burst between
   *        second dark states)
   * @param blinks2 Average number of blinks into the second dark state
   * @param useGeometricBlinkingDistribution Set to true to use the geometric distribution (default
   *        is Poisson)
   * @param rng the random generator
   */
  public StandardFluorophoreSequenceModel(double averageActivationTime, int id, double[] xyz,
      double onTime, double offTime, double offTime2, double blinks1, double blinks2,
      boolean useGeometricBlinkingDistribution, UniformRandomProvider rng) {
    super(id, xyz);
    init(SamplerUtils.createExponentialSampler(rng, averageActivationTime).sample(), onTime,
        offTime, offTime2, blinks1, blinks2, useGeometricBlinkingDistribution, rng);
  }

  /**
   * Construct a new flourophore.
   *
   * @param id The identifier
   * @param xyz The [x,y,z] coordinates
   * @param startT The activation time
   * @param onTime Average on-state time
   * @param offTime Average off-state time for the first dark state
   * @param offTime2 Average off-state time for the second dark state
   * @param blinks1 Average number of blinks int the first dark state (used for each burst between
   *        second dark states)
   * @param blinks2 Average number of blinks into the second dark state
   * @param useGeometricBlinkingDistribution Set to true to use the geometric distribution (default
   *        is Poisson)
   * @param rng the random generator
   */
  public StandardFluorophoreSequenceModel(int id, double[] xyz, double startT, double onTime,
      double offTime, double offTime2, double blinks1, double blinks2,
      boolean useGeometricBlinkingDistribution, UniformRandomProvider rng) {
    super(id, xyz);
    init(startT, onTime, offTime, offTime2, blinks1, blinks2, useGeometricBlinkingDistribution,
        rng);
  }

  private void init(double startT, double onTime, double offTime, double offTime2, double blinks1,
      double blinks2, boolean useGeometricBlinkingDistribution, UniformRandomProvider rand) {
    // Model two dark states: short and long. The second offTime and blinks1 is for the long dark
    // state:
    //
    // ++-+-+++-+.................+-+--++-+................................+--+++-+
    //
    // + = on
    // - = Short dark state
    // . = Long dark state

    // Note: 1+blinks1 is the number of on-states

    final DoubleArrayList sequence = new DoubleArrayList();

    // The exponential distribution is just scaled by the mean
    final ContinuousSampler sampler = SamplerUtils.createExponentialSampler(rand);

    // Perform a set number of long blinks
    final int nLongBlinks = getBlinks(useGeometricBlinkingDistribution, rand, blinks2);
    double time = startT;
    for (int n = 0; n <= nLongBlinks; n++) {
      // For each burst between long blinks perform a number of short blinks
      final int nShortBlinks = getBlinks(useGeometricBlinkingDistribution, rand, blinks1);

      // Starts on the current time
      sequence.add(time);
      // Stops after the on-time
      time += sampler.sample() * onTime;
      sequence.add(time);

      // Remaining bursts
      for (int i = 0; i < nShortBlinks; i++) {
        // Next burst starts after the short off-time
        time += sampler.sample() * offTime;
        sequence.add(time);
        // Stops after the on-time
        time += sampler.sample() * onTime;
        sequence.add(time);
      }

      // Add the long dark state if there are more bursts.
      time += sampler.sample() * offTime2;
    }

    // Convert the sequence to the burst sequence array
    setBurstSequence(sequence.toDoubleArray());
  }

  /**
   * Get the number of blinks using the specified random data generator using a Poisson or Geometric
   * distribution.
   *
   * @param useGeometricBlinkingDistribution Set to true to use the geometric distribution (default
   *        is Poisson)
   * @param rand the random generator
   * @param mean the mean
   * @return The number of blinks
   */
  public static int getBlinks(boolean useGeometricBlinkingDistribution, UniformRandomProvider rand,
      double mean) {
    if (mean > 0) {
      return (useGeometricBlinkingDistribution)
          ? SamplerUtils.createGeometricSamplerFromMean(rand, mean).sample()
          : PoissonSamplers.nextPoissonSample(rand, mean);
    }
    return 0;
  }
}
