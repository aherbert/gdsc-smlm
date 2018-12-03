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
package uk.ac.sussex.gdsc.smlm.model;

import org.apache.commons.math3.random.RandomDataGenerator;

import gnu.trove.list.array.TDoubleArrayList;

/**
 * Contains a continuous-time model for a blinking fluorophore. Assumes a constant activation laser
 * and a simple exponential model for the average activation time. <p> Based on the work of Coltharp
 * et al (2012) Accurate Construction of photoactivated localization microscopy images for
 * quantitative measurements. PLOS One 7, Issue 12, pp 1-15
 */
public class StandardFluorophoreSequenceModel extends FluorophoreSequenceModel {
  /**
   * Construct a new flourophore.
   *
   * @param tAct Average time for activation
   * @param id The identifier
   * @param xyz The [x,y,z] coordinates
   * @param tOn Average on-state time
   * @param tOff Average off-state time for the first dark state
   * @param tOff2 Average off-state time for the second dark state
   * @param nBlinks Average number of blinks int the first dark state (used for each burst between
   *        second dark states)
   * @param nBlinks2 Average number of blinks into the second dark state
   * @param useGeometricBlinkingDistribution Set to true to use the geometric distribution (default
   *        is Poisson)
   */
  public StandardFluorophoreSequenceModel(double tAct, int id, double[] xyz, double tOn,
      double tOff, double tOff2, double nBlinks, double nBlinks2,
      boolean useGeometricBlinkingDistribution) {
    this(tAct, id, xyz, tOn, tOff, tOff2, nBlinks, nBlinks2, useGeometricBlinkingDistribution,
        new RandomDataGenerator());
  }

  /**
   * Construct a new flourophore.
   *
   * @param tAct Average time for activation
   * @param id The identifier
   * @param xyz The [x,y,z] coordinates
   * @param tOn Average on-state time
   * @param tOff Average off-state time for the first dark state
   * @param tOff2 Average off-state time for the second dark state
   * @param nBlinks Average number of blinks int the first dark state (used for each burst between
   *        second dark states)
   * @param nBlinks2 Average number of blinks into the second dark state
   * @param useGeometricBlinkingDistribution Set to true to use the geometric distribution (default
   *        is Poisson)
   * @param randomGenerator the random generator
   */
  public StandardFluorophoreSequenceModel(double tAct, int id, double[] xyz, double tOn,
      double tOff, double tOff2, double nBlinks, double nBlinks2,
      boolean useGeometricBlinkingDistribution, RandomDataGenerator randomGenerator) {
    super(id, xyz);
    init(randomGenerator.nextExponential(tAct), tOn, tOff, tOff2, nBlinks, nBlinks2,
        useGeometricBlinkingDistribution, randomGenerator);
  }

  /**
   * Construct a new flourophore.
   *
   * @param id The identifier
   * @param xyz The [x,y,z] coordinates
   * @param startT The activation time
   * @param tOn Average on-state time
   * @param tOff Average off-state time for the first dark state
   * @param tOff2 Average off-state time for the second dark state
   * @param nBlinks Average number of blinks int the first dark state (used for each burst between
   *        second dark states)
   * @param nBlinks2 Average number of blinks into the second dark state
   * @param useGeometricBlinkingDistribution Set to true to use the geometric distribution (default
   *        is Poisson)
   */
  public StandardFluorophoreSequenceModel(int id, double[] xyz, double startT, double tOn,
      double tOff, double tOff2, double nBlinks, double nBlinks2,
      boolean useGeometricBlinkingDistribution) {
    super(id, xyz);
    init(startT, tOn, tOff, tOff2, nBlinks, nBlinks2, useGeometricBlinkingDistribution,
        new RandomDataGenerator());
  }

  /**
   * Construct a new flourophore.
   *
   * @param id The identifier
   * @param xyz The [x,y,z] coordinates
   * @param startT The activation time
   * @param tOn Average on-state time
   * @param tOff Average off-state time for the first dark state
   * @param tOff2 Average off-state time for the second dark state
   * @param nBlinks Average number of blinks int the first dark state (used for each burst between
   *        second dark states)
   * @param nBlinks2 Average number of blinks into the second dark state
   * @param useGeometricBlinkingDistribution Set to true to use the geometric distribution (default
   *        is Poisson)
   * @param randomGenerator the random generator
   */
  public StandardFluorophoreSequenceModel(int id, double[] xyz, double startT, double tOn,
      double tOff, double tOff2, double nBlinks, double nBlinks2,
      boolean useGeometricBlinkingDistribution, RandomDataGenerator randomGenerator) {
    super(id, xyz);
    init(startT, tOn, tOff, tOff2, nBlinks, nBlinks2, useGeometricBlinkingDistribution,
        randomGenerator);
  }

  private void init(double t, double tOn, double tOff, double tOff2, double nBlinks,
      double nBlinks2, boolean useGeometricBlinkingDistribution, RandomDataGenerator rand) {
    // Model two dark states: short and long. The second tOff and nBlinks is for the long dark
    // state:
    //
    // ++-+-+++-+.................+-+--++-+................................+--+++-+
    //
    // + = on
    // - = Short dark state
    // . = Long dark state

    // Note: 1+nBlinks is the number of on-states

    final TDoubleArrayList sequence = new TDoubleArrayList();

    // Perform a set number of long blinks
    final int nLongBlinks = getBlinks(useGeometricBlinkingDistribution, rand, nBlinks2);
    for (int n = 0; n <= nLongBlinks; n++) {
      // For each burst between long blinks perform a number of short blinks
      final int nShortBlinks = getBlinks(useGeometricBlinkingDistribution, rand, nBlinks);

      // Starts on the current time
      sequence.add(t);
      // Stops after the on-time
      t += rand.nextExponential(tOn);
      sequence.add(t);

      // Remaining bursts
      for (int i = 0; i < nShortBlinks; i++) {
        // Next burst starts after the short off-time
        t += rand.nextExponential(tOff);
        sequence.add(t);
        // Stops after the on-time
        t += rand.nextExponential(tOn);
        sequence.add(t);
      }

      // Add the long dark state if there are more bursts.
      t += rand.nextExponential(tOff2);
    }

    // Convert the sequence to the burst sequence array
    setBurstSequence(sequence.toArray());
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
  public static int getBlinks(boolean useGeometricBlinkingDistribution, RandomDataGenerator rand,
      double mean) {
    if (mean > 0) {
      return (useGeometricBlinkingDistribution) ? nextGeometric(rand, mean)
          : (int) rand.nextPoisson(mean);
    }
    return 0;
  }

  private static int nextGeometric(RandomDataGenerator rand, double mean) {
    // Use a geometric distribution by sampling the floor from the exponential.
    // Geometric distribution where k { 0, 1, 2, ... }
    // See: http://en.wikipedia.org/wiki/Geometric_distribution#Related_distributions
    final double p = 1 / (1 + mean);
    return (int) Math.floor(Math.log(rand.nextUniform(0, 1, true)) / Math.log(1 - p));
  }
}
