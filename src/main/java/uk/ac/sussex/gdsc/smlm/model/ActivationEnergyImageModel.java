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

package uk.ac.sussex.gdsc.smlm.model;

import org.apache.commons.rng.UniformRandomProvider;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;

/**
 * Contains a model for an image of blinking fluorophores under pulsed activation illumination.
 * Activation energy is sampled from an exponential distribution. Fluorophores are created when the
 * activation energy has been achieved under the given illumination.
 *
 * <p>Based on the work of Coltharp et al (2012) Accurate Construction of photoactivated
 * localization microscopy images for quantitative measurements. PLOS One 7, Issue 12, pp 1-15
 */
public class ActivationEnergyImageModel extends ImageModel {
  private double activationEnergy;
  private SpatialIllumination illumination;

  /**
   * Construct a new image model.
   *
   * @param activationEnergy Average energy for activation
   * @param illumination The illumination model
   * @param onTime Average on-state time
   * @param offTime1 Average off-state time for the first dark state
   * @param offTime2 Average off-state time for the second dark state
   * @param blinks1 Average number of blinks int the first dark state (used for each burst between
   *        second dark states)
   * @param blinks2 Average number of blinks into the second dark state
   * @param rng the random generator for creating the image
   */
  public ActivationEnergyImageModel(double activationEnergy, SpatialIllumination illumination,
      double onTime, double offTime1, double offTime2, double blinks1, double blinks2,
      UniformRandomProvider rng) {
    super(onTime, offTime1, offTime2, blinks1, blinks2, rng);
    init(activationEnergy, illumination);
  }

  private void init(double activationEnergy, SpatialIllumination illumination) {
    checkParameter("activationEnergy", activationEnergy);
    if (illumination == null) {
      throw new IllegalArgumentException("SpatialIllumination is null");
    }
    this.activationEnergy = activationEnergy;
    this.illumination = illumination;
  }

  /**
   * Gets the activation energy.
   *
   * @return the average energy for activation.
   */
  public double getActivationEnergy() {
    return activationEnergy;
  }

  @Override
  protected double createActivationTime(double[] xyz) {
    return getimeActivationTime(xyz, frameLimit);
  }

  @Override
  protected FluorophoreSequenceModel createFluorophore(int id, double[] xyz, double timeAct) {
    return new StandardFluorophoreSequenceModel(id, xyz, timeAct, onTime, offTime1, offTime2,
        blinks1, blinks2, isUseGeometricDistribution(), getRandom());
  }

  private double getimeActivationTime(double[] xyz, int frames) {
    final double activation =
        SamplerUtils.createExponentialSampler(getRandom(), activationEnergy).sample();
    double energy = 0;
    for (int t = 0; t < frames; t++) {
      // Q. Should the molecule be moving during the activation phase?
      final double[] photons = illumination.getPulsedPhotons(xyz, t + 1);

      energy += photons[0]; // pulse energy
      if (energy > activation) {
        return t;
      }

      energy += photons[1]; // during energy
      if (energy > activation) {
        // Interpolate
        return t + 1 - (energy - activation) / photons[1];
      }
    }
    return frames; // default to the number of frames.
  }
}
