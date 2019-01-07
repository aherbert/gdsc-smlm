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

package uk.ac.sussex.gdsc.smlm.model;

/**
 * Contains a model for an image of blinking fluorophores under constant activation illumination.
 *
 *
 * <p>Based on the work of Coltharp et al (2012) Accurate Construction of photoactivated
 * localization microscopy images for quantitative measurements. PLOS One 7, Issue 12, pp 1-15
 */
public class ActivationTimeImageModel extends ImageModel {
  private double activationTime;

  /**
   * Instantiates a new activation time image model.
   *
   * @param activationTime Average time for activation
   * @param onTime Average on-state time
   * @param offTime1 Average off-state time for the first dark state
   * @param offTime2 Average off-state time for the second dark state
   * @param blinks1 Average number of blinks in the first dark state (used for each burst between
   *        second dark states)
   * @param blinks2 Average number of blinks into the second dark state
   */
  public ActivationTimeImageModel(double activationTime, double onTime, double offTime1,
      double offTime2, double blinks1, double blinks2) {
    super(onTime, offTime1, offTime2, blinks1, blinks2);
    init(activationTime);
  }

  private void init(double activationTime) {
    checkParameter("activationTime", activationTime);
    this.activationTime = activationTime;
  }

  /**
   * Gets the average time for activation.
   *
   * @return the activationTime.
   */
  public double getActivationTime() {
    return activationTime;
  }

  @Override
  protected double createActivationTime(double[] xyz) {
    return getRandom().nextExponential(activationTime);
  }

  @Override
  protected FluorophoreSequenceModel createFluorophore(int id, double[] xyz,
      double activationTime) {
    return new StandardFluorophoreSequenceModel(id, xyz, activationTime, onTime, offTime1, offTime2,
        blinks1, blinks2, isUseGeometricDistribution(), getRandom());
  }
}
