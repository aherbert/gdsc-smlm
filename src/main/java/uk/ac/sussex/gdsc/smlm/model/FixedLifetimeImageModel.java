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

import org.apache.commons.rng.UniformRandomProvider;

/**
 * Contains a model for an image of fixed lifetime fluorophores. All fluorphores will have the same
 * on time. The activation time will be the incremented by the on-time plus the dark time between
 * fluorophores.
 */
public class FixedLifetimeImageModel extends ImageModel {
  private double next;

  /**
   * Construct a new image model.
   *
   * @param onTime Fixed on-state time
   * @param offTime1 Dark time between successive fluorophores
   * @param rng the random generator for creating the image
   */
  public FixedLifetimeImageModel(double onTime, double offTime1, UniformRandomProvider rng) {
    super(onTime, offTime1, 0, 0, 0, rng);
  }

  @Override
  protected double createActivationTime(double[] xyz) {
    final double timeAct = next + getRandom().nextDouble();
    // @formatter:off
    // Ensure at least offTime1 full dark frames between lifetimes:
    // Frames:    |      |      |      |
    //         ------|
    //              end             |--------
    //                              start
    //                   offTime1
    // @formatter:on
    final int unit = (int) Math.ceil(offTime1);
    final int endT = (int) Math.ceil((timeAct + onTime) / unit);
    next = (endT + 1) * unit;
    return timeAct;
  }

  @Override
  protected FluorophoreSequenceModel createFluorophore(int id, double[] xyz, double timeAct) {
    return new SimpleFluorophoreSequenceModel(id, xyz, timeAct, onTime);
  }
}
