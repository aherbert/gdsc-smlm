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
package uk.ac.sussex.gdsc.smlm.model.camera;

import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;

import java.awt.Rectangle;
import java.util.Arrays;

/**
 * Base class for the camera model.
 *
 * @author Alex Herbert
 */
public abstract class BaseCameraModel implements CameraModel, Cloneable {
  /**
   * Check bias is finite.
   *
   * @param bias the bias
   */
  public void checkBias(float bias) {
    if (!Double.isFinite(bias)) {
      throw new IllegalArgumentException("Bias must be a finite number");
    }
  }

  /**
   * Check gain is strictly positive.
   *
   * @param gain the gain
   */
  public void checkGain(float gain) {
    if (!(gain <= Double.MAX_VALUE && gain > 0)) {
      throw new IllegalArgumentException("Gain must be strictly positive");
    }
  }

  /**
   * Check variance is positive.
   *
   * @param variance the variance
   */
  public void checkVariance(float variance) {
    if (!(variance <= Double.MAX_VALUE && variance >= 0)) {
      throw new IllegalArgumentException("Variance must be positive");
    }
  }

  /**
   * Create a new array.
   *
   * @param bounds the bounds
   * @param value the value
   * @return the float[]
   */
  protected static float[] newArray(Rectangle bounds, float value) {
    if (bounds == null || bounds.width <= 0 || bounds.height <= 0) {
      return new float[0];
    }
    final float[] data = new float[bounds.width * bounds.height];
    Arrays.fill(data, value);
    return data;
  }

  /**
   * Convert the variance to weights (1/variance). Any value of the variance that is not strictly
   * positive is set to the minimum variance above zero.
   *
   * @param variance the variance
   * @return the weights
   */
  public static float[] toWeights(float[] variance) {
    final float[] w = SimpleArrayUtils.ensureStrictlyPositive(variance);
    // If all the weights are zero then the first item will be zero
    // and there are no weights
    if (w[0] == 0) {
      Arrays.fill(w, 1f);
      return w;
    }
    for (int i = 0; i < w.length; i++) {
      w[i] = (float) (1.0 / w[i]);
    }
    return w;
  }
}
