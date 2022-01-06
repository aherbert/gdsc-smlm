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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 * Set the display values of an image to render it as a binary mask, all values above zero are
 * white.
 *
 * <p>This does not support 32-bit images with negative values.</p>
 */
public class BinaryDisplay implements PlugInFilter {
  private ImagePlus imp;

  @Override
  public int setup(String arg, ImagePlus imp) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (imp == null) {
      return DONE;
    }

    if ("reset".equals(arg)) {
      final ImageProcessor ip = imp.getProcessor();
      ip.reset();
      if (ip instanceof ByteProcessor) {
        // Reset to the entire range
        ip.resetMinAndMax();
      } else {
        // Short and FloatProcessor store the min and max in the snapshot
        // so restore it from the reset values.
        ip.setMinAndMax(ip.getMin(), ip.getMax());
      }
      imp.updateAndDraw();
      return DONE;
    }

    this.imp = imp;
    return DOES_8G | DOES_16 | DOES_32;
  }

  @Override
  public void run(ImageProcessor ip) {
    double max;
    if (ip instanceof ByteProcessor) {
      max = 1;
    } else {
      // Short and FloatProcessor store the current min and max in the snapshot
      ip.snapshot();
      if (ip instanceof FloatProcessor) {
        // Compute the lowest value above zero
        float min = Float.POSITIVE_INFINITY;
        for (int i = 0; i < ip.getPixelCount(); i++) {
          final float value = ip.getf(i);
          if (value > 0 && value < min) {
            min = value;
          }
        }
        max = min;
      } else {
        max = 1;
      }
    }

    ip.setMinAndMax(0, max);
    imp.updateAndDraw();
  }
}
