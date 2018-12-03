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
package uk.ac.sussex.gdsc.smlm.ij.plugins;

import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 * Set the display values of an image to render it as a binary mask, all non-zero values are white.
 */
public class BinaryDisplay implements PlugInFilter {
  private ImagePlus imp;

  /** {@inheritDoc} */
  @Override
  public int setup(String arg, ImagePlus imp) {
    SMLMUsageTracker.recordPlugin(this.getClass(), arg);

    if (imp == null) {
      return DONE;
    }

    if (arg.equals("reset")) {
      final ImageProcessor ip = imp.getProcessor();
      ip.reset();
      imp.setProcessor(ip);
      imp.resetDisplayRange();
      imp.updateAndDraw();
      return DONE;
    }

    this.imp = imp;
    return DOES_ALL;
  }

  /** {@inheritDoc} */
  @Override
  public void run(ImageProcessor ip) {
    // float min = Float.POSITIVE_INFINITY;
    // for (int i=0; i<ip.getPixelCount(); i++)
    // {
    // final float value = ip.getf(i);
    // if (value == 0)
    // continue;
    // if (value < min)
    // min = value;
    // }
    // ip.setMinAndMax(0, min);
    // imp.updateAndDraw();

    final FloatProcessor fp = new FloatProcessor(ip.getWidth(), ip.getHeight());
    final float[] data = (float[]) fp.getPixels();
    for (int i = 0; i < ip.getPixelCount(); i++) {
      final float value = ip.getf(i);
      if (value == 0) {
        continue;
      }
      data[i] = 1;
    }

    ip.snapshot();
    ip.setPixels(0, fp);
    ip.setMinAndMax(0, 1);
    imp.updateAndDraw();
  }
}
