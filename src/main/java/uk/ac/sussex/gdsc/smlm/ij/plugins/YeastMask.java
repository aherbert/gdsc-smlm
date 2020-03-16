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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.process.Blitter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import java.util.concurrent.atomic.AtomicReference;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;

/**
 * This plugin creates a mask image of a Yeast cell for use in diffusion simulations.
 */
public class YeastMask implements PlugIn {
  private static final String TITLE = "Yeast Mask";

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    double length = 8;
    double radius = 1.5;
    double nucleus = 0.9;
    double nmPerPixel = 107;
    double nmPerSlice = 20;
    boolean excludeNucleus = true;
    boolean squareOutput = true;
    int border = 3;
    boolean is2D;

    Settings() {
      // Set defaults
      length = 8;
      radius = 1.5;
      nucleus = 0.9;
      nmPerPixel = 107;
      nmPerSlice = 20;
      excludeNucleus = true;
      squareOutput = true;
      border = 3;
    }

    Settings(Settings source) {
      length = source.length;
      radius = source.radius;
      nucleus = source.nucleus;
      nmPerPixel = source.nmPerPixel;
      nmPerSlice = source.nmPerSlice;
      excludeNucleus = source.excludeNucleus;
      squareOutput = source.squareOutput;
      border = source.border;
      is2D = source.is2D;
    }

    Settings copy() {
      return new Settings(this);
    }

    /**
     * Load a copy of the settings.
     *
     * @return the settings
     */
    static Settings load() {
      return lastSettings.get().copy();
    }

    /**
     * Save the settings.
     */
    void save() {
      lastSettings.set(this);
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (!showDialog()) {
      return;
    }

    createMask();
  }

  private boolean showDialog() {
    settings = Settings.load();

    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("yeast-mask"));

    gd.addMessage("Create a mask of a yeast cell as a tube plus end-caps");
    gd.addSlider("Tube_length (um)", 10, 20, settings.length);
    gd.addSlider("Radius (um)", 0.5, 5, settings.radius);
    gd.addCheckbox("Exclude_nucleus", settings.excludeNucleus);
    gd.addSlider("Nucleus (fraction)", 0.5, 1, settings.nucleus);
    gd.addNumericField("Pixel_pitch", settings.nmPerPixel, 1, 6, "nm");
    gd.addNumericField("Pixel_depth", settings.nmPerSlice, 1, 6, "nm");
    gd.addCheckbox("Square_output", settings.squareOutput);
    gd.addSlider("Border", 0, 10, settings.border);
    gd.addCheckbox("2D", settings.is2D);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.length = gd.getNextNumber();
    settings.radius = gd.getNextNumber();
    settings.excludeNucleus = gd.getNextBoolean();
    settings.nucleus = gd.getNextNumber();
    settings.nmPerPixel = gd.getNextNumber();
    settings.nmPerSlice = gd.getNextNumber();
    settings.squareOutput = gd.getNextBoolean();
    settings.border = (int) gd.getNextNumber();
    settings.is2D = gd.getNextBoolean();

    if (settings.radius < 0.5) {
      settings.radius = 0.5;
    }
    if (settings.length < 0) {
      settings.length = 0;
    }
    if (settings.nmPerPixel < 1) {
      settings.nmPerPixel = 1;
    }
    if (settings.nmPerSlice < 1) {
      settings.nmPerSlice = 1;
    }

    settings.save();

    return true;
  }

  private void createMask() {
    // Create the dimensions
    final int hw = (int) Math.ceil(settings.radius * 1000 / settings.nmPerPixel);
    final int hd = (int) Math.ceil(settings.radius * 1000 / settings.nmPerSlice);

    final int width = 2 * hw + 1;
    final int depth = 2 * hd + 1;

    ImageStack stack = createHemiSphere(width, depth);

    // Extend the centre circle of the sphere into a tube of the required length
    final int h = (int) Math.ceil(settings.length * 1000 / settings.nmPerPixel);
    if (h > 0) {
      final ImageStack newStack = new ImageStack(width, stack.getHeight() + h, stack.getSize());
      for (int slice = 1; slice <= stack.getSize(); slice++) {
        final byte[] pixels = (byte[]) stack.getPixels(slice);
        final byte[] newPixels = new byte[width * newStack.getHeight()];
        newStack.setPixels(newPixels, slice);
        System.arraycopy(pixels, 0, newPixels, 0, pixels.length);
        // Get the final strip to be extended
        final int offset = pixels.length - width;
        int target = pixels.length;
        for (int i = 0; i < h; i++) {
          System.arraycopy(pixels, offset, newPixels, target, width);
          target += width;
        }
      }
      stack = newStack;
    }

    // Copy the hemi-sphere onto the end
    final ImageStack newStack = new ImageStack(width, stack.getHeight() + hw, stack.getSize());
    for (int slice = 1; slice <= stack.getSize(); slice++) {
      final byte[] pixels = (byte[]) stack.getPixels(slice);
      final byte[] newPixels = new byte[width * newStack.getHeight()];
      newStack.setPixels(newPixels, slice);
      System.arraycopy(pixels, 0, newPixels, 0, pixels.length);
      // Copy the hemi-sphere
      int source = 0;
      int target = newPixels.length - width;
      for (int i = 0; i < hw; i++) {
        System.arraycopy(pixels, source, newPixels, target, width);
        target -= width;
        source += width;
      }
    }
    stack = newStack;

    if (settings.excludeNucleus) {
      final ImageStack stack2 = createNucleusSphere(width, depth);
      final int xloc = (stack.getWidth() - stack2.getWidth()) / 2;
      final int yloc = (stack.getHeight() - stack2.getHeight()) / 2;
      final int offset = (stack.getSize() - stack2.getSize()) / 2;
      for (int slice = 1; slice <= stack2.getSize(); slice++) {
        final ImageProcessor ip = stack.getProcessor(slice + offset);
        final ImageProcessor ip2 = stack2.getProcessor(slice);
        ip.copyBits(ip2, xloc, yloc, Blitter.SUBTRACT);
      }
    }

    if (settings.squareOutput && stack.getWidth() != stack.getHeight()) {
      final ImageStack stack2 = new ImageStack(stack.getHeight(), stack.getHeight());
      final int end = stack.getHeight() - stack.getWidth();
      for (int slice = 1; slice <= stack.getSize(); slice++) {
        final ImageProcessor ip = stack.getProcessor(slice);
        final ImageProcessor ip2 = new ByteProcessor(stack2.getWidth(), stack2.getHeight());
        stack2.addSlice(ip2);
        for (int xloc = 0; xloc <= end; xloc += stack.getWidth()) {
          ip2.insert(ip, xloc, 0);
        }
      }
      stack = stack2;
    }

    if (settings.border > 0) {
      final ImageStack stack2 = new ImageStack(stack.getWidth() + 2 * settings.border,
          stack.getHeight() + 2 * settings.border);
      for (int slice = 1; slice <= stack.getSize(); slice++) {
        final ImageProcessor ip = stack.getProcessor(slice);
        final ImageProcessor ip2 = new ByteProcessor(stack2.getWidth(), stack2.getHeight());
        stack2.addSlice(ip2);
        ip2.insert(ip, settings.border, settings.border);
      }
      stack = stack2;
    }

    ImagePlus imp;
    if (settings.is2D) {
      // TODO - Remove this laziness since we should really just do a 2D image
      final int centre = stack.getSize() / 2;
      imp = ImageJUtils.display(TITLE, stack.getProcessor(centre));
    } else {
      imp = ImageJUtils.display(TITLE, stack);
    }

    // Calibrate
    final Calibration cal = new Calibration();
    cal.setUnit("um");
    cal.pixelWidth = cal.pixelHeight = settings.nmPerPixel / 1000;
    cal.pixelDepth = settings.nmPerSlice / 1000;
    imp.setCalibration(cal);
  }

  /**
   * Create a sphere using the given pixel width and stack depth using a fraction of the original
   * cell radius.
   *
   * @param width the width
   * @param depth the depth
   * @return A sphere
   */
  private ImageStack createNucleusSphere(int width, int depth) {
    // Create a sphere. This could be done exploiting symmetry to be more efficient
    // but is left as a simple implementation
    final double centreX = width * 0.5;
    final double centreZ = depth * 0.5;

    // Precompute squares for the width
    final double[] s = new double[width];
    for (int iy = 0; iy < width; iy++) {
      final double y = (centreX - (iy + 0.5)) * settings.nmPerPixel;
      s[iy] = y * y;
    }

    final ImageStack stack = new ImageStack(width, width, depth);
    final byte on = (byte) 255;
    final double r = settings.radius * 1000 * settings.nucleus;
    final double r2 = r * r;
    for (int iz = 0; iz < depth; iz++) {
      final double z = (centreZ - (iz + 0.5)) * settings.nmPerSlice;
      final double z2 = z * z;
      final byte[] mask = new byte[width * width];
      for (int iy = 0, i = 0; iy < width; iy++) {
        final double y2z2 = s[iy] + z2;
        for (int ix = 0; ix < width; ix++, i++) {
          final double d2 = s[ix] + y2z2;
          if (d2 < r2) {
            mask[i] = on;
          }
        }
      }
      stack.setPixels(mask, iz + 1);
    }

    return stack;
  }

  /**
   * Create a hemi-sphere using the given pixel width and stack depth using the original cell
   * radius.
   *
   * @param width the width
   * @param depth the depth
   * @return A hemi-sphere
   */
  private ImageStack createHemiSphere(int width, int depth) {
    // Create a sphere. This could be done exploiting symmetry to be more efficient
    // but is left as a simple implementation
    final double centreX = width * 0.5;
    final double centreZ = depth * 0.5;

    // Precompute squares for the width
    final double[] s = new double[width];
    for (int iy = 0; iy < width; iy++) {
      final double y = (centreX - (iy + 0.5)) * settings.nmPerPixel;
      s[iy] = y * y;
    }

    final int halfHeight = 1 + width / 2;
    final ImageStack stack = new ImageStack(width, halfHeight, depth);
    final byte on = (byte) 255;
    final double r = settings.radius * 1000;
    final double r2 = r * r;
    for (int iz = 0; iz < depth; iz++) {
      final double z = (centreZ - (iz + 0.5)) * settings.nmPerSlice;
      final double z2 = z * z;
      final byte[] mask = new byte[width * halfHeight];
      for (int iy = 0, i = 0; iy < halfHeight; iy++) {
        final double y2z2 = s[iy] + z2;
        for (int ix = 0; ix < width; ix++, i++) {
          final double d2 = s[ix] + y2z2;
          if (d2 < r2) {
            mask[i] = on;
          }
        }
      }
      stack.setPixels(mask, iz + 1);
    }

    return stack;
  }
}
