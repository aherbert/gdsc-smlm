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

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import java.util.concurrent.atomic.AtomicReference;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;

/**
 * This plugin creates a mask image stack using an XY and XZ mask image.
 */
public class DepthMask implements PlugIn {
  /** The on-value for a byte mask. */
  private static final byte ON = (byte) 255;
  private static final String TITLE = "Depth Mask";

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());
    String titleXy;
    String titleXz;
    String titleYz;

    Settings() {
      // Set defaults
      titleXy = "";
      titleXz = "";
      titleYz = "";
    }

    Settings(Settings source) {
      titleXy = source.titleXy;
      titleXz = source.titleXz;
      titleYz = source.titleYz;
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
      return INSTANCE.get().copy();
    }

    /**
     * Save the settings.
     */
    void save() {
      INSTANCE.set(this);
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
    gd.addHelp(HelpUrls.getUrl("depth-mask"));

    gd.addMessage("Create a mask stack using XY, XZ and YZ mask images");

    final String[] maskList = ImageJUtils.getImageList(ImageJUtils.SINGLE);
    gd.addChoice("Mask_XY", maskList, settings.titleXy);
    gd.addChoice("Mask_XZ", maskList, settings.titleXz);
    gd.addChoice("Mask_YZ", maskList, settings.titleYz);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.titleXy = gd.getNextChoice();
    settings.titleXz = gd.getNextChoice();
    settings.titleYz = gd.getNextChoice();
    settings.save();

    return true;
  }

  private void createMask() {
    final ImagePlus impXy = WindowManager.getImage(settings.titleXy);
    if (impXy == null) {
      IJ.error(TITLE, "No XY mask");
      return;
    }
    final ImagePlus impXz = WindowManager.getImage(settings.titleXz);
    if (impXz == null) {
      IJ.error(TITLE, "No XZ mask");
      return;
    }
    final ImagePlus impYz = WindowManager.getImage(settings.titleYz);
    if (impYz == null) {
      IJ.error(TITLE, "No YZ mask");
      return;
    }
    if (impXy.getWidth() != impXz.getWidth()) {
      IJ.error(TITLE, "XY mask width does not match XZ mask width");
      return;
    }
    if (impXy.getHeight() != impYz.getWidth()) {
      IJ.error(TITLE, "XY mask height does not match YZ mask width");
      return;
    }
    if (impXz.getHeight() != impYz.getHeight()) {
      IJ.error(TITLE, "XZ mask height does not match YZ mask height");
      return;
    }

    final int maxx = impXy.getWidth();
    final int maxy = impXy.getHeight();
    final int maxz = impXz.getHeight();
    final ImageStack stack = new ImageStack(maxx, maxy, maxz);
    final byte[] maskXy = getMask(impXy);
    final byte[] maskXz = getMask(impXz);
    final byte[] maskYz = getMask(impYz);
    for (int z = 0; z < maxz; z++) {
      final byte[] mask = maskXy.clone();

      //// Simple method
      // for (int y = 0, i = 0; y < maxy; y++, i++)
      // for (int x = 0; x < maxx; x++, i++)
      // {
      // if (maskXZ[z * maxx + x] == 0)
      // mask[i] = 0;
      // else if (maskYZ[z * maxy + y] == 0)
      // mask[i] = 0;
      // }

      for (int x = 0, i = maxx * z; x < maxx; x++, i++) {
        if (maskXz[i] == 0) {
          // Blank all the (x,y) for this X
          for (int y = 0, xy = x; y < maxy; y++, xy += maxx) {
            mask[xy] = 0;
          }
        }
      }

      for (int y = 0, i = maxy * z; y < maxy; y++, i++) {
        if (maskYz[i] == 0) {
          // Blank all the (x,y) for this Y
          for (int x = 0, xy = y * maxx; x < maxx; x++, xy++) {
            mask[xy] = 0;
          }
        }
      }

      stack.setPixels(mask, z + 1);
    }
    ImageJUtils.display(TITLE, stack);
  }

  private static byte[] getMask(ImagePlus impXy) {
    final byte[] mask = (byte[]) impXy.getProcessor().convertToByte(false).getPixels();
    // Make binary
    for (int i = 0; i < mask.length; i++) {
      if (mask[i] != 0) {
        mask[i] = ON;
      }
    }
    return mask;
  }
}
