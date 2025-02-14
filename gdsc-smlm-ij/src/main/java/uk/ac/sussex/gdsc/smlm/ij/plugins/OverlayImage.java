/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Package
 *
 * Software for single molecule localisation microscopy (SMLM) in ImageJ
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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
import ij.Undo;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.ImageRoi;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import java.awt.Rectangle;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Function;

/**
 * This plugin is extracted from ij.plugins.OverlayCommands to allow an image to be added with a
 * transparent background. An option to overlay a stack has been added.
 */
public class OverlayImage implements PlugIn {
  private static final String TITLE = "Overlay image";

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());

    int opacity;
    boolean transparent;
    boolean replace;
    boolean stack;
    String title;

    Settings() {
      // Set defaults
      opacity = 100;
      transparent = true;
      replace = true;
      stack = true;
      title = "";
    }

    Settings(Settings source) {
      opacity = source.opacity;
      transparent = source.transparent;
      replace = source.replace;
      stack = source.stack;
      title = source.title;
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

    addImage();
  }

  /**
   * Adapted from ij.plugins.OverlayCommands#addImage(boolean) with: the additional option for
   * setting the zero pixels to transparent; support for overlay of a stack.
   */
  void addImage() {
    final int[] wList = WindowManager.getIDList();
    if (wList == null || wList.length < 2) {
      IJ.error(TITLE, "The command requires at least two open images.");
      return;
    }
    final ImagePlus imp = IJ.getImage();
    String[] titles = new String[wList.length];
    int count = 0;
    for (final int id : wList) {
      final ImagePlus imp2 = WindowManager.getImage(id);
      if (imp2 != null && imp2.getID() != imp.getID() && imp.getWidth() >= imp2.getWidth()
          && imp.getHeight() >= imp2.getHeight()) {
        titles[count++] = imp2.getTitle();
      }
    }
    if (count < 1) {
      IJ.error(TITLE, "The command requires at least one valid overlay image.");
      return;
    }
    titles = Arrays.copyOf(titles, count);

    int x = 0;
    int y = 0;
    final Roi roi = imp.getRoi();
    if (roi != null && roi.isArea()) {
      final Rectangle r = roi.getBounds();
      x = r.x;
      y = r.y;
    }

    final Settings settings = Settings.load();
    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addChoice("Image to add:", titles, settings.title);
    gd.addNumericField("X location:", x, 0);
    gd.addNumericField("Y location:", y, 0);
    gd.addNumericField("Opacity (0-100%):", settings.opacity, 0);
    gd.addCheckbox("Transparent background", settings.transparent);
    gd.addCheckbox("Replace overlay", settings.replace);
    gd.addCheckbox("Use stack", settings.stack);

    gd.addHelp(HelpUrls.getUrl("overlay-image"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }

    settings.title = gd.getNextChoice();
    int xx = (int) gd.getNextNumber();
    int yy = (int) gd.getNextNumber();
    settings.opacity = (int) gd.getNextNumber();
    settings.transparent = gd.getNextBoolean();
    settings.replace = gd.getNextBoolean();
    settings.stack = gd.getNextBoolean();
    settings.save();

    final ImagePlus imp2 = WindowManager.getImage(settings.title);
    if (imp2.getID() == imp.getID()) {
      IJ.error(TITLE, "Image to be added cannot be the same as\n\"" + imp.getTitle() + "\".");
      return;
    }
    if (imp2.getWidth() > imp.getWidth() && imp2.getHeight() > imp.getHeight()) {
      IJ.error(TITLE, "Image to be added cannnot be larger than\n\"" + imp.getTitle() + "\".");
      return;
    }

    Overlay overlayList = imp.getOverlay();
    if (overlayList == null || settings.replace) {
      overlayList = new Overlay();
    }

    final Function<ImageProcessor, ImageRoi> createRoi = ip -> {
      final ImageRoi imageRoi = new ImageRoi(xx, yy, ip);
      imageRoi.setZeroTransparent(settings.transparent);
      imageRoi.setName(imp2.getShortTitle());
      if (settings.opacity != 100) {
        imageRoi.setOpacity(settings.opacity / 100.0);
      }
      return imageRoi;
    };

    // Support overlay of stacks
    final int nc = imp2.getNChannels();
    final int nz = imp2.getNSlices();
    final int nt = imp2.getNFrames();
    if (settings.stack && imp.getNChannels() == nc && imp.getNSlices() == nz
        && imp.getNFrames() == nt) {
      final ImageStack stack2 = imp2.getImageStack();
      for (int c = 1; c <= nc; c++) {
        for (int z = 1; z <= nz; z++) {
          for (int t = 1; t <= nt; t++) {
            // 1-based stack index
            // See ImagePlus.getStackIndex(c, z, t);
            final int i = (t - 1) * nc * nz + (z - 1) * nc + c;
            final ImageRoi r = createRoi.apply(stack2.getProcessor(i));
            // See Roi.setPosition(ImagePlus):
            if (imp2.isHyperStack()) {
              r.setPosition(c, z, t);
            } else {
              r.setPosition(i);
            }
            overlayList.add(r);
          }
        }
      }
    } else {
      // Single frame
      overlayList.add(createRoi.apply(imp2.getProcessor()));
    }

    imp.setOverlay(overlayList);
    Undo.setup(Undo.OVERLAY_ADDITION, imp);
    imp.updateAndDraw();
  }
}
