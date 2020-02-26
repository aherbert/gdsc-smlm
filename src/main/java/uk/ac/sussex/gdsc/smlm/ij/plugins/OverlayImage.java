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

import ij.IJ;
import ij.ImagePlus;
import ij.Undo;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.ImageRoi;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import java.awt.Rectangle;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicReference;

/**
 * This plugin is extracted from ij.plugins.OverlayCommands to allow an image to be added with a
 * transparent background.
 */
public class OverlayImage implements PlugIn {
  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    int opacity;
    boolean transparent;
    boolean replace;
    String title;

    Settings() {
      // Set defaults
      opacity = 100;
      transparent = true;
      replace = true;
      title = "";
    }

    Settings(Settings source) {
      opacity = source.opacity;
      transparent = source.transparent;
      replace = source.replace;
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

    addImage();
  }

  /**
   * Adapted from ij.plugins.OverlayCommands#addImage(boolean) with the additional option for
   * setting the zero pixels to transparent.
   */
  void addImage() {
    final ImagePlus imp = IJ.getImage();
    final int[] wList = WindowManager.getIDList();
    if (wList == null || wList.length < 2) {
      IJ.error("Add Image...", "The command requires at least two open images.");
      return;
    }
    String[] titles = new String[wList.length];
    int count = 0;
    for (int i = 0; i < wList.length; i++) {
      final ImagePlus imp2 = WindowManager.getImage(wList[i]);
      if (imp2 != null && imp2 != imp && imp.getWidth() >= imp2.getWidth()
          && imp.getHeight() >= imp2.getHeight()) {
        titles[count++] = imp2.getTitle();
      }
    }
    if (count < 1) {
      IJ.error("Add Image...", "The command requires at least one valid overlay image.");
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

    settings = Settings.load();
    final GenericDialog gd = new GenericDialog("Add Image...");
    gd.addChoice("Image to add:", titles, settings.title);
    gd.addNumericField("X location:", x, 0);
    gd.addNumericField("Y location:", y, 0);
    gd.addNumericField("Opacity (0-100%):", settings.opacity, 0);
    gd.addCheckbox("Transparent background", settings.transparent);
    gd.addCheckbox("Replace overlay", settings.replace);

    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }

    settings.title = gd.getNextChoice();
    x = (int) gd.getNextNumber();
    y = (int) gd.getNextNumber();
    settings.opacity = (int) gd.getNextNumber();
    settings.transparent = gd.getNextBoolean();
    settings.replace = gd.getNextBoolean();
    settings.save();

    final ImagePlus overlay = WindowManager.getImage(settings.title);
    if (overlay == imp) {
      IJ.error("Add Image...",
          "Image to be added cannot be the same as\n\"" + imp.getTitle() + "\".");
      return;
    }
    if (overlay.getWidth() > imp.getWidth() && overlay.getHeight() > imp.getHeight()) {
      IJ.error("Add Image...",
          "Image to be added cannnot be larger than\n\"" + imp.getTitle() + "\".");
      return;
    }

    final ImageRoi roi2 = new ImageRoi(x, y, overlay.getProcessor());
    roi2.setZeroTransparent(settings.transparent);
    roi2.setName(overlay.getShortTitle());
    if (settings.opacity != 100) {
      roi2.setOpacity(settings.opacity / 100.0);
    }

    Overlay overlayList = imp.getOverlay();
    if (overlayList == null || settings.replace) {
      overlayList = new Overlay();
    }
    overlayList.add(roi2);
    imp.setOverlay(overlayList);
    Undo.setup(Undo.OVERLAY_ADDITION, imp);
  }
}
