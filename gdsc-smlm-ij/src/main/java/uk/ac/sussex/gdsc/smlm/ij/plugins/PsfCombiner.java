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
import ij.Macro;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.plugin.frame.Recorder;
import ij.process.Blitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.MultiDialog;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.ImagePSF;
import uk.ac.sussex.gdsc.smlm.ij.settings.ImagePsfHelper;
import uk.ac.sussex.gdsc.smlm.ij.utils.ImageJImageConverter;

/**
 * Produces an average PSF image from multiple PSF images.
 *
 * <p>The input images must be a z-stack of a PSF. These can be produced using the PSFCreator
 * plugin.
 */
public class PsfCombiner implements PlugIn {
  private static final String TITLE = "PSF Combiner";

  private static AtomicReference<List<String>> lastSelected = new AtomicReference<>();
  private final List<Psf> input = new LinkedList<>();

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    // Build a list of suitable images
    final List<String> titles = createImageList();
    if (titles.size() < 2) {
      IJ.error(TITLE, "No suitable PSF images to combine");
      return;
    }

    final MultiDialog md = new MultiDialog("Select PSFs", titles);
    md.setSelected(lastSelected.get());

    md.setHelpUrl(HelpUrls.getUrl("psf-combiner"));
    md.showDialog();

    if (md.wasCancelled()) {
      return;
    }

    final List<String> selected = md.getSelectedResults();
    if (selected.size() < 2) {
      IJ.error(TITLE, "Require at least 2 PSF images to combine");
      return;
    }

    lastSelected.set(selected);

    for (final String title : selected) {
      input.add(new Psf(title));
    }

    if (input.isEmpty()) {
      return;
    }

    if (input.size() < 2) {
      return;
    }

    // GenericDialog gd = new GenericDialog(TITLE);
    // gd.addMessage("Set the maximum z-depth +/- from the PSF centre");
    // gd.addSlider("Z-depth", 20, 200, zDepth);
    // gd.showDialog();
    // if (gd.wasCanceled())
    // return;
    // zDepth = Math.abs((int) gd.getNextNumber());
    //
    // for (PSF psf : input)
    // psf.crop(zDepth);

    combineImages();
  }

  /**
   * Creates the image list. Images must be greyscale, square and a stack of a single channel.
   *
   * @return the list
   */
  public static List<String> createImageList() {
    final List<String> titles = new LocalList<>();
    final int[] ids = WindowManager.getIDList();
    if (ids != null) {
      for (final int id : ids) {
        final ImagePlus imp = WindowManager.getImage(id);
        if (imp != null
            // Image must be greyscale
            && (imp.getType() == ImagePlus.GRAY8 || imp.getType() == ImagePlus.GRAY16
                || imp.getType() == ImagePlus.GRAY32)
            // Image must be square and a stack of a single channel
            && (imp.getWidth() == imp.getHeight() && imp.getNChannels() == 1)
            // Check if these are PSF images created by the SMLM plugins
            && containsPsf(imp)) {
          titles.add(imp.getTitle());
        }
      }
    }
    return titles;
  }

  private static boolean containsPsf(ImagePlus imp) {
    return PsfDrift.getPsfSettings(imp) != null;
  }

  private void combineImages() {
    final double nmPerPixel = getNmPerPixel();
    if (nmPerPixel <= 0) {
      return;
    }
    final double nmPerSlice = getNmPerSlice();
    if (nmPerSlice <= 0) {
      return;
    }

    // Find the lowest & highest dimensions
    int minStart = Integer.MAX_VALUE;
    int maxStart = Integer.MIN_VALUE;
    int minEnd = Integer.MAX_VALUE;
    int maxEnd = Integer.MIN_VALUE;
    int minSize = Integer.MAX_VALUE;
    int maxSize = 0;
    for (final Psf psf : input) {
      if (minStart > psf.start) {
        minStart = psf.start;
      }
      if (maxStart < psf.start) {
        maxStart = psf.start;
      }
      if (maxEnd < psf.getEnd()) {
        maxEnd = psf.getEnd();
      }
      if (minEnd > psf.getEnd()) {
        minEnd = psf.getEnd();
      }
      if (maxSize < psf.getSize()) {
        maxSize = psf.getSize();
      }
      if (minSize > psf.getSize()) {
        minSize = psf.getSize();
      }
    }

    int size = maxSize;
    int shift = -minStart;
    int depth = maxEnd - minStart + 1;

    // Option to crop. Do this before processing as it will make the plugin faster
    if (minStart < maxStart || minEnd < maxEnd || minSize < maxSize) {
      boolean crop;
      if (ImageJUtils.isMacro()) {
        final String options = Macro.getOptions();
        crop = options.contains(" crop");
      } else {
        final GenericDialog gd = new GenericDialog(TITLE);
        ImageJUtils.addMessage(gd,
            "The range of the PSFs is different:\nStart %d to %d\nEnd %d to %d\n"
                + "Size %d to %d\n \nCrop to the smallest?",
            minStart, maxStart, minEnd, maxEnd, minSize, maxSize);
        gd.enableYesNoCancel();
        gd.addHelp(HelpUrls.getUrl("psf-combiner"));
        gd.showDialog();
        if (gd.wasCanceled()) {
          return;
        }
        crop = gd.wasOKed();
      }

      if (crop) {
        Recorder.recordOption("crop");

        for (final Psf psf : input) {
          psf.crop(maxStart, minEnd, minSize);
        }

        size = minSize;
        shift = -maxStart;
        depth = minEnd - maxStart + 1;
      }
    }

    // Shift all stacks
    int totalImages = 0;
    for (final Psf psf : input) {
      psf.start += shift;
      totalImages += psf.psfSettings.getImageCount();
    }

    // Create a stack to hold all the images
    // Create a stack to hold the sum of the weights
    final ImageStack stack = new ImageStack(size, size, depth);
    final ImageStack stackW = new ImageStack(size, size, depth);
    for (int n = 1; n <= depth; n++) {
      stack.setPixels(new float[size * size], n);
      stackW.setPixels(new float[size * size], n);
    }

    // Insert all the PSFs
    IJ.showStatus("Creating combined image ...");
    int imageNo = 0;
    final double fraction = 1.0 / input.size();
    for (final Psf psf : input) {
      double progress = imageNo * fraction;
      final ImageStack psfStack = psf.psfStack;
      final int w = psf.getSize();
      final int offsetXy = (size - w) / 2;
      final int offsetZ = psf.start;
      final double weight = (1.0 * psf.psfSettings.getImageCount()) / totalImages;
      final FloatProcessor wp = new FloatProcessor(w, w, SimpleArrayUtils
          .newFloatArray(psfStack.getWidth() * psfStack.getHeight(), (float) weight));
      final double increment = fraction / psfStack.getSize();
      for (int n = 1; n <= psfStack.getSize(); n++) {
        progress += increment;
        IJ.showProgress(progress);

        // Get the data and adjust using the weight
        final float[] psfData = ImageJImageConverter.getData(psfStack.getProcessor(n));
        for (int i = 0; i < psfData.length; i++) {
          psfData[i] *= weight;
        }

        // Insert into the combined PSF
        final int slice = n + offsetZ;
        ImageProcessor ip = stack.getProcessor(slice);
        ip.copyBits(new FloatProcessor(w, w, psfData), offsetXy, offsetXy, Blitter.ADD);

        // Insert the weights
        ip = stackW.getProcessor(slice);
        ip.copyBits(wp, offsetXy, offsetXy, Blitter.ADD);
      }
      imageNo++;
    }

    // Normalise
    for (int n = 1; n <= depth; n++) {
      stack.getProcessor(n).copyBits(stackW.getProcessor(n), 0, 0, Blitter.DIVIDE);
    }

    ImageJUtils.finished();

    final ImagePlus imp = ImageJUtils.display("Combined PSF", stack);
    imp.setSlice(1 + shift);
    imp.resetDisplayRange();
    imp.updateAndDraw();

    final double fwhm = getFwhm();
    imp.setProperty("Info", ImagePsfHelper.toString(
        ImagePsfHelper.create(imp.getSlice(), nmPerPixel, nmPerSlice, totalImages, fwhm)));

    ImageJUtils.log("%s : z-centre = %d, nm/Pixel = %s, nm/Slice = %s, %d images, FWHM = %s\n",
        imp.getTitle(), imp.getSlice(), MathUtils.rounded(nmPerPixel),
        MathUtils.rounded(nmPerSlice), totalImages, MathUtils.rounded(fwhm));
  }

  private double getNmPerPixel() {
    final double nmPerPixel = input.get(0).psfSettings.getPixelSize();
    for (final Psf psf : input) {
      if (psf.psfSettings.getPixelSize() != nmPerPixel) {
        IJ.error(TITLE, "Different pixel size resolutions for the input PSFs");
        return -1;
      }
    }
    return nmPerPixel;
  }

  private double getNmPerSlice() {
    final double nmPerSlice = input.get(0).psfSettings.getPixelDepth();
    for (final Psf psf : input) {
      if (psf.psfSettings.getPixelDepth() != nmPerSlice) {
        IJ.error(TITLE, "Different pixel depth resolutions for the input PSFs");
        return -1;
      }
    }
    return nmPerSlice;
  }

  private double getFwhm() {
    double fwhm = 0;
    for (final Psf psf : input) {
      fwhm += psf.psfSettings.getFwhm();
    }
    return fwhm / input.size();
  }

  /**
   * Represent a PSF with an ImageStack and ImagePSF settings.
   */
  private static class Psf {
    ImagePSF psfSettings;
    int start;
    ImageStack psfStack;

    Psf(String title) {
      final ImagePlus imp = WindowManager.getImage(title);
      if (imp == null) {
        throw new IllegalArgumentException("No image with title: " + title);
      }

      this.psfSettings = PsfDrift.getPsfSettings(imp);
      if (psfSettings == null) {
        throw new IllegalArgumentException("Unknown PSF settings for image: " + title);
      }
      final int zCentre = psfSettings.getCentreImage();
      if (zCentre < 1 || zCentre > imp.getStackSize()) {
        throw new IllegalArgumentException(
            "z-centre must be within the stack size: " + imp.getStackSize());
      }
      start = 1 - zCentre;
      psfStack = imp.getImageStack();
    }

    void crop(int maxStart, int minEnd, int size) {
      final int removeBorder = (getSize() - size) / 2;

      final int removeStart = maxStart - start; // Should be positive
      final int removeEnd = getEnd() - minEnd; // Should be positive

      if (removeStart < 0 || removeEnd < 0 || removeBorder < 0) {
        throw new IllegalArgumentException("Bad PSF crop");
      }

      if (removeStart > 0 || removeEnd > 0 || removeBorder > 0) {
        final int depth = psfStack.getSize() - removeStart - removeEnd;
        psfStack = psfStack.crop(removeBorder, removeBorder, removeStart, size, size, depth);

        // Update range
        start += removeStart;
      }
    }

    int getEnd() {
      return start + psfStack.getSize();
    }

    int getSize() {
      return psfStack.getWidth();
    }
  }
}
