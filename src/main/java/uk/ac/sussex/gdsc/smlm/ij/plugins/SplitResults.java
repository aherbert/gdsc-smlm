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

import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.utils.ObjectAnalyzer;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;
import uk.ac.sussex.gdsc.smlm.results.procedures.XYRResultProcedure;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

/**
 * Splits PeakFit results into separate datasets using an input mask of objects.
 */
public class SplitResults implements PlugIn {
  private static final String TITLE = "Split Results";
  private static String inputOption = "";
  private static String objectMask = "";
  private static boolean showObjectMask;
  private static boolean nonMaskDataset;

  /** {@inheritDoc} */
  @Override
  public void run(String arg) {
    SMLMUsageTracker.recordPlugin(this.getClass(), arg);

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "There are no fitting results in memory");
      return;
    }
    final String[] items = ImageJUtils.getImageList(ImageJUtils.GREY_8_16);
    if (items.length == 0) {
      IJ.error(TITLE, "There are no suitable mask images");
      return;
    }

    // Show a dialog allowing the results set to be filtered
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Select a dataset to split");
    ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);
    gd.addChoice("Object_mask", items, objectMask);
    gd.addCheckbox("Show_object_mask", showObjectMask);
    gd.addCheckbox("Non_mask_dataset", nonMaskDataset);
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }

    inputOption = ResultsManager.getInputSource(gd);
    objectMask = gd.getNextChoice();
    showObjectMask = gd.getNextBoolean();
    nonMaskDataset = gd.getNextBoolean();

    final MemoryPeakResults results =
        ResultsManager.loadInputResults(inputOption, false, null, null);
    if (results == null || results.size() == 0) {
      IJ.error(TITLE, "No results could be loaded");
      return;
    }

    final ImagePlus imp = WindowManager.getImage(objectMask);
    if (imp == null) {
      IJ.error(TITLE, "No object mask could be found");
      return;
    }

    splitResults(results, imp.getProcessor());
  }

  private static void splitResults(MemoryPeakResults results, ImageProcessor ip) {
    IJ.showStatus("Splitting " + TextUtils.pleural(results.size(), "result"));

    // Create an object mask
    final ObjectAnalyzer objectAnalyzer = new ObjectAnalyzer(ip, false);

    final int maxx = ip.getWidth();
    final int maxy = ip.getHeight();

    final float scaleX = (float) results.getBounds().width / maxx;
    final float scaleY = (float) results.getBounds().height / maxy;

    // Create a results set for each object
    final int maxObject = objectAnalyzer.getMaxObject();
    final MemoryPeakResults[] resultsSet = new MemoryPeakResults[maxObject + 1];
    for (int object = 0; object <= maxObject; object++) {
      final MemoryPeakResults newResults = new MemoryPeakResults();
      newResults.copySettings(results);
      newResults.setName(results.getName() + " " + object);
      resultsSet[object] = newResults;
    }

    final int[] mask = objectAnalyzer.getObjectMask();

    if (showObjectMask) {
      final ImageProcessor objectIp =
          (maxObject <= 255) ? new ByteProcessor(maxx, maxy) : new ShortProcessor(maxx, maxy);
      for (int i = 0; i < mask.length; i++) {
        objectIp.set(i, mask[i]);
      }
      final ImagePlus imp = ImageJUtils.display(objectMask + " Objects", objectIp);
      imp.setDisplayRange(0, maxObject);
      imp.updateAndDraw();
    }

    // Process the results mapping them to their objects
    final Counter i = new Counter();
    final int size = results.size();
    final int step = ImageJUtils.getProgressInterval(size);
    results.forEach(DistanceUnit.PIXEL, new XYRResultProcedure() {
      @Override
      public void executeXYR(float xx, float yy, PeakResult result) {
        if (i.incrementAndGet() % step == 0) {
          IJ.showProgress(i.getCount(), size);
        }

        // Map to the mask objects
        final int object;
        final int x = (int) (xx / scaleX);
        final int y = (int) (yy / scaleY);
        if (x < 0 || x >= maxx || y < 0 || y >= maxy) {
          object = 0;
        } else {
          final int index = y * maxx + x;
          // Q. Is this bounds check needed given the above check shows that x,y
          // is within the bounds of the image processor?
          if (index < 0 || index >= mask.length) {
            object = 0;
          } else {
            object = mask[index];
          }
        }
        resultsSet[object].add(result);
      }
    });
    IJ.showProgress(1);

    // Add the new results sets to memory
    i.reset();
    for (int object = (nonMaskDataset) ? 0 : 1; object <= maxObject; object++) {
      if (resultsSet[object].isNotEmpty()) {
        MemoryPeakResults.addResults(resultsSet[object]);
        i.increment();
      }
    }

    IJ.showStatus("Split " + TextUtils.pleural(results.size(), "result") + " into "
        + TextUtils.pleural(i.getCount(), "set"));
  }
}
