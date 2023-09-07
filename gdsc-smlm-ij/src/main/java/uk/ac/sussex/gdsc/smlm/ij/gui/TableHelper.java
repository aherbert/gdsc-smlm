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

package uk.ac.sussex.gdsc.smlm.ij.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import java.io.File;
import java.util.Arrays;
import java.util.function.Consumer;
import java.util.function.Supplier;
import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.OffsetPointRoi;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.XmlUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationHelper;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.IJImageSource;
import uk.ac.sussex.gdsc.smlm.ij.SeriesImageSource;
import uk.ac.sussex.gdsc.smlm.ij.plugins.TiffSeriesViewer.TiffSeriesVirtualStack;
import uk.ac.sussex.gdsc.smlm.results.ImageSource;
import uk.ac.sussex.gdsc.smlm.results.ImageSource.ReadHint;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.sort.FramePeakResultComparator;

/**
 * A helper class for the results tables.
 */
final class TableHelper {
  /**
   * No public constructor.
   */
  private TableHelper() {}

  /**
   * Save results to memory.
   *
   * @param title the title of the data
   * @param name the proposed name of the dataset (can be null)
   * @param resultsSupplier the results supplier
   * @return the result name if saved, else the original name argument
   */
  static String saveResults(String title, String name,
      Supplier<MemoryPeakResults> resultsSupplier) {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog("Save Results");
    String saveName = TextUtils.isNullOrEmpty(name) ? title : name;
    gd.addStringField("Results_set_name", saveName, 30);
    gd.showDialog();
    if (gd.wasCanceled()) {
      return name;
    }
    saveName = gd.getNextString();
    if (TextUtils.isNullOrEmpty(saveName)) {
      IJ.error("No results set name");
      return name;
    }
    final MemoryPeakResults results = resultsSupplier.get();
    results.setName(saveName);
    MemoryPeakResults.addResults(results);
    return saveName;
  }

  /**
   * Update the source from available 2D/3D images.
   *
   * @param action the consumer for the new source
   */
  static void updateSource(Consumer<ImageSource> action) {
    // Assumes 3D stack (no channel/time)
    final String[] list = ImageJUtils.getImageList(imp -> imp.getNDimensions() <= 3);
    if (list.length == 0) {
      IJ.log("No suitable source images (require a 3D stack)");
      return;
    }
    final ExtendedGenericDialog gd = new ExtendedGenericDialog("Update source image");
    gd.addChoice("Image", list, 0);
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    final ImagePlus imp = WindowManager.getImage(gd.getNextChoice());
    if (imp == null) {
      return;
    }
    action.accept(new IJImageSource(imp));
  }

  /**
   * Show the source info in the ImageJ log window.
   *
   * @param title the title
   * @param source the source
   */
  static void showInfo(String title, ImageSource source) {
    final StringBuilder sb = new StringBuilder(title).append(" source");
    if (source == null) {
      sb.append(" = NA");
    } else {
      sb.append('\n').append(XmlUtils.prettyPrintXml(source.toXml()));
    }
    // Note:
    // If a raw path is printed to the ImageJ log double-clicking it will open the image.
    // We could separate these onto multiple lines:
    // <path>/path/to/image.tif</path>
    // <string>/path/to/image.tif</string>
    IJ.log(sb.toString());
  }

  /**
   * Show the source image.
   *
   * @param source the source
   */
  static void showImage(ImageSource source) {
    if (source == null) {
      return;
    }
    // Check if already open
    final ImagePlus imp = WindowManager.getImage(source.getName());
    if (imp != null) {
      imp.getWindow().toFront();
      return;
    }

    // Check if an ImageJ image source
    if (source instanceof IJImageSource) {
      final IJImageSource imageSource = (IJImageSource) source;
      final String path = imageSource.getPath();
      if (path != null && new File(path).exists()) {
        IJ.showStatus("Opening image ...");
        IJ.open(path);
        IJ.showStatus("");
      } else {
        IJ.log("Cannot find the image source: " + path);
      }
      return;
    }
    // Open a SeriesImageSource.
    if (source instanceof SeriesImageSource) {
      final SeriesImageSource imageSource = (SeriesImageSource) source;
      imageSource.setBufferLimit(0); // No memory buffer
      imageSource.setReadHint(ReadHint.NONSEQUENTIAL);
      if (!source.open()) {
        IJ.log("Cannot open the series image source");
        return;
      }
      new TiffSeriesVirtualStack(imageSource).show();
    }
  }

  /**
   * Overlay the results on the image.
   *
   * @param source the source
   * @param calibration the calibration
   * @param results the results
   */
  static void showOverlay(ImageSource source, Calibration calibration, PeakResult[] results) {
    if (source == null || results.length == 0) {
      return;
    }
    final String title = source.getOriginal().getName();
    final ImagePlus imp = WindowManager.getImage(title);
    // Assumes 3D stack (no channel/time)
    if ((imp == null) || (imp.getNDimensions() > 3)) {
      return;
    }
    try {
      final TypeConverter<DistanceUnit> converter =
          CalibrationHelper.getDistanceConverter(calibration, DistanceUnit.PIXEL);
      final Overlay o = new Overlay();

      if (results.length == 1) {
        final PeakResult p = results[0];
        final PointRoi roi = new OffsetPointRoi(converter.convert(p.getXPosition()),
            converter.convert(p.getYPosition()));
        roi.setPointType(3);
        roi.setPosition(p.getFrame());
        o.add(roi);
      } else {
        Arrays.sort(results, FramePeakResultComparator.INSTANCE);
        final float[] ox = new float[results.length];
        final float[] oy = new float[results.length];
        int size = 0;
        int frame = results[0].getFrame() - 1;
        for (final PeakResult r : results) {
          if (frame != r.getFrame()) {
            if (size != 0) {
              // Arrays are copied in PolygonRoi
              final PointRoi roi = new OffsetPointRoi(ox, oy, size);
              roi.setPointType(3);
              roi.setPosition(frame);
              size = 0;
              o.add(roi);
            }
            frame = r.getFrame();
          }
          ox[size] = converter.convert(r.getXPosition());
          oy[size] = converter.convert(r.getYPosition());
          size++;
        }
        if (size > 0) {
          final PointRoi roi = new OffsetPointRoi(ox, oy, size);
          roi.setPointType(3);
          roi.setPosition(frame);
          o.add(roi);
        }
      }
      imp.setOverlay(o);
      final PeakResult p = results[0];
      imp.setSlice(p.getFrame());
      ImageJUtils.adjustSourceRect(imp, 0, (int) converter.convert(p.getXPosition()),
          (int) converter.convert(p.getYPosition()));
      imp.getWindow().toFront();
    } catch (final ConversionException ignored) {
      // Ignore
    }
  }
}
