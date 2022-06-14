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
import ij.WindowManager;
import ij.plugin.PlugIn;
import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.util.function.Consumer;
import java.util.function.Predicate;
import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.data.utils.IdentityTypeConverter;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog.OptionListener;
import uk.ac.sussex.gdsc.core.ij.roi.CoordinatePredicate;
import uk.ac.sussex.gdsc.core.ij.roi.CoordinatePredicateUtils;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.CropResultsSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResultValueParameter;
import uk.ac.sussex.gdsc.smlm.results.predicates.MinMaxPeakResultPredicate;
import uk.ac.sussex.gdsc.smlm.results.predicates.PassPeakResultPredicate;
import uk.ac.sussex.gdsc.smlm.results.procedures.MinMaxResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyrResultProcedure;

/**
 * Filters PeakFit results that are stored in memory using various fit criteria.
 */
public class CropResults implements PlugIn {
  private static final String TITLE = "Crop Results";

  /** The text description of the options for the output cropped results name. */
  static final String[] NAME_OPTIONS = {"Name", "Suffix", "Sequence"};

  /** The option to specify the entire name for the output cropped results name. */
  static final int NAME_OPTION_NAME = 0;
  /** The option to specify a suffix for the output cropped results name. */
  static final int NAME_OPTION_SUFFIX = 1;
  /**
   * The option to specify a sequence for the output cropped results name. The output name will a
   * suffix plus the current value of the name counter.
   */
  static final int NAME_OPTION_SEQUENCE = 2;

  private LocalList<String> titles;
  /** The settings. */
  CropResultsSettings.Builder settings;
  private boolean myUseRoi;
  /** The results. */
  MemoryPeakResults results;
  private String outputName;
  private MinMaxResultProcedure minMax;
  private TypeConverter<DistanceUnit> converter;
  private boolean myLimitZ;

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "There are no fitting results in memory");
      return;
    }

    // Build a list of all images with a region ROI
    titles = new LocalList<>(WindowManager.getWindowCount());
    for (final int imageId : ImageJUtils.getIdList()) {
      final ImagePlus imp = WindowManager.getImage(imageId);
      if (imp != null && imp.getRoi() != null && imp.getRoi().isArea()) {
        titles.add(imp.getTitle());
      }
    }
    final boolean roiMode = "roi".equals(arg);
    if (roiMode && titles.isEmpty()) {
      IJ.error(TITLE, "No images with an ROI");
      return;
    }

    settings = SettingsManager.readCropResultsSettings(0).toBuilder();

    // Show a dialog allowing the results set to be filtered
    final ExtendedGenericDialog gd = createDialog(roiMode);
    gd.addMessage("Select a dataset to crop");
    ResultsManager.addInput(gd, settings.getInputOption(), InputSource.MEMORY);
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    settings.setInputOption(ResultsManager.getInputSource(gd));
    results = ResultsManager.loadInputResults(settings.getInputOption(), false, null, null);
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(TITLE, "No results could be loaded");
      IJ.showStatus("");
      return;
    }

    // Allow z-filtering
    if (results.is3D()) {
      minMax = new MinMaxResultProcedure(results, new PeakResultValueParameter(PeakResult.Z));
    }

    if (roiMode) {
      runRoiCrop();
    } else {
      runCrop();
    }

    SettingsManager.writeSettings(settings);
  }

  private void runCrop() {
    if (!showCropDialog()) {
      return;
    }
    cropResults();
  }

  private boolean showCropDialog() {
    final ExtendedGenericDialog gd = createDialog(false);

    final Rectangle bounds = results.getBounds(true);
    results.is3D();

    gd.addMessage(String.format("Current bounds: x=%d,y=%d,w=%d,h=%d", bounds.x, bounds.y,
        bounds.width, bounds.height));
    gd.addNumericField("Border", settings.getBorder(), 2);
    gd.addCheckbox("Select_region", settings.getSelectRegion());
    gd.addNumericField("X", settings.getX(), 2);
    gd.addNumericField("Y", settings.getY(), 2);
    gd.addNumericField("Width", settings.getWidth(), 2);
    gd.addNumericField("Height", settings.getHeight(), 2);
    if (!titles.isEmpty()) {
      gd.addCheckbox("Use_ROI", settings.getUseRoi());
      final String[] items = titles.toArray(new String[0]);
      gd.addChoice("Image", items, settings.getRoiImage());
    }
    addStandardFields(gd);
    gd.addCheckbox("Reset_origin", settings.getResetOrigin());

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.setBorder(Math.max(0, gd.getNextNumber()));
    settings.setSelectRegion(gd.getNextBoolean());
    settings.setX(gd.getNextNumber());
    settings.setY(gd.getNextNumber());
    settings.setWidth(Math.max(0, gd.getNextNumber()));
    settings.setHeight(Math.max(0, gd.getNextNumber()));
    if (!titles.isEmpty()) {
      myUseRoi = gd.getNextBoolean();
      settings.setUseRoi(myUseRoi);
      settings.setRoiImage(gd.getNextChoice());
    }
    readStandardFields(gd);
    settings.setResetOrigin(gd.getNextBoolean());

    gd.collectOptions();

    return validateOutputName();
  }

  private void addStandardFields(final ExtendedGenericDialog gd) {
    if (minMax != null) {
      // 3D crop options
      double min = minMax.getMinimum();
      double max = minMax.getMaximum();

      final double maxz = Math.min(settings.getMaxZ(), max);
      final double minz = Math.max(settings.getMinZ(), min);

      // Display in nm
      converter = new IdentityTypeConverter<>(null);
      String unit = "";

      final DistanceUnit nativeUnit = results.getDistanceUnit();
      if (nativeUnit != null) {
        unit = UnitHelper.getShortName(nativeUnit);
        try {
          converter =
              CalibrationHelper.getDistanceConverter(results.getCalibration(), DistanceUnit.NM);
          unit = UnitHelper.getShortName(DistanceUnit.NM);
        } catch (final ConversionException ignored) {
          // No native units
        }
      }
      min = converter.convert(min);
      max = converter.convert(max);

      final String msg = String.format("%.2f <= z <= %.2f (%s)", min, max, unit);

      min = Math.floor(min);
      max = Math.ceil(max);

      gd.addMessage(msg);
      gd.addCheckbox("Limit Z-depth", settings.getLimitZ());
      gd.addSlider("minZ", min, max, converter.convert(minz));
      gd.addSlider("maxZ", min, max, converter.convert(maxz));
    }

    gd.addChoice("Name_option", NAME_OPTIONS, settings.getNameOption(),
        new OptionListener<Integer>() {
          @Override
          public boolean collectOptions(Integer value) {
            settings.setNameOption(value);
            return collectOptions(false);
          }

          @Override
          public boolean collectOptions() {
            return collectOptions(true);
          }

          private boolean collectOptions(boolean silent) {
            final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
            if (settings.getNameOption() == NAME_OPTION_NAME) {
              final String name = (TextUtils.isNullOrEmpty(settings.getOutputName()))
                  ? (results.getName() + " Cropped")
                  : settings.getOutputName();
              egd.addStringField("Output_name", name, MathUtils.clip(60, 120, name.length()));
            } else if (settings.getNameOption() == NAME_OPTION_SUFFIX) {
              final String name = (TextUtils.isNullOrEmpty(settings.getNameSuffix())) ? " Cropped"
                  : settings.getNameSuffix();
              egd.addStringField("Name_suffix", name, MathUtils.clip(20, 60, name.length()));
            } else if (settings.getNameOption() == NAME_OPTION_SEQUENCE) {
              final String name = settings.getNameSuffix();
              egd.addStringField("Name_suffix", name, MathUtils.clip(20, 60, name.length()));
              int counter = settings.getNameCounter();
              if (counter < 1) {
                counter = 1;
              }
              egd.addNumericField("Name_counter", counter, 0);
            } else {
              throw new IllegalStateException("Unknown name option: " + settings.getNameOption());
            }
            egd.setSilent(silent);
            egd.showDialog(true, gd);
            if (egd.wasCanceled()) {
              return false;
            }
            if (settings.getNameOption() == NAME_OPTION_NAME) {
              settings.setOutputName(egd.getNextString());
            } else if (settings.getNameOption() == NAME_OPTION_SUFFIX) {
              settings.setNameSuffix(egd.getNextString());
            } else if (settings.getNameOption() == NAME_OPTION_SEQUENCE) {
              settings.setNameSuffix(egd.getNextString());
              settings.setNameCounter(Math.max(1, (int) egd.getNextNumber()));
            }

            return true;
          }
        });
    gd.addCheckbox("Preserve_bounds", settings.getPreserveBounds());
  }

  private void readStandardFields(ExtendedGenericDialog gd) {
    if (minMax != null) {
      // 3D crop options
      myLimitZ = gd.getNextBoolean();
      settings.setLimitZ(myLimitZ);
      settings.setMinZ(converter.convertBack(gd.getNextNumber()));
      settings.setMaxZ(converter.convertBack(gd.getNextNumber()));
    }

    settings.setNameOption(gd.getNextChoiceIndex());
    settings.setPreserveBounds(gd.getNextBoolean());
  }

  private boolean validateOutputName() {
    if (settings.getNameOption() == NAME_OPTION_NAME) {
      outputName = settings.getOutputName();
      if (TextUtils.isNullOrEmpty(outputName)) {
        IJ.error(TITLE, "No output name");
        return false;
      }
    } else if (settings.getNameOption() == NAME_OPTION_SUFFIX) {
      final String suffix = settings.getNameSuffix();
      if (TextUtils.isNullOrEmpty(suffix)) {
        IJ.error(TITLE, "No output suffix");
        return false;
      }
      outputName = results.getName() + suffix;
    } else if (settings.getNameOption() == NAME_OPTION_SEQUENCE) {
      outputName = results.getName();
      final String suffix = settings.getNameSuffix();
      if (!TextUtils.isNullOrEmpty(suffix)) {
        outputName += suffix;
      }
      final int count = settings.getNameCounter();
      outputName += count;
      settings.setNameCounter(count + 1); // Increment for next time
    }
    return true;
  }

  /**
   * Apply the filters to the data.
   */
  private void cropResults() {
    final MemoryPeakResults newResults = createNewResults();

    // These bounds are integer. But this is because the results are meant to come from an image.
    final Rectangle integerBounds = results.getBounds(true);

    // The crop bounds can be floating point...

    // Border
    final double border = settings.getBorder();
    final double xx = integerBounds.x + border;
    final double yy = integerBounds.y + border;
    final double w = Math.max(0, integerBounds.width - 2 * border);
    final double h = Math.max(0, integerBounds.height - 2 * border);
    Rectangle2D pixelBounds = new Rectangle2D.Double(xx, yy, w, h);

    // Bounding box
    if (settings.getSelectRegion()) {
      final Rectangle2D boxBounds = new Rectangle2D.Double(settings.getX(), settings.getY(),
          settings.getWidth(), settings.getHeight());
      pixelBounds = pixelBounds.createIntersection(boxBounds);
    }

    // If an ROI was chosen from an image, scale the roi to the bounds of this dataset
    // and create another intersection
    if (myUseRoi) {
      final ImagePlus imp = WindowManager.getImage(settings.getRoiImage());
      if (imp != null && imp.getRoi() != null) {
        final Rectangle roi = imp.getRoi().getBounds();
        final int roiImageWidth = imp.getWidth();
        final int roiImageHeight = imp.getHeight();

        final double xscale = (double) roiImageWidth / integerBounds.width;
        final double yscale = (double) roiImageHeight / integerBounds.height;

        final Rectangle2D roiBounds = new Rectangle2D.Double(roi.x / xscale, roi.y / yscale,
            roi.width / xscale, roi.height / yscale);
        pixelBounds = pixelBounds.createIntersection(roiBounds);
      }
    }

    final Rectangle2D bounds = pixelBounds;

    final Predicate<PeakResult> testZ = getZFilter();
    // Copy the results if the origin is reset
    final Consumer<PeakResult> consumer =
        settings.getResetOrigin() ? r -> newResults.add(r.copy()) : newResults::add;

    if (bounds.getWidth() > 0 && bounds.getHeight() > 0) {
      results.forEach(DistanceUnit.PIXEL, (XyrResultProcedure) (x, y, result) -> {
        if (bounds.contains(x, y) && testZ.test(result)) {
          consumer.accept(result);
        }
      });
    }

    if (settings.getPreserveBounds()) {
      newResults.setBounds(integerBounds);
    } else {
      // Get the lower and upper integer limits
      final int x = (int) Math.floor(bounds.getX());
      final int y = (int) Math.floor(bounds.getY());
      final int ux = (int) Math.ceil(bounds.getX() + bounds.getWidth());
      final int uy = (int) Math.ceil(bounds.getY() + bounds.getHeight());
      // Ensure the width and height are at least 1
      newResults.setBounds(new Rectangle(x, y, Math.max(1, ux - x), Math.max(1, uy - y)));

      if (settings.getResetOrigin()) {
        newResults.translate(-x, -y);
      }
    }

    IJ.showStatus(newResults.size() + " Cropped localisations");
  }

  /**
   * Gets the z filter if cropping the z coordinate range.
   *
   * @return the z filter
   */
  private Predicate<PeakResult> getZFilter() {
    if (myLimitZ) {
      return new MinMaxPeakResultPredicate((float) settings.getMinZ(), (float) settings.getMaxZ(),
          new PeakResultValueParameter(PeakResult.Z));
    }
    return PassPeakResultPredicate.INSTANCE;
  }

  private MemoryPeakResults createNewResults() {
    final MemoryPeakResults newResults = new MemoryPeakResults();
    newResults.copySettings(results);
    newResults.setName(outputName);
    MemoryPeakResults.addResults(newResults);
    return newResults;
  }

  private void runRoiCrop() {
    if (!showRoiCropDialog()) {
      return;
    }
    roiCropResults();
  }

  private boolean showRoiCropDialog() {
    final ExtendedGenericDialog gd = createDialog(true);

    final String[] items = titles.toArray(new String[0]);
    gd.addMessage("Use ROI from ...");
    gd.addChoice("Image", items, settings.getRoiImage());
    addStandardFields(gd);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.setRoiImage(gd.getNextChoice());
    readStandardFields(gd);

    gd.collectOptions();

    return validateOutputName();
  }

  private static ExtendedGenericDialog createDialog(boolean roiMode) {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl(roiMode ? "roi-crop-results" : "crop-results"));
    return gd;
  }

  private void roiCropResults() {
    final ImagePlus imp = WindowManager.getImage(settings.getRoiImage());
    if (imp == null) {
      IJ.error(TITLE, "No ROI image: " + settings.getRoiImage());
      return;
    }
    final CoordinatePredicate roiTest =
        CoordinatePredicateUtils.createContainsPredicate(imp.getRoi());
    if (roiTest == null) {
      IJ.error(TITLE, "Not an area ROI");
      return;
    }

    // These bounds are integer. But this is because the results are meant to come from an image.
    final Rectangle integerBounds = results.getBounds(true);

    // Scale the results to the size of the image with the ROI
    final int roiImageWidth = imp.getWidth();
    final int roiImageHeight = imp.getHeight();
    final double ox = integerBounds.getX();
    final double oy = integerBounds.getY();
    final double xscale = roiImageWidth / integerBounds.getWidth();
    final double yscale = roiImageHeight / integerBounds.getHeight();

    final Predicate<PeakResult> testZ = getZFilter();
    final MemoryPeakResults newResults = createNewResults();

    results.forEach(DistanceUnit.PIXEL, (XyrResultProcedure) (x, y, result) -> {
      if (roiTest.test((x - ox) * xscale, (y - oy) * yscale) && testZ.test(result)) {
        newResults.add(result);
      }
    });

    if (settings.getPreserveBounds()) {
      newResults.setBounds(integerBounds);
    } else {
      newResults.setBounds(null);
      newResults.getBounds(true);
    }

    IJ.showStatus(newResults.size() + " Cropped localisations");
  }
}
