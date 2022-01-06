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

import gnu.trove.map.hash.TIntObjectHashMap;
import ij.IJ;
import ij.ImageListener;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.ImageRoi;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.Blitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.LUT;
import java.awt.AWTEvent;
import java.awt.Choice;
import java.awt.Color;
import java.awt.Label;
import java.awt.Rectangle;
import java.awt.Scrollbar;
import java.awt.TextField;
import java.awt.event.ItemEvent;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.atomic.AtomicReference;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.ImageAdapter;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.NonBlockingExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.OffsetPointRoi;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper.LutColour;
import uk.ac.sussex.gdsc.core.match.AucCalculator;
import uk.ac.sussex.gdsc.core.match.BasePoint;
import uk.ac.sussex.gdsc.core.match.ClassificationResult;
import uk.ac.sussex.gdsc.core.match.Coordinate;
import uk.ac.sussex.gdsc.core.match.FractionalAssignment;
import uk.ac.sussex.gdsc.core.match.ImmutableFractionalAssignment;
import uk.ac.sussex.gdsc.core.match.RankedScoreCalculator;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.RampedScore;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.StoredData;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.DataFilterMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.DataFilterType;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitEngineSettings;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PsfProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.TemplateProtos.TemplateSettings;
import uk.ac.sussex.gdsc.smlm.engine.FitConfiguration;
import uk.ac.sussex.gdsc.smlm.engine.FitEngineConfiguration;
import uk.ac.sussex.gdsc.smlm.filters.MaximaSpotFilter;
import uk.ac.sussex.gdsc.smlm.filters.Spot;
import uk.ac.sussex.gdsc.smlm.filters.SpotFilterHelper;
import uk.ac.sussex.gdsc.smlm.ij.IJImageSource;
import uk.ac.sussex.gdsc.smlm.ij.plugins.PeakFit.FitConfigurationProvider;
import uk.ac.sussex.gdsc.smlm.ij.plugins.PeakFit.RelativeParameterProvider;
import uk.ac.sussex.gdsc.smlm.ij.plugins.benchmark.CreateData;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.model.camera.CameraModel;
import uk.ac.sussex.gdsc.smlm.model.camera.FakePerPixelCameraModel;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;

/**
 * Runs the candidate maxima identification on the image and provides a preview using an overlay.
 */
public class SpotFinderPreview implements ExtendedPlugInFilter {
  private static final String TITLE = "Spot Finder Preview";

  private static final int FLAGS = DOES_16 | DOES_8G | DOES_32 | NO_CHANGES;
  private FitEngineConfiguration config;
  private FitConfiguration fitConfig;
  private Overlay overlay;
  private ImagePlus imp;
  private boolean preview;
  private Label label;
  private TIntObjectHashMap<List<Coordinate>> actualCoordinates;

  private int currentSlice;
  private MaximaSpotFilter filter;

  // All the fields that will be updated when reloading the configuration file
  private Choice textCameraModelName;
  private Choice textPsf;
  private Choice textDataFilterType;
  private Choice textDataFilterMethod;
  private TextField textSmooth;
  private Choice textDataFilterMethod2;
  private TextField textSmooth2;
  private TextField textSearch;
  private TextField textBorder;

  // For adjusting the selction sliders
  private Scrollbar topNScrollBar;
  private Scrollbar selectScrollBar;

  private boolean refreshing;
  private NonBlockingExtendedGenericDialog gd;

  private final SpotFilterHelper spotFilterHelper = new SpotFilterHelper();

  private Calibration lastCalibration;
  private FitEngineSettings lastFitEngineSettings;
  private PSF lastPsf;

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    DataFilterMethod defaultDataFilterMethod;
    double defaultSmooth;
    double distance;
    double lowerDistance;
    boolean multipleMatches;
    boolean showTP;
    boolean showFP;
    int topN;
    int select;
    int neighbourRadius;

    Settings() {
      // Set defaults
      final FitEngineConfiguration c = new FitEngineConfiguration();
      defaultDataFilterMethod = c.getDataFilterMethod(0);
      defaultSmooth = c.getDataFilterParameterValue(0);
      distance = 1.5;
      lowerDistance = 50;
      showTP = true;
      showFP = true;
      topN = 100;
      select = 1;
      neighbourRadius = 4;
    }

    Settings(Settings source) {
      defaultDataFilterMethod = source.defaultDataFilterMethod;
      defaultSmooth = source.defaultSmooth;
      distance = source.distance;
      lowerDistance = source.lowerDistance;
      multipleMatches = source.multipleMatches;
      showTP = source.showTP;
      showFP = source.showFP;
      topN = source.topN;
      select = source.select;
      neighbourRadius = source.neighbourRadius;
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
  public int setup(String arg, ImagePlus imp) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (imp == null) {
      IJ.noImage();
      return DONE;
    }

    final Roi roi = imp.getRoi();
    if (roi != null && roi.getType() != Roi.RECTANGLE) {
      IJ.error("Rectangular ROI required");
      return DONE;
    }

    return FLAGS;
  }

  @Override
  public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
    settings = Settings.load();
    this.overlay = imp.getOverlay();
    this.imp = imp;

    // Saved by reference so do it once now
    settings.save();

    // The image is locked by the PlugInFilterRunner so unlock it to allow scroll.
    // This should be OK as the image data is not modified and only the overlay is
    // adjusted. If another plugin changes the image then the preview should update
    // the overlay and it will be obvious to the user to turn this plugin off.
    imp.unlock();

    config = SettingsManager.readFitEngineConfiguration(0);
    fitConfig = config.getFitConfiguration();

    gd = new NonBlockingExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("spot-finder-preview"));
    gd.addMessage("Preview candidate maxima");

    final String[] templates = ConfigurationTemplate.getTemplateNames(true);
    gd.addChoice("Template", templates, templates[0]);

    final String[] models = CameraModelManager.listCameraModels(true);
    gd.addChoice("Camera_model_name", models, fitConfig.getCameraModelName());

    PeakFit.addPsfOptions(gd, (FitConfigurationProvider) () -> fitConfig);
    final PeakFit.SimpleFitEngineConfigurationProvider provider =
        new PeakFit.SimpleFitEngineConfigurationProvider(config);
    PeakFit.addDataFilterOptions(gd, provider);
    gd.addChoice("Spot_filter_2", SettingsManager.getDataFilterMethodNames(),
        config.getDataFilterMethod(1, settings.defaultDataFilterMethod).ordinal());
    PeakFit.addRelativeParameterOptions(gd,
        new RelativeParameterProvider(2.5, 4.5, "Smoothing_2", provider) {
          @Override
          void setAbsolute(boolean absolute) {
            final FitEngineConfiguration c =
                fitEngineConfigurationProvider.getFitEngineConfiguration();
            final DataFilterMethod m = c.getDataFilterMethod(1, settings.defaultDataFilterMethod);
            final double smooth = c.getDataFilterParameterValue(1, settings.defaultSmooth);
            c.setDataFilter(m, smooth, absolute, 1);
          }

          @Override
          boolean isAbsolute() {
            return fitEngineConfigurationProvider.getFitEngineConfiguration()
                .getDataFilterParameterAbsolute(1, false);
          }

          @Override
          double getValue() {
            return fitEngineConfigurationProvider.getFitEngineConfiguration()
                .getDataFilterParameterValue(1, settings.defaultSmooth);
          }
        });

    PeakFit.addSearchOptions(gd, provider);
    PeakFit.addBorderOptions(gd, provider);

    // Find if this image was created with ground truth data
    if (imp.getID() == CreateData.getImageId()) {
      final MemoryPeakResults results = CreateData.getResults();
      if (results != null) {
        gd.addSlider("Match_distance", 0, 2.5, settings.distance);
        gd.addSlider("Lower_match_distance (%)", 0, 100, settings.lowerDistance);
        gd.addCheckbox("Multiple_matches", settings.multipleMatches);
        gd.addCheckbox("Show_TP", settings.showTP);
        gd.addCheckbox("Show_FP", settings.showFP);
        gd.addMessage("");
        label = (Label) gd.getMessage();
        final boolean integerCoords = false;
        actualCoordinates = ResultsMatchCalculator.getCoordinates(results, integerCoords);
      }
    }
    if (label == null) {
      // If no ground truth data add options to show the spots by their rank
      // and number of neighbours
      gd.addSlider("Top_N", 0, 100, settings.topN);
      topNScrollBar = gd.getLastScrollbar();
      gd.addSlider("Select", 0, 100, settings.select);
      selectScrollBar = gd.getLastScrollbar();
      gd.addSlider("Neigbour_radius", 0, 10, settings.neighbourRadius);
    }

    ImageListener imageListener = null;

    if (ImageJUtils.isShowGenericDialog()) {
      // Listen for changes in the dialog options
      gd.addOptionCollectedListener(event -> {
        // Just run on the current processor
        if (preview) {
          run(imp.getProcessor());
        }
      });
      // Listen for changes to an image
      imageListener = new ImageAdapter() {
        @Override
        public void imageUpdated(ImagePlus imp) {
          if (SpotFinderPreview.this.imp.getID() == imp.getID() && preview
              && imp.getCurrentSlice() != currentSlice && filter != null) {
            run(imp.getProcessor(), filter);
          }
        }
      };
      ImagePlus.addImageListener(imageListener);

      // Support template settings
      final Vector<TextField> numerics = gd.getNumericFields();
      final Vector<Choice> choices = gd.getChoices();

      final Iterator<TextField> nu = numerics.iterator();
      final Iterator<Choice> ch = choices.iterator();

      final Choice textTemplate = ch.next();
      textTemplate.removeItemListener(gd);
      textTemplate.removeKeyListener(gd);
      textTemplate.addItemListener(this::itemStateChanged);

      textCameraModelName = ch.next();
      textPsf = ch.next();
      textDataFilterType = ch.next();
      textDataFilterMethod = ch.next();
      textSmooth = nu.next();
      textDataFilterMethod2 = ch.next();
      textSmooth2 = nu.next();
      textSearch = nu.next();
      textBorder = nu.next();
    }

    gd.addPreviewCheckbox(pfr);
    gd.addDialogListener(this::dialogItemChanged);
    gd.setOKLabel("Save");
    gd.setCancelLabel("Close");
    gd.showDialog();

    if (imageListener != null) {
      ImagePlus.removeImageListener(imageListener);
    }

    if (!gd.wasCanceled() && !SettingsManager.writeSettings(config, SettingsManager.FLAG_SILENT)) {
      IJ.error(TITLE, "Failed to save settings");
    }

    // Reset
    imp.setOverlay(overlay);

    return DONE;
  }

  private boolean dialogItemChanged(GenericDialog gd, @SuppressWarnings("unused") AWTEvent event) {
    if (refreshing) {
      return false;
    }

    gd.getNextChoice(); // Ignore template

    // Set a camera model
    fitConfig.setCameraModelName(gd.getNextChoice());

    fitConfig.setPsfType(PeakFit.getPsfTypeValues()[gd.getNextChoiceIndex()]);

    config.setDataFilterType(gd.getNextChoiceIndex());
    config.setDataFilter(gd.getNextChoiceIndex(), Math.abs(gd.getNextNumber()), 0);
    config.setDataFilter(gd.getNextChoiceIndex(), Math.abs(gd.getNextNumber()), 1);
    config.setSearch(gd.getNextNumber());
    config.setBorder(gd.getNextNumber());
    if (label == null) {
      settings.topN = (int) gd.getNextNumber();
      settings.select = (int) gd.getNextNumber();
      settings.neighbourRadius = (int) gd.getNextNumber();
    } else {
      settings.distance = gd.getNextNumber();
      settings.lowerDistance = gd.getNextNumber();
      settings.multipleMatches = gd.getNextBoolean();
      settings.showTP = gd.getNextBoolean();
      settings.showFP = gd.getNextBoolean();
    }
    preview = gd.getNextBoolean();

    ((ExtendedGenericDialog) gd).collectOptions();

    final boolean result = !gd.invalidNumber();
    if (!preview) {
      setLabel("");
      this.imp.setOverlay(overlay);
    }
    // For astigmatism PSF.
    // TODO - See if this is slowing the preview down. If so only do if the PSF type changes.
    if (!PeakFit.configurePsfModel(config, PeakFit.FLAG_NO_SAVE)) {
      return false;
    }
    return result;
  }

  private void setLabel(String message) {
    if (label != null) {
      label.setText(message);
    }
  }

  @Override
  public void run(ImageProcessor ip) {
    if (refreshing) {
      return;
    }

    final Rectangle bounds = ip.getRoi();

    // Only do this if the settings changed
    final Calibration calibration = fitConfig.getCalibration();
    final FitEngineSettings fitEngineSettings = config.getFitEngineSettings();
    final PSF psf = fitConfig.getPsf();

    boolean newCameraModel = filter == null;
    if (!calibration.equals(lastCalibration)) {
      newCameraModel = true;
      // Set a camera model.
      // We have to set the camera type too to avoid configuration errors.
      CameraModel cameraModel = CameraModelManager.load(fitConfig.getCameraModelName());
      if (cameraModel == null) {
        cameraModel = new FakePerPixelCameraModel(0, 1, 1);
        fitConfig.setCameraType(CameraType.EMCCD);
      } else {
        fitConfig.setCameraType(CameraType.SCMOS);

        // Support cropped origin selection.
        final Rectangle sourceBounds = IJImageSource.getBounds(imp);
        cameraModel = PeakFit.cropCameraModel(cameraModel, sourceBounds, null, true);
        if (cameraModel == null) {
          gd.getPreviewCheckbox().setState(false);
          return;
        }
      }
      fitConfig.setCameraModel(cameraModel);
    }

    if (newCameraModel || !fitEngineSettings.equals(lastFitEngineSettings)
        || !psf.equals(lastPsf)) {
      // Configure a jury filter
      if (config.getDataFilterType() == DataFilterType.JURY
          && !PeakFit.configureDataFilter(config, PeakFit.FLAG_NO_SAVE)) {
        gd.getPreviewCheckbox().setState(false);
        return;
      }

      try {
        filter = config.createSpotFilter();
      } catch (final Exception ex) {
        filter = null;
        this.imp.setOverlay(overlay);
        // Required for ImageJ to disable the preview
        throw new IllegalStateException("Unable to create spot filter", ex);
      }
      ImageJUtils.log(filter.getDescription());
    }

    lastCalibration = calibration;
    lastFitEngineSettings = fitEngineSettings;
    lastPsf = psf;

    // This code can probably be removed since the crop is done above.
    if (fitConfig.getCameraTypeValue() == CameraType.SCMOS_VALUE) {
      // Instead just warn if the roi cannot be extracted from the selected model
      // or there is a mismatch
      final Rectangle modelBounds = fitConfig.getCameraModel().getBounds();
      if (modelBounds != null) {
        if (!modelBounds.contains(bounds)) {
          //@formatter:off
          ImageJUtils.log(
              "WARNING: Camera model bounds [x=%d,y=%d,width=%d,height=%d]"
                  + " does not contain image target bounds [x=%d,y=%d,width=%d,height=%d]",
              modelBounds.x, modelBounds.y, modelBounds.width, modelBounds.height,
              bounds.x, bounds.y, bounds.width, bounds.height
          );
          //@formatter:on

          // Warn if the model bounds are mismatched than the image as this may be an incorrect
          // selection for the camera model
        } else if (modelBounds.x != 0 || modelBounds.y != 0 || modelBounds.width > ip.getWidth()
            || modelBounds.height > ip.getHeight()) {
          //@formatter:off
          ImageJUtils.log(
              "WARNING: Probably an incorrect camera model!\n"
                  + "Model bounds [x=%d,y=%d,width=%d,height=%d]\n"
                  + "do not match the image target bounds [width=%d,height=%d].",
              modelBounds.x, modelBounds.y, modelBounds.width, modelBounds.height,
              ip.getWidth(),  ip.getHeight()
          );
          //@formatter:on
        }
      }
    }

    run(ip, filter);
  }

  private void run(ImageProcessor ip, MaximaSpotFilter filter) {
    if (refreshing) {
      return;
    }

    currentSlice = imp.getCurrentSlice();

    final Rectangle bounds = ip.getRoi();

    // Crop to the ROI
    FloatProcessor fp = ip.crop().toFloat(0, null);

    float[] data = (float[]) fp.getPixels();

    final int width = fp.getWidth();
    final int height = fp.getHeight();

    // Store the mean bias and gain of the region data.
    // This is used to correctly overlay the filtered data on the original image.
    double bias = 0;
    double gain = 1;
    boolean adjust = false;

    // Set weights
    final CameraModel cameraModel = fitConfig.getCameraModel();
    if (!(cameraModel instanceof FakePerPixelCameraModel)) {
      // This should be done on the normalised data
      final float[] w = cameraModel.getNormalisedWeights(bounds);
      filter.setWeights(w, width, height);
      data = data.clone();
      if (data.length < ip.getPixelCount()) {
        adjust = true;
        bias = MathUtils.sum(cameraModel.getBias(bounds)) / data.length;
        gain = MathUtils.sum(cameraModel.getGain(bounds)) / data.length;
      }
      cameraModel.removeBiasAndGain(bounds, data);
    }

    final Spot[] spots = filter.rank(data, width, height);
    data = filter.getPreprocessedData();

    final int size = spots.length;
    if (topNScrollBar != null) {
      topNScrollBar.setMaximum(size);
      selectScrollBar.setMaximum(size);
    }

    fp = new FloatProcessor(width, height, data);
    final FloatProcessor out = new FloatProcessor(ip.getWidth(), ip.getHeight());
    out.copyBits(ip, 0, 0, Blitter.COPY);
    if (adjust) {
      fp.multiply(gain);
      fp.add(bias);
    }
    out.insert(fp, bounds.x, bounds.y);
    final double min = fp.getMin();
    final double max = fp.getMax();
    out.setMinAndMax(min, max);

    final Overlay o = new Overlay();
    o.add(new ImageRoi(0, 0, out));

    if (label != null) {
      // Get results for frame
      final Coordinate[] actual =
          ResultsMatchCalculator.getCoordinates(actualCoordinates, imp.getCurrentSlice());

      final Coordinate[] predicted = new Coordinate[size];
      for (int i = 0; i < size; i++) {
        predicted[i] = new BasePoint(spots[i].x + bounds.x, spots[i].y + bounds.y);
      }

      // Compute assignments
      final LocalList<FractionalAssignment> fractionalAssignments =
          new LocalList<>(3 * predicted.length);
      final double matchDistance = settings.distance * fitConfig.getInitialPeakStdDev();
      final RampedScore score =
          RampedScore.of(matchDistance, matchDistance * settings.lowerDistance / 100, false);
      final double dmin = matchDistance * matchDistance;
      final int nActual = actual.length;
      final int nPredicted = predicted.length;
      for (int j = 0; j < nPredicted; j++) {
        // Centre in the middle of the pixel
        final float x = predicted[j].getX() + 0.5f;
        final float y = predicted[j].getY() + 0.5f;
        // Any spots that match
        for (int i = 0; i < nActual; i++) {
          final double dx = (x - actual[i].getX());
          final double dy = (y - actual[i].getY());
          final double d2 = dx * dx + dy * dy;
          if (d2 <= dmin) {
            final double d = Math.sqrt(d2);
            final double s = score.score(d);

            if (s == 0) {
              continue;
            }

            double distance = 1 - s;
            if (distance == 0) {
              // In the case of a match below the distance thresholds
              // the distance will be 0. To distinguish between candidates all below
              // the thresholds just take the closest.
              // We know d2 is below dmin so we subtract the delta.
              distance -= (dmin - d2);
            }

            // Store the match
            fractionalAssignments.add(new ImmutableFractionalAssignment(i, j, distance, s));
          }
        }
      }

      final FractionalAssignment[] assignments =
          fractionalAssignments.toArray(new FractionalAssignment[0]);

      // Compute matches
      final RankedScoreCalculator calc =
          RankedScoreCalculator.create(assignments, nActual - 1, nPredicted - 1);
      final boolean save = settings.showTP || settings.showFP;
      final double[] calcScore = calc.score(nPredicted, settings.multipleMatches, save);
      final ClassificationResult result =
          RankedScoreCalculator.toClassificationResult(calcScore, nActual);

      // Compute AUC and max jaccard (and plot)
      final double[][] curve =
          RankedScoreCalculator.getPrecisionRecallCurve(assignments, nActual, nPredicted);
      final double[] precision = curve[0];
      final double[] recall = curve[1];
      final double[] jaccard = curve[2];
      final double auc = AucCalculator.auc(precision, recall);

      // Show scores
      final String scoreLabel = String.format("Slice=%d, AUC=%s, R=%s, Max J=%s",
          imp.getCurrentSlice(), MathUtils.rounded(auc), MathUtils.rounded(result.getRecall()),
          MathUtils.rounded(MathUtils.maxDefault(0, jaccard)));
      setLabel(scoreLabel);

      // Plot
      String title = TITLE + " Performance";
      Plot plot = new Plot(title, "Spot Rank", "");
      final double[] rank = SimpleArrayUtils.newArray(precision.length, 0, 1.0);
      plot.setLimits(0, nPredicted, 0, 1.05);
      plot.setColor(Color.blue);
      plot.addPoints(rank, precision, Plot.LINE);
      plot.setColor(Color.red);
      plot.addPoints(rank, recall, Plot.LINE);
      plot.setColor(Color.black);
      plot.addPoints(rank, jaccard, Plot.LINE);
      plot.setColor(Color.black);
      plot.addLabel(0, 0, scoreLabel);

      final WindowOrganiser windowOrganiser = new WindowOrganiser();
      ImageJUtils.display(title, plot, 0, windowOrganiser);

      title = TITLE + " Precision-Recall";
      plot = new Plot(title, "Recall", "Precision");
      plot.setLimits(0, 1, 0, 1.05);
      plot.setColor(Color.red);
      plot.addPoints(recall, precision, Plot.LINE);
      plot.drawLine(recall[recall.length - 1], precision[recall.length - 1],
          recall[recall.length - 1], 0);
      plot.setColor(Color.black);
      plot.addLabel(0, 0, scoreLabel);
      ImageJUtils.display(title, plot, 0, windowOrganiser);

      windowOrganiser.tile();

      // Create Rois for TP and FP
      if (save) {
        final double[] matchScore =
            RankedScoreCalculator.getMatchScore(calc.getScoredAssignments(), nPredicted);
        int matches = 0;
        for (int i = 0; i < matchScore.length; i++) {
          if (matchScore[i] != 0) {
            matches++;
          }
        }
        if (settings.showTP) {
          final float[] x = new float[matches];
          final float[] y = new float[x.length];
          int count = 0;
          for (int i = 0; i < matchScore.length; i++) {
            if (matchScore[i] != 0) {
              final BasePoint p = (BasePoint) predicted[i];
              x[count] = p.getX() + 0.5f;
              y[count] = p.getY() + 0.5f;
              count++;
            }
          }
          addRoi(0, o, x, y, count, Color.green);
        }
        if (settings.showFP) {
          final float[] x = new float[nPredicted - matches];
          final float[] y = new float[x.length];
          int count = 0;
          for (int i = 0; i < matchScore.length; i++) {
            if (matchScore[i] == 0) {
              final BasePoint p = (BasePoint) predicted[i];
              x[count] = p.getX() + 0.5f;
              y[count] = p.getY() + 0.5f;
              count++;
            }
          }
          addRoi(0, o, x, y, count, Color.red);
        }
      }
    } else {
      final WindowOrganiser wo = new WindowOrganiser();

      // Option to show the number of neighbours within a set pixel box radius
      final int[] count =
          spotFilterHelper.countNeighbours(spots, width, height, settings.neighbourRadius);

      // Show as histogram the totals...
      new HistogramPlotBuilder(TITLE, StoredData.create(count), "Neighbours").setIntegerBins(true)
          .setPlotLabel("Radius = " + settings.neighbourRadius).show(wo);

      // TODO - Draw n=0, n=1 on the image overlay

      final LUT lut = LutHelper.createLut(LutColour.FIRE_LIGHT);
      // These are copied by the ROI
      final float[] x = new float[1];
      final float[] y = new float[1];
      // Plot the intensity
      final double[] intensity = new double[size];
      final double[] rank = SimpleArrayUtils.newArray(size, 1, 1.0);
      final int top = (settings.topN > 0) ? settings.topN : size;
      final int size_1 = size - 1;
      for (int i = 0; i < size; i++) {
        intensity[i] = spots[i].intensity;
        if (i < top) {
          x[0] = spots[i].x + bounds.x + 0.5f;
          y[0] = spots[i].y + bounds.y + 0.5f;
          final Color c = LutHelper.getColour(lut, size_1 - i, size);
          addRoi(0, o, x, y, 1, c, 2, 1);
        }
      }

      final String title = TITLE + " Intensity";
      final Plot plot = new Plot(title, "Rank", "Intensity");
      plot.setColor(Color.blue);
      plot.addPoints(rank, intensity, Plot.LINE);
      if (settings.topN > 0 && settings.topN < size) {
        plot.setColor(Color.magenta);
        plot.drawLine(settings.topN, 0, settings.topN, intensity[settings.topN - 1]);
      }
      if (settings.select > 0 && settings.select < size) {
        plot.setColor(Color.yellow);
        final int index = settings.select - 1;
        final double in = intensity[index];
        plot.drawLine(settings.select, 0, settings.select, in);
        x[0] = spots[index].x + bounds.x + 0.5f;
        y[0] = spots[index].y + bounds.y + 0.5f;
        final Color c = LutHelper.getColour(lut, size_1 - settings.select, size);
        addRoi(0, o, x, y, 1, c, 3, 3);
        plot.setColor(Color.black);
        plot.addLabel(0, 0, "Selected spot intensity = " + MathUtils.rounded(in));
      }

      ImageJUtils.display(title, plot, 0, wo);

      wo.tile();
    }

    imp.setOverlay(o);
  }

  /**
   * Adds the roi to the overlay as a circle.
   *
   * @param frame the frame
   * @param overlay the overlay
   * @param x the x coordinates
   * @param y the y coordinates
   * @param npoints the number of points
   * @param colour the colour
   * @see PointRoi#setPosition(int)
   * @see PointRoi#setFillColor(Color)
   * @see PointRoi#setStrokeColor(Color)
   */
  public static void addRoi(int frame, Overlay overlay, float[] x, float[] y, int npoints,
      Color colour) {
    // Add as a circle
    addRoi(frame, overlay, x, y, npoints, colour, 3);
  }

  /**
   * Adds the roi to the overlay.
   *
   * @param frame the frame
   * @param overlay the overlay
   * @param x the x coordinates
   * @param y the y coordinates
   * @param npoints the number of points
   * @param colour the colour
   * @param pointType the point type
   * @see PointRoi#setPosition(int)
   * @see PointRoi#setPointType(int)
   * @see PointRoi#setFillColor(Color)
   * @see PointRoi#setStrokeColor(Color)
   */
  public static void addRoi(int frame, Overlay overlay, float[] x, float[] y, int npoints,
      Color colour, int pointType) {
    addRoi(frame, overlay, x, y, npoints, colour, pointType, 0);
  }

  /**
   * Adds the roi to the overlay.
   *
   * @param frame the frame
   * @param overlay the overlay
   * @param x the x coordinates
   * @param y the y coordinates
   * @param npoints the number of points
   * @param colour the colour
   * @param pointType the point type
   * @param size the size
   * @see PointRoi#setPosition(int)
   * @see PointRoi#setPointType(int)
   * @see PointRoi#setFillColor(Color)
   * @see PointRoi#setStrokeColor(Color)
   * @see PointRoi#setSize(int)
   */
  public static void addRoi(int frame, Overlay overlay, float[] x, float[] y, int npoints,
      Color colour, int pointType, int size) {
    if (npoints == 0) {
      return;
    }
    final PointRoi roi = new OffsetPointRoi(x, y, npoints);
    roi.setPointType(pointType);
    roi.setFillColor(colour);
    roi.setStrokeColor(colour);
    if (frame != 0) {
      roi.setPosition(frame);
    }
    if (size != 0) {
      roi.setSize(size);
    }
    overlay.add(roi);
  }

  @Override
  public void setNPasses(int passes) {
    // Nothing to do
  }

  private void itemStateChanged(ItemEvent event) {
    if (event.getSource() instanceof Choice) {
      // Update the settings from the template
      final Choice choice = (Choice) event.getSource();
      final String templateName = choice.getSelectedItem();

      // Get the configuration template
      final TemplateSettings template = ConfigurationTemplate.getTemplate(templateName);

      if (template != null) {
        refreshing = true;

        IJ.log("Applying template: " + templateName);

        for (final String note : template.getNotesList()) {
          IJ.log(note);
        }

        final boolean custom = ConfigurationTemplate.isCustomTemplate(templateName);
        if (template.hasPsf()) {
          refreshSettings(template.getPsf(), custom);
        }
        if (template.hasFitEngineSettings()) {
          refreshSettings(template.getFitEngineSettings(), custom);
        }

        refreshing = false;
      }
    }
  }

  private void refreshSettings(PSF psf, boolean isCustomTemplate) {
    if (!isCustomTemplate || psf == null) {
      return;
    }

    // Do not use set() as we support merging a partial PSF
    fitConfig.mergePsf(psf);

    textPsf.select(PsfProtosHelper.getName(fitConfig.getPsfType()));
  }

  /**
   * Refresh settings.
   *
   * <p>If this is a custom template then use all the settings. If a default template then leave
   * some existing spot settings untouched as the user may have updated them (e.g. PSF width).
   *
   * @param fitEngineSettings the config
   * @param isCustomTemplate True if a custom template.
   */
  private void refreshSettings(FitEngineSettings fitEngineSettings, boolean isCustomTemplate) {
    // Set the configuration
    // This will clear everything and merge the configuration so
    // remove the fit settings (as we do not care about those).

    this.config.setFitEngineSettings(fitEngineSettings.toBuilder().clearFitSettings().build());
    fitConfig = this.config.getFitConfiguration();

    textCameraModelName.select(fitConfig.getCameraModelName());
    textDataFilterType
        .select(SettingsManager.getDataFilterTypeNames()[config.getDataFilterType().ordinal()]);
    textDataFilterMethod.select(
        SettingsManager.getDataFilterMethodNames()[config.getDataFilterMethod(0).ordinal()]);
    textSmooth.setText("" + config.getDataFilterParameterValue(0));
    if (config.getDataFiltersCount() > 1) {
      textDataFilterMethod2.select(SettingsManager.getDataFilterMethodNames()[config
          .getDataFilterMethod(1, settings.defaultDataFilterMethod).ordinal()]);
      textSmooth2.setText("" + config.getDataFilterParameterValue(1, settings.defaultSmooth));
      // XXX - What about the Absolute/Relative flag?
    }
    textSearch.setText("" + config.getSearch());
    textBorder.setText("" + config.getBorder());
  }
}
