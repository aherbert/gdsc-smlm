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
import ij.gui.Line;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import java.awt.Label;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicReference;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.NonBlockingExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.SoftLock;
import uk.ac.sussex.gdsc.core.utils.SortUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.StoredData;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.NamedObject;
import uk.ac.sussex.gdsc.smlm.data.config.ConfigurationException;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;
import uk.ac.sussex.gdsc.smlm.fitting.JumpDistanceAnalysis;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.PrecisionResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.sort.IdFramePeakResultComparator;

/**
 * Analyses the track lengths of traced data.
 */
public class FilterMolecules implements PlugIn {
  private static final String TITLE = "Filter Molecules";

  /** The plugin settings. */
  Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());

    String inputOption;
    FilterMode filterMode;
    boolean removeSingles;
    String outputName;
    boolean outputSuffix;

    // FilterDiffusionCoefficient
    double lowerDThreshold;
    double upperDThreshold;
    int minimumLength;

    Settings() {
      // Set defaults
      inputOption = "";
      filterMode = FilterMode.D;
      removeSingles = true;
      outputName = "Filtered";
      outputSuffix = true;
    }

    Settings(Settings source) {
      inputOption = source.inputOption;
      filterMode = source.filterMode;
      removeSingles = source.removeSingles;
      outputName = source.outputName;
      outputSuffix = source.outputSuffix;
      lowerDThreshold = source.lowerDThreshold;
      upperDThreshold = source.upperDThreshold;
      minimumLength = source.minimumLength;
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

  /**
   * The filter mode.
   */
  private enum FilterMode implements NamedObject {
    D("Diffusion Coefficient", "D");

    private final String name;
    private final String shortName;

    FilterMode(String name) {
      this(name, name);
    }

    FilterMode(String name, String shortName) {
      this.name = name;
      this.shortName = shortName;
    }

    @Override
    public String getName() {
      return name;
    }

    @Override
    public String getShortName() {
      return shortName;
    }

    /**
     * Get the filter mode for the given ordinal. The input ordinal is clipped to the allowed range
     * of the enum.
     *
     * @param ordinal the ordinal
     * @return the filter mode
     */
    static FilterMode forOrdinal(int ordinal) {
      final FilterMode[] values = values();
      return values[MathUtils.clip(0, values.length - 1, ordinal)];
    }
  }

  /**
   * Class to filter using the diffusion coefficient.
   */
  private class FilterDiffusionCoefficient {
    private TypeConverter<DistanceUnit> distanceConverter;
    private TypeConverter<TimeUnit> timeConverter;
    private final SoftLock lock = new SoftLock();
    private double lastLowerMsdThreshold = -1;
    private double lastUpperMsdThreshold = -1;
    private int lastMinimumLength = -1;
    private double error;
    private double[] msds; // MSD of trace
    private int[] lengths; // Length of trace
    private int[] ids; // trace id

    private IntOpenHashSet filtered;
    private Label filteredLabel;
    private double maxDiffusionCoefficientHistogram;
    private double maxLengthHistogram;

    private int lastid = -1;
    private float lastx;
    private float lasty;
    private int startFrame;
    private int lastFrame;
    private int totalJump;
    private double sumSquared;
    private final DoubleArrayList msdList = new DoubleArrayList();
    private final IntArrayList lengthList = new IntArrayList();
    private final IntArrayList idList = new IntArrayList();

    /**
     * Filter the molecules using the diffusion coefficient.
     *
     * @param results the results
     */
    void run(MemoryPeakResults results) {
      try {
        distanceConverter = results.getDistanceConverter(DistanceUnit.UM);
        timeConverter = results.getTimeConverter(TimeUnit.SECOND);
      } catch (final ConversionException | ConfigurationException ex) {
        IJ.error(TITLE, "Cannot convert units to um or seconds: " + ex.getMessage());
        return;
      }

      // Get the localisation error (4s^2) in raw units^2
      double precision = 0;
      try {
        final PrecisionResultProcedure p = new PrecisionResultProcedure(results);
        p.getPrecision();

        // Precision in nm using the median
        precision = new Percentile().evaluate(p.precisions, 50);
        // Convert from nm to um to raw units
        final double rawPrecision = distanceConverter.convertBack(precision / 1e3);
        // Get the localisation error (4s^2) in units^2
        error = 4 * rawPrecision * rawPrecision;
      } catch (final MathIllegalArgumentException | DataException ex) {
        ImageJUtils.log(TITLE + " - Unable to compute precision: " + ex.getMessage());
      }

      // Analyse the track lengths
      results.sort(IdFramePeakResultComparator.INSTANCE);
      // Ensure the first result triggers an id change
      lastid = results.getFirst().getId() - 1;
      results.forEach(this::processTrackLength);
      store(); // For the final track

      // Create histogram of D and lengths
      // Add sliders for the thresholds.
      // Dynamically put ROI onto the histograms for the thresholds.

      // When the thresholds are adjusted recompute the size of the output dataset
      // and show in the dialog.
      // When click save extract all the results that match the ID filter.

      msds = msdList.toDoubleArray();
      lengths = lengthList.toIntArray();
      ids = idList.toIntArray();

      // Sort by MSD
      final int[] indices = SimpleArrayUtils.natural(msds.length);
      SortUtils.sortIndices(indices, msds, false);
      final double[] msds2 = msds.clone();
      final int[] lengths2 = lengths.clone();
      final int[] ids2 = ids.clone();
      for (int i = 0; i < indices.length; i++) {
        msds[i] = msds2[indices[i]];
        lengths[i] = lengths2[indices[i]];
        ids[i] = ids2[indices[i]];
      }

      final WindowOrganiser wo = new WindowOrganiser();

      // Histogram the diffusion coefficients
      final Statistics s = Statistics.create(msds);
      final double diffusionCoefficientMean = s.getMean();
      final double diffusionCoefficientMax = MathUtils.max(msds);
      final String msg = String.format("Average D per track = %s um^2/s. Max = %s um^2/s",
          MathUtils.rounded(diffusionCoefficientMean), MathUtils.rounded(diffusionCoefficientMax));
      final HistogramPlot diffusionCoefficientHistogramPlot =
          new HistogramPlotBuilder("Trace Diffusion Coefficient", StoredData.create(msds),
              "D (um^2/s)").setRemoveOutliersOption(1).setPlotLabel(msg).build();
      diffusionCoefficientHistogramPlot.show(wo);
      final Plot diffusionCoefficientPlot = diffusionCoefficientHistogramPlot.getPlot();

      // Histogram the lengths
      final HistogramPlot lengthHistogramPlot =
          new HistogramPlotBuilder("Trace Jumps", StoredData.create(lengths), "Length").build();
      lengthHistogramPlot.show(wo);
      final Plot lengthPlot = lengthHistogramPlot.getPlot();

      // Interactive analysis
      final NonBlockingExtendedGenericDialog gd = new NonBlockingExtendedGenericDialog(TITLE);
      ImageJUtils.addMessage(gd,
          "Filter %s using the track diffusion coefficient (D).\n"
              + "Localisation error has been subtracted from jumps (%s nm).",
          TextUtils.pleural(ids.length, "molecule"), MathUtils.rounded(precision));
      gd.addMessage(msg);

      // Build a nice slider range from the histogram limits
      double[] xvalues = diffusionCoefficientHistogramPlot.getPlotXValues();
      double min = xvalues[0];
      double max = xvalues[xvalues.length - 1];
      maxDiffusionCoefficientHistogram = diffusionCoefficientHistogramPlot.getPlotMaxY();
      if (max - min < 5) {
        // Because sliders are used when the range is <5 and floating point
        gd.addSlider("Lower_D_threshold", min, max, settings.lowerDThreshold);
        gd.addSlider("Upper_D_threshold", min, max, settings.upperDThreshold);
      } else {
        gd.addNumericField("Lower_D_threshold", settings.lowerDThreshold, 2, 6, "um^2/s");
        gd.addNumericField("Upper_D_threshold", settings.upperDThreshold, 2, 6, "um^2/s");
      }
      xvalues = lengthHistogramPlot.getPlotXValues();
      min = xvalues[0];
      max = xvalues[xvalues.length - 1];
      maxLengthHistogram = lengthHistogramPlot.getPlotMaxY();
      gd.addSlider("Minimum_length", min, max, settings.minimumLength, 1);

      // Compute initial filtered set.
      filtered = new IntOpenHashSet();
      // Add label to dialog containing the set size.
      gd.addMessage("");
      filteredLabel = (Label) gd.getMessage();
      // Add label with the output results set name.
      gd.addMessage("Output name: " + createResultsName(results, settings));

      filter();
      draw(diffusionCoefficientPlot, lengthPlot);
      wo.tile();

      // Add an interactive dialog listener to dynamically recompute the filtered set
      Runnable callback = () -> {
        /* do nothing */ };
      if (ImageJUtils.isShowGenericDialog()) {
        final ExecutorService es = Executors.newSingleThreadExecutor();
        callback = () -> es.shutdown();
        gd.addDialogListener((gd1, event) -> {
          settings.lowerDThreshold = gd1.getNextNumber();
          settings.upperDThreshold = gd1.getNextNumber();
          if (settings.upperDThreshold < settings.lowerDThreshold) {
            settings.upperDThreshold = diffusionCoefficientMax;
          }
          settings.minimumLength = (int) gd1.getNextNumber();
          update(es, diffusionCoefficientPlot, lengthPlot);
          return true;
        });
      }
      gd.setOKLabel("Save dataset");
      gd.setCancelLabel("Close");
      gd.addHelp(HelpUrls.getUrl("filter-molecules"));
      gd.showDialog();

      callback.run();
      if (gd.wasCanceled()) {
        return;
      }

      createResults(results, settings);
    }

    private void processTrackLength(PeakResult peakResult) {
      final int id = peakResult.getId();
      final int frame = peakResult.getFrame();
      final float x = peakResult.getXPosition();
      final float y = peakResult.getYPosition();
      if (lastid == id) {
        // Compute the jump
        final int jump = frame - lastFrame;
        // Get the raw distance but subtract the expected localisation error
        final double d2 = Math.max(0, MathUtils.distance2(lastx, lasty, x, y) - error);
        // We expect the Mean Squared Distance (MSD) to scale linearly
        // with time so just weight each jump by the time gap.
        // However we apply a correction factor for diffusion with frames.
        sumSquared += JumpDistanceAnalysis.convertObservedToActual(d2, jump);
        totalJump += jump;
      } else {
        // New Id
        store();
        lastid = id;
        startFrame = frame;
        sumSquared = 0;
        totalJump = 0;
      }
      lastFrame = frame;
      lastx = x;
      lasty = y;
    }

    private void store() {
      if (lastid == 0 || totalJump == 0) {
        return;
      }
      final int len = lastFrame - startFrame;
      final double msd = sumSquared / totalJump;
      // 4D = MSD => D = MSD / 4
      // Mean squared distance is in raw units squared.
      // Convert twice since this is a squared distance then from frames to seconds.
      // Use the convertBack() since this is a divide in the final units: um^2/s
      msdList
          .add(timeConverter.convertBack(distanceConverter.convert(distanceConverter.convert(msd)))
              / 4.0);
      lengthList.add(len);
      idList.add(lastid);
    }

    private void filter() {
      lastLowerMsdThreshold = settings.lowerDThreshold;
      lastUpperMsdThreshold = settings.upperDThreshold;
      lastMinimumLength = settings.minimumLength;
      filtered.clear();
      // Find the index in the MSD array
      int lowerIndex = Arrays.binarySearch(msds, lastLowerMsdThreshold);
      if (lowerIndex < 0) {
        lowerIndex = -(lowerIndex + 1);
      } else {
        // Ensure the index is at the first point (inclusive lower)
        while (lowerIndex > 0 && msds[lowerIndex - 1] >= lastLowerMsdThreshold) {
          lowerIndex--;
        }
      }
      int upperIndex = Arrays.binarySearch(msds, lastUpperMsdThreshold);
      if (upperIndex < 0) {
        upperIndex = -(upperIndex + 1);
      } else {
        // Ensure the index is beyond the last point (exclusive upper)
        while (upperIndex < msds.length && msds[upperIndex] <= lastUpperMsdThreshold) {
          upperIndex++;
        }
      }
      // Filter by length
      final int minLength = lastMinimumLength;
      for (int i = lowerIndex; i < upperIndex; i++) {
        if (lengths[i] >= minLength) {
          filtered.add(ids[i]);
        }
      }
    }

    private void draw(Plot diffusionCoefficientPlot, Plot lengthPlot) {
      filteredLabel.setText(TextUtils.pleural(filtered.size(), "molecule"));
      // Map the thresholds to the plot area.
      double lx = diffusionCoefficientPlot.scaleXtoPxl(lastLowerMsdThreshold);
      final double ux = diffusionCoefficientPlot.scaleXtoPxl(lastUpperMsdThreshold);
      double uy = diffusionCoefficientPlot.scaleYtoPxl(maxDiffusionCoefficientHistogram);
      double ly = diffusionCoefficientPlot.scaleYtoPxl(0);
      diffusionCoefficientPlot.getImagePlus().setRoi(new Roi(lx, uy, ux - lx, ly - uy));
      lx = lengthPlot.scaleXtoPxl(lastMinimumLength);
      uy = lengthPlot.scaleYtoPxl(maxLengthHistogram);
      ly = lengthPlot.scaleYtoPxl(0);
      lengthPlot.getImagePlus().setRoi(new Line(lx, uy, lx, ly));
    }

    void update(ExecutorService es, Plot diffusionCoefficientPlot, Plot lengthPlot) {
      if (lock.acquire()) {
        // Run in a new thread to allow the GUI to continue updating
        es.submit(() -> {
          try {
            // Continue while the parameter is changing
            while (lastLowerMsdThreshold != settings.lowerDThreshold
                || lastUpperMsdThreshold != settings.upperDThreshold
                || lastMinimumLength != settings.minimumLength) {
              filter();
              draw(diffusionCoefficientPlot, lengthPlot);
            }
          } finally {
            // Ensure the running flag is reset
            lock.release();
          }
        });
      }
    }

    private void createResults(MemoryPeakResults results, Settings settings) {
      final MemoryPeakResults out = FilterMolecules.createResults(results, settings);
      results.forEach((PeakResultProcedure) r -> {
        if (filtered.contains(r.getId())) {
          out.add(r);
        }
      });
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "No localisations in memory");
      return;
    }

    if (!showDialog()) {
      return;
    }

    // Load the results
    MemoryPeakResults results =
        ResultsManager.loadInputResults(settings.inputOption, false, null, null);
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(TITLE, "No results could be loaded");
      return;
    }

    // Allow reordering when filtering
    results = results.copy();

    if (settings.removeSingles) {
      results.removeIf(p -> p.getId() <= 0);
      final Int2IntOpenHashMap count = new Int2IntOpenHashMap(results.size());
      results.forEach((PeakResultProcedure) r -> count.addTo(r.getId(), 1));
      results.removeIf(p -> count.get(p.getId()) == 1);
      if (results.isEmpty()) {
        IJ.error(TITLE, "No results after filtering singles");
        return;
      }
    }

    if (settings.filterMode == FilterMode.D) {
      new FilterDiffusionCoefficient().run(results);
    } else {
      IJ.error(TITLE, "Unknown filter mode: " + settings.filterMode);
    }
  }

  private boolean showDialog() {
    settings = Settings.load();
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Filter molecules");
    ResultsManager.addInput(gd, "Input", settings.inputOption, InputSource.MEMORY_CLUSTERED);
    gd.addChoice("Filter_mode", SettingsManager.getNames((Object[]) FilterMode.values()),
        settings.filterMode.ordinal());
    gd.addCheckbox("Remove_singles", settings.removeSingles);
    gd.addStringField("Output_name", settings.outputName, 50);
    gd.addCheckbox("Output_suffix", settings.outputSuffix);
    gd.addHelp(HelpUrls.getUrl("filter-molecules"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    settings.inputOption = ResultsManager.getInputSource(gd);
    settings.filterMode = FilterMode.forOrdinal(gd.getNextChoiceIndex());
    settings.removeSingles = gd.getNextBoolean();
    settings.outputName = gd.getNextString();
    settings.outputSuffix = gd.getNextBoolean();
    settings.save();
    return true;
  }

  /**
   * Creates the results and adds them to memory.
   *
   * @param sourceResults the source results
   * @param settings the settings
   * @return the memory peak results
   */
  static MemoryPeakResults createResults(MemoryPeakResults sourceResults, Settings settings) {
    final MemoryPeakResults out = new MemoryPeakResults();
    out.copySettings(sourceResults);
    out.setName(createResultsName(sourceResults, settings));
    MemoryPeakResults.addResults(out);
    return out;
  }

  // TODO - improve this to have suffix or full name options.
  // Update the docs for the plugin.

  /**
   * Creates the results name.
   *
   * @param sourceResults the source results
   * @param settings the settings
   * @return the string
   */
  static String createResultsName(MemoryPeakResults sourceResults, Settings settings) {
    if (settings.outputSuffix) {
      return sourceResults.getName() + " " + settings.outputName.trim();
    }
    return settings.outputName;
  }
}
