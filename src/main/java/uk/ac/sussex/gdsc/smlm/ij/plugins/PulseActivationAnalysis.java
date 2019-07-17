/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

import uk.ac.sussex.gdsc.core.clustering.DensityCounter;
import uk.ac.sussex.gdsc.core.clustering.DensityCounter.Molecule;
import uk.ac.sussex.gdsc.core.data.NotImplementedException;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.NonBlockingExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.Plot2;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.TurboList;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrencyUtils;
import uk.ac.sussex.gdsc.core.utils.rng.RandomUtils;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.core.utils.rng.SplitMix;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationHelper;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageMode;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageType;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsSettings;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.results.ImageJImagePeakResults;
import uk.ac.sussex.gdsc.smlm.ij.results.ImagePeakResultsFactory;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.Cluster.CentroidMethod;
import uk.ac.sussex.gdsc.smlm.results.IdPeakResult;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResultStoreList;
import uk.ac.sussex.gdsc.smlm.results.PeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResultsList;
import uk.ac.sussex.gdsc.smlm.results.Trace;
import uk.ac.sussex.gdsc.smlm.results.TraceManager;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyrResultProcedure;

import ij.CompositeImage;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.gui.ShapeRoi;
import ij.plugin.PlugIn;
import ij.plugin.frame.Recorder;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.UnitSphereSampler;
import org.apache.commons.rng.sampling.distribution.ContinuousUniformSampler;
import org.apache.commons.rng.sampling.distribution.DiscreteSampler;
import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import org.apache.commons.rng.sampling.distribution.ZigguratNormalizedGaussianSampler;
import org.apache.commons.rng.simple.RandomSource;

import java.awt.AWTEvent;
import java.awt.Checkbox;
import java.awt.Choice;
import java.awt.Color;
import java.awt.Rectangle;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Perform multi-channel super-resolution imaging by means of photo-switchable probes and pulsed
 * light activation.
 *
 * <p>This plugin is based on the methods described in: Mark Bates, Bo Huang, Graham T. Dempsey,
 * Xiaowei Zhuang (2007). Multicolor Super-Resolution Imaging with Photo-Switchable Fluorescent
 * Probes. Science 317, 1749. DOI: 10.1126/science.1146598.
 */
public class PulseActivationAnalysis implements PlugIn {
  private String title = "Activation Analysis";

  private static final Correction[] specificCorrection;
  private static final Correction[] nonSpecificCorrection;

  static {
    final EnumSet<Correction> correction = EnumSet.allOf(Correction.class);
    specificCorrection = correction.toArray(new Correction[correction.size()]);
    correction.remove(Correction.SUBTRACTION);
    nonSpecificCorrection = correction.toArray(new Correction[correction.size()]);
  }

  private static String inputOption = "";
  // Note: Set defaults to work with the 3-channel simulation
  private static int channels = 3;
  private static final int MAX_CHANNELS = 3;

  private static int repeatInterval = 30;
  private static final int[] startFrame = {1, 11, 21};
  /** The crosstalk between channels. */
  private static final double[] ct = new double[6];
  private static final String[] ctNames = {"21", "31", "12", "32", "13", "23"};
  private static final int C21 = 0;
  private static final int C31 = 1;
  private static final int C12 = 2;
  private static final int C32 = 3;
  private static final int C13 = 4;
  private static final int C23 = 5;
  private static int darkFramesForNewActivation = 1;

  private static int targetChannel = 1;

  private static double densityRadius = 35;
  private static int minNeighbours = 5;
  private static int specificCorrectionIndex = Correction.SUBTRACTION.ordinal();
  private static double[] specificCorrectionCutoff = {50, 50, 50};
  private static int nonSpecificCorrectionIndex;
  private static double nonSpecificCorrectionCutoff = 50;

  // Simulation settings
  private UniformRandomProvider rng;

  private UniformRandomProvider getUniformRandomProvider() {
    if (rng == null) {
      rng = RandomSource.create(RandomSource.XOR_SHIFT_1024_S);
    }
    return rng;
  }

  private static int[] sim_numberOfMolecules = {1000, 1000, 1000};
  private static SimulationDistribution[] sim_distribution =
      {SimulationDistribution.CIRCLE, SimulationDistribution.LINE, SimulationDistribution.POINT};
  private static double[] sim_precision = {15, 15, 15}; // nm
  private static int sim_cycles = 1000;
  private static int sim_size = 256;
  private static double sim_nmPerPixel = 100;
  private static double sim_activationDensity = 0.1; // molecules/micrometer
  private static double sim_nonSpecificFrequency = 0.01;

  // private ResultsSettings resultsSettings;
  private ResultsSettings.Builder resultsSettingsBuilder;
  private MemoryPeakResults results;
  private Trace[] traces;

  // The output. Used for the loop functionality
  private PeakResultsList[] output;
  private static final Color[] colors = new Color[] {Color.RED, Color.GREEN, Color.BLUE};
  private static final String[] magnifications;

  static {
    final ArrayList<String> list = new ArrayList<>();
    for (int i = 1; i <= 256; i *= 2) {
      list.add(Integer.toString(i));
    }
    magnifications = list.toArray(new String[list.size()]);
  }

  private static String magnification = magnifications[1];
  private Choice magnificationChoice;
  private Checkbox previewCheckBox;

  private Activation[] specificActivations;
  private Activation[] nonSpecificActivations;
  private int[] counts;

  private int nextPeakResultId;

  private enum Correction {
    //@formatter:off
    NONE{ @Override
    public String getName() { return "None"; }},
    SUBTRACTION{ @Override
    public String getName() { return "Subtraction"; }},
    MOST_LIKELY{ @Override
    public String getName() { return "Most likely"; }},
    WEIGHTED_RANDOM{ @Override
    public String getName() { return "Weighted random"; }};
    //@formatter:on

    @Override
    public String toString() {
      return getName();
    }

    /**
     * Gets the name.
     *
     * @return the name
     */
    public abstract String getName();
  }

  private enum SimulationDistribution {
    //@formatter:off
    POINT{ @Override
    public String getName() { return "Point"; }},
    LINE{ @Override
    public String getName() { return "Line"; }},
    CIRCLE{ @Override
    public String getName() { return "Circle"; }};
    //@formatter:on

    @Override
    public String toString() {
      return getName();
    }

    /**
     * Gets the name.
     *
     * @return the name
     */
    public abstract String getName();
  }

  private abstract class Shape {
    float x;
    float y;

    Shape(float x, float y) {
      this.x = x;
      this.y = y;
    }

    boolean canSample() {
      return true;
    }

    float[] getPosition() {
      return new float[] {x, y};
    }

    abstract float[] sample(UniformRandomProvider rng);
  }

  private class Point extends Shape {
    float[] xy;

    Point(float x, float y) {
      super(x, y);
      xy = super.getPosition();
    }

    @Override
    boolean canSample() {
      return false;
    }

    @Override
    float[] getPosition() {
      return xy;
    }

    @Override
    float[] sample(UniformRandomProvider rng) {
      throw new NotImplementedException();
    }
  }

  private class Line extends Shape {
    double radius;
    double length;
    float sina;
    float cosa;

    /**
     * Instantiates a new line.
     *
     * @param x the x
     * @param y the y
     * @param angle the angle (in radians)
     * @param radius the radius
     */
    Line(float x, float y, double angle, double radius) {
      super(x, y);
      this.radius = radius;
      length = 2 * radius;
      sina = (float) Math.sin(angle);
      cosa = (float) Math.cos(angle);
    }

    @Override
    float[] sample(UniformRandomProvider rng) {
      final float p = (float) (-radius + rng.nextDouble() * length);
      return new float[] {sina * p + x, cosa * p + y};
    }
  }

  private class Circle extends Shape {
    double radius;

    Circle(float x, float y, double radius) {
      super(x, y);
      this.radius = radius;
    }

    @Override
    float[] sample(UniformRandomProvider rng) {
      final double[] v = new UnitSphereSampler(2, rng).nextVector();
      return new float[] {(float) (v[0] * radius + x), (float) (v[1] * radius + y)};
    }
  }

  private static class Activation implements Molecule {
    final Trace trace;
    float x;
    float y;
    final int channel;
    int currentChannel;

    Activation(Trace trace, int channel) {
      this.trace = trace;
      final float[] centroid = trace.getCentroid(CentroidMethod.SIGNAL_WEIGHTED);
      x = centroid[0];
      y = centroid[1];
      this.channel = channel;
      currentChannel = channel;
    }

    boolean hasChannel() {
      return channel != 0;
    }

    int getChannel() {
      return channel - 1;
    }

    boolean hasCurrentChannel() {
      return currentChannel != 0;
    }

    int getCurrentChannel() {
      return currentChannel - 1;
    }

    @Override
    public float getX() {
      return x;
    }

    @Override
    public float getY() {
      return y;
    }

    @Override
    public int getId() {
      // Allow the ID to be updated from the original channel by using a current channel field
      // Note: the ID must be zero or above
      return currentChannel;
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    switch (isSimulation()) {
      case -1:
        // Cancelled
        return;
      case 1:
        // OK'd
        runSimulation();
        return;
      case 0:
      default:
        // Most common to not run the simulation
        break;
    }

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(title, "No localisations in memory");
      return;
    }

    final boolean crosstalkMode = "crosstalk".equals(arg);

    if (!showDialog(crosstalkMode)) {
      return;
    }

    // Load the results
    results = ResultsManager.loadInputResults(inputOption, false, DistanceUnit.PIXEL, null);
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(title, "No results could be loaded");
      return;
    }

    if (!results.isCalibrated()) {
      IJ.error(title, "Results must have basic calibration (pixel pitch and gain)");
      return;
    }

    // Get the traces
    traces = TraceManager.convert(results);
    if (traces == null || traces.length == 0) {
      IJ.error(title, "No traces could be loaded");
      return;
    }

    if (!showPulseCycleDialog()) {
      return;
    }

    createActivations();

    if (crosstalkMode) {
      runCrosstalkAnalysis();
    } else {
      runPulseAnalysis();
    }
  }

  private boolean showDialog(boolean crosstalkMode) {
    title = ((crosstalkMode) ? "Crosstalk " : "Pulse ") + title;

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(title);

    if (crosstalkMode) {
      gd.addMessage("Analyse crosstalk activation rate");
    } else {
      gd.addMessage("Count & plot molecules activated after a pulse");
    }

    ResultsManager.addInput(gd, "Input", inputOption, InputSource.MEMORY_CLUSTERED);
    final int min = (crosstalkMode) ? 2 : 1;
    gd.addSlider("Channels", min, MAX_CHANNELS, channels);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    inputOption = ResultsManager.getInputSource(gd);
    channels = (int) gd.getNextNumber();
    if (channels < min || channels > MAX_CHANNELS) {
      IJ.error(title, "Channels must be between " + min + " and " + MAX_CHANNELS);
      return false;
    }

    return true;
  }

  private boolean showPulseCycleDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(title);

    gd.addMessage("Specify the pulse cycle");

    gd.addNumericField("Repeat_interval", repeatInterval, 0);
    gd.addNumericField("Dark_frames_for_new_activation", darkFramesForNewActivation, 0);
    for (int c = 1; c <= channels; c++) {
      gd.addNumericField("Activation_frame_C" + c, startFrame[c - 1], 0);
    }

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    repeatInterval = (int) gd.getNextNumber();
    if (repeatInterval < channels) {
      IJ.error(title, "Repeat interval must be greater than the number of channels: " + channels);
      return false;
    }
    darkFramesForNewActivation = Math.max(1, (int) gd.getNextNumber());
    for (int c = 1; c <= channels; c++) {
      final int frame = (int) gd.getNextNumber();
      if (frame < 1 || frame > repeatInterval) {
        IJ.error(title, "Channel " + c + " activation frame must within the repeat interval");
        return false;
      }
      startFrame[c - 1] = frame;
    }

    // Check all start frames are unique
    for (int i = 0; i < channels; i++) {
      for (int j = i + 1; j < channels; j++) {
        if (startFrame[i] == startFrame[j]) {
          IJ.error(title, "Start frames must be unique for each channel");
          return false;
        }
      }
    }

    return true;
  }

  /**
   * Creates the activations. This splits the input traces into continuous chains of localisations.
   * Each chain is an activation. A new activation is created if there are more than the configured
   * number of dark frames since the last localisation. The start frame for the activation defines
   * the channel the activation is assigned to (this may be channel 0 if the start frame is not in a
   * pulse start frame).
   */
  @SuppressWarnings("null")
  private void createActivations() {
    final TurboList<Activation> activations = new TurboList<>(traces.length);

    // Activations are only counted if there are at least
    // n frames between localisations.
    final int n = darkFramesForNewActivation + 1;

    for (final Trace trace : traces) {
      trace.sort(); // Time-order

      final PeakResultStoreList points = trace.getPoints();

      // Define the frame for a new activation
      int nextActivationStartFrame = Integer.MIN_VALUE;
      Trace current = null;
      int channel = 0;
      for (int j = 0; j < points.size(); j++) {
        final PeakResult p = points.get(j);
        // Check if this is an activation
        if (p.getFrame() >= nextActivationStartFrame) {
          if (current != null) {
            // Store the last
            activations.add(new Activation(current, channel));
          }

          // Create a new activation
          current = new Trace(p);
          channel = getChannel(p);
        } else {
          // This is the same chain of localisations
          current.add(p);
        }
        nextActivationStartFrame = p.getEndFrame() + n;
      }

      if (current != null) {
        activations.add(new Activation(current, channel));
      }
    }

    save(activations);
  }

  private void save(TurboList<Activation> list) {
    // Count the activations per channel
    // Note: Channels are 0-indexed in the activations
    counts = new int[channels];
    for (int i = list.size(); i-- > 0;) {
      final Activation result = list.getf(i);
      if (result.hasChannel()) {
        counts[result.getChannel()]++;
      }
    }

    // Store specific activations
    final int sum = (int) MathUtils.sum(counts);
    specificActivations = new Activation[sum];
    final int nonSpecificActivationsSize = list.size() - sum;
    nonSpecificActivations = new Activation[nonSpecificActivationsSize];
    for (int i = list.size(), c1 = 0, c2 = 0; i-- > 0;) {
      final Activation result = list.getf(i);
      if (result.hasChannel()) {
        specificActivations[c1++] = result;
      } else {
        nonSpecificActivations[c2++] = result;
      }
    }

    // Output activation rates
    final int[] frameCount = new int[channels + 1];
    int firstFrame = results.getMinFrame();
    int lastFrame = results.getMaxFrame();

    // if (false)
    // {
    // for (int t = firstFrame; t <= lastFrame; t++)
    // {
    // frameCount[getChannel(t)]++;
    // }
    // }
    // else
    // {
    // Move the ends to the repeat interval
    while (firstFrame % repeatInterval != 1) {
      frameCount[getChannel(firstFrame++)]++;
    }
    while (lastFrame % repeatInterval != 0) {
      frameCount[getChannel(lastFrame--)]++;
    }
    final int total = lastFrame - firstFrame + 1;
    final int cycles = total / repeatInterval;
    for (int c = 1; c <= channels; c++) {
      frameCount[c] += cycles;
    }
    final int remaining = (total - channels * cycles);
    frameCount[0] += remaining;
    // }

    printRate("Background", nonSpecificActivationsSize, frameCount[0]);
    for (int c = 1; c <= channels; c++) {
      printRate("Channel " + c, counts[c - 1], frameCount[c]);
    }
  }

  private static void printRate(String title, int count, int numberOfFrames) {
    ImageJUtils.log("Activation rate : %s = %d/%d = %s per frame", title, count, numberOfFrames,
        MathUtils.rounded((double) count / numberOfFrames));
  }

  private static int getChannel(PeakResult result) {
    return getChannel(result.getFrame());
  }

  private static int getChannel(int frame) {
    // Classify if within a channel activation start frame
    final int mod = frame % repeatInterval;
    for (int i = 0; i < channels; i++) {
      if (mod == startFrame[i]) {
        return i + 1;
      }
    }
    return 0;
  }

  private void runCrosstalkAnalysis() {
    // Determine the cross talk ratio.
    // This is done by imaging only a single photo-switchable probe with the
    // same activation pulse imaging routine used for multi-colour imaging.
    // Concept:
    // A probe is meant to turn on in a frame following a pulse from a specific wavelength.
    // Multi-wavelengths can be used with probes responding to each wavelength. However
    // each probe may be activated by the 'wrong' wavelength. This is crosstalk.
    // The idea is to understand how many times the probe will turn on in a
    // frame following a pulse from the other lasers.

    // To determine the crosstalk ratio we must have a single probe imaged with the full
    // multi-wavelength pulse cycle. We then count how many times a probe activated by
    // the correct wavelength is activated by the others.
    // Crosstalk for each wavelength is then the fraction of times molecules were activated
    // by the 'wrong' wavelength.

    if (!showCrossTalkAnalysisDialog()) {
      return;
    }

    final double[] crosstalk = computeCrosstalk(counts, targetChannel - 1);

    // Store the cross talk.
    // Crosstalk from M into N is defined as the number of times the molecule that should be
    // activated by a pulse from channel M is activated by a pulse from channel N.
    // targetChannel = M
    // activationChannel = N
    int index1;
    int index2 = -1;
    if (channels == 2) {
      if (targetChannel == 1) {
        index1 = setCrosstalk(C12, crosstalk[1]);
      } else {
        index1 = setCrosstalk(C21, crosstalk[0]);
      }

      // 3-channel
    } else if (targetChannel == 1) {
      index1 = setCrosstalk(C12, crosstalk[1]);
      index2 = setCrosstalk(C13, crosstalk[2]);
    } else if (targetChannel == 2) {
      index1 = setCrosstalk(C21, crosstalk[0]);
      index2 = setCrosstalk(C23, crosstalk[2]);
    } else {
      index1 = setCrosstalk(C31, crosstalk[0]);
      index2 = setCrosstalk(C32, crosstalk[1]);
    }

    // Show fraction activations histogram. So we have to set the sum to 1
    final double sum = MathUtils.sum(crosstalk);
    for (int i = 0; i < crosstalk.length; i++) {
      crosstalk[i] /= sum;
    }

    // Plot a histogram
    final double[] x = SimpleArrayUtils.newArray(channels, 0.5, 1);
    final double[] y = crosstalk;
    final Plot2 plot = new Plot2(title, "Channel", "Fraction activations");
    plot.setLimits(0, channels + 1, 0, 1);
    plot.setXMinorTicks(false);
    plot.addPoints(x, y, Plot2.BAR);
    String label = String.format("Crosstalk %s = %s", ctNames[index1], MathUtils.round(ct[index1]));
    if (index2 > -1) {
      label += String.format(", %s = %s", ctNames[index2], MathUtils.round(ct[index2]));
    }
    plot.addLabel(0, 0, label);
    ImageJUtils.display(title, plot);
  }

  /**
   * Compute crosstalk.
   *
   * <p>"The crosstalk ratios can be calculated from the ratios of incorrectly to correctly colored
   * localizations."
   *
   * @param count the count
   * @param target the target
   * @return the double[]
   */
  private static double[] computeCrosstalk(int[] count, int target) {
    final double[] crosstalk = new double[count.length];
    for (int c = 0; c < count.length; c++) {
      crosstalk[c] = (double) count[c] / count[target];
    }
    return crosstalk;
  }

  private static int setCrosstalk(int index, double value) {
    ct[index] = value;
    return index;
  }

  private boolean showCrossTalkAnalysisDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(title);

    gd.addMessage(TextUtils
        .wrap("Crosstalk analysis requires a sample singly labelled with only one photo-switchable"
            + " probe and imaged with the full pulse lifecycle. The probe should be activated by"
            + " the pulse in the target channel. Activations from the pulse in other channels"
            + " is crosstalk.", 80));

    final String[] ch = new String[channels];
    for (int i = 0; i < ch.length; i++) {
      ch[i] = "Channel " + (i + 1);
    }

    gd.addChoice("Target", ch, "Channel " + targetChannel);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    targetChannel = gd.getNextChoiceIndex() + 1;

    return true;
  }

  /**
   * Unmix the observed local densities into the actual densities for 2-channels.
   *
   * <p>Crosstalk from M into N is defined as the number of times the molecule that should be
   * activated by a pulse from channel M is activated by a pulse from channel N. A value less than 1
   * is expected (otherwise the fluorophore is not being specifically activated by channel M).
   *
   * @param od1 the observed density in channel 1
   * @param od2 the observed density in channel 2
   * @param c21 the crosstalk from channel 2 into channel 1
   * @param c12 the crosstalk from channel 1 into channel 2
   * @return the actual densities [d1, d2]
   */
  public static double[] unmix(double od1, double od2, double c21, double c12) {
    // Solve the equations:
    // od1 = d1 + c21 * d2
    // od2 = d2 + c12 * d1
    // This is done by direct substitution
    final double d1 = (od1 - c21 * od2) / (1 - c12 * c21);
    final double d2 = od2 - c12 * d1;
    // Assuming od1 and od2 are positive and c12 and c21 are
    // between 0 and 1 then we do not need to check the bounds.
    // d1 = Maths.clip(0, od1, d1);
    // d2 = Maths.clip(0, od2, d2);
    return new double[] {d1, d2};
  }

  /**
   * Unmix the observed local densities into the actual densities for 3-channels.
   *
   * <p>Crosstalk from M into N is defined as the number of times the molecule that should be
   * activated by a pulse from channel M is activated by a pulse from channel N. A value less than 1
   * is expected (otherwise the fluorophore is not being specifically activated by channel M).
   *
   * @param od1 the observed density in channel 1
   * @param od2 the observed density in channel 2
   * @param od3 the observed density in channel 3
   * @param c21 the crosstalk from channel 2 into channel 1
   * @param c31 the crosstalk from channel 3 into channel 1
   * @param c12 the crosstalk from channel 1 into channel 2
   * @param c32 the crosstalk from channel 3 into channel 2
   * @param c13 the crosstalk from channel 1 into channel 3
   * @param c23 the crosstalk from channel 2 into channel 3
   * @return the actual densities [d1, d2, d3]
   */
  public static double[] unmix(double od1, double od2, double od3, double c21, double c31,
      double c12, double c32, double c13, double c23) {
    // Solve the linear equations: A * X = B
    // od1 = d1 + c21 * d2 + c31 * d3
    // od2 = d2 + c12 * d1 + c32 * d3
    // od3 = d3 + c13 * d1 + c23 * d2

    // Use matrix inversion so that: X = A^-1 * B
    // @CHECKSTYLE.OFF: LocalVariableName
    double a = 1;
    double b = c21;
    double c = c31;
    double d = c12;
    double e = 1;
    double f = c32;
    double g = c13;
    double h = c23;
    double i = 1;
    // @CHECKSTYLE.ON: LocalVariableName

    final double A = (e * i - f * h);
    final double B = -(d * i - f * g);
    final double C = (d * h - e * g);

    final double det = a * A + b * B + c * C;

    final double det_recip = 1.0 / det;

    if (!Double.isFinite(det_recip)) {
      // Failed so reset to the observed densities
      return new double[] {od1, od2, od3};
    }

    final double D = -(b * i - c * h);
    final double E = (a * i - c * g);
    final double F = -(a * h - b * g);
    final double G = (b * f - c * e);
    final double H = -(a * f - c * d);
    final double I = (a * e - b * d);

    a = det_recip * A;
    b = det_recip * D;
    c = det_recip * G;
    d = det_recip * B;
    e = det_recip * E;
    f = det_recip * H;
    g = det_recip * C;
    h = det_recip * F;
    i = det_recip * I;

    final double[] x = new double[3];
    x[0] = a * od1 + b * od2 + c * od3;
    x[1] = d * od1 + e * od2 + f * od3;
    x[2] = g * od1 + h * od2 + i * od3;

    // Use matrix decomposition
    // // Note: The linear solver uses LU decomposition.
    // // We cannot use a faster method as the matrix A is not symmetric.
    // final LinearSolver<DenseMatrix64F> linearSolver = LinearSolverFactory.linear(3);
    // final DenseMatrix64F A = new DenseMatrix64F(3, 3);
    // A.set(0, 1);
    // A.set(1, c21);
    // A.set(2, c31);
    // A.set(3, c12);
    // A.set(4, 1);
    // A.set(5, c32);
    // A.set(6, c13);
    // A.set(7, c23);
    // A.set(8, 1);
    // final DenseMatrix64F B = new DenseMatrix64F(3, 1);
    // B.set(0, od1);
    // B.set(1, od2);
    // B.set(2, od3);
    //
    // if (!linearSolver.setA(A))
    // {
    // // Failed so reset to the observed densities
    // return new double[] { od1, od2, od3 };
    // }
    //
    // // Input B is not modified to so we can re-use for output X
    // linearSolver.solve(B, B);
    //
    // final double[] x = B.getData();

    // Due to floating-point error in the decomposition we check the bounds
    x[0] = MathUtils.clip(0, od1, x[0]);
    x[1] = MathUtils.clip(0, od2, x[1]);
    x[2] = MathUtils.clip(0, od3, x[2]);

    return x;
  }

  private class RunSettings {
    double densityRadius;
    int minNeighbours;
    Correction specificCorrection = Correction.NONE;
    double[] specificCorrectionCutoff;
    Correction nonSpecificCorrection = Correction.NONE;
    double nonSpecificCorrectionCutoff;

    final ResultsSettings resultsSettings;

    RunSettings(ResultsSettings resultsSettings) {
      // Copy the current settings required
      if (channels > 1) {
        this.densityRadius = PulseActivationAnalysis.densityRadius / results.getNmPerPixel();
        this.minNeighbours = PulseActivationAnalysis.minNeighbours;

        specificCorrection = getCorrection(PulseActivationAnalysis.specificCorrection,
            PulseActivationAnalysis.specificCorrectionIndex);
        this.specificCorrectionCutoff = new double[channels];
        for (int i = channels; i-- > 0;) {
          // Convert from percentage to a probability
          this.specificCorrectionCutoff[i] =
              PulseActivationAnalysis.specificCorrectionCutoff[i] / 100.0;
        }

        nonSpecificCorrection = getCorrection(PulseActivationAnalysis.nonSpecificCorrection,
            PulseActivationAnalysis.nonSpecificCorrectionIndex);
        this.nonSpecificCorrectionCutoff =
            PulseActivationAnalysis.nonSpecificCorrectionCutoff / 100.0;
      }
      this.resultsSettings = resultsSettings;
    }

    Correction getCorrection(Correction[] correction, int index) {
      if (index >= 0 && index < correction.length) {
        return correction[index];
      }
      return Correction.NONE;
    }

    public boolean newUnmixSettings(RunSettings lastRunSettings) {
      if (lastRunSettings == null) {
        return true;
      }
      if (lastRunSettings.densityRadius != densityRadius) {
        return true;
      }
      if (lastRunSettings.minNeighbours != minNeighbours) {
        return true;
      }
      if (lastRunSettings.specificCorrection != specificCorrection) {
        return true;
      }
      if (specificCorrection != Correction.NONE) {
        for (int i = channels; i-- > 0;) {
          if (lastRunSettings.specificCorrectionCutoff[i] != specificCorrectionCutoff[i]) {
            return true;
          }
        }
      }
      return false;
    }

    public boolean newNonSpecificCorrectionSettings(RunSettings lastRunSettings) {
      if (lastRunSettings == null) {
        return true;
      }
      if (lastRunSettings.densityRadius != densityRadius) {
        return true;
      }
      if (lastRunSettings.minNeighbours != minNeighbours) {
        return true;
      }
      if (lastRunSettings.nonSpecificCorrection != nonSpecificCorrection) {
        return true;
      }
      if (nonSpecificCorrection != Correction.NONE) {
        if (lastRunSettings.nonSpecificCorrectionCutoff != nonSpecificCorrectionCutoff) {
          return true;
        }
      }
      return false;
    }
  }

  // Here we use a simple workflow with only one worker since the results are
  // written straight back to this class' objects
  private Workflow<RunSettings, Object> workflow;

  private void runPulseAnalysis() {
    // Use a simple workflow with one worker
    workflow = new Workflow<>();
    workflow.add(new WorkflowWorker<RunSettings, Object>() {
      @Override
      public boolean equalSettings(RunSettings current, RunSettings previous) {
        return false;
      }

      @Override
      public boolean equalResults(Object current, Object previous) {
        return false;
      }

      @Override
      public Pair<RunSettings, Object> doWork(Pair<RunSettings, Object> work) {
        synchronized (PulseActivationAnalysis.this.analysisLock) {
          PulseActivationAnalysis.this.runAnalysis(work.getKey());
        }
        return work;
      }
    });

    workflow.start();

    final boolean cancelled = !showPulseAnalysisDialog();

    workflow.shutdown(cancelled);

    if (executor != null) {
      if (cancelled) {
        // Stop immediately
        executor.shutdownNow();
      } else {
        executor.shutdown();
      }
    }
  }

  private boolean showPulseAnalysisDialog() {
    final NonBlockingExtendedGenericDialog gd = new NonBlockingExtendedGenericDialog(title);

    gd.addMessage("Plot molecules activated after a pulse");
    String[] correctionNames = null;
    String[] assignmentNames = null;

    if (channels > 1) {
      if (channels == 2) {
        gd.addNumericField("Crosstalk_21", ct[C21], 3);
        gd.addNumericField("Crosstalk_12", ct[C12], 3);
      } else {
        for (int i = 0; i < ctNames.length; i++) {
          gd.addNumericField("Crosstalk_" + ctNames[i], ct[i], 3);
        }
      }

      gd.addNumericField("Local_density_radius", densityRadius, 0, 6, "nm");
      gd.addSlider("Min_neighbours", 0, 15, minNeighbours);
      correctionNames = SettingsManager.getNames((Object[]) specificCorrection);
      gd.addChoice("Crosstalk_correction", correctionNames,
          correctionNames[specificCorrectionIndex]);
      for (int c = 1; c <= channels; c++) {
        gd.addSlider("Crosstalk_correction_cutoff_C" + c + "(%)", 0, 100,
            specificCorrectionCutoff[c - 1]);
      }
      assignmentNames = SettingsManager.getNames((Object[]) nonSpecificCorrection);
      gd.addChoice("Nonspecific_assigment", assignmentNames,
          assignmentNames[nonSpecificCorrectionIndex]);
      gd.addSlider("Nonspecific_assignment_cutoff (%)", 0, 100, nonSpecificCorrectionCutoff);
    }

    resultsSettingsBuilder = SettingsManager.readResultsSettings(0).toBuilder();
    ResultsManager.addImageResultsOptions(gd, resultsSettingsBuilder, 0);

    // ResultsImageSettings s = resultsSettings.getResultsImageSettings();
    // gd.addChoice("Image", SettingsManager.getResultsImageTypeNames(),
    // SettingsManager.getResultsImageTypeNames()[s.getImageTypeValue()]);
    // gd.addCheckbox("Weighted", s.getWeighted());
    // gd.addCheckbox("Equalised", s.getEqualised());
    // gd.addSlider("Image_Precision (nm)", 5, 30, s.getAveragePrecision());
    // gd.addSlider("Image_Scale", 1, 15, s.getScale());

    previewCheckBox = gd.addAndGetCheckbox("Preview", false);

    final String buttonLabel = "Draw loop";
    gd.addMessage("Click '" + buttonLabel + "' to draw the current ROIs in a loop view");
    gd.addAndGetButton(buttonLabel, this::actionPerformed);
    magnificationChoice = gd.addAndGetChoice("Magnification", magnifications, magnification);

    gd.addDialogListener(this::dialogItemChanged);
    gd.addOptionCollectedListener(event -> addWork(previewCheckBox.getState()));

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    // The dialog was OK'd so run if work was staged in the workflow.
    if (workflow.isStaged()) {
      workflow.runStaged();
    }

    // Record options for a macro since the NonBlockingDialog does not
    if (Recorder.record) {
      if (channels > 1) {
        // Suppress null warnings
        if (correctionNames == null || assignmentNames == null) {
          throw new IllegalStateException();
        }

        if (channels == 2) {
          Recorder.recordOption("Crosstalk_21", Double.toString(ct[C21]));
          Recorder.recordOption("Crosstalk_12", Double.toString(ct[C12]));
        } else {
          for (int i = 0; i < ctNames.length; i++) {
            Recorder.recordOption("Crosstalk_" + ctNames[i], Double.toString(ct[i]));
          }
        }

        Recorder.recordOption("Local_density_radius", Double.toString(densityRadius));
        Recorder.recordOption("Min_neighbours", Integer.toString(minNeighbours));

        Recorder.recordOption("Crosstalk_correction", correctionNames[specificCorrectionIndex]);
        for (int c = 1; c <= channels; c++) {
          Recorder.recordOption("Crosstalk_correction_cutoff_C" + c,
              Double.toString(specificCorrectionCutoff[c - 1]));
        }

        Recorder.recordOption("Nonspecific_assigment", assignmentNames[nonSpecificCorrectionIndex]);
        Recorder.recordOption("Nonspecific_assignment_cutoff (%)",
            Double.toString(nonSpecificCorrectionCutoff));
      }

      final ResultsImageSettings s = resultsSettingsBuilder.getResultsImageSettings();
      Recorder.recordOption("Image",
          SettingsManager.getResultsImageTypeNames()[s.getImageTypeValue()]);
      if (s.getWeighted()) {
        Recorder.recordOption("Weighted");
      }
      if (s.getEqualised()) {
        Recorder.recordOption("Equalised");
      }
      Recorder.recordOption("Image_Precision", Double.toString(s.getAveragePrecision()));
      Recorder.recordOption("Image_Scale", Double.toString(s.getScale()));
    }

    SettingsManager.writeSettings(resultsSettingsBuilder);

    return true;
  }

  private boolean dialogItemChanged(GenericDialog gd, AWTEvent event) {
    // The event is null when the NonBlockingExtendedGenericDialog is first shown.
    // Do not ignore this if a macro.
    if (event == null && !ImageJUtils.isMacro()) {
      return true;
    }

    // Check arguments
    try {
      if (channels > 1) {
        if (channels == 2) {
          ct[C21] = gd.getNextNumber();
          ct[C12] = gd.getNextNumber();
          validateCrosstalk(C21);
          validateCrosstalk(C12);
        } else {
          ct[C21] = gd.getNextNumber();
          ct[C31] = gd.getNextNumber();
          ct[C12] = gd.getNextNumber();
          ct[C32] = gd.getNextNumber();
          ct[C13] = gd.getNextNumber();
          ct[C23] = gd.getNextNumber();
          for (int i = 0; i < ct.length; i += 2) {
            validateCrosstalk(i, i + 1);
          }
        }

        densityRadius = Math.abs(gd.getNextNumber());
        minNeighbours = Math.abs((int) gd.getNextNumber());
        specificCorrectionIndex = gd.getNextChoiceIndex();
        for (int c = 1; c <= channels; c++) {
          specificCorrectionCutoff[c - 1] = (int) gd.getNextNumber();
          validatePercentage("Crosstalk_correction_cutoff_C" + c, specificCorrectionCutoff[c - 1]);
        }
        nonSpecificCorrectionIndex = gd.getNextChoiceIndex();
        nonSpecificCorrectionCutoff = gd.getNextNumber();
        validatePercentage("Nonspecific_assignment_cutoff", nonSpecificCorrectionCutoff);
      }
    } catch (final IllegalArgumentException ex) {
      IJ.error(title, ex.getMessage());
      return false;
    }

    resultsSettingsBuilder.getResultsImageSettingsBuilder()
        .setImageTypeValue(gd.getNextChoiceIndex());

    // ResultsImageSettings.Builder s = resultsSettingsBuilder.getResultsImageSettingsBuilder();
    // s.setImageTypeValue(gd.getNextChoiceIndex());
    // s.setWeighted(gd.getNextBoolean());
    // s.setEqualised(gd.getNextBoolean());
    // s.setAveragePrecision(gd.getNextNumber());
    // s.setScale(gd.getNextNumber());

    final boolean preview = gd.getNextBoolean();

    if (gd.invalidNumber()) {
      return false;
    }

    // This should not fail as this plugin built the dialog
    ((NonBlockingExtendedGenericDialog) gd).collectOptions();

    addWork(preview);

    return true;
  }

  private void addWork(boolean preview) {
    final RunSettings settings = new RunSettings(resultsSettingsBuilder.build());
    if (preview) {
      // Run the settings
      workflow.run(settings);
      workflow.startPreview();
    } else {
      workflow.stopPreview();
      // Stage the work but do not run
      workflow.stage(settings);
    }
  }

  private static void validateCrosstalk(int index) {
    final String name = "Crosstalk " + ctNames[index];
    ParameterUtils.isPositive(name, ct[index]);
    ParameterUtils.isBelow(name, ct[index], 0.5);
  }

  private static void validateCrosstalk(int index1, int index2) {
    validateCrosstalk(index1);
    validateCrosstalk(index2);
    ParameterUtils.isBelow("Crosstalk " + ctNames[index1] + " + " + ctNames[index2],
        ct[index1] + ct[index2], 0.5);
  }

  private static void validatePercentage(String name, double percentage) {
    ParameterUtils.isPositive(name, percentage);
    ParameterUtils.isEqualOrBelow(name, percentage, 100);
  }

  /**
   * The analysis lock. All access to state modified by the analysis must be synchronized on this.
   */
  private final Object analysisLock = new Object();
  private DensityCounter dc;
  private int[][] density;
  private int numberOfThreads;
  private ExecutorService executor;
  private TurboList<Future<?>> futures;
  /** The last run settings. All access to the must be synchronized. */
  private RunSettings lastRunSettings;

  /**
   * Run the analysis. This modifies state and so should be synchronized.
   *
   * <p>Note: We check against the last settings and only repeat what is necessary.
   *
   * @param runSettings the run settings
   */
  private void runAnalysis(RunSettings runSettings) {
    if (runSettings == null) {
      lastRunSettings = null;
      return;
    }

    IJ.showStatus("Analysing ...");

    // Assign all activations to a channel.
    // This is only necessary when we have more than 1 channel. If we have 1 channel then
    // no correction method is specified.
    boolean changed = false;
    if (runSettings.newUnmixSettings(lastRunSettings)) {
      changed = true;

      // Reset
      for (int i = specificActivations.length; i-- > 0;) {
        final Activation result = specificActivations[i];
        result.currentChannel = result.channel;
      }

      if (runSettings.specificCorrection != Correction.NONE) {
        // Use a density counter that can put all the activations on a grid.
        // It has a method to count the number of activations within a radius that
        // belong to each channel.

        // Add only those with specific activations. Non-specific activations are ignored.
        createDensityCounter((float) runSettings.densityRadius);

        // Do this all together: it uses a faster algorithm and we can cache the results
        if (density == null) {
          IJ.showStatus("Computing observed density");
          density = dc.countAll(channels);
        }

        final SplitMix sm = new SplitMix(System.currentTimeMillis());

        // -=-=-=--=-=-
        // Unmix the specific activations to their correct channel.
        // -=-=-=--=-=-
        IJ.showStatus("Unmixing");
        createThreadPool();

        final int[] newChannel = new int[specificActivations.length];

        final int nPerThread =
            (int) Math.ceil((double) specificActivations.length / numberOfThreads);
        for (int from = 0; from < specificActivations.length;) {
          final int to = Math.min(from + nPerThread, specificActivations.length);
          futures.add(executor.submit(new SpecificUnmixWorker(runSettings, density, newChannel,
              from, to, sm.copyAndJump())));
          from = to;
        }
        waitToFinish();

        // Update the channel assignment
        for (int i = specificActivations.length; i-- > 0;) {
          specificActivations[i].currentChannel = newChannel[i];
        }
      }
    }

    // -=-=-=--=-=-
    // Assign non-specific activations
    // -=-=-=--=-=-
    if (changed || runSettings.newNonSpecificCorrectionSettings(lastRunSettings)) {
      // Reset
      for (int i = nonSpecificActivations.length; i-- > 0;) {
        final Activation result = nonSpecificActivations[i];
        result.currentChannel = result.channel;
      }

      if (runSettings.nonSpecificCorrection != Correction.NONE) {
        createDensityCounter((float) runSettings.densityRadius);

        final SplitMix sm = new SplitMix(System.currentTimeMillis());

        IJ.showStatus("Non-specific assignment");
        createThreadPool();

        final int[] newChannel = new int[nonSpecificActivations.length];

        final int nPerThread =
            (int) Math.ceil((double) nonSpecificActivations.length / numberOfThreads);
        for (int from = 0; from < nonSpecificActivations.length;) {
          final int to = Math.min(from + nPerThread, nonSpecificActivations.length);
          futures.add(executor.submit(
              new NonSpecificUnmixWorker(runSettings, dc, newChannel, from, to, sm.copyAndJump())));
          from = to;
        }
        waitToFinish();

        // Update the channel assignment
        for (int i = nonSpecificActivations.length; i-- > 0;) {
          nonSpecificActivations[i].currentChannel = newChannel[i];
        }
      }
    }

    // Set-up outputs for each channel
    IJ.showStatus("Creating outputs");
    output = new PeakResultsList[channels];
    for (int c = 0; c < channels; c++) {
      output[c] = createOutput(c + 1);
    }

    // Create a results set with only those molecules assigned to a channel
    int count = write(output, specificActivations, 0);
    count = write(output, nonSpecificActivations, count);

    int size = 0;
    for (int c = 0; c < channels; c++) {
      output[c].end();
      size += output[c].size();
    }

    // Collate image into a stack
    if (channels > 1 && runSettings.resultsSettings.getResultsImageSettings()
        .getImageType() != ResultsImageType.DRAW_NONE) {
      final ImageProcessor[] images = new ImageProcessor[channels];
      for (int c = 0; c < channels; c++) {
        images[c] = getImage(output[c]);
      }
      displayComposite(images, results.getName() + " " + title);
    }

    lastRunSettings = runSettings;

    IJ.showStatus(String.format("%d/%s, %d/%s", count, TextUtils.pleural(traces.length, "Trace"),
        size, TextUtils.pleural(results.size(), "Result")));
  }

  private static void displayComposite(ImageProcessor[] images, String name) {
    ImageStack stack = null; // We do not yet know the size
    for (int i = 0; i < images.length; i++) {
      final ImageProcessor ip = images[i];
      if (stack == null) {
        stack = new ImageStack(ip.getWidth(), ip.getHeight());
      }
      ip.setColorModel(null);
      stack.addSlice("C" + (i + 1), ip);
    }

    // Create a composite
    ImagePlus imp = new ImagePlus(name, stack);
    imp.setDimensions(images.length, 1, 1);
    final CompositeImage ci = new CompositeImage(imp, IJ.COMPOSITE);

    // Make it easier to see
    // ij.plugin.ContrastEnhancerce = new ij.plugin.ContrastEnhancer();
    // double saturated = 0.35;
    // ce.stretchHistogram(ci, saturated);

    autoAdjust(ci, ci.getProcessor());

    imp = WindowManager.getImage(name);
    if (imp != null && imp.isComposite()) {
      ci.setMode(imp.getCompositeMode());
      imp.setImage(ci);
      imp.getWindow().toFront();
    } else {
      ci.show();
      imp = ci;
    }

    if (WindowManager.getWindow("Channels") == null) {
      IJ.run("Channels Tool...");
      final Window w = WindowManager.getWindow("Channels");
      if (w == null) {
        return;
      }
      final Window w2 = imp.getWindow();
      if (w2 == null) {
        return;
      }
      final java.awt.Point p = w2.getLocation();
      p.x += w2.getWidth();
      w.setLocation(p);
    }
  }

  /**
   * Auto adjust. Copied from {@link ij.plugin.frame.ContrastAdjuster }.
   *
   * <p>Although the ContrastAdjuster records its actions as 'run("Enhance Contrast",
   * "saturated=0.35");' it actually does something else which makes the image easier to see than
   * the afore mentioned command.
   *
   * @param imp the image
   * @param ip the image
   */
  private static void autoAdjust(ImagePlus imp, ImageProcessor ip) {
    final ij.measure.Calibration cal = imp.getCalibration();
    imp.setCalibration(null);
    final ImageStatistics stats = imp.getStatistics(); // get uncalibrated stats
    imp.setCalibration(cal);
    final int limit = stats.pixelCount / 10;
    final int[] histogram = stats.histogram;
    int autoThreshold = 0;
    if (autoThreshold < 10) {
      autoThreshold = 5000;
    } else {
      autoThreshold /= 2;
    }
    final int threshold = stats.pixelCount / autoThreshold;
    int index = -1;
    boolean found = false;
    int count;
    do {
      index++;
      count = histogram[index];
      if (count > limit) {
        count = 0;
      }
      found = count > threshold;
    } while (!found && index < 255);
    final int hmin = index;
    index = 256;
    do {
      index--;
      count = histogram[index];
      if (count > limit) {
        count = 0;
      }
      found = count > threshold;
    } while (!found && index > 0);
    final int hmax = index;
    if (hmax >= hmin) {
      double min = stats.histMin + hmin * stats.binSize;
      double max = stats.histMin + hmax * stats.binSize;
      if (Double.compare(min, max) == 0) {
        min = stats.min;
        max = stats.max;
      }
      imp.setDisplayRange(min, max);
    } else {
      reset(imp);
    }
  }

  private static void reset(ImagePlus imp) {
    final int bitDepth = imp.getBitDepth();
    double defaultMin;
    double defaultMax;
    if (bitDepth == 16 || bitDepth == 32) {
      imp.resetDisplayRange();
      defaultMin = imp.getDisplayRangeMin();
      defaultMax = imp.getDisplayRangeMax();
    } else {
      defaultMin = 0;
      defaultMax = 255;
    }
    imp.setDisplayRange(defaultMin, defaultMax);
  }

  private void createThreadPool() {
    if (executor == null) {
      numberOfThreads = Prefs.getThreads();
      executor = Executors.newFixedThreadPool(numberOfThreads);
      futures = new TurboList<>(numberOfThreads);
    }
  }

  private void createDensityCounter(float densityRadius) {
    if (dc == null || dc.getRadius() != densityRadius) {
      dc = new DensityCounter(specificActivations, densityRadius, false);
      // Clear cache of density
      density = null;
    }
  }

  private void waitToFinish() {
    ConcurrencyUtils.waitForCompletionUnchecked(futures, ex -> Logger
        .getLogger(getClass().getName()).log(Level.SEVERE, "Failed to finish analysis", ex));
    futures.clear();
  }

  private abstract class UnmixWorker {
    final int[] newChannel;
    final int from;
    final int to;
    UniformRandomProvider rng;
    int[] assignedChannel = new int[channels];
    double[] probability = new double[channels];

    public UnmixWorker(int[] newChannel, int from, int to, UniformRandomProvider rng) {
      this.newChannel = newChannel;
      this.from = from;
      this.to = to;
      this.rng = rng;
    }

    int weightedRandomSelection(double cutoff) {
      double sum = 0;
      for (int j = channels; j-- > 0;) {
        if (probability[j] > cutoff) {
          sum += probability[j];
        } else {
          probability[j] = 0;
        }
      }
      if (sum == 0) {
        return 0;
      }

      final double sum2 = sum * rng.nextDouble();
      sum = 0;
      for (int j = channels; j-- > 0;) {
        sum += probability[j];
        if (sum >= sum2) {
          return j + 1;
        }
      }
      // This should not happen
      return 0;
    }

    int mostLikelySelection(double cutoff) {
      double max = cutoff;
      int size = 0;
      for (int j = channels; j-- > 0;) {
        if (probability[j] > max) {
          size = 1;
          max = probability[j];
          assignedChannel[0] = j;
        } else if (probability[j] == max) {
          // Equal so store all for a random pick
          assignedChannel[size++] = j;
        }
      }

      if (size == 0) {
        return 0;
      }

      return (size > 1) ? assignedChannel[rng.nextInt(size)] + 1 : assignedChannel[0] + 1;
    }
  }

  /**
   * For processing the unmixing of specific channel activations.
   */
  private class SpecificUnmixWorker extends UnmixWorker implements Runnable {
    final RunSettings runSettings;
    final int[][] density;

    public SpecificUnmixWorker(RunSettings runSettings, int[][] density, int[] newChannel, int from,
        int to, UniformRandomProvider rng) {
      super(newChannel, from, to, rng);
      this.runSettings = runSettings;
      this.density = density;
    }

    @Override
    public void run() {
      for (int i = from; i < to; i++) {
        // Observed density
        final int[] obsDen = density[i];

        // Compute the number of neighbours.
        int neighbours = 0;
        for (int j = 1; j <= channels; j++) {
          neighbours += obsDen[j];
        }

        // Do not unmix if there are not enough neighbours.
        // Note this will count the target activation so
        // use <= to ensure the neighbours is above the min.
        if (neighbours <= runSettings.minNeighbours) {
          newChannel[i] = 0;
          continue;
        }

        // Current channel (1-indexed)
        int ch = specificActivations[i].channel;

        // Compute the true local densities
        double[] den;
        if (channels == 2) {
          den = unmix(obsDen[1], obsDen[2], ct[C12], ct[C12]);
        } else {
          den = unmix(obsDen[1], obsDen[2], obsDen[3], ct[C21], ct[C31], ct[C12], ct[C32], ct[C13],
              ct[C23]);
        }

        // Apply crosstalk correction
        if (runSettings.specificCorrection == Correction.SUBTRACTION) {
          // Compute the probability it is correct:
          // This is a measure of how much crosstalk effected the observed density.
          // (This is taken from Bates et al, 2007)
          final double pc = den[ch - 1] / obsDen[ch];

          // Remove it if below the subtraction threshold
          if (pc < runSettings.specificCorrectionCutoff[ch - 1]) {
            ch = 0;
          }
        } else {
          // Compute the probability of each channel as:
          // p(i) = di / (d1 + d2 + ... + dn)
          // Note this is different from computing the probability of the channel being correct.
          // That probability is an indication of how much crosstalk has effected the observed
          // density.
          // This value is a simple probability using the local density in each channel.
          double sum = 0;
          for (int j = channels; j-- > 0;) {
            sum += den[j];
          }
          // Note that since this is a specific activation we can assume the molecule will be
          // self-counted within the radius and that d will never be zero in every channel
          for (int j = channels; j-- > 0;) {
            probability[j] = den[j] / sum;
          }

          if (runSettings.specificCorrection == Correction.WEIGHTED_RANDOM) {
            ch = weightedRandomSelection(runSettings.specificCorrectionCutoff[ch - 1]);
          } else {
            ch = mostLikelySelection(runSettings.specificCorrectionCutoff[ch - 1]);
          }
        }

        newChannel[i] = ch;
      }
    }
  }

  /**
   * For processing the unmixing of specific channel activations.
   */
  private class NonSpecificUnmixWorker extends UnmixWorker implements Runnable {
    final RunSettings runSettings;
    final DensityCounter dc;

    public NonSpecificUnmixWorker(RunSettings runSettings, DensityCounter dc, int[] newChannel,
        int from, int to, UniformRandomProvider rng) {
      super(newChannel, from, to, rng);
      this.runSettings = runSettings;
      this.dc = dc;
    }

    @Override
    public void run() {
      // TODO - We could do other non-specific assignments.
      // e.g. Compute probability for each channel and assign
      // using a weighted random selection

      for (int i = from; i < to; i++) {
        int ch = 0;

        // Assume the observed density is the true local density
        // (i.e. cross talk correction of specific activations is perfect)
        final int[] d = dc.count(nonSpecificActivations[i], channels);

        // Compute the probability of each channel as:
        // p(i) = di / (d1 + d2 + ... + dn)
        double sum = 0;
        for (int j = 1; j <= channels; j++) {
          sum += d[j];
        }

        // Do not unmix if there are not enough neighbours.
        if (sum >= runSettings.minNeighbours) {
          for (int j = channels; j-- > 0;) {
            probability[j] = d[j + 1] / sum;
          }

          if (runSettings.nonSpecificCorrection == Correction.WEIGHTED_RANDOM) {
            ch = weightedRandomSelection(runSettings.nonSpecificCorrectionCutoff);
          } else {
            ch = mostLikelySelection(runSettings.nonSpecificCorrectionCutoff);
          }
        }

        newChannel[i] = ch;
      }
    }
  }

  private PeakResultsList createOutput(int channel) {
    final PeakResultsList outputList = new PeakResultsList();
    outputList.copySettings(results);
    if (channels > 1) {
      outputList.setName(results.getName() + " " + title + " C" + channel);
    } else {
      outputList.setName(results.getName() + " " + title);
    }

    // Store the set in memory
    final MemoryPeakResults memoryResults = new MemoryPeakResults(this.results.size());
    outputList.addOutput(memoryResults);
    MemoryPeakResults.addResults(memoryResults);

    // Draw the super-resolution image
    final Rectangle bounds = results.getBounds(true);
    addImageResults(outputList, results.getName(), bounds, results.getNmPerPixel(),
        results.getGain(), resultsSettingsBuilder.getResultsImageSettings());

    outputList.begin();

    return outputList;
  }

  private static void addImageResults(PeakResultsList resultsList, String title, Rectangle bounds,
      double nmPerPixel, double gain, ResultsImageSettings imageSettings) {
    if (imageSettings.getImageType() != ResultsImageType.DRAW_NONE) {
      final ImageJImagePeakResults image = ImagePeakResultsFactory.createPeakResultsImage(
          imageSettings.getImageType(), imageSettings.getWeighted(), imageSettings.getEqualised(),
          title, bounds, nmPerPixel, gain, imageSettings.getScale(),
          imageSettings.getAveragePrecision(), ResultsImageMode.IMAGE_ADD);
      image.setLiveImage(false);
      image.setDisplayImage(channels == 1);
      resultsList.addOutput(image);
    }
  }

  private static int write(PeakResultsList[] output, Activation[] activations, int count) {
    for (int i = activations.length; i-- > 0;) {
      final Activation result = activations[i];
      if (result.hasCurrentChannel()) {
        count++;
        output[result.getCurrentChannel()].addAll(result.trace.getPoints());
      }
    }
    return count;
  }

  private static ImageProcessor getImage(PeakResultsList peakResultsList) {
    final PeakResults[] list = peakResultsList.toArray();
    final ImageJImagePeakResults image = (ImageJImagePeakResults) list[1];
    return image.getImagePlus().getProcessor();
  }

  private int isSimulation() {
    if (ImageJUtils.isExtraOptions()) {
      final GenericDialog gd = new GenericDialog(title);
      gd.addMessage("Perform a crosstalk simulation?");
      gd.enableYesNoCancel();
      gd.showDialog();
      if (gd.wasOKed()) {
        return 1;
      }
      if (gd.wasCanceled()) {
        return -1;
      }
    }
    return 0;
  }

  private void runSimulation() {
    title += " Simulation";

    if (!showSimulationDialog()) {
      return;
    }

    final long start = System.currentTimeMillis();
    final UniformRandomProvider rng = getUniformRandomProvider();

    // Draw the molecule positions
    ImageJUtils.showStatus("Simulating molecules ...");
    final float[][][] molecules = new float[3][][];
    final MemoryPeakResults[] results = new MemoryPeakResults[3];
    final Calibration calibration = CalibrationHelper.create(sim_nmPerPixel, 1, 100);
    final Rectangle bounds = new Rectangle(0, 0, sim_size, sim_size);
    for (int c = 0; c < 3; c++) {
      molecules[c] = simulateMolecules(rng, c);

      // Create a dataset to store the activations
      final MemoryPeakResults r = new MemoryPeakResults();
      r.setCalibration(calibration);
      r.setBounds(bounds);
      r.setName(title + " C" + (c + 1));
      results[c] = r;
    }

    // Simulate activation
    ImageJUtils.showStatus("Simulating activations ...");
    for (int c = 0; c < 3; c++) {
      simulateActivations(rng, molecules, c, results);
    }

    // Combine
    ImageJUtils.showStatus("Producing simulation output ...");
    final MemoryPeakResults r = new MemoryPeakResults();
    r.setCalibration(calibration);
    r.setBounds((Rectangle) bounds.clone());
    r.setName(title);

    final ImageProcessor[] images = new ImageProcessor[3];
    for (int c = 0; c < 3; c++) {
      final PeakResult[] list = results[c].toArray();
      r.addAll(list);

      // Draw the unmixed activations
      final ImageJImagePeakResults image = ImagePeakResultsFactory.createPeakResultsImage(
          ResultsImageType.DRAW_LOCALISATIONS, true, true, title, bounds, sim_nmPerPixel, 1,
          1024.0 / sim_size, 0, ResultsImageMode.IMAGE_ADD);
      image.setCalibration(calibration);
      image.setLiveImage(false);
      image.setDisplayImage(false);
      image.begin();
      image.addAll(list);
      image.end();
      images[c] = image.getImagePlus().getProcessor();
    }
    displayComposite(images, title);

    // Add to memory. Set the composite dataset first.
    MemoryPeakResults.addResults(r);
    for (int c = 0; c < 3; c++) {
      MemoryPeakResults.addResults(results[c]);
    }

    // TODO:
    // Show an image of what it looks like with no unmixing, i.e. colours allocated
    // from the frame

    ImageJUtils.showStatus(
        "Simulation complete: " + TextUtils.millisToString(System.currentTimeMillis() - start));
  }

  private float[][] simulateMolecules(UniformRandomProvider rng, int channel) {
    final int n = sim_numberOfMolecules[channel];
    final float[][] molecules = new float[n][];
    if (n == 0) {
      return molecules;
    }

    // Draw the shapes
    final Shape[] shapes = createShapes(rng, channel);

    // Sample positions from within the shapes
    final boolean canSample = shapes[0].canSample();
    int count = 0;
    while (count < n) {
      float[] coords;
      if (canSample) {
        final int next = rng.nextInt(shapes.length);
        coords = shapes[next].sample(rng);
      } else {
        coords = shapes[count % shapes.length].getPosition();
      }

      // Avoid out-of-bounds positions
      if (outOfBounds(coords[0]) || outOfBounds(coords[1])) {
        continue;
      }
      molecules[count++] = coords;
    }
    return molecules;
  }

  private Shape[] createShapes(UniformRandomProvider rng, int channel) {
    Shape[] shapes;
    final double min = sim_size / 20.0;
    final double max = sim_size / 10.0;
    final double range = max - min;
    switch (sim_distribution[channel]) {
      case CIRCLE:
        shapes = new Shape[10];
        for (int i = 0; i < shapes.length; i++) {
          final float x = nextCoordinate(rng);
          final float y = nextCoordinate(rng);
          final double radius = rng.nextDouble() * range + min;
          shapes[i] = new Circle(x, y, radius);
        }
        break;

      case LINE:
        shapes = new Shape[10];
        for (int i = 0; i < shapes.length; i++) {
          final float x = nextCoordinate(rng);
          final float y = nextCoordinate(rng);
          final double angle = rng.nextDouble() * Math.PI;
          final double radius = rng.nextDouble() * range + min;
          shapes[i] = new Line(x, y, angle, radius);
        }

        break;

      case POINT:
      default:
        shapes = new Shape[sim_numberOfMolecules[channel]];
        for (int i = 0; i < shapes.length; i++) {
          final float x = nextCoordinate(rng);
          final float y = nextCoordinate(rng);
          shapes[i] = new Point(x, y);
        }
    }
    return shapes;
  }

  private static float nextCoordinate(UniformRandomProvider rng) {
    return (float) rng.nextDouble() * sim_size;
  }

  private static boolean outOfBounds(float value) {
    return value < 0 || value > sim_size;
  }

  private void simulateActivations(UniformRandomProvider rng, float[][][] molecules, int channel,
      MemoryPeakResults[] results) {
    final int n = molecules[channel].length;
    if (n == 0) {
      return;
    }

    // Compute desired number per frame
    final double umPerPixel = sim_nmPerPixel / 1000;
    final double um2PerPixel = umPerPixel * umPerPixel;
    final double area = sim_size * sim_size * um2PerPixel;
    final double nPerFrame = area * sim_activationDensity;

    // Compute the activation probability (but set an upper limit so not all are on in every frame)
    final double p = Math.min(0.5, nPerFrame / n);

    // Determine the other channels activation probability using crosstalk
    final double[] p0 = {p, p, p};
    int index1;
    int index2;
    int c1;
    int c2;
    switch (channel) {
      case 0:
        index1 = C12;
        index2 = C13;
        c1 = 1;
        c2 = 2;
        break;
      case 1:
        index1 = C21;
        index2 = C23;
        c1 = 0;
        c2 = 2;
        break;
      case 2:
      default:
        index1 = C31;
        index2 = C32;
        c1 = 0;
        c2 = 1;
        break;
    }
    p0[c1] *= ct[index1];
    p0[c2] *= ct[index2];

    // Assume 10 frames after each channel pulse => 30 frames per cycle
    final double precision = sim_precision[channel] / sim_nmPerPixel;

    final DiscreteSampler[] bd = new DiscreteSampler[4];
    for (int i = 0; i < 3; i++) {
      bd[i] = createBinomialDistribution(rng, n, p0[i]);
    }

    final int[] frames = new int[27];
    for (int i = 1, j = 0; i <= 30; i++) {
      if (i % 10 == 1) {
        // Skip specific activation frames
        continue;
      }
      frames[j++] = i;
    }
    bd[3] = createBinomialDistribution(rng, n, p * sim_nonSpecificFrequency);

    // Count the actual cross talk
    final int[] count = new int[3];

    for (int i = 0, t = 1; i < sim_cycles; i++, t += 30) {
      count[0] +=
          simulateActivations(rng, bd[0], molecules[channel], results[channel], t, precision);
      count[1] +=
          simulateActivations(rng, bd[1], molecules[channel], results[channel], t + 10, precision);
      count[2] +=
          simulateActivations(rng, bd[2], molecules[channel], results[channel], t + 20, precision);
      // Add non-specific activations
      if (bd[3] != null) {
        for (final int t2 : frames) {
          simulateActivations(rng, bd[3], molecules[channel], results[channel], t2, precision);
        }
      }
    }

    // Report simulated cross talk
    final double[] crosstalk = computeCrosstalk(count, channel);
    ImageJUtils.log("Simulated crosstalk C%s  %s=>%s, C%s  %s=>%s", ctNames[index1],
        MathUtils.rounded(ct[index1]), MathUtils.rounded(crosstalk[c1]), ctNames[index2],
        MathUtils.rounded(ct[index2]), MathUtils.rounded(crosstalk[c2]));
  }

  private int simulateActivations(UniformRandomProvider rng, DiscreteSampler bd,
      float[][] molecules, MemoryPeakResults results, int time, double precision) {
    if (bd == null) {
      return 0;
    }
    final int size = molecules.length;
    final int numberOfSamples = bd.sample();
    // Sample
    final NormalizedGaussianSampler gauss = new ZigguratNormalizedGaussianSampler(rng);
    final int[] samples = RandomUtils.sample(numberOfSamples, size, rng);
    for (final int index : samples) {
      final float[] xy = molecules[index];
      float x;
      float y;
      do {
        x = (float) (xy[0] + gauss.sample() * precision);
      } while (outOfBounds(x));
      do {
        y = (float) (xy[1] + gauss.sample() * precision);
      } while (outOfBounds(y));

      results.add(createResult(time, x, y));
    }
    return samples.length;
  }

  private static DiscreteSampler createBinomialDistribution(UniformRandomProvider rand, int trials,
      double pvalue) {
    if (pvalue == 0) {
      return null;
    }
    return SamplerUtils.createBinomialSampler(rand, trials, pvalue);
  }

  private IdPeakResult createResult(int time, float x, float y) {
    // We add them as if tracing is perfect. So each peak result has a new ID.
    // This allows the output of the simulation to be used directly by the pulse analysis code.
    final IdPeakResult r = new IdPeakResult(time, x, y, 1, ++nextPeakResultId);
    r.setNoise(1); // So it appears calibrated
    return r;
  }

  private boolean showSimulationDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(title);

    final SimulationDistribution[] distributionValues = SimulationDistribution.values();
    final String[] distribution = SettingsManager.getNames((Object[]) distributionValues);

    // Random crosstalk if not set
    if (MathUtils.max(ct) == 0) {
      // Have some crosstalk
      final ContinuousUniformSampler sampler =
          new ContinuousUniformSampler(getUniformRandomProvider(), 0.05, 0.15);
      for (int i = 0; i < ct.length; i++) {
        ct[i] = sampler.sample();
      }
    }

    // Three channel
    for (int c = 0; c < 3; c++) {
      final String ch = "_C" + (c + 1);
      gd.addNumericField("Molcules" + ch, sim_numberOfMolecules[c], 0);
      gd.addChoice("Distribution" + ch, distribution, distribution[sim_distribution[c].ordinal()]);
      gd.addNumericField("Precision_" + ch, sim_precision[c], 3);
      gd.addNumericField("Crosstalk_" + ctNames[2 * c], ct[2 * c], 3);
      gd.addNumericField("Crosstalk_" + ctNames[2 * c + 1], ct[2 * c + 1], 3);
    }
    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    int count = 0;
    for (int c = 0; c < 3; c++) {
      sim_numberOfMolecules[c] = (int) Math.abs(gd.getNextNumber());
      if (sim_numberOfMolecules[c] > 0) {
        count++;
      }
      sim_distribution[c] = distributionValues[gd.getNextChoiceIndex()];
      sim_precision[c] = Math.abs(gd.getNextNumber());
      ct[2 * c] = Math.abs(gd.getNextNumber());
      ct[2 * c + 1] = Math.abs(gd.getNextNumber());
    }

    if (gd.invalidNumber()) {
      return false;
    }
    if (count < 2) {
      IJ.error(title, "Simulation requires at least 2 channels");
      return false;
    }

    try {
      for (int i = 0; i < ct.length; i += 2) {
        if (sim_numberOfMolecules[i / 2] > 0) {
          validateCrosstalk(i, i + 1);
        }
      }
    } catch (final IllegalArgumentException ex) {
      IJ.error(title, ex.getMessage());
      return false;
    }

    return true;
  }

  private void actionPerformed(@SuppressWarnings("unused") ActionEvent event) {
    final ImagePlus imp = WindowManager.getImage(results.getName() + " " + title);
    if (imp == null || output == null) {
      return;
    }

    // List the ROIs
    final Roi imageRoi = imp.getRoi();
    if (imageRoi == null || !imageRoi.isArea()) {
      return;
    }
    Roi[] rois;
    if (imageRoi instanceof ShapeRoi) {
      rois = ((ShapeRoi) imageRoi).getRois();
    } else {
      rois = new Roi[] {imageRoi};
    }

    for (int i = 0; i < rois.length; i++) {
      drawLoop(imp, rois[i], i + 1);
    }
  }

  private void drawLoop(ImagePlus imp, Roi roi, int number) {
    if (!roi.isArea()) {
      return;
      // System.out.println(roi);
    }

    // Map the ROI to a crop of the results set
    final Rectangle roiBounds = roi.getBounds();
    final Rectangle resultsBounds = results.getBounds(true);

    //@formatter:off
    final Rectangle2D.Double r = new Rectangle2D.Double(
        resultsBounds.width * (double)roiBounds.x / imp.getWidth(),
        resultsBounds.height * (double)roiBounds.y / imp.getHeight(),
        resultsBounds.width * (double)roiBounds.width / imp.getWidth(),
        resultsBounds.height * (double)roiBounds.height / imp.getHeight());
    //@formatter:on
    // System.out.println(r);

    final int x = (int) r.getX();
    final int y = (int) r.getY();

    final int magnification = getMagnification();

    // For each result set crop out the localisation and construct an overlay
    final Overlay o = new Overlay();
    for (int i = 0; i < output.length; i++) {
      final Color color = colors[i];

      // The first result is the memory results
      final MemoryPeakResults results = (MemoryPeakResults) output[i].getOutput(0);
      results.forEach(DistanceUnit.PIXEL, new XyrResultProcedure() {
        @Override
        public void executeXyr(float xx, float yy, PeakResult result) {
          if (r.contains(xx, yy)) {
            add(o, (xx - x) * magnification, (yy - y) * magnification, color);
          }
        }
      });
    }

    // This results in a change of shape depending on where the roi is positioned
    int width = (int) Math.ceil(r.getWidth());
    int height = (int) Math.ceil(r.getHeight());
    width *= magnification;
    height *= magnification;
    final ImageProcessor ip = new ByteProcessor(width, height);

    final String title = imp.getTitle() + " Loop " + number;
    imp = WindowManager.getImage(title);
    if (imp == null) {
      imp = new ImagePlus(title, ip);
      imp.show();
      // ImageCanvas ic = imp.getWindow().getCanvas();
      // for (int i = 10; i-- > 0;)
      // ic.zoomIn(imp.getWidth() / 2, imp.getHeight() / 2);
      // ic.setMagnification(32);
    } else {
      imp.setProcessor(ip);
    }
    imp.setOverlay(o);
  }

  private int getMagnification() {
    magnification = magnificationChoice.getSelectedItem();
    try {
      return Integer.parseInt(magnification);
    } catch (final NumberFormatException ex) {
      return 1;
    }
  }

  private static void add(Overlay overlay, float x, float y, Color color) {
    final PointRoi roi = new PointRoi(x, y);
    roi.setStrokeColor(color);
    roi.setFillColor(color);
    roi.setPointType(1); // PointRoi.CROSSHAIR);
    roi.setSize(1); // PointRoi.TINY);
    overlay.add(roi);
  }
}
