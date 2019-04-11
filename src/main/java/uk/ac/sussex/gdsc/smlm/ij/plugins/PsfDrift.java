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

import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.NonBlockingExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.Plot2;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.TurboList;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrencyUtils;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitEngineSettings;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.ImagePSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.Offset;
import uk.ac.sussex.gdsc.smlm.data.config.PsfProtosHelper;
import uk.ac.sussex.gdsc.smlm.engine.FitConfiguration;
import uk.ac.sussex.gdsc.smlm.engine.FitEngineConfiguration;
import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.fitting.FunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.Gaussian2DFitter;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.ij.IJImageSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.ImagePsfHelper;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.model.ImagePsfModel;

import gnu.trove.list.array.TDoubleArrayList;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Line;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.plugin.PlugIn;

import org.apache.commons.lang3.concurrent.ConcurrentRuntimeException;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.rng.sampling.distribution.PoissonSampler;
import org.apache.commons.rng.simple.RandomSource;

import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Label;
import java.awt.TextField;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicReference;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Produces an drift curve for a PSF image using fitting.
 *
 * <p>The input images must be a z-stack of a PSF. These can be produced using the PSFCreator
 * plugin.
 */
public class PsfDrift implements PlugIn {
  private static final String TITLE = "PSF Drift";

  private ImagePlus imp;
  private ImagePSF psfSettings;
  private static FitConfiguration fitConfig = new FitConfiguration();

  private int centrePixel;
  private int total;
  private double[][] results;

  private final WindowOrganiser windowOrganiser = new WindowOrganiser();

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    String title;
    boolean useOffset;
    double scale;
    double zDepth;
    int gridSize;
    double recallLimit;
    int regionSize;
    boolean backgroundFitting;
    boolean offsetFitting;
    double startOffset;
    boolean comFitting;
    boolean useSampling;
    double photons;
    double photonLimit;
    int positionsToAverage;
    double smoothing;
    boolean updateCentre;
    boolean updateHwhm;

    Settings() {
      // Set defaults
      title = "";
      scale = 10;
      zDepth = 1000;
      gridSize = 10;
      recallLimit = 0.25;
      regionSize = 5;
      offsetFitting = true;
      startOffset = 0.5;
      comFitting = true;
      photons = 1000;
      photonLimit = 0.25;
      positionsToAverage = 5;
      smoothing = 0.1;
      updateCentre = true;
      updateHwhm = true;
    }

    Settings(Settings source) {
      title = source.title;
      useOffset = source.useOffset;
      scale = source.scale;
      zDepth = source.zDepth;
      gridSize = source.gridSize;
      recallLimit = source.recallLimit;
      regionSize = source.regionSize;
      backgroundFitting = source.backgroundFitting;
      offsetFitting = source.offsetFitting;
      startOffset = source.startOffset;
      comFitting = source.comFitting;
      useSampling = source.useSampling;
      photons = source.photons;
      photonLimit = source.photonLimit;
      positionsToAverage = source.positionsToAverage;
      smoothing = source.smoothing;
      updateCentre = source.updateCentre;
      updateHwhm = source.updateHwhm;
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
     * Save the settings. This can be called only once as it saves via a reference.
     */
    void save() {
      lastSettings.set(this);
    }
  }

  private static class Job {
    final int z;
    final double cx;
    final double cy;
    final int index;

    public Job(int z, double cx, double cy, int index) {
      this.z = z;
      this.cx = cx;
      this.cy = cy;
      this.index = index;
    }

    public Job() {
      this(0, 0, 0, -1);
    }

    @Override
    public String toString() {
      return String.format("z=%d, cx=%.2f, cy=%.2f", z, cx, cy);
    }
  }

  /**
   * Used to allow multi-threading of the fitting method.
   */
  private class Worker implements Runnable {
    volatile boolean finished;
    final ImagePsfModel psf;
    final BlockingQueue<Job> jobs;
    final FitConfiguration fitConfig2;
    final Ticker ticker;
    final double sx;
    final double sy;
    final double pixelPitch;
    final double[][] xy;
    final int width;
    final int w2;
    PoissonSampler poissonSampler;

    private double[] lb;
    private double[] ub;
    private double[] lc;
    private double[] uc;

    public Worker(BlockingQueue<Job> jobs, ImagePsfModel psf, int width, FitConfiguration fitConfig,
        Ticker ticker) {
      this.jobs = jobs;
      this.psf = psf.copy(null);
      this.fitConfig2 = fitConfig.createCopy();
      this.ticker = ticker;
      sx = fitConfig.getInitialXSd();
      sy = fitConfig.getInitialYSd();
      pixelPitch = psfSettings.getPixelSize() * settings.scale;
      xy = getStartPoints();
      this.width = width;
      w2 = width * width;
      if (settings.useSampling) {
        // TOOD: This could be updated to use an input RNG
        poissonSampler =
            new PoissonSampler(RandomSource.create(RandomSource.SPLIT_MIX_64), settings.photons);
      }

      createBounds();
    }

    @Override
    public void run() {
      try {
        while (true) {
          final Job job = jobs.take();
          if (job == null || job.index < 0) {
            break;
          }
          if (!finished) {
            // Only run if not finished to allow queue to be emptied
            run(job);
            ticker.tick();
          }
        }
      } catch (final InterruptedException ex) {
        ConcurrencyUtils.interruptAndThrowUncheckedIf(!finished, ex);
      } finally {
        finished = true;
      }
    }

    private void run(Job job) {
      if (IJ.escapePressed()) {
        finished = true;
        return;
      }

      final double cx = centrePixel + job.cx;
      final double cy = centrePixel + job.cy;

      // Draw the PSF
      final double[] data = new double[w2];
      if (poissonSampler != null) {
        final int p = poissonSampler.sample();
        psf.sample3D(data, width, width, p, cx, cy, job.z);
      } else {
        psf.create3D(data, width, width, settings.photons, cx, cy, job.z, false);
      }

      // Fit the PSF. Do this from different start positions.

      // Get the background and signal estimate
      final double background =
          (settings.backgroundFitting) ? Gaussian2DFitter.getBackground(data, width, width, 1) : 0;
      final double signal = BenchmarkFit.getSignal(data, background);

      if (settings.comFitting) {
        // Get centre-of-mass estimate, then subtract the centre that will be added later
        BenchmarkFit.getCentreOfMass(data, width, width, xy[xy.length - 1]);
        xy[xy.length - 1][0] -= cx;
        xy[xy.length - 1][1] -= cy;
      }

      final double[] initialParams = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
      initialParams[Gaussian2DFunction.BACKGROUND] = background;
      initialParams[Gaussian2DFunction.SIGNAL] = signal;
      initialParams[Gaussian2DFunction.X_SD] = sx;
      initialParams[Gaussian2DFunction.Y_SD] = sy;

      int resultPosition = job.index;
      for (final double[] centre : xy) {
        // Do fitting
        final double[] params = initialParams.clone();
        params[Gaussian2DFunction.X_POSITION] = cx + centre[0];
        params[Gaussian2DFunction.Y_POSITION] = cy + centre[1];
        fitConfig2.initialise(1, width, width, params);
        final FunctionSolver solver = fitConfig2.getFunctionSolver();
        if (solver.isBounded()) {
          setBounds(solver);
        } else if (solver.isConstrained()) {
          setConstraints(solver);
        }
        final FitStatus status = solver.fit(data, null, params, null);
        // Account for 0.5 pixel offset during fitting
        params[Gaussian2DFunction.X_POSITION] += 0.5;
        params[Gaussian2DFunction.Y_POSITION] += 0.5;
        if (isValid(status, params, width)) {
          // XXX Decide what results are needed for analysis
          // Store all the results for later analysis
          // results[resultPosition] = params;
          // Store only the drift
          results[resultPosition] =
              new double[] {pixelPitch * (params[Gaussian2DFunction.X_POSITION] - cx),
                  pixelPitch * (params[Gaussian2DFunction.Y_POSITION] - cy), job.z};
          // System.out.printf("Fit " + job + ". %f,%f\n", results[resultPosition][0],
          // results[resultPosition][1]);
        } else {
          // System.out.println("Failed to fit " + job + ". " + status);
        }
        resultPosition += total;
      }
    }

    private void setBounds(FunctionSolver solver) {
      solver.setBounds(lb, ub);
    }

    private void createBounds() {
      if (ub == null) {
        ub = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
        lb = new double[ub.length];

        // Background could be zero so always have an upper limit
        ub[Gaussian2DFunction.BACKGROUND] = 1;
        lb[Gaussian2DFunction.SIGNAL] = settings.photons * settings.photonLimit;
        ub[Gaussian2DFunction.SIGNAL] = settings.photons * 2;
        ub[Gaussian2DFunction.X_POSITION] = width;
        ub[Gaussian2DFunction.Y_POSITION] = width;
        lb[Gaussian2DFunction.ANGLE] = -Math.PI;
        ub[Gaussian2DFunction.ANGLE] = Math.PI;
        lb[Gaussian2DFunction.Z_POSITION] = Double.NEGATIVE_INFINITY;
        ub[Gaussian2DFunction.Z_POSITION] = Double.POSITIVE_INFINITY;
        final double wf = 1.5;
        lb[Gaussian2DFunction.X_SD] = sx / wf;
        ub[Gaussian2DFunction.X_SD] = sx * 5;
        lb[Gaussian2DFunction.Y_SD] = sy / wf;
        ub[Gaussian2DFunction.Y_SD] = sy * 5;
      }
    }

    private void setConstraints(FunctionSolver solver) {
      if (uc == null) {
        lc = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
        uc = new double[lc.length];
        Arrays.fill(lc, Float.NEGATIVE_INFINITY);
        Arrays.fill(uc, Float.POSITIVE_INFINITY);
        lc[Gaussian2DFunction.BACKGROUND] = 0;
        lc[Gaussian2DFunction.SIGNAL] = 0;
      }
      solver.setConstraints(lc, uc);
    }

    private boolean isValid(FitStatus status, double[] params, int size) {
      if (status != FitStatus.OK) {
        return false;
      }

      // Reject fits that are outside the bounds of the data
      if (params[Gaussian2DFunction.X_POSITION] < 0 || params[Gaussian2DFunction.Y_POSITION] < 0
          || params[Gaussian2DFunction.X_POSITION] > size
          || params[Gaussian2DFunction.Y_POSITION] > size) {
        return false;
      }

      // Reject fits that do not correctly estimate the signal
      if (params[Gaussian2DFunction.SIGNAL] < lb[Gaussian2DFunction.SIGNAL]
          || params[Gaussian2DFunction.SIGNAL] > ub[Gaussian2DFunction.SIGNAL]) {
        return false;
      }

      // Reject fits that have a background too far from zero
      // TODO - configure this better
      // if (params[Gaussian2DFunction.BACKGROUND] < -10 || params[Gaussian2DFunction.BACKGROUND] >
      // 10)
      // {
      // return false;
      // }

      // Q. Should we do width bounds checking?
      if (fitConfig2.isXSdFitting()
          && (params[Gaussian2DFunction.X_SD] < lb[Gaussian2DFunction.X_SD]
              || params[Gaussian2DFunction.X_SD] > ub[Gaussian2DFunction.X_SD])) {
        return false;
      }
      return !(fitConfig2.isYSdFitting()
          && (params[Gaussian2DFunction.Y_SD] < lb[Gaussian2DFunction.Y_SD]
              || params[Gaussian2DFunction.Y_SD] > ub[Gaussian2DFunction.Y_SD]));
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if ("hwhm".equals(arg)) {
      showHwhm();
      return;
    }

    // Build a list of suitable images
    final List<String> titles = createImageList(true);

    if (titles.isEmpty()) {
      IJ.error(TITLE, "No suitable PSF images");
      return;
    }

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    settings = Settings.load();
    gd.addMessage("Select the input PSF image");
    gd.addChoice("PSF", titles.toArray(new String[titles.size()]), settings.title);
    gd.addCheckbox("Use_offset", settings.useOffset);
    gd.addNumericField("Scale", settings.scale, 2);
    gd.addNumericField("z_depth", settings.zDepth, 2, 6, "nm");
    gd.addNumericField("Grid_size", settings.gridSize, 0);
    gd.addSlider("Recall_limit", 0.01, 1, settings.recallLimit);

    gd.addSlider("Region_size", 2, 20, settings.regionSize);
    gd.addCheckbox("Background_fitting", settings.backgroundFitting);
    PeakFit.addPsfOptions(gd, fitConfig);
    gd.addChoice("Fit_solver", SettingsManager.getFitSolverNames(),
        fitConfig.getFitSolver().ordinal());
    // We need these to set bounds for any bounded fitters
    gd.addSlider("Min_width_factor", 0, 0.99, fitConfig.getMinWidthFactor());
    gd.addSlider("Width_factor", 1, 4.5, fitConfig.getMaxWidthFactor());
    gd.addCheckbox("Offset_fit", settings.offsetFitting);
    gd.addNumericField("Start_offset", settings.startOffset, 3);
    gd.addCheckbox("Include_CoM_fit", settings.comFitting);
    gd.addCheckbox("Use_sampling", settings.useSampling);
    gd.addNumericField("Photons", settings.photons, 0);
    gd.addSlider("Photon_limit", 0, 1, settings.photonLimit);
    gd.addSlider("Smoothing", 0, 0.5, settings.smoothing);

    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }

    settings.title = gd.getNextChoice();
    settings.useOffset = gd.getNextBoolean();
    settings.scale = gd.getNextNumber();
    settings.zDepth = gd.getNextNumber();
    settings.gridSize = (int) gd.getNextNumber();
    settings.recallLimit = gd.getNextNumber();
    settings.regionSize = (int) Math.abs(gd.getNextNumber());
    settings.backgroundFitting = gd.getNextBoolean();
    fitConfig.setPsfType(PeakFit.getPsfTypeValues()[gd.getNextChoiceIndex()]);
    fitConfig.setFitSolver(gd.getNextChoiceIndex());
    fitConfig.setMinWidthFactor(gd.getNextNumber());
    fitConfig.setMaxWidthFactor(gd.getNextNumber());
    settings.offsetFitting = gd.getNextBoolean();
    settings.startOffset = Math.abs(gd.getNextNumber());
    settings.comFitting = gd.getNextBoolean();
    settings.useSampling = gd.getNextBoolean();
    settings.photons = Math.abs(gd.getNextNumber());
    settings.photonLimit = Math.abs(gd.getNextNumber());
    settings.smoothing = Math.abs(gd.getNextNumber());
    settings.save();

    gd.collectOptions();

    if (!settings.comFitting && !settings.offsetFitting) {
      IJ.error(TITLE, "No initial fitting positions");
      return;
    }

    if (settings.regionSize < 1) {
      settings.regionSize = 1;
    }

    if (gd.invalidNumber()) {
      return;
    }

    imp = WindowManager.getImage(settings.title);
    if (imp == null) {
      IJ.error(TITLE, "No PSF image for image: " + settings.title);
      return;
    }
    psfSettings = getPsfSettings(imp);
    if (psfSettings == null) {
      IJ.error(TITLE, "No PSF settings for image: " + settings.title);
      return;
    }

    // Configure the fit solver. We must wrap the settings with a
    // FitEngineConfiguration to pass to the PeakFit method
    final FitEngineSettings fitEngineSettings = FitProtosHelper.defaultFitEngineSettings;
    final FitEngineConfiguration config = new FitEngineConfiguration(fitEngineSettings,
        SettingsManager.readCalibration(0), PsfProtosHelper.defaultOneAxisGaussian2DPSF);
    config.getFitConfiguration().setFitSettings(fitConfig.getFitSettings());
    if (!PeakFit.configurePsfModel(config)) {
      return;
    }
    if (!PeakFit.configureFitSolver(config, IJImageSource.getBounds(imp), null,
        PeakFit.FLAG_NO_SAVE)) {
      return;
    }
    fitConfig = config.getFitConfiguration();

    computeDrift();
  }

  private void computeDrift() {
    // Create a grid of XY offset positions between 0-1 for PSF insert
    final double[] grid = new double[settings.gridSize];
    for (int i = 0; i < grid.length; i++) {
      grid[i] = (double) i / settings.gridSize;
    }

    // Configure fitting region
    final int w = 2 * settings.regionSize + 1;
    centrePixel = w / 2;

    // Check region size using the image PSF
    final double newPsfWidth = imp.getWidth() / settings.scale;
    if (Math.ceil(newPsfWidth) > w) {
      ImageJUtils.log(TITLE + ": Fitted region size (%d) is smaller than the scaled PSF (%.1f)", w,
          newPsfWidth);
    }

    // Create robust PSF fitting settings
    final double a = psfSettings.getPixelSize() * settings.scale;
    final double sa = PsfCalculator.squarePixelAdjustment(
        psfSettings.getPixelSize() * (psfSettings.getFwhm() / Gaussian2DFunction.SD_TO_FWHM_FACTOR),
        a);
    fitConfig.setInitialPeakStdDev(sa / a);
    fitConfig.setBackgroundFitting(settings.backgroundFitting);
    fitConfig.setNotSignalFitting(false);
    fitConfig.setComputeDeviations(false);
    fitConfig.setDisableSimpleFilter(true);

    // Create the PSF over the desired z-depth
    final int depth = (int) Math.round(settings.zDepth / psfSettings.getPixelDepth());
    int startSlice = psfSettings.getCentreImage() - depth;
    int endSlice = psfSettings.getCentreImage() + depth;
    final int nSlices = imp.getStackSize();
    startSlice = MathUtils.clip(1, nSlices, startSlice);
    endSlice = MathUtils.clip(1, nSlices, endSlice);

    final ImagePsfModel psf = createImagePsf(startSlice, endSlice, settings.scale);

    final int minz = startSlice - psfSettings.getCentreImage();
    final int maxz = endSlice - psfSettings.getCentreImage();

    final int nZ = maxz - minz + 1;
    final int gridSize2 = grid.length * grid.length;
    total = nZ * gridSize2;

    // Store all the fitting results
    final int nStartPoints = getNumberOfStartPoints();
    results = new double[total * nStartPoints][];

    // TODO - Add ability to iterate this, adjusting the current offset in the PSF
    // each iteration

    // Create a pool of workers
    final int threadCount = Prefs.getThreads();
    final Ticker ticker = ImageJUtils.createTicker(total, threadCount, "Fitting...");
    final BlockingQueue<Job> jobs = new ArrayBlockingQueue<>(threadCount * 2);
    final List<Thread> threads = new LinkedList<>();
    for (int i = 0; i < threadCount; i++) {
      final Worker worker = new Worker(jobs, psf, w, fitConfig, ticker);
      final Thread t = new Thread(worker);
      threads.add(t);
      t.start();
    }

    // Fit
    outer: for (int z = minz, i = 0; z <= maxz; z++) {
      for (int x = 0; x < grid.length; x++) {
        for (int y = 0; y < grid.length; y++, i++) {
          if (IJ.escapePressed()) {
            break outer;
          }
          put(jobs, new Job(z, grid[x], grid[y], i));
        }
      }
    }

    // If escaped pressed then do not need to stop the workers, just return
    if (ImageJUtils.isInterrupted()) {
      ImageJUtils.finished();
      return;
    }

    // Finish all the worker threads by passing in a null job
    for (int i = 0; i < threads.size(); i++) {
      put(jobs, new Job());
    }

    // Wait for all to finish
    for (int i = 0; i < threads.size(); i++) {
      try {
        threads.get(i).join();
      } catch (final InterruptedException ex) {
        Thread.currentThread().interrupt();
        throw new ConcurrentRuntimeException("Unexpected interrupt", ex);
      }
    }
    threads.clear();

    ImageJUtils.finished();

    // Plot the average and SE for the drift curve
    // Plot the recall
    final double[] zPosition = new double[nZ];
    final double[] avX = new double[nZ];
    final double[] seX = new double[nZ];
    final double[] avY = new double[nZ];
    final double[] seY = new double[nZ];
    final double[] recall = new double[nZ];
    for (int z = minz, i = 0; z <= maxz; z++, i++) {
      final Statistics statsX = new Statistics();
      final Statistics statsY = new Statistics();
      for (int s = 0; s < nStartPoints; s++) {
        int resultPosition = i * gridSize2 + s * total;
        final int endResultPosition = resultPosition + gridSize2;
        while (resultPosition < endResultPosition) {
          if (results[resultPosition] != null) {
            statsX.add(results[resultPosition][0]);
            statsY.add(results[resultPosition][1]);
          }
          resultPosition++;
        }
      }
      zPosition[i] = z * psfSettings.getPixelDepth();
      avX[i] = statsX.getMean();
      seX[i] = statsX.getStandardError();
      avY[i] = statsY.getMean();
      seY[i] = statsY.getStandardError();
      recall[i] = (double) statsX.getN() / (nStartPoints * gridSize2);
    }

    // Find the range from the z-centre above the recall limit
    int centre = 0;
    for (int slice = startSlice, i = 0; slice <= endSlice; slice++, i++) {
      if (slice == psfSettings.getCentreImage()) {
        centre = i;
        break;
      }
    }
    if (recall[centre] < settings.recallLimit) {
      return;
    }
    int start = centre;
    int end = centre;
    for (int i = centre; i-- > 0;) {
      if (recall[i] < settings.recallLimit) {
        break;
      }
      start = i;
    }
    for (int i = centre; ++i < recall.length;) {
      if (recall[i] < settings.recallLimit) {
        break;
      }
      end = i;
    }

    final int iterations = 1;
    LoessInterpolator loess = null;
    if (settings.smoothing > 0) {
      loess = new LoessInterpolator(settings.smoothing, iterations);
    }

    final double[][] smoothx =
        displayPlot("Drift X", "X (nm)", zPosition, avX, seX, loess, start, end);
    final double[][] smoothy =
        displayPlot("Drift Y", "Y (nm)", zPosition, avY, seY, loess, start, end);
    displayPlot("Recall", "Recall", zPosition, recall, null, null, start, end);

    windowOrganiser.tile();

    // Ask the user if they would like to store them in the image
    final GenericDialog gd = new GenericDialog(TITLE);
    gd.enableYesNoCancel();
    gd.hideCancelButton();
    startSlice = psfSettings.getCentreImage() - (centre - start);
    endSlice = psfSettings.getCentreImage() + (end - centre);
    ImageJUtils.addMessage(gd,
        "Save the drift to the PSF?\n \nSlices %d (%s nm) - %d (%s nm) above recall limit",
        startSlice, MathUtils.rounded(zPosition[start]), endSlice,
        MathUtils.rounded(zPosition[end]));
    gd.addMessage("Optionally average the end points to set drift outside the limits.\n"
        + "(Select zero to ignore)");
    gd.addSlider("Number_of_points", 0, 10, settings.positionsToAverage);
    gd.showDialog();
    if (gd.wasOKed()) {
      settings.positionsToAverage = Math.abs((int) gd.getNextNumber());
      final Map<Integer, Offset> oldOffset = psfSettings.getOffsetsMap();
      final boolean useOldOffset = settings.useOffset && !oldOffset.isEmpty();
      final TurboList<double[]> offset = new TurboList<>();
      final double pitch = psfSettings.getPixelSize();
      int index = 0;
      for (int i = start, slice = startSlice; i <= end; slice++, i++) {
        index = findCentre(zPosition[i], smoothx, index);
        if (index == -1) {
          ImageJUtils.log("Failed to find the offset for depth %.2f", zPosition[i]);
          continue;
        }
        // The offset should store the difference to the centre in pixels so divide by the pixel
        // pitch
        double cx = smoothx[1][index] / pitch;
        double cy = smoothy[1][index] / pitch;
        if (useOldOffset) {
          final Offset o = oldOffset.get(slice);
          if (o != null) {
            cx += o.getCx();
            cy += o.getCy();
          }
        }
        offset.add(new double[] {slice, cx, cy});
      }
      addMissingOffsets(startSlice, endSlice, nSlices, offset);
      final Offset.Builder offsetBuilder = Offset.newBuilder();
      final ImagePSF.Builder imagePsfBuilder = psfSettings.toBuilder();
      for (final double[] o : offset) {
        final int slice = (int) o[0];
        offsetBuilder.setCx(o[1]);
        offsetBuilder.setCy(o[2]);
        imagePsfBuilder.putOffsets(slice, offsetBuilder.build());
      }
      imagePsfBuilder.putNotes(TITLE, String.format("Solver=%s, Region=%d",
          PeakFit.getSolverName(fitConfig), settings.regionSize));
      imp.setProperty("Info", ImagePsfHelper.toString(imagePsfBuilder));
    }
  }

  private static int findCentre(double depth, double[][] smoothx, int index) {
    while (index < smoothx[0].length) {
      if (smoothx[0][index] == depth) {
        return index;
      }
      index++;
    }
    return -1;
  }

  private void addMissingOffsets(int startSlice, int endSlice, int slices,
      TurboList<double[]> offset) {
    // Add an offset for the remaining slices
    if (settings.positionsToAverage > 0) {
      double cx = 0;
      double cy = 0;
      int count = 0;
      for (int i = 0; count < settings.positionsToAverage && i < offset.size(); i++, count++) {
        cx += offset.get(i)[1];
        cy += offset.get(i)[2];
      }
      cx /= count;
      cy /= count;
      double cx2 = 0;
      double cy2 = 0;
      double n2 = 0;
      for (int i = offset.size(); n2 < settings.positionsToAverage && i-- > 0;) {
        cx2 += offset.get(i)[1];
        cy2 += offset.get(i)[2];
        n2++;
      }
      cx2 /= n2;
      cy2 /= n2;

      for (int slice = 1; slice < startSlice; slice++) {
        offset.add(new double[] {slice, cx, cy});
      }
      for (int slice = endSlice + 1; slice <= slices; slice++) {
        offset.add(new double[] {slice, cx2, cy2});
      }
      Collections.sort(offset, (a1, a2) -> Double.compare(a1[0], a2[0]));
    }
  }

  private double[][] displayPlot(String title, String yLabel, double[] x, double[] y, double[] se,
      LoessInterpolator loess, int start, int end) {
    // Extract non NaN numbers
    double[] newX = new double[x.length];
    double[] newY = new double[x.length];
    int count = 0;
    for (int i = 0; i < x.length; i++) {
      if (!Double.isNaN(y[i])) {
        newX[count] = x[i];
        newY[count] = y[i];
        count++;
      }
    }
    newX = Arrays.copyOf(newX, count);
    newY = Arrays.copyOf(newY, count);

    title = TITLE + " " + title;
    final Plot2 plot = new Plot2(title, "z (nm)", yLabel);
    final double[] limitsx = MathUtils.limits(x);
    double[] limitsy = new double[2];
    if (se != null) {
      if (count > 0) {
        limitsy = new double[] {newY[0] - se[0], newY[0] + se[0]};
        for (int i = 1; i < newY.length; i++) {
          limitsy[0] = MathUtils.min(limitsy[0], newY[i] - se[i]);
          limitsy[1] = MathUtils.max(limitsy[1], newY[i] + se[i]);
        }
      }
    } else if (count > 0) {
      limitsy = MathUtils.limits(newY);
    }
    final double rangex = Math.max(0.05 * (limitsx[1] - limitsx[0]), 0.1);
    final double rangey = Math.max(0.05 * (limitsy[1] - limitsy[0]), 0.1);
    plot.setLimits(limitsx[0] - rangex, limitsx[1] + rangex, limitsy[0] - rangey,
        limitsy[1] + rangey);

    if (loess == null) {
      addPoints(plot, Plot.LINE, newX, newY, x[start], x[end]);
    } else {
      addPoints(plot, Plot.DOT, newX, newY, x[start], x[end]);
      newY = loess.smooth(newX, newY);
      addPoints(plot, Plot.LINE, newX, newY, x[start], x[end]);
    }
    if (se != null) {
      plot.setColor(Color.magenta);
      for (int i = 0; i < x.length; i++) {
        if (!Double.isNaN(y[i])) {
          plot.drawLine(x[i], y[i] - se[i], x[i], y[i] + se[i]);
        }
      }

      // Draw the start and end lines for the valid range
      plot.setColor(Color.green);
      plot.drawLine(x[start], limitsy[0], x[start], limitsy[1]);
      plot.drawLine(x[end], limitsy[0], x[end], limitsy[1]);
    } else {
      // draw a line for the recall limit
      plot.setColor(Color.magenta);
      plot.drawLine(limitsx[0] - rangex, settings.recallLimit, limitsx[1] + rangex,
          settings.recallLimit);
    }
    ImageJUtils.display(title, plot, windowOrganiser);

    return new double[][] {newX, newY};
  }

  private static void addPoints(Plot plot, int shape, double[] x, double[] y, double lower,
      double upper) {
    if (x.length == 0) {
      return;
      // Split the line into three:
      // 1. All points up to and including lower
      // 2. All points between lower and upper inclusive
      // 3. All point from upper upwards
    }

    // Plot the main curve first
    addPoints(plot, shape, x, y, lower, upper, Color.blue);
    // Then plot the others
    addPoints(plot, shape, x, y, x[0], lower, Color.red);
    addPoints(plot, shape, x, y, upper, x[x.length - 1], Color.red);
  }

  private static void addPoints(Plot plot, int shape, double[] x, double[] y, double lower,
      double upper, Color color) {
    double[] x2 = new double[x.length];
    double[] y2 = new double[y.length];
    int count = 0;
    for (int i = 0; i < x.length; i++) {
      if (x[i] >= lower && x[i] <= upper) {
        x2[count] = x[i];
        y2[count] = y[i];
        count++;
      }
    }
    if (count == 0) {
      return;
    }
    x2 = Arrays.copyOf(x2, count);
    y2 = Arrays.copyOf(y2, count);
    plot.setColor(color);
    plot.addPoints(x2, y2, shape);
  }

  private ImagePsfModel createImagePsf(int lower, int upper, double scale) {
    final int zCentre = psfSettings.getCentreImage();

    final double unitsPerPixel = 1.0 / scale;
    final double unitsPerSlice = 1; // So we can move from -depth to depth

    // Extract data uses index not slice number as arguments so subtract 1
    final double noiseFraction = 1e-3;
    final float[][] image = CreateData.extractImageStack(imp, lower - 1, upper - 1);
    final ImagePsfModel model =
        new ImagePsfModel(image, zCentre - lower, unitsPerPixel, unitsPerSlice, noiseFraction);

    // Add the calibrated centres
    final Map<Integer, Offset> oldOffset = psfSettings.getOffsetsMap();
    if (settings.useOffset && !oldOffset.isEmpty()) {
      final int sliceOffset = lower;
      for (final Entry<Integer, Offset> entry : oldOffset.entrySet()) {
        model.setRelativeCentre(entry.getKey() - sliceOffset, entry.getValue().getCx(),
            entry.getValue().getCy());
      }
    } else {
      // Use the CoM if present
      final double cx = psfSettings.getXCentre();
      final double cy = psfSettings.getYCentre();
      if (cx != 0 || cy != 0) {
        for (int slice = 0; slice < image.length; slice++) {
          model.setCentre(slice, cx, cy);
        }
      }
    }

    return model;
  }

  private static void put(BlockingQueue<Job> jobs, Job job) {
    try {
      jobs.put(job);
    } catch (final InterruptedException ex) {
      Logger.getLogger(PsfDrift.class.getName()).log(Level.SEVERE, "Unexpected interruption", ex);
      Thread.currentThread().interrupt();
    }
  }

  /**
   * Gets the starting points for the fitting.
   *
   * @return The starting points for the fitting
   */
  private double[][] getStartPoints() {
    final double[][] xy = new double[getNumberOfStartPoints()][];
    int ii = 0;

    if (settings.offsetFitting) {
      if (settings.startOffset == 0) {
        xy[ii++] = new double[] {0, 0};
      } else {
        // Fit using region surrounding the point. Use -1,-1 : -1:1 : 1,-1 : 1,1 directions at
        // startOffset pixels total distance
        final double distance = Math.sqrt(settings.startOffset * settings.startOffset * 0.5);

        for (int x = -1; x <= 1; x += 2) {
          for (int y = -1; y <= 1; y += 2) {
            xy[ii++] = new double[] {x * distance, y * distance};
          }
        }
      }
    }
    // Add space for centre-of-mass at the end of the array
    if (settings.comFitting) {
      xy[ii] = new double[2];
    }
    return xy;
  }

  private int getNumberOfStartPoints() {
    int points = (settings.offsetFitting) ? 1 : 0;
    if (settings.startOffset > 0) {
      points *= 4;
    }
    return (settings.comFitting) ? points + 1 : points;
  }

  private static List<String> createImageList(boolean requireFwhm) {
    final List<String> titles = new LinkedList<>();
    final int[] ids = WindowManager.getIDList();
    if (ids != null) {
      for (final int id : ids) {
        final ImagePlus imp = WindowManager.getImage(id);
        if (imp != null) {
          // Image must be greyscale
          if (imp.getType() == ImagePlus.GRAY8 || imp.getType() == ImagePlus.GRAY16
              || imp.getType() == ImagePlus.GRAY32) {
            // Image must be square and a stack of a single channel
            if (imp.getWidth() == imp.getHeight() && imp.getNChannels() == 1) {
              // Check if these are PSF images created by the SMLM plugins
              final ImagePSF psfSettings = getPsfSettings(imp);
              if (psfSettings != null) {
                if (psfSettings.getCentreImage() <= 0) {
                  ImageJUtils
                      .log(TITLE + ": Unknown PSF z-centre setting for image: " + imp.getTitle());
                  continue;
                }
                if (psfSettings.getPixelSize() <= 0) {
                  ImageJUtils
                      .log(TITLE + ": Unknown PSF nm/pixel setting for image: " + imp.getTitle());
                  continue;
                }
                if (psfSettings.getPixelDepth() <= 0) {
                  ImageJUtils
                      .log(TITLE + ": Unknown PSF nm/slice setting for image: " + imp.getTitle());
                  continue;
                }
                if (requireFwhm && psfSettings.getFwhm() <= 0) {
                  ImageJUtils
                      .log(TITLE + ": Unknown PSF FWHM setting for image: " + imp.getTitle());
                  continue;
                }

                titles.add(imp.getTitle());
              }
            }
          }
        }
      }
    }
    return titles;
  }

  /**
   * Gets the PSF settings from the Info property.
   *
   * @param imp the imp
   * @return the PSF settings
   * @see ImagePlus#getProperty(String)
   */
  public static ImagePSF getPsfSettings(ImagePlus imp) {
    final Object info = imp.getProperty("Info");
    if (info != null) {
      return ImagePsfHelper.fromString(info.toString());
    }
    return null;
  }

  private void showHwhm() {
    // Build a list of suitable images
    final List<String> titles = createImageList(false);

    if (titles.isEmpty()) {
      IJ.error(TITLE, "No suitable PSF images");
      return;
    }

    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addMessage("Approximate the volume of the PSF as a Gaussian and\n"
        + "compute the equivalent Gaussian width.");
    settings = Settings.load();
    gd.addChoice("PSF", titles.toArray(new String[titles.size()]), settings.title);
    gd.addCheckbox("Use_offset", settings.useOffset);
    gd.addSlider("Smoothing", 0, 0.5, settings.smoothing);

    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }

    settings.title = gd.getNextChoice();
    settings.useOffset = gd.getNextBoolean();
    settings.smoothing = gd.getNextNumber();
    settings.save();

    imp = WindowManager.getImage(settings.title);
    if (imp == null) {
      IJ.error(TITLE, "No PSF image for image: " + settings.title);
      return;
    }
    psfSettings = getPsfSettings(imp);
    if (psfSettings == null) {
      IJ.error(TITLE, "No PSF settings for image: " + settings.title);
      return;
    }

    final int size = imp.getStackSize();
    final ImagePsfModel psf = createImagePsf(1, size, 1);

    final double[] w0 = psf.getAllHwhm0();
    final double[] w1 = psf.getAllHwhm1();

    // Get current centre
    final int centre = psfSettings.getCentreImage();

    // Extract valid values (some can be NaN)
    double[] slice0;
    double[] slice1;
    double[] sw0 = new double[w0.length];
    double[] sw1 = new double[w1.length];
    {
      final TDoubleArrayList s0 = new TDoubleArrayList(w0.length);
      final TDoubleArrayList s1 = new TDoubleArrayList(w0.length);
      int c0 = 0;
      int c1 = 0;
      for (int i = 0; i < w0.length; i++) {
        if (Double.isFinite(w0[i])) {
          s0.add(i + 1);
          sw0[c0++] = w0[i];
        }
        if (Double.isFinite(w1[i])) {
          s1.add(i + 1);
          sw1[c1++] = w1[i];
        }
      }
      if (c0 == 0 && c1 == 0) {
        IJ.error(TITLE, "No computed HWHM for image: " + settings.title);
        return;
      }
      slice0 = s0.toArray();
      sw0 = Arrays.copyOf(sw0, c0);
      slice1 = s1.toArray();
      sw1 = Arrays.copyOf(sw1, c1);
    }

    // Smooth
    if (settings.smoothing > 0) {
      final LoessInterpolator loess = new LoessInterpolator(settings.smoothing, 1);
      sw0 = loess.smooth(slice0, sw0);
      sw1 = loess.smooth(slice1, sw1);
    }

    // int newCentre = 0;
    // double minW = Double.POSITIVE_INFINITY;
    final TDoubleArrayList minWx = new TDoubleArrayList();
    final TDoubleArrayList minWy = new TDoubleArrayList();
    for (int i = 0; i < w0.length; i++) {
      double weight = 0;
      if (Double.isFinite(w0[i])) {
        if (Double.isFinite(w1[i])) {
          weight = w0[i] * w1[i];
        } else {
          weight = w0[i] * w0[i];
        }
      } else if (Double.isFinite(w1[i])) {
        weight = w1[i] * w1[i];
      }

      if (weight != 0) {
        minWx.add(i + 1);
        minWy.add(Math.sqrt(weight));
      }
    }

    // Smooth the combined line
    final double[] cx = minWx.toArray();
    double[] cy = minWy.toArray();
    if (settings.smoothing > 0) {
      final LoessInterpolator loess = new LoessInterpolator(settings.smoothing, 1);
      cy = loess.smooth(cx, cy);
    }
    final int newCentre = SimpleArrayUtils.findMinIndex(cy);

    // Convert to FWHM
    final double fwhm = psfSettings.getFwhm();

    // Widths are in pixels
    final String title = TITLE + " HWHM";
    final Plot plot = new Plot(title, "Slice", "HWHM (px)");
    double[] limits = MathUtils.limits(sw0);
    limits = MathUtils.limits(limits, sw1);
    final double maxY = limits[1] * 1.05;
    plot.setLimits(1, size, 0, maxY);
    plot.setColor(Color.red);
    plot.addPoints(slice0, sw0, Plot.LINE);
    plot.setColor(Color.blue);
    plot.addPoints(slice1, sw1, Plot.LINE);
    plot.setColor(Color.magenta);
    plot.addPoints(cx, cy, Plot.LINE);
    plot.setColor(Color.black);
    plot.addLabel(0, 0, "X=red; Y=blue, Combined=Magenta");
    final PlotWindow pw = ImageJUtils.display(title, plot);

    // Show a non-blocking dialog to allow the centre to be updated ...
    // Add a label and dynamically update when the centre is moved.
    final NonBlockingExtendedGenericDialog gd2 = new NonBlockingExtendedGenericDialog(TITLE);
    final double scale = psfSettings.getPixelSize();
    //@formatter:off
    ImageJUtils.addMessage(gd2,
        "Update the PSF information?\n \n" +
        "Current z-centre = %d, FHWM = %s px (%s nm)\n",
        centre, MathUtils.rounded(fwhm), MathUtils.rounded(fwhm * scale));
    //@formatter:on
    gd2.addSlider("z-centre", cx[0], cx[cx.length - 1], newCentre);
    final TextField tf = gd2.getLastTextField();
    gd2.addMessage("");
    gd2.addAndGetButton("Reset", event -> tf.setText(Integer.toString(newCentre)));
    final Label label = gd2.getLastLabel();
    gd2.addCheckbox("Update_centre", settings.updateCentre);
    gd2.addCheckbox("Update_HWHM", settings.updateHwhm);
    gd2.enableYesNoCancel();
    gd2.hideCancelButton();
    final UpdateDialogListener dl =
        new UpdateDialogListener(cx, cy, maxY, newCentre, scale, pw, label);
    gd2.addDialogListener(dl);
    gd2.showDialog();
    if (gd2.wasOKed() && (settings.updateCentre || settings.updateHwhm)) {
      final ImagePSF.Builder b = psfSettings.toBuilder();
      if (settings.updateCentre) {
        b.setCentreImage(dl.centre);
      }
      if (settings.updateHwhm) {
        b.setFwhm(dl.getFwhm());
      }
      imp.setProperty("Info", ImagePsfHelper.toString(b));
    }
  }

  private class UpdateDialogListener implements DialogListener {
    int offset;
    double[] cy;
    double maxY;
    int centre;
    double scale;
    PlotWindow pw;
    Label label;
    boolean drawing;

    UpdateDialogListener(double[] cx, double[] cy, double maxY, int centre, double scale,
        PlotWindow pw, Label label) {
      offset = (int) cx[0];
      this.cy = cy;

      // Interpolate missing values
      final int upper = cx.length - 1;
      if (cx[upper] - cx[0] != upper) {
        final LinearInterpolator in = new LinearInterpolator();
        final PolynomialSplineFunction f = in.interpolate(cx, cy);
        cx = SimpleArrayUtils.newArray(upper + 1, cx[0], 1.0);
        this.cy = new double[cx.length];
        for (int i = 0; i < cx.length; i++) {
          this.cy[i] = f.value(cx[i]);
        }
      }

      this.maxY = maxY;
      this.centre = centre;
      this.scale = scale;
      this.pw = pw;
      this.label = label;
      drawing = ImageJUtils.isShowGenericDialog();
      if (drawing) {
        update();
      }
    }

    private void update() {
      final double fwhm = getFwhm();
      label.setText(String.format("FWHM = %s px (%s nm)", MathUtils.rounded(fwhm),
          MathUtils.rounded(fwhm * scale)));

      final Plot plot = pw.getPlot();

      final double x = plot.scaleXtoPxl(centre);
      final double min = plot.scaleYtoPxl(0);
      final double max = plot.scaleYtoPxl(maxY);

      pw.getImagePlus().setRoi(new Line(x, min, x, max));

      imp.setSlice(centre);
      imp.resetDisplayRange();
      imp.updateAndDraw();
    }

    double getFwhm() {
      return 2 * cy[centre - offset];
    }

    @Override
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent event) {
      centre = (int) gd.getNextNumber();
      settings.updateCentre = gd.getNextBoolean();
      settings.updateHwhm = gd.getNextBoolean();
      update();
      return true;
    }
  }
}
