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
import ij.Menus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GUI;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.plugin.ZProjector;
import ij.plugin.filter.GaussianBlur;
import ij.plugin.frame.PlugInFrame;
import ij.process.ImageProcessor;
import ij.text.TextWindow;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.Choice;
import java.awt.Color;
import java.awt.Component;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Label;
import java.awt.Panel;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.WindowEvent;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicReference;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.DefaultListModel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import uk.ac.sussex.gdsc.core.ij.ImageAdapter;
import uk.ac.sussex.gdsc.core.ij.ImageJTrackProgress;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.OffsetPointRoi;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrencyUtils;
import uk.ac.sussex.gdsc.smlm.data.config.PsfProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.engine.FitConfiguration;
import uk.ac.sussex.gdsc.smlm.fitting.FitResult;
import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.fitting.Gaussian2DFitter;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.ij.IJImageSource;
import uk.ac.sussex.gdsc.smlm.ij.utils.ImageRoiPainter;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.Trace;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyrResultProcedure;

/**
 * Allows analysis of the signal and on/off times for fixed fluorophore spots in an image stack.
 */
public class SpotAnalysis extends PlugInFrame
    implements ActionListener, ItemListener, ListSelectionListener {
  private static final long serialVersionUID = 20190108L;

  private static final String PLUGIN_TITLE = "Spot Analysis";

  // Image titles
  private static final String RAW_MEAN_TITLE = PLUGIN_TITLE + " Raw mean";
  private static final String RAW_SD_TITLE = PLUGIN_TITLE + " Raw SD";
  private static final String RAW_SPLOT_TITLE = PLUGIN_TITLE + " Raw spot";
  private static final String BLUR_SPOT_TITLE = PLUGIN_TITLE + " Blur spot";
  private static final String AVG_SPOT_TITLE = PLUGIN_TITLE + " Average spot";
  private static final String[] RESULTS_TITLES =
      {RAW_MEAN_TITLE, RAW_SD_TITLE, RAW_SPLOT_TITLE, BLUR_SPOT_TITLE, AVG_SPOT_TITLE};

  private static final AtomicReference<Frame> INSTANCE = new AtomicReference<>();
  private static final AtomicReference<TextWindow> RESULTS_WINDOW = new AtomicReference<>();

  private TextWindow resultsWindow;

  private Choice inputChoice;
  private TextField widthTextField;
  private TextField blurTextField;
  private TextField gainTextField;
  private TextField exposureTextField;
  private TextField smoothingTextField;
  private Button profileButton;
  private Button addButton;
  private Button deleteButton;
  private Button saveButton;
  private Button saveTracesButton;
  private Button helpButton;
  private Label currentLabel;
  private Label rawFittedLabel;
  private Label blurFittedLabel;
  @SuppressWarnings("rawtypes")
  private DefaultListModel listModel;
  @SuppressWarnings("rawtypes")
  private JList onFramesList;

  private static final String OPT_LOCATION = "CT.location";

  private int runMode;
  private transient ImagePlus imp;
  /** The raw imp. */
  transient ImagePlus rawImp;
  private transient ImagePlus blurImp;

  private double gain;
  private double msPerFrame;
  private double[] xValues;
  private double[] rawMean;
  private double[] smoothMean;
  private double[] rawSd;
  private double[] smoothSd;
  private int area;
  private int currentSlice;
  private Rectangle areaBounds;

  private final Set<Spot> onFrames = new TreeSet<>();
  private final IntArrayList candidateFrames = new IntArrayList();

  private final Int2ObjectOpenHashMap<Trace> traces = new Int2ObjectOpenHashMap<>();
  private int id;
  private boolean updated;

  private final transient Object runLock = new Object();

  // Stores the list of images last used in the selection options
  private List<String> imageList = new ArrayList<>();

  /**
   * Represent a spot in a frame with a signal intensity.
   */
  private class Spot implements Comparable<Spot> {
    int frame;
    double signal;

    Spot(int frame, double signal) {
      this.frame = frame;
      this.signal = signal;
    }

    @Override
    public String toString() {
      return String.format("%d : %.2f", frame, getSignal(frame));
    }

    @Override
    public int compareTo(Spot other) {
      return Integer.compare(frame, other.frame);
    }

    @Override
    public boolean equals(Object obj) {
      if (!(obj instanceof Spot)) {
        return false;
      }
      final Spot o = (Spot) obj;
      return frame == o.frame;
    }

    @Override
    public int hashCode() {
      return Integer.hashCode(frame);
    }
  }

  /**
   * Represent a spot and a trace.
   */
  private static class TraceResult {
    Spot spot;
    Trace trace;

    TraceResult(Spot spot, Trace trace) {
      this.spot = spot;
      this.trace = trace;
    }
  }

  /**
   * Use a runnable for the image generation to allow multi-threaded operation. Input parameters
   * that are manipulated should have synchronized methods.
   */
  private class BlurWorker implements Runnable {
    Ticker ticker;
    ImageStack inputStack;
    ImageStack outputStack;
    int slice;
    int slices;
    Rectangle bounds;
    double blur;

    BlurWorker(Ticker ticker, ImageStack inputStack, int slice, int slices, Rectangle bounds,
        double blur, ImageStack outputStack) {
      this.ticker = ticker;
      this.inputStack = inputStack;
      this.slice = slice;
      this.slices = slices;
      this.bounds = bounds;
      this.blur = blur;
      this.outputStack = outputStack;
    }

    @Override
    public void run() {
      final GaussianBlur gb = new GaussianBlur();
      for (int i = 0; i < slices && slice <= inputStack.getSize(); i++, slice++) {
        ImageJUtils.showStatus(
            () -> String.format("Calculating blur ... %.1f%%", 100.0 * ticker.getProgress()));
        final ImageProcessor ip = inputStack.getProcessor(slice).duplicate();
        ip.setRoi(bounds);
        ip.snapshot();
        gb.blurGaussian(ip, blur, blur, 0.002);
        outputStack.setPixels(ip.getPixels(), slice);
        ticker.tick();
      }
    }
  }

  /**
   * Instantiates a new spot analysis.
   */
  public SpotAnalysis() {
    super(PLUGIN_TITLE);
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (WindowManager.getImageCount() == 0) {
      IJ.showMessage(PLUGIN_TITLE, "No images opened.");
      return;
    }

    final Frame frame = INSTANCE.get();

    if ("add".equals(arg)) {
      if (frame != null) {
        ((SpotAnalysis) frame).addFrame();
      }
      return;
    }

    if (frame != null) {
      frame.toFront();
      return;
    }

    INSTANCE.set(this);
    IJ.register(SpotAnalysis.class);
    WindowManager.addWindow(this);
    ImagePlus.addImageListener(new ImageAdapter() {
      @Override
      public void imageUpdated(ImagePlus imp) {
        SpotAnalysis.this.imageUpdated(imp);
      }
    });

    createFrame();
    setup();

    addKeyListener(new KeyAdapter() {
      @Override
      public void keyTyped(KeyEvent event) {
        switch (event.getKeyChar()) {
          case 'a':
            addFrame();
            break;
          case ',':
            rawImp.setSlice(rawImp.getSlice() + 1);
            break;
          case '.':
            rawImp.setSlice(rawImp.getSlice() - 1);
            break;
          default:
            break;
        }
      }
    });

    pack();
    final Point loc = Prefs.getLocation(OPT_LOCATION);
    if (loc != null) {
      setLocation(loc);
    } else {
      GUI.center(this);
    }
    if (IJ.isMacOSX()) {
      setResizable(false);
    }
    setVisible(true);

    // Install shortcut for adding a slice
    final String command = "Spot Analysis (Add)";
    final String shortcut = "6";
    final String plugin = "ij.plugin.Hotkeys(" + "\"" + command + "\")";
    Menus.installPlugin(plugin, Menus.SHORTCUTS_MENU, "*" + command, shortcut, IJ.getInstance());
  }

  private void setup() {
    if (WindowManager.getCurrentImage() == null) {
      return;
    }
    fillImagesList();
  }

  private void fillImagesList() {
    // Find the currently open images
    final ArrayList<String> newImageList = new ArrayList<>();

    for (final int imageId : ImageJUtils.getIdList()) {
      final ImagePlus image = WindowManager.getImage(imageId);

      // Image must be greyscale stacks
      if (image != null && (image.getType() == ImagePlus.GRAY8
          || image.getType() == ImagePlus.GRAY16 || image.getType() == ImagePlus.GRAY32)
          && image.getStackSize() > 2) {
        // Exclude previous results
        if (previousResult(image.getTitle())) {
          continue;
        }

        newImageList.add(image.getTitle());
      }
    }

    // Check if the image list has changed
    if (imageList.equals(newImageList)) {
      return;
    }

    imageList = newImageList;

    // Re-populate the image lists
    final String oldChoice = inputChoice.getSelectedItem();
    inputChoice.removeAll();

    for (final String imageTitle : newImageList) {
      inputChoice.add(imageTitle);
    }

    // Ensure the drop-downs are resized
    pack();

    // Restore previous selection
    inputChoice.select(oldChoice);
  }

  private static boolean previousResult(String title) {
    for (final String resultTitle : RESULTS_TITLES) {
      if (title.startsWith(resultTitle)) {
        return true;
      }
    }
    return false;
  }

  @Override
  public synchronized void actionPerformed(ActionEvent event) {
    final Object actioner = event.getSource();

    if (actioner == null) {
      return;
    }

    // Note: runLock and runMode are only used in this method

    synchronized (runLock) {
      if (runMode > 0) {
        return;
      }

      if (((Button) actioner == profileButton) && (parametersReady())) {
        runMode = 1;
      } else if ((Button) actioner == addButton) {
        runMode = 2;
      } else if ((Button) actioner == deleteButton) {
        runMode = 3;
      } else if ((Button) actioner == saveButton) {
        runMode = 4;
      } else if ((Button) actioner == saveTracesButton) {
        runMode = 5;
      } else if ((Button) actioner == helpButton) {
        ImageJUtils.showUrl(HelpUrls.getUrl("spot-analysis"));
        runMode = 0;
      }
    }

    if (runMode > 0) {
      final Thread thread = new Thread(() -> {
        synchronized (runLock) {
          try {
            switch (runMode) {
              case 1:
                createProfile();
                break;
              case 2:
                addFrame();
                break;
              case 3:
                deleteFrames();
                break;
              case 4:
                saveSpot();
                break;
              case 5:
                saveTraces();
                break;
              default:
                // Do nothing
                break;
            }
          } finally {
            runMode = 0;
          }
        }

        // Alert anything waiting for the actionPerformed method
        synchronized (SpotAnalysis.this) {
          super.notifyAll();
        }
      }, PLUGIN_TITLE);
      thread.start();
    }

    // Alert anything waiting for the actionPerformed method
    super.notifyAll();
  }

  @Override
  public void itemStateChanged(ItemEvent event) {
    // Ignore
  }

  @Override
  public void windowClosing(WindowEvent event) {
    Prefs.saveLocation(OPT_LOCATION, getLocation());
    close();
  }

  @Override
  public void close() {
    INSTANCE.set(null);
    super.close();
  }

  @Override
  public void windowActivated(WindowEvent event) {
    fillImagesList();

    super.windowActivated(event);
    WindowManager.setWindow(this);
  }

  private void createProfile() {
    if (!parametersReady()) {
      return;
    }

    double psfWidth;
    double blur;

    // Read settings
    try {
      psfWidth = Double.parseDouble(widthTextField.getText());
      blur = Double.parseDouble(blurTextField.getText());
      gain = Double.parseDouble(gainTextField.getText());
      msPerFrame = Double.parseDouble(exposureTextField.getText());
    } catch (final NumberFormatException ex) {
      IJ.error(PLUGIN_TITLE, "Invalid numbers in the input parameters");
      return;
    }

    final ImagePlus image = WindowManager.getImage(inputChoice.getSelectedItem());

    // This should not be a problem but leave it in for now
    if (image == null || (image.getType() != ImagePlus.GRAY8 && image.getType() != ImagePlus.GRAY16
        && image.getType() != ImagePlus.GRAY32)) {
      IJ.showMessage(PLUGIN_TITLE, "Images must be grayscale.");
      return;
    }

    final Roi roi = image.getRoi();
    if (roi == null || !roi.isArea()) {
      IJ.showMessage(PLUGIN_TITLE, "Image must have an area ROI");
      return;
    }

    final int recommendedSize = (int) Math.ceil(8 * psfWidth);
    final Rectangle bounds = roi.getBounds();
    if (bounds.width < recommendedSize || bounds.height < recommendedSize) {
      IJ.showMessage(PLUGIN_TITLE,
          String.format("Recommend using an ROI of at least %d x %d for the PSF width",
              recommendedSize, recommendedSize));
      return;
    }

    // Check no existing spots are within the ROI
    if (resultsWithinBounds(bounds)) {
      final GenericDialog gd = new GenericDialog(PLUGIN_TITLE);
      gd.enableYesNoCancel();
      gd.hideCancelButton();
      gd.addMessage("The results list contains a spot within the selected bounds\n \n"
          + "Do you want to continue?");
      gd.showDialog();
      if (!gd.wasOKed()) {
        return;
      }
    }

    runCreateProfile(image, bounds, psfWidth, blur);
  }

  private boolean parametersReady() {
    if (inputChoice.getItemCount() == 0) {
      IJ.showMessage(PLUGIN_TITLE,
          "No available images. Images must be 8-bit or 16-bit grayscale.");
      return false;
    }

    if (!onFrames.isEmpty() && updated) {
      final GenericDialog gd = new GenericDialog(PLUGIN_TITLE);
      gd.enableYesNoCancel();
      gd.hideCancelButton();
      gd.addMessage(
          "The list contains unsaved selected frames. Creating a new profile will erase them.\n \n"
              + "Do you want to continue?");
      gd.showDialog();
      if (!gd.wasOKed()) {
        return false;
      }
      clearSelectedFrames();
    }

    return (inputChoice.getSelectedIndex() != -1);
  }

  private boolean resultsWithinBounds(Rectangle bounds) {
    if (ImageJUtils.isShowing(resultsWindow)) {
      final float minx = bounds.x;
      final float maxx = minx + bounds.width;
      final float miny = bounds.y;
      final float maxy = miny + bounds.height;

      for (int i = 0; i < resultsWindow.getTextPanel().getLineCount(); i++) {
        final String line = resultsWindow.getTextPanel().getLine(i);
        try (Scanner s = new Scanner(line)) {
          s.useDelimiter("\t");
          s.nextInt();
          final float cx = s.nextFloat(); // cx
          final float cy = s.nextFloat(); // cy
          if (cx >= minx && cx <= maxx && cy >= miny && cy <= maxy) {
            return true;
          }
        }
      }
    }
    return false;
  }

  private void runCreateProfile(ImagePlus imp, Rectangle bounds, double psfWidth, double blur) {
    areaBounds = bounds;
    this.imp = imp;
    area = bounds.width * bounds.height;
    clearSelectedFrames();

    // Get a profile through the images
    IJ.showStatus("Calculating raw profile");

    final int nSlices = imp.getStackSize();

    final ImageStack rawSpot = new ImageStack(bounds.width, bounds.height, nSlices);
    final double[][] profile = extractSpotProfile(imp, bounds, rawSpot);

    // Retain the existing display range
    double min = 0;
    double max = Double.POSITIVE_INFINITY;
    if (rawImp != null) {
      min = rawImp.getDisplayRangeMin();
      max = rawImp.getDisplayRangeMax();
    }
    rawImp = showSpot(RAW_SPLOT_TITLE, rawSpot);
    if (max != Double.POSITIVE_INFINITY) {
      rawImp.setDisplayRange(min, max);
    }

    rawMean = profile[0];
    rawSd = profile[1];

    // Check if there are fitted results in memory
    addCandidateFrames(imp.getTitle());

    updateProfilePlots();

    if (blur > 0) {
      IJ.showStatus("Calculating blur ...");

      final ImageStack stack = imp.getImageStack();
      final ImageStack newStack =
          new ImageStack(stack.getWidth(), stack.getHeight(), stack.getSize());
      // Multi-thread the blur stage
      final ExecutorService threadPool = Executors.newFixedThreadPool(Prefs.getThreads());
      final List<Future<?>> futures = new LinkedList<>();

      final Ticker ticker = Ticker.create(new ImageJTrackProgress(true), nSlices, true);
      final int slices = 5;
      ImageJUtils.showSlowProgress(0, nSlices);
      for (int n = 1; n <= nSlices; n += slices) {
        futures.add(threadPool
            .submit(new BlurWorker(ticker, stack, n, slices, bounds, blur * psfWidth, newStack)));
      }

      IJ.showStatus("Calculating blur ... Finishing");
      threadPool.shutdown();
      ConcurrencyUtils.waitForCompletionUnchecked(futures);
      ImageJUtils.clearSlowProgress();
      IJ.showStatus("Calculating blur ... Drawing");

      final ImageStack blurSpot = new ImageStack(bounds.width, bounds.height, nSlices);
      extractSpotProfile(new ImagePlus("Blur", newStack), bounds, blurSpot);
      // Retain the existing display range
      max = Double.POSITIVE_INFINITY;
      if (blurImp != null) {
        min = blurImp.getDisplayRangeMin();
        max = blurImp.getDisplayRangeMax();
      }
      blurImp = showSpot(BLUR_SPOT_TITLE, blurSpot);
      if (max != Double.POSITIVE_INFINITY) {
        blurImp.setDisplayRange(min, max);
      }
      IJ.showStatus("");
    } else {
      blurImp = null;
    }

    // Add a z-projection of the blur/original image
    final ZProjector project = new ZProjector((blurImp == null) ? rawImp : blurImp);
    project.setMethod(ZProjector.AVG_METHOD);
    project.doProjection();
    showSpot(AVG_SPOT_TITLE, project.getProjection().getImageStack());

    if (!candidateFrames.isEmpty()) {
      // Set the first candidate frame
      rawImp.setSlice(candidateFrames.elements()[0]);
    } else {
      updateCurrentSlice(rawImp.getCurrentSlice());
    }

    IJ.showStatus("");
  }

  private void clearSelectedFrames() {
    currentSlice = -1;
    onFrames.clear();
    listModel.clear();
    candidateFrames.clear();
    updated = false;
  }

  private double[][] extractSpotProfile(ImagePlus imp, Rectangle bounds, ImageStack rawSpot) {
    final int nSlices = imp.getStackSize();
    final IJImageSource rawSource = new IJImageSource(imp);
    rawSource.open();

    final double[][] profile = new double[2][nSlices];
    for (int n = 0; n < nSlices; n++) {
      IJ.showProgress(n, nSlices);
      final float[] data = rawSource.next(bounds);
      rawSpot.setPixels(data, n + 1);
      final Statistics stats = Statistics.create(data);
      profile[0][n] = stats.getMean() / gain;
      profile[1][n] = stats.getStandardDeviation() / gain;
    }

    return profile;
  }

  private static ImagePlus showSpot(String title, ImageStack spot) {
    final WindowOrganiser windowOrganiser = new WindowOrganiser();
    final ImagePlus imp = ImageJUtils.display(title, spot, windowOrganiser);
    if (windowOrganiser.isNotEmpty() || imp.getWindow().getCanvas().getMagnification() == 1) {
      for (int i = 9; i-- > 0;) {
        imp.getWindow().getCanvas().zoomIn(imp.getWidth() / 2, imp.getHeight() / 2);
      }
    }
    return imp;
  }

  private void addCandidateFrames(String title) {
    for (final MemoryPeakResults r : MemoryPeakResults.getAllResults()) {
      if (r.getSource() instanceof IJImageSource && r.getSource().getName().equals(title)) {
        final float minx = areaBounds.x;
        final float maxx = minx + areaBounds.width;
        final float miny = areaBounds.y;
        final float maxy = miny + areaBounds.height;

        r.forEach(DistanceUnit.PIXEL, (XyrResultProcedure) (x, y, result) -> {
          if (result.getXPosition() >= minx && result.getXPosition() <= maxx
              && result.getYPosition() >= miny && result.getYPosition() <= maxy) {
            candidateFrames.add(result.getFrame());
          }
        });
      }
    }
  }

  private void updateProfilePlots() {
    xValues = SimpleArrayUtils.newArray(rawMean.length, 1.0, 1.0);
    smoothMean = interpolate(xValues, rawMean);
    smoothSd = interpolate(xValues, rawSd);

    drawProfiles();
  }

  private double[] interpolate(double[] xValues, double[] yValues) {
    // Smooth the values not in the current on-frames
    double[] newX = Arrays.copyOf(xValues, xValues.length);
    double[] newY = Arrays.copyOf(yValues, yValues.length);

    for (final Spot s : onFrames) {
      newX[s.frame - 1] = -1;
    }
    int count = 0;
    for (int i = 0; i < newX.length; i++) {
      if (newX[i] == -1) {
        continue;
      }
      newX[count] = newX[i];
      newY[count] = newY[i];
      count++;
    }
    newX = Arrays.copyOf(newX, count);
    newY = Arrays.copyOf(newY, count);
    double smoothing = 0.25;
    try {
      smoothing = Double.parseDouble(smoothingTextField.getText());
      if (smoothing < 0.01 || smoothing > 0.9) {
        smoothing = 0.25;
      }
    } catch (final NumberFormatException ignored) {
      // Ignore
    }

    final LoessInterpolator loess = new LoessInterpolator(smoothing, 1);
    final PolynomialSplineFunction f = loess.interpolate(newX, newY);

    // Interpolate
    final double[] plotSmooth = new double[xValues.length];
    for (int i = 0; i < xValues.length; i++) {
      // Cannot interpolate outside the bounds of the input data
      if (xValues[i] < newX[0]) {
        plotSmooth[i] = newY[0];
      } else if (xValues[i] > newX[newX.length - 1]) {
        plotSmooth[i] = newY[newX.length - 1];
      } else {
        plotSmooth[i] = f.value(xValues[i]);
      }
    }

    return plotSmooth;
  }

  private void drawProfiles() {
    showProfile(RAW_MEAN_TITLE, "Mean", xValues, rawMean, smoothMean);
    showProfile(RAW_SD_TITLE, "SD", xValues, rawSd, smoothSd);
  }

  private void showProfile(String title, String yTitle, double[] xValues, double[] yValues,
      double[] yValues2) {
    final Plot plot = new Plot(title, "Frame", yTitle);
    plot.addPoints(xValues, yValues, Plot.LINE);
    final double[] limits = MathUtils.limits(yValues);
    plot.setLimits(xValues[0], xValues[xValues.length - 1], limits[0], limits[1]);
    plot.draw();

    plot.setColor(Color.red);
    plot.addPoints(xValues, yValues2, Plot.LINE);

    plot.setColor(Color.magenta);

    // Add the on-frames
    if (!onFrames.isEmpty()) {
      final double[] onx = new double[onFrames.size()];
      final double[] ony = new double[onx.length];
      int count = 0;
      for (final Spot s : onFrames) {
        onx[count] = s.frame;
        ony[count] = yValues[s.frame - 1];
        count++;
      }
      plot.addPoints(onx, ony, Plot.CIRCLE);
    }

    // Add the candidate frames
    if (!candidateFrames.isEmpty()) {
      plot.setColor(Color.cyan);
      final double[] onx = new double[candidateFrames.size()];
      final double[] ony = new double[onx.length];
      int count = 0;
      final int[] e = candidateFrames.elements();
      for (int i = 0; i < candidateFrames.size(); i++) {
        final int frame = e[i];
        onx[count] = frame;
        ony[count] = yValues[frame - 1];
        count++;
      }
      plot.addPoints(onx, ony, Plot.BOX);
      plot.setColor(Color.magenta);
    }

    // Overlay current position
    plot.addPoints(new double[] {rawImp.getCurrentSlice(), rawImp.getCurrentSlice()}, limits,
        Plot.LINE);

    plot.setColor(Color.blue);
    ImageJUtils.display(title, plot, ImageJUtils.NO_TO_FRONT);
  }

  /**
   * Adds the frame.
   */
  @SuppressWarnings("unchecked")
  void addFrame() {
    if (rawImp != null) {
      final int slice = rawImp.getCurrentSlice();
      final double signal = getSignal(slice);
      final Spot s = new Spot(slice, signal);
      if (onFrames.add(s)) {
        onFramesList.clearSelection();
        // Find the location to insert in order
        int index = 0;
        while (index < listModel.size()) {
          final Spot s2 = (Spot) listModel.get(index);
          if (s.compareTo(s2) < 0) {
            break;
          }
          index++;
        }
        listModel.add(index, s);
        updateProfilePlots();
        updated = true;
      }
    }
  }

  /**
   * Gets the signal.
   *
   * @param slice the slice
   * @return the signal
   */
  double getSignal(int slice) {
    return (rawMean[slice - 1] - smoothMean[slice - 1]) * area;
  }

  private void deleteFrames() {
    final int[] indices = onFramesList.getSelectedIndices();
    if (indices.length > 0) {
      updated = true;
    }

    onFramesList.clearSelection();
    for (int i = indices.length; i-- > 0;) {
      final Spot removed = (Spot) listModel.get(indices[i]);
      onFrames.remove(removed);
    }
    if (onFrames.isEmpty()) {
      listModel.clear();
    } else {
      for (int i = indices.length; i-- > 0;) {
        listModel.remove(indices[i]);
      }
    }
  }

  private void saveSpot() {
    if (onFrames.isEmpty() || !updated) {
      return;
    }

    createResultsWindow();

    id++;
    double signal = 0;
    Trace trace = null;
    final float psfWidth = Float.parseFloat(widthTextField.getText());
    final float cx = areaBounds.x + areaBounds.width / 2.0f;
    final float cy = areaBounds.y + areaBounds.height / 2.0f;
    for (final Spot s : onFrames) {
      // Get the signal again since the background may have changed.
      final double spotSignal = getSignal(s.frame);
      signal += spotSignal;
      final float[] params = Gaussian2DPeakResultHelper.createOneAxisParams(0, (float) (spotSignal),
          cx, cy, 0, psfWidth);
      final PeakResult result =
          new PeakResult(s.frame, (int) cx, (int) cy, 0, 0, 0, 0, params, null);
      if (trace == null) {
        trace = new Trace(result);
      } else {
        trace.add(result);
      }
    }

    if (trace == null) {
      return;
    }

    final Statistics tOn = Statistics.create(trace.getOnTimes());
    final Statistics tOff = Statistics.create(trace.getOffTimes());
    resultsWindow.append(String.format("%d\t%.1f\t%.1f\t%s\t%s\t%s\t%d\t%s\t%s\t%s", id, cx, cy,
        MathUtils.rounded(signal, 4), MathUtils.rounded(tOn.getSum() * msPerFrame, 3),
        MathUtils.rounded(tOff.getSum() * msPerFrame, 3), trace.getBlinks() - 1,
        MathUtils.rounded(tOn.getMean() * msPerFrame, 3),
        MathUtils.rounded(tOff.getMean() * msPerFrame, 3), imp.getTitle()));

    // Save the individual on/off times for use in creating a histogram
    traces.put(id, trace);

    updated = false;
  }

  private void createResultsWindow() {
    resultsWindow = ImageJUtils.refresh(RESULTS_WINDOW,
        () -> new TextWindow(PLUGIN_TITLE + " Results",
            "Id\tcx\tcy\tSignal\tt-On (ms)\tt-Off (ms)\tBlinks\tAv.t-On\tAv.t-Off\tSource", "", 600,
            200));
  }

  private void saveTraces() {
    if (!onFrames.isEmpty() && updated) {
      final GenericDialog gd = new GenericDialog(PLUGIN_TITLE);
      gd.enableYesNoCancel();
      gd.hideCancelButton();
      gd.addMessage("The list contains unsaved selected frames.\n \nDo you want to continue?");
      gd.showDialog();
      if (!gd.wasOKed()) {
        return;
      }
    }

    // For all spots in the results window, get the ID and then save the traces to memory
    if (!ImageJUtils.isShowing(resultsWindow)) {
      return;
    }

    // Create a results set in memory
    final MemoryPeakResults results = new MemoryPeakResults();
    results.setName(PLUGIN_TITLE);
    results.begin();
    MemoryPeakResults.addResults(results);

    final ArrayList<TraceResult> traceResults =
        new ArrayList<>(resultsWindow.getTextPanel().getLineCount());
    for (int i = 0; i < resultsWindow.getTextPanel().getLineCount(); i++) {
      final String line = resultsWindow.getTextPanel().getLine(i);
      try (Scanner s = new Scanner(line)) {
        s.useDelimiter("\t");
        int id = -1;
        double signal = -1;
        // Be careful as the text panel may not contain what we expect, i.e. empty lines, etc
        if (s.hasNextInt()) {
          id = s.nextInt();
          try {
            s.nextDouble(); // cx
            s.nextDouble(); // cy
            signal = s.nextDouble();
          } catch (final NoSuchElementException ignored) {
            // Ignore
          }
        }

        if (id != -1 && signal != -1) {
          final Trace trace = traces.get(id);
          if (trace != null) {
            results.addAll(trace.getPoints());
            traceResults.add(new TraceResult(new Spot(id, signal), trace));
          }
        }
      }
    }

    results.end();

    saveTracesToFile(traceResults);

    IJ.showStatus("Saved traces");
  }

  private void saveTracesToFile(List<TraceResult> traceResults) {
    final String resultsDirectory = IJ.getDirectory(PLUGIN_TITLE);
    if (resultsDirectory == null) {
      return;
    }

    // Save the traces to a single file.
    // Also save the blinks and on/off times into data files for histogram analysis

    try (
        BufferedWriter files0 = openBufferedWriter(resultsDirectory + "traces.txt",
            String.format("#ms/frame = %s%n#Id\tcx\tcy\tsignal\tn-Blinks\tStart\tStop\t...",
                MathUtils.rounded(msPerFrame, 3)));
        BufferedWriter files1 = openBufferedWriter(resultsDirectory + "tOn.txt", "");
        BufferedWriter files2 = openBufferedWriter(resultsDirectory + "tOff.txt", "");
        BufferedWriter files3 = openBufferedWriter(resultsDirectory + "blinks.txt", "");
        BufferedWriter files4 = openBufferedWriter(resultsDirectory + "signal.txt", "")) {
      for (final TraceResult traceResult : traceResults) {
        final StringBuilder sb = new StringBuilder();
        sb.append(traceResult.spot.frame).append('\t')
            .append(traceResult.trace.getHead().getXPosition()).append('\t')
            .append(traceResult.trace.getHead().getYPosition()).append('\t')
            .append(traceResult.spot.signal).append('\t');
        final int nBlinks = traceResult.trace.getBlinks() - 1;
        sb.append(nBlinks);

        final int[] on = traceResult.trace.getOnTimes();
        final int[] off = traceResult.trace.getOffTimes();
        int time = traceResult.trace.getHead().getFrame();
        for (int i = 0; i < on.length; i++) {
          writeLine(files1, Double.toString(msPerFrame * on[i]));

          sb.append('\t').append(time).append('\t').append(time + on[i] - 1);
          if (off != null && i < off.length) {
            writeLine(files2, Double.toString(msPerFrame * off[i]));
            time += on[i] + off[i];
          }
        }
        writeLine(files0, sb.toString());
        writeLine(files3, Integer.toString(nBlinks));
        writeLine(files4, String.format("# Id=%d, Blinks=%d, Signal=%f", traceResult.spot.frame,
            nBlinks, traceResult.spot.signal));
        for (int k = 0; k < traceResult.trace.size(); k++) {
          final PeakResult r = traceResult.trace.get(k);
          writeLine(files4, String.format("%d %f", r.getFrame(), r.getIntensity()));
        }
      }
    } catch (final IOException ex) {
      // Q. Add better handling of errors?
      Logger.getLogger(getClass().getName()).log(Level.WARNING, ex,
          () -> "Failed to save traces to results directory: " + resultsDirectory);
    }
  }

  private static BufferedWriter openBufferedWriter(String filename, String header)
      throws IOException {
    final BufferedWriter tracesFile = Files.newBufferedWriter(Paths.get(filename));
    if (header != null && header.length() > 0) {
      writeLine(tracesFile, header);
    }
    return tracesFile;
  }

  private static void writeLine(BufferedWriter tracesFile, String line) throws IOException {
    tracesFile.write(line);
    tracesFile.newLine();
  }

  @SuppressWarnings({"rawtypes", "unchecked"})
  private void createFrame() {
    final Panel mainPanel = new Panel();
    add(mainPanel);

    inputChoice = new Choice();
    mainPanel.add(createChoicePanel(inputChoice, ""));

    widthTextField = new TextField();
    mainPanel.add(createTextPanel(widthTextField, "PSF width", "1.2"));

    blurTextField = new TextField();
    mainPanel.add(createTextPanel(blurTextField, "Blur (relative to width)", "1"));

    gainTextField = new TextField();
    mainPanel.add(createTextPanel(gainTextField, "Gain", "37.7"));

    exposureTextField = new TextField();
    mainPanel.add(createTextPanel(exposureTextField, "ms/Frame", "20"));

    smoothingTextField = new TextField();
    mainPanel.add(createTextPanel(smoothingTextField, "Smoothing", "0.25"));

    profileButton = new Button("Profile");
    profileButton.addActionListener(this);
    addButton = new Button("Add");
    addButton.addActionListener(this);
    deleteButton = new Button("Remove");
    deleteButton.addActionListener(this);
    saveButton = new Button("Save");
    saveButton.addActionListener(this);
    saveTracesButton = new Button("Save Traces");
    saveTracesButton.addActionListener(this);
    helpButton = new Button("Help");
    helpButton.addActionListener(this);

    currentLabel = new Label();
    mainPanel.add(createLabelPanel(currentLabel, "", ""));

    rawFittedLabel = new Label();
    mainPanel.add(createLabelPanel(rawFittedLabel, "", ""));

    blurFittedLabel = new Label();
    mainPanel.add(createLabelPanel(blurFittedLabel, "", ""));

    final JPanel buttonPanel = new JPanel();
    final FlowLayout l = new FlowLayout();
    l.setVgap(0);
    buttonPanel.setLayout(l);
    buttonPanel.add(profileButton, BorderLayout.CENTER);
    buttonPanel.add(addButton, BorderLayout.CENTER);
    buttonPanel.add(deleteButton, BorderLayout.CENTER);
    buttonPanel.add(saveButton, BorderLayout.CENTER);
    buttonPanel.add(saveTracesButton, BorderLayout.CENTER);
    buttonPanel.add(helpButton, BorderLayout.CENTER);

    mainPanel.add(buttonPanel);

    listModel = new DefaultListModel();
    onFramesList = new JList(listModel);
    onFramesList.setVisibleRowCount(15);
    onFramesList.addListSelectionListener(this);

    final JScrollPane scrollPane = new JScrollPane(onFramesList);
    mainPanel.add(scrollPane);

    final GridBagLayout mainGrid = new GridBagLayout();
    int y = 0;
    final GridBagConstraints c = new GridBagConstraints();
    c.gridx = 0;
    c.fill = GridBagConstraints.BOTH;
    c.anchor = GridBagConstraints.WEST;
    c.gridwidth = 1;
    c.insets = new Insets(2, 2, 2, 2);

    for (final Component comp : mainPanel.getComponents()) {
      c.gridy = y++;
      mainGrid.setConstraints(comp, c);
    }

    mainPanel.setLayout(mainGrid);
  }

  private static Panel createChoicePanel(Choice list, String label) {
    final Panel panel = new Panel();
    panel.setLayout(new BorderLayout());
    final Label listLabel = new Label(label, 0);
    panel.add(listLabel, BorderLayout.WEST);
    panel.add(list, BorderLayout.CENTER);
    return panel;
  }

  private static Panel createTextPanel(TextField textField, String label, String value) {
    final Panel panel = new Panel();
    panel.setLayout(new BorderLayout());
    final Label listLabel = new Label(label, 0);
    textField.setText(value);
    panel.add(listLabel, BorderLayout.WEST);
    panel.add(textField, BorderLayout.CENTER);
    return panel;
  }

  private static Panel createLabelPanel(Label field, String label, String value) {
    final Panel panel = new Panel();
    panel.setLayout(new BorderLayout());
    final Label listLabel = new Label(label, 0);
    field.setText(value);
    panel.add(listLabel, BorderLayout.WEST);
    panel.add(field, BorderLayout.CENTER);
    return panel;
  }

  /**
   * Image updated.
   *
   * @param imp the imp
   */
  void imageUpdated(ImagePlus imp) {
    ImagePlus from = null;
    ImagePlus to = null;
    if (imp == rawImp) {
      from = rawImp;
      to = blurImp;
    } else if (imp == blurImp) {
      from = blurImp;
      to = rawImp;
    }

    if (from != null) {
      final int slice = from.getCurrentSlice();

      updateCurrentSlice(slice);

      if (to != null && to.getCurrentSlice() != slice) {
        to.setSlice(slice);
      }
    }
  }

  private void updateCurrentSlice(int slice) {
    if (slice != currentSlice) {
      currentSlice = slice;
      final double signal = getSignal(slice);
      final double noise = smoothSd[slice - 1];
      currentLabel.setText(String.format("Frame %d: Signal = %s, SNR = %s", slice,
          MathUtils.rounded(signal, 4), MathUtils.rounded(signal / noise, 3)));

      drawProfiles();

      // Fit the PSF using a Gaussian
      final FitConfiguration fitConfiguration = FitConfiguration.create();
      fitConfiguration.setPsf(PsfProtosHelper.DefaultOneAxisGaussian2dPsf.INSTANCE);
      fitConfiguration.setFixedPsf(true);
      fitConfiguration.setBackgroundFitting(true);
      fitConfiguration.setSignalStrength(0);
      fitConfiguration.setCoordinateShift(rawImp.getWidth() / 4.0f);
      fitConfiguration.setComputeResiduals(false);
      fitConfiguration.setComputeDeviations(false);
      final Gaussian2DFitter gf = new Gaussian2DFitter(fitConfiguration);
      double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
      final double psfWidth = Double.parseDouble(widthTextField.getText());
      params[Gaussian2DFunction.BACKGROUND] = smoothMean[slice - 1];
      params[Gaussian2DFunction.SIGNAL] = (gain * signal);
      params[Gaussian2DFunction.X_POSITION] = rawImp.getWidth() / 2.0f;
      params[Gaussian2DFunction.Y_POSITION] = rawImp.getHeight() / 2.0f;
      params[Gaussian2DFunction.X_SD] = params[Gaussian2DFunction.Y_SD] = psfWidth;
      float[] data = (float[]) rawImp.getImageStack().getProcessor(slice).getPixels();
      FitResult fitResult = gf.fit(SimpleArrayUtils.toDouble(data), rawImp.getWidth(),
          rawImp.getHeight(), 1, params, new boolean[1]);
      if (fitResult.getStatus() == FitStatus.OK) {
        params = fitResult.getParameters();
        final double spotSignal = params[Gaussian2DFunction.SIGNAL] / gain;
        rawFittedLabel.setText(String.format("Raw fit: Signal = %s, SNR = %s",
            MathUtils.rounded(spotSignal, 4), MathUtils.rounded(spotSignal / noise, 3)));
        ImageRoiPainter.addRoi(rawImp, slice, new OffsetPointRoi(
            params[Gaussian2DFunction.X_POSITION], params[Gaussian2DFunction.Y_POSITION]));
      } else {
        rawFittedLabel.setText("");
        rawImp.setOverlay(null);
      }

      // Fit the PSF using a Gaussian
      if (blurImp == null) {
        return;
      }

      params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
      params[Gaussian2DFunction.BACKGROUND] = (float) smoothMean[slice - 1];
      params[Gaussian2DFunction.SIGNAL] = (float) (gain * signal);
      params[Gaussian2DFunction.X_POSITION] = rawImp.getWidth() / 2.0f;
      params[Gaussian2DFunction.Y_POSITION] = rawImp.getHeight() / 2.0f;
      params[Gaussian2DFunction.X_SD] = params[Gaussian2DFunction.Y_SD] = psfWidth;
      data = (float[]) blurImp.getImageStack().getProcessor(slice).getPixels();
      fitResult = gf.fit(SimpleArrayUtils.toDouble(data), rawImp.getWidth(), rawImp.getHeight(), 1,
          params, new boolean[1]);
      if (fitResult.getStatus() == FitStatus.OK) {
        params = fitResult.getParameters();
        final double spotSignal = params[Gaussian2DFunction.SIGNAL] / gain;
        blurFittedLabel.setText(String.format("Blur fit: Signal = %s, SNR = %s",
            MathUtils.rounded(spotSignal, 4), MathUtils.rounded(spotSignal / noise, 3)));
        ImageRoiPainter.addRoi(blurImp, slice, new OffsetPointRoi(
            params[Gaussian2DFunction.X_POSITION], params[Gaussian2DFunction.Y_POSITION]));
      } else {
        blurFittedLabel.setText("");
        blurImp.setOverlay(null);
      }
    }
  }

  @Override
  public void valueChanged(ListSelectionEvent event) {
    if (!event.getValueIsAdjusting()) {
      final int index = onFramesList.getSelectedIndex();
      if (index >= 0 && index < listModel.size()) {
        final Spot spot = (Spot) listModel.get(index);
        rawImp.setSlice(spot.frame);
      }
    }
  }
}
