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

package uk.ac.sussex.gdsc.smlm.ij.results;

import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.data.utils.Converter;
import uk.ac.sussex.gdsc.core.data.utils.IdentityTypeConverter;
import uk.ac.sussex.gdsc.core.data.utils.Rounder;
import uk.ac.sussex.gdsc.core.data.utils.RounderUtils;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.ConfigurationException;
import uk.ac.sussex.gdsc.smlm.data.config.UnitConverterUtils;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.ij.utils.ImageRoiPainter;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResultConversionHelper;
import uk.ac.sussex.gdsc.smlm.utils.CoordinateProvider;

import gnu.trove.list.array.TIntArrayList;

import ij.WindowManager;
import ij.text.TextPanel;
import ij.text.TextWindow;

import java.awt.Frame;
import java.awt.event.MouseListener;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Saves the fit results to an ImageJ results table.
 *
 * <p>The table supports mouse click events to draw the selected coordinates on the original source
 * image using the ImageROIPainter.
 */
public class ImageJTablePeakResults extends ImageJAbstractPeakResults
    implements CoordinateProvider {
  /**
   * Converter to change the distances to pixels. It is created in {@link #begin()} and may be an
   * identity converter.
   */
  protected TypeConverter<DistanceUnit> toPixelConverter;

  /** The calculator used to compute precision. */
  private Gaussian2DPeakResultCalculator calculator;

  private boolean canComputePrecision;

  private PeakResultConversionHelper helper;
  private Converter[] converters;
  private Converter ic;
  private int[] outIndices;

  private DistanceUnit distanceUnit;
  private IntensityUnit intensityUnit;
  private AngleUnit angleUnit;
  private boolean showPrecision;
  private boolean computePrecision;
  private int roundingPrecision;

  // Store the ROI painters that have been attached to TextPanels so they can be updated
  // with a new image source
  private static final Map<TextPanel, ImageRoiPainter> map = new ConcurrentHashMap<>();

  private boolean showDeviations;
  private boolean showEndFrame;
  private boolean showId;
  private boolean showFittingData;
  private boolean showNoiseData;
  private boolean showZ;
  private boolean clearAtStart;
  private boolean hideSourceText;
  private String frameColumnName = "T";
  private String source;
  private String sourceText;
  private String tableTitle = "Fit Results";
  private boolean newWindow;
  private TextWindow resultsWindow;
  private TextPanel tp;
  private ImageRoiPainter roiPainter;
  private boolean addCounter;

  /** Set to true if the table is active. */
  protected boolean tableActive;
  private int nextRepaintSize;
  private double repaintInterval = 0.1;
  private Rounder rounder;

  private int indexT = -1;
  private int indexX = -1;
  private int indexY = -1;

  private int size;

  /**
   * Instantiates a new IJ table peak results.
   *
   * @param showDeviations Set to true to show deviations
   */
  public ImageJTablePeakResults(boolean showDeviations) {
    this.showDeviations = showDeviations;
  }

  /**
   * Instantiates a new IJ table peak results.
   *
   * @param showDeviations Set to true to show deviations
   * @param source the source
   */
  public ImageJTablePeakResults(boolean showDeviations, String source) {
    this.showDeviations = showDeviations;
    this.source = source;
  }

  /**
   * Instantiates a new IJ table peak results.
   *
   * @param showDeviations Set to true to show deviations
   * @param source the source
   * @param clearAtStart Set to true to clear table contents in {@link #begin()}
   */
  public ImageJTablePeakResults(boolean showDeviations, String source, boolean clearAtStart) {
    this.showDeviations = showDeviations;
    this.source = source;
    this.clearAtStart = clearAtStart;
  }

  @Override
  public void begin() {
    tableActive = false;

    // Set-up unit processing that requires the calibration
    toPixelConverter = new IdentityTypeConverter<>(null);
    calculator = null;
    canComputePrecision = false;
    rounder = RounderUtils.create(roundingPrecision);

    // We must correctly convert all the PSF parameter types
    helper = new PeakResultConversionHelper(getCalibration(), getPsf());
    helper.setIntensityUnit(intensityUnit);
    helper.setDistanceUnit(distanceUnit);
    helper.setAngleUnit(angleUnit);
    converters = helper.getConverters();

    if (hasCalibration()) {
      if (showPrecision) {
        if (computePrecision) {
          try {
            calculator = Gaussian2DPeakResultHelper.create(getPsf(), getCalibrationReader(),
                Gaussian2DPeakResultHelper.LSE_PRECISION);
            canComputePrecision = true;
          } catch (final ConfigurationException ex) {
            // Cannot compute precision
          } catch (final ConversionException ex) {
            // Cannot compute precision
          }
        }
      }

      try {
        if (helper.hasDistanceConverter()) {
          toPixelConverter = UnitConverterUtils.createConverter(distanceUnit, DistanceUnit.PIXEL,
              getCalibrationReader().getNmPerPixel());
        }
      } catch (final ConversionException ex) {
        // Gracefully fail so ignore this
      }
    }

    ic = converters[PeakResult.INTENSITY];
    outIndices = SimpleArrayUtils.natural(converters.length);
    if (!showZ) {
      final TIntArrayList list = new TIntArrayList(outIndices);
      list.remove(PeakResult.Z);
      outIndices = list.toArray();
    }
    // Update the calibration if converters were created
    if (helper.isCalibrationChanged()) {
      setCalibration(helper.getCalibration());
    }

    createSourceText();
    createResultsWindow();
    if (clearAtStart) {
      tp.clear();
    }
    size = 0;
    // Let some results appear before drawing.
    // ImageJ will auto-layout columns if it has less than 10 rows
    nextRepaintSize = 9;
    tableActive = true;
  }

  /**
   * Clear the table contents.
   */
  public void clear() {
    tp.clear();
    size = 0;
    // Let some results appear before drawing.
    // ImageJ will auto-layout columns if it has less than 10 rows
    nextRepaintSize = 9;
  }

  /**
   * Create the result window (if it is not available).
   */
  private void createResultsWindow() {
    final String header = createResultsHeader();

    roiPainter = null;
    for (final Frame f : WindowManager.getNonImageWindows()) {
      if (f != null && tableTitle.equals(f.getTitle()) && f instanceof TextWindow) {
        resultsWindow = (TextWindow) f;

        // Check if the existing table matches the desired header
        final String currentHeader = resultsWindow.getTextPanel().getColumnHeadings();
        if (!currentHeader.startsWith(header)) {
          resultsWindow = null;
          continue;
        }

        roiPainter = map.get(resultsWindow.getTextPanel());
        break;
      }
    }

    newWindow = false;
    if (!ImageJUtils.isShowing(resultsWindow)) {
      newWindow = true;
      resultsWindow = new TextWindow(tableTitle, header, "", 800, 300);
      roiPainter = new ImageRoiPainter(resultsWindow.getTextPanel(), "", this);

      // The ROI painter adds itself to the TextPanel as a mouse listener. However
      // the TextPanel addMouseListener() adds to the private TextCanvas object so it
      // cannot be retrieved. Store the painter in a global lookup table.
      map.put(resultsWindow.getTextPanel(), roiPainter);
    }

    tp = resultsWindow.getTextPanel();

    if (roiPainter != null && getSource() != null) {
      roiPainter.setTitle(getSource().getOriginal().getName());

      // Update the coordinate provider (avoids memory leaks with old objects lying around)
      roiPainter.setCoordProvider(this);

      // Get the headings for extracting the coordinates
      final String[] headings = tp.getColumnHeadings().split("\t");
      indexT = indexX = indexY = -1;
      for (int i = 0; i < headings.length; i++) {
        if (headings[i].equals(frameColumnName)) {
          indexT = i;
          continue;
        }
        // Allow for units
        if (headings[i].equals("X") || headings[i].startsWith("X (")) {
          indexX = i;
          continue;
        }
        if (headings[i].equals("Y") || headings[i].startsWith("Y (")) {
          indexY = i;
          continue;
        }
      }
    }
  }

  /**
   * Checks if is new window.
   *
   * @return true, if is new window
   */
  public boolean isNewWindow() {
    return newWindow;
  }

  private String createResultsHeader() {
    final String[] names = helper.getNames();
    final String[] unitNames = helper.getUnitNames();

    final StringBuilder sb = new StringBuilder();
    if (addCounter) {
      sb.append("#\t");
    }
    if (sourceText != null) {
      sb.append("Source\t");
    }
    sb.append(frameColumnName);
    if (showEndFrame) {
      sb.append("\tEnd ").append(frameColumnName);
    }
    if (showId) {
      sb.append("\tId");
    }
    if (showFittingData) {
      sb.append("\torigX");
      sb.append("\torigY");
      sb.append("\torigValue");
      sb.append("\tError");
    }
    if (showNoiseData) {
      sb.append("\tNoise");
      if (!TextUtils.isNullOrEmpty(unitNames[PeakResult.INTENSITY])) {
        sb.append(" (").append(unitNames[PeakResult.INTENSITY]).append(')');
      }
      sb.append("\tMean");
      if (!TextUtils.isNullOrEmpty(unitNames[PeakResult.INTENSITY])) {
        sb.append(" (").append(unitNames[PeakResult.INTENSITY]).append(')');
      }
      sb.append("\tSNR");
    }

    for (int i = 0; i < outIndices.length; i++) {
      sb.append('\t').append(names[outIndices[i]]);
      if (!TextUtils.isNullOrEmpty(unitNames[outIndices[i]])) {
        sb.append(" (").append(unitNames[outIndices[i]]).append(')');
      }
      addDeviation(sb);
    }
    if (showPrecision) {
      sb.append("\tPrecision (nm)");
    }
    return sb.toString();
  }

  private void createSourceText() {
    if (hideSourceText) {
      sourceText = null;
      return;
    }
    final StringBuilder sb = new StringBuilder();
    if (source != null) {
      sb.append(source);
    } else if (getSource() != null) {
      sb.append(getSource().getName());
    }
    if (getBounds() != null) {
      if (sb.length() > 0) {
        sb.append(": ");
      }
      sb.append(getBoundsString());
    }
    if (sb.length() > 0) {
      sb.append('\t');
      sourceText = sb.toString();
    }
  }

  private void addDeviation(StringBuilder sb) {
    if (showDeviations) {
      sb.append("\t+/-");
    }
  }

  @Override
  public void add(int frame, int origX, int origY, float origValue, double error, float noise,
      float meanIntensity, float[] params, float[] paramsDev) {
    addPeak(frame, frame, 0, origX, origY, origValue, error, noise, meanIntensity, params,
        paramsDev, -1);
  }

  @Override
  public void add(PeakResult result) {
    addPeak(result.getFrame(), result.getEndFrame(), result.getId(), result.getOrigX(),
        result.getOrigY(), result.getOrigValue(), result.getError(), result.getNoise(),
        result.getMeanIntensity(), result.getParameters(), result.getParameterDeviations(),
        result.getPrecision());
  }

  private void addPeak(int frame, int endFrame, int id, int origX, int origY, float origValue,
      double error, float noise, float meanIntensity, float[] params, float[] paramsStdDev,
      double precision) {
    if (!tableActive) {
      return;
    }

    final StringBuilder sb =
        addStandardData(frame, endFrame, id, origX, origY, origValue, error, noise, meanIntensity);
    if (isShowDeviations()) {
      if (paramsStdDev != null) {
        for (int i = 0; i < outIndices.length; i++) {
          addFloat(sb, converters[outIndices[i]].convert(params[outIndices[i]]));
          addFloat(sb, converters[outIndices[i]].convert(paramsStdDev[outIndices[i]]));
        }
      } else {
        for (int i = 0; i < outIndices.length; i++) {
          addFloat(sb, converters[outIndices[i]].convert(params[outIndices[i]]));
          sb.append("\t0");
        }
      }
    } else {
      for (int i = 0; i < outIndices.length; i++) {
        addFloat(sb, converters[outIndices[i]].convert(params[outIndices[i]]));
      }
    }
    if (isShowPrecision()) {
      // The default precision in a peak result is NaN so this compare will be false
      if (precision >= 0) {
        addPrecision(sb, precision, false);
      } else if (canComputePrecision) {
        addPrecision(sb, calculator.getLsePrecision(params, noise), true);
      } else {
        sb.append("\t0");
      }
    }
    append(sb.toString());
  }

  private StringBuilder addStandardData(int frame, int endFrame, int id, int origX, int origY,
      float origValue, double error, float noise, float meanIntensity) {
    final StringBuilder sb = new StringBuilder();
    if (addCounter) {
      sb.append(size + 1).append('\t');
    }
    if (sourceText != null) {
      sb.append(sourceText);
    }
    // Do not calibrate the original values
    // if (showCalibratedValues)
    // sb.append(frame).append(String.format("\t%g", origX)).append(String.format("\t%g", origY));
    // else
    sb.append(frame);
    if (showEndFrame) {
      sb.append('\t').append(endFrame);
    }
    if (showId) {
      sb.append('\t').append(id);
    }
    if (showFittingData) {
      sb.append('\t').append(origX).append('\t').append(origY);
      addFloat(sb, origValue);
      addDouble(sb, error);
    }
    if (showNoiseData) {
      // These should be converted
      addFloat(sb, ic.convert(noise));
      addFloat(sb, ic.convert(meanIntensity));
      addFloat(sb, meanIntensity / noise);
    }
    return sb;
  }

  private void addFloat(StringBuilder sb, float value) {
    sb.append('\t').append(rounder.toString(value));
  }

  private void addDouble(StringBuilder sb, double value) {
    sb.append('\t').append(rounder.toString(value));
  }

  private void addPrecision(StringBuilder sb, double value, boolean computed) {
    addDouble(sb, value);
    if (computed) {
      sb.append('*');
    }
  }

  private void append(String result) {
    // Support for periodic refresh
    synchronized (tp) {
      addResult(result);
    }
    updateTable();
  }

  private void addResult(String result) {
    size++;
    tp.appendWithoutUpdate(result);
  }

  private void updateTable() {
    if (size < nextRepaintSize) {
      return;
    }

    if (!resultsWindow.isShowing()) {
      tableActive = false;
      return;
    }

    drawTable();
  }

  private void drawTable() {
    synchronized (tp) {
      nextRepaintSize = (int) (size + size * repaintInterval);
      tp.updateDisplay();
    }
  }

  @Override
  public void addAll(PeakResult[] results) {
    if (!tableActive) {
      return;
    }
    int counter = 0;
    for (final PeakResult result : results) {
      addPeak(result.getFrame(), result.getEndFrame(), result.getId(), result.getOrigX(),
          result.getOrigY(), result.getOrigValue(), result.getError(), result.getNoise(),
          result.getMeanIntensity(), result.getParameters(), result.getParameterDeviations(),
          result.getPrecision());
      if (counter++ > 31) {
        if (!tableActive) {
          return;
        }
        counter = 0;
      }
    }
  }

  @Override
  public int size() {
    return size;
  }

  @Override
  public void end() {
    tableActive = false;
    drawTable();
  }

  /**
   * Forces the table to be updated with the current contents.
   */
  public void flush() {
    drawTable();
  }

  /**
   * Gets the frame column name.
   *
   * @return the name of the frame column.
   */
  public String getFrameColumnName() {
    return frameColumnName;
  }

  /**
   * Sets the frame column name.
   *
   * @param frameColumnName the name of the frame column
   */
  public void setFrameColumnName(String frameColumnName) {
    this.frameColumnName = frameColumnName;
  }

  @Override
  public boolean isActive() {
    return tableActive;
  }

  /**
   * Gets the table title.
   *
   * @return the table title.
   */
  public String getTableTitle() {
    return tableTitle;
  }

  /**
   * Use to set the title of the table. If an existing table exists with the same title then it will
   * be appended, otherwise a new table is created.
   *
   * @param tableTitle the table title
   */
  public void setTableTitle(String tableTitle) {
    if (tableTitle != null && tableTitle.length() > 0) {
      this.tableTitle = tableTitle;
    }
  }

  /**
   * Check if the deviations of the parameters should be shown.
   *
   * @return True if the deviations of the parameters should be shown.
   */
  public boolean isShowDeviations() {
    return showDeviations;
  }

  /**
   * Set if the deviations of the parameters should be shown.
   *
   * @param showDeviations True if the deviations of the parameters should be shown
   */
  public void setShowDeviations(boolean showDeviations) {
    this.showDeviations = showDeviations;
  }

  /**
   * Check if the table should be cleared in {@link #begin()}.
   *
   * @return True if the table should be cleared in {@link #begin()}
   */
  public boolean isClearAtStart() {
    return clearAtStart;
  }

  /**
   * Set if the table should be cleared in {@link #begin()}.
   *
   * @param clearAtStart True if the table should be cleared in {@link #begin()}
   */
  public void setClearAtStart(boolean clearAtStart) {
    this.clearAtStart = clearAtStart;
  }

  /**
   * Checks if a counter will be displayed in the table.
   *
   * @return true, if a counter will be displayed in the table.
   */
  public boolean isAddCounter() {
    return addCounter;
  }

  /**
   * Set if a counter will be displayed in the table.
   *
   * @param addCounter True if a counter will be displayed in the table.
   */
  public void setAddCounter(boolean addCounter) {
    this.addCounter = addCounter;
  }

  /**
   * Checks if the source text will be added to each entry.
   *
   * @return true, if hiding the source text
   */
  public boolean isHideSourceText() {
    return hideSourceText;
  }

  /**
   * Sets the hide source text flag.
   *
   * @param hideSourceText the new hide source text flag
   */
  public void setHideSourceText(boolean hideSourceText) {
    this.hideSourceText = hideSourceText;
  }

  /**
   * Gets the results window.
   *
   * @return the resultsWindow.
   */
  public TextWindow getResultsWindow() {
    return resultsWindow;
  }

  /**
   * Gets the coordinates.
   *
   * @param line the line
   * @return the coordinates
   */
  @Override
  public double[] getCoordinates(String line) {
    // Extract the startT and x,y coordinates from the PeakResult line
    final String[] fields = line.split("\t");
    try {
      final int startT = Integer.valueOf(fields[indexT]);
      final double x = Double.valueOf(fields[indexX]);
      final double y = Double.valueOf(fields[indexY]);
      return new double[] {startT, toPixelConverter.convert(x), toPixelConverter.convert(y)};
    } catch (final ArrayIndexOutOfBoundsException ex) {
      // Will happen if any index is still at the default of -1 or if there are not enough fields
    } catch (final NumberFormatException ex) {
      // In case any field is not a number
    }
    return null;
  }

  /**
   * Check if showing the results end frame in the table.
   *
   * @return True if showing the results end frame in the table.
   */
  public boolean isShowEndFrame() {
    return showEndFrame;
  }

  /**
   * Set whether to show the results end frame in the table.
   *
   * @param showEndFrame True to show the results end frame in the table
   */
  public void setShowEndFrame(boolean showEndFrame) {
    this.showEndFrame = showEndFrame;
  }

  /**
   * Check if showing the results Id in the table.
   *
   * @return True if showing the results Id in the table.
   */
  public boolean isShowId() {
    return showId;
  }

  /**
   * Set whether to show the results Id in the table.
   *
   * @param showId True to show the results Id in the table
   */
  public void setShowId(boolean showId) {
    this.showId = showId;
  }

  /**
   * Check whether to show the fitting data (original pixel data and fit error) in the table.
   *
   * @return True to show the fitting data (original pixel data and fit error) in the table.
   */
  public boolean isShowFittingData() {
    return showFittingData;
  }

  /**
   * Set whether to show the fitting data (original pixel data and fit error) in the table.
   *
   * @param showFittingData If true then show the fitting data (original pixel data and fit error)
   *        in the table
   */
  public void setShowFittingData(boolean showFittingData) {
    this.showFittingData = showFittingData;
  }

  /**
   * Check if showing the noise and SNR in the table.
   *
   * @return True to show the noise and SNR in the table.
   */
  public boolean isShowNoiseData() {
    return showNoiseData;
  }

  /**
   * Set whether to show the noise and SNR in the table.
   *
   * @param showNoiseData If true then show the noise and SNR in the table
   */
  public void setShowNoiseData(boolean showNoiseData) {
    this.showNoiseData = showNoiseData;
  }

  /**
   * Checks if showing the Z column.
   *
   * @return true, if is show Z
   */
  public boolean isShowZ() {
    return showZ;
  }

  /**
   * Set to true to show the Z column.
   *
   * @param showZ the new show Z
   */
  public void setShowZ(boolean showZ) {
    this.showZ = showZ;
  }

  /**
   * Image will be repainted when a fraction of new results have been added.
   *
   * @param repaintInterval the repaint interval to set (range 0.001-1)
   */
  public void setRepaintInterval(double repaintInterval) {
    this.repaintInterval = MathUtils.clip(0.001, 1.0, repaintInterval);
  }

  /**
   * Gets the repaint interval.
   *
   * @return the repaint interval
   */
  public double getRepaintInterval() {
    return repaintInterval;
  }

  /**
   * Gets the distance unit.
   *
   * @return the distance unit
   */
  public DistanceUnit getDistanceUnit() {
    return distanceUnit;
  }

  /**
   * Sets the distance unit.
   *
   * @param distanceUnit the new distance unit
   */
  public void setDistanceUnit(DistanceUnit distanceUnit) {
    this.distanceUnit = distanceUnit;
  }

  /**
   * Gets the intensity unit.
   *
   * @return the intensity unit
   */
  public IntensityUnit getIntensityUnit() {
    return intensityUnit;
  }

  /**
   * Sets the intensity unit.
   *
   * @param intensityUnit the new intensity unit
   */
  public void setIntensityUnit(IntensityUnit intensityUnit) {
    this.intensityUnit = intensityUnit;
  }

  /**
   * Gets the angle unit.
   *
   * @return the angle unit
   */
  public AngleUnit getAngleUnit() {
    return angleUnit;
  }

  /**
   * Sets the angle unit.
   *
   * @param angleUnit the new angle unit
   */
  public void setAngleUnit(AngleUnit angleUnit) {
    this.angleUnit = angleUnit;
  }

  /**
   * Checks if the precision will be computed.
   *
   * @return true, if the precision will be computed
   */
  public boolean isShowPrecision() {
    return showPrecision;
  }

  /**
   * Sets the show precision flag.
   *
   * @param showPrecision set to true to show the precision and write to the output
   */
  public void setShowPrecision(boolean showPrecision) {
    this.showPrecision = showPrecision;
  }

  /**
   * Checks if the precision will be computed if needed. This is only relevant if show precision is
   * true (see {@link #isShowPrecision()}).
   *
   * @return true, if the precision will be computed
   */
  public boolean isComputePrecision() {
    return computePrecision;
  }

  /**
   * Sets the compute precision flag. This is only relevant if show precision is true (see
   * {@link #isShowPrecision()}).
   *
   * @param computePrecision set to true to compute the precision if needed
   */
  public void setComputePrecision(boolean computePrecision) {
    this.computePrecision = computePrecision;
  }

  /**
   * Gets the rounding precision.
   *
   * @return the rounding precision
   */
  public int getRoundingPrecision() {
    return roundingPrecision;
  }

  /**
   * Sets the rounding precision.
   *
   * @param roundingPrecision the new rounding precision
   */
  public void setRoundingPrecision(int roundingPrecision) {
    this.roundingPrecision = roundingPrecision;
  }

  /**
   * Select an index from the text panel.
   *
   * @param selectedIndex the selected index
   */
  public void select(int selectedIndex) {
    if (selectedIndex < 0 || selectedIndex >= tp.getLineCount()) {
      return;
    }
    tp.setSelection(selectedIndex, selectedIndex);
    if (roiPainter != null) {
      roiPainter.selected(selectedIndex);
    }
  }

  /**
   * Select a range of indices from the text panel.
   *
   * @param selectionStart the selection start
   * @param selectionEnd the selection end
   */
  public void select(int selectionStart, int selectionEnd) {
    if (selectionStart < 0 || selectionStart >= tp.getLineCount()) {
      return;
    }
    if (selectionEnd < selectionStart || selectionEnd >= tp.getLineCount()) {
      return;
    }
    tp.setSelection(selectionStart, selectionEnd);
    if (roiPainter != null) {
      roiPainter.selected(selectionStart, selectionEnd);
    }
  }

  /**
   * Adds the mouse listener to the text panel.
   *
   * @param listener the listener
   */
  public void addTextPanelMouseListener(MouseListener listener) {
    tp.addMouseListener(listener);
  }

  /**
   * Removes the mouse listener to the text panel.
   *
   * @param listener the listener
   */
  public void removeTextPanelMouseListener(MouseListener listener) {
    tp.removeMouseListener(listener);
  }

  /**
   * Gets the text panel.
   *
   * @return the text panel
   */
  public TextPanel getTextPanel() {
    return tp;
  }

  // TODO - Extend the selection functionality so that selectionListeners can be added
  // to receive notification whan items from the table are selected.
  // public void addSelectionListener()
  // {
  //
  // }
  //
  // public void removeSelectionListener()
  // {
  //
  // }
}
