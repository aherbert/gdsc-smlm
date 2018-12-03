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
package uk.ac.sussex.gdsc.smlm.results;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.InputMismatchException;
import java.util.NoSuchElementException;
import java.util.Scanner;

import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.data.utils.Converter;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.TurboList;
import uk.ac.sussex.gdsc.smlm.data.config.ConfigurationException;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

/**
 * Saves the fit results to file.
 */
public class TextFilePeakResults extends SMLMFilePeakResults {
  private Gaussian2DPeakResultCalculator calculator;

  private PeakResultConversionHelper helper;
  private Converter[] converters;

  private DistanceUnit distanceUnit = null;
  private IntensityUnit intensityUnit = null;
  private AngleUnit angleUnit = null;
  private boolean computePrecision = false;
  private boolean canComputePrecision = false;

  private OutputStreamWriter out;

  /**
   * Instantiates a new text file peak results.
   *
   * @param filename the filename
   */
  public TextFilePeakResults(String filename) {
    super(filename);
  }

  /**
   * Instantiates a new text file peak results.
   *
   * @param filename the filename
   * @param showDeviations Set to true to show deviations
   */
  public TextFilePeakResults(String filename, boolean showDeviations) {
    super(filename, showDeviations);
  }

  /**
   * Instantiates a new text file peak results.
   *
   * @param filename the filename
   * @param showDeviations Set to true to show deviations
   * @param showEndFrame Set to true to show the end frame
   */
  public TextFilePeakResults(String filename, boolean showDeviations, boolean showEndFrame) {
    super(filename, showDeviations, showEndFrame);
  }

  /**
   * Instantiates a new text file peak results.
   *
   * @param filename the filename
   * @param showDeviations Set to true to show deviations
   * @param showEndFrame Set to true to show the end frame
   * @param showId Set to true to show the id
   */
  public TextFilePeakResults(String filename, boolean showDeviations, boolean showEndFrame,
      boolean showId) {
    super(filename, showDeviations, showEndFrame, showId);
  }

  /**
   * Instantiates a new text file peak results.
   *
   * @param filename the filename
   * @param showDeviations Set to true to show deviations
   * @param showEndFrame Set to true to show the end frame
   * @param showId Set to true to show the id
   * @param showPrecision Set to true to show the precision
   */
  public TextFilePeakResults(String filename, boolean showDeviations, boolean showEndFrame,
      boolean showId, boolean showPrecision) {
    super(filename, showDeviations, showEndFrame, showId, showPrecision);
  }

  @Override
  protected void openOutput() {
    try {
      out = new OutputStreamWriter(fos, "UTF-8");
    } catch (final UnsupportedEncodingException ex) {
      throw new RuntimeException(ex);
    }
  }

  @Override
  protected void write(String data) {
    try {
      out.write(data);
    } catch (final IOException ex) {
      closeOutput();
    }
  }

  @Override
  protected void closeOutput() {
    if (fos == null) {
      return;
    }

    try {
      // Make sure we close the writer since it may be buffered
      out.close();
    } catch (final Exception ex) {
      // Ignore exception
    } finally {
      fos = null;
    }
  }

  /** {@inheritDoc} */
  @Override
  public void begin() {
    calculator = null;
    canComputePrecision = false;

    if (isShowPrecision() && hasCalibration()) {
      // Determine if we can compute the precision using the current settings
      if (computePrecision) {
        try {
          calculator = Gaussian2DPeakResultHelper.create(getPSF(), getCalibrationReader(),
              Gaussian2DPeakResultHelper.LSE_PRECISION);
          canComputePrecision = true;
        } catch (final ConfigurationException ex) {
          // Ignore
        } catch (final ConversionException ex) {
          // Ignore
        }
      }
    }

    // We must correctly convert all the PSF parameter types
    helper = new PeakResultConversionHelper(getCalibration(), getPSF());
    helper.setIntensityUnit(intensityUnit);
    helper.setDistanceUnit(distanceUnit);
    helper.setAngleUnit(angleUnit);
    converters = helper.getConverters();
    // Update the calibration if converters were created
    if (helper.isCalibrationChanged()) {
      setCalibration(helper.getCalibration());
    }

    super.begin();
  }

  /**
   * @return The names of the fields in each record. Will be the last comment of the header
   */
  @Override
  protected String[] getFieldNames() {
    final String[] unitNames = helper.getUnitNames();

    final ArrayList<String> names = new ArrayList<>(20);
    if (isShowId()) {
      names.add("Id");
    }
    names.add(peakIdColumnName);
    if (isShowEndFrame()) {
      names.add("End " + peakIdColumnName);
    }
    names.add("origX");
    names.add("origY");
    names.add("origValue");
    names.add("Error");
    String noiseField = "Noise";
    if (!TextUtils.isNullOrEmpty(unitNames[PeakResult.INTENSITY])) {
      noiseField += " (" + (unitNames[PeakResult.INTENSITY] + ")");
    }
    names.add(noiseField);
    String meanIntensityField = "Mean";
    if (!TextUtils.isNullOrEmpty(unitNames[PeakResult.INTENSITY])) {
      meanIntensityField += " (" + (unitNames[PeakResult.INTENSITY] + ")");
    }
    names.add(meanIntensityField);

    final String[] fields = helper.getNames();

    for (int i = 0; i < fields.length; i++) {
      String f = fields[i];
      // Add units
      if (!TextUtils.isNullOrEmpty(unitNames[i])) {
        f += " (" + unitNames[i] + ")";
      }
      names.add(f);
      if (isShowDeviations()) {
        names.add("+/-");
      }
    }
    if (isShowPrecision()) {
      names.add("Precision (nm)");
    }
    return names.toArray(new String[names.size()]);
  }

  /** {@inheritDoc} */
  @Override
  public void add(int peak, int origX, int origY, float origValue, double error, float noise,
      float meanIntensity, float[] params, float[] paramsStdDev) {
    if (fos == null) {
      return;
    }

    final StringBuilder sb = new StringBuilder();
    addStandardData(sb, 0, peak, peak, origX, origY, origValue, error, noise, meanIntensity);

    // Add the parameters
    if (isShowDeviations()) {
      if (paramsStdDev != null) {
        checkSize(converters.length, params, paramsStdDev);
        for (int i = 0; i < converters.length; i++) {
          add(sb, converters[i].convert(params[i]));
          add(sb, converters[i].convert(paramsStdDev[i]));
        }
      } else {
        checkSize(converters.length, params);
        for (int i = 0; i < converters.length; i++) {
          add(sb, converters[i].convert(params[i]));
          sb.append("\t0");
        }
      }
    } else {
      checkSize(converters.length, params);
      for (int i = 0; i < converters.length; i++) {
        add(sb, converters[i].convert(params[i]));
      }
    }
    if (isShowPrecision()) {
      if (canComputePrecision) {
        addPrecision(sb, calculator.getLSEPrecision(params, noise), true);
      } else {
        sb.append("\t0");
      }
    }
    sb.append('\n');
    writeResult(1, sb.toString());
  }

  private void addStandardData(StringBuilder sb, final int id, final int peak, final int endPeak,
      final int origX, final int origY, final float origValue, final double error,
      final float noise, float meanIntensity) {
    if (isShowId()) {
      sb.append(id);
      sb.append('\t');
    }
    sb.append(peak);
    if (isShowEndFrame()) {
      sb.append('\t').append(endPeak);
    }
    sb.append('\t').append(origX).append('\t').append(origY);
    add(sb, origValue);
    add(sb, error);
    add(sb, converters[PeakResult.INTENSITY].convert(noise));
    add(sb, converters[PeakResult.INTENSITY].convert(meanIntensity));
  }

  private static void add(StringBuilder sb, float value) {
    sb.append('\t').append(value);
  }

  private static void add(StringBuilder sb, double value) {
    sb.append('\t').append(value);
  }

  private static void addPrecision(StringBuilder sb, double value, boolean computed) {
    // Cast to a float as the precision is probably limited in significant figures
    sb.append('\t').append((float) value);
    if (computed) {
      sb.append('*');
    }
  }

  @Override
  public void add(PeakResult result) {
    if (fos == null) {
      return;
    }

    final StringBuilder sb = new StringBuilder();
    add(sb, result);
    writeResult(1, sb.toString());
  }

  private void add(StringBuilder sb, PeakResult result) {
    addStandardData(sb, result.getId(), result.getFrame(), result.getEndFrame(), result.getOrigX(),
        result.getOrigY(), result.getOrigValue(), result.getError(), result.getNoise(),
        result.getMeanIntensity());

    // Add the parameters
    final float[] params = result.getParameters();
    if (isShowDeviations()) {
      final float[] paramsStdDev = result.getParameterDeviations();
      if (paramsStdDev != null) {
        checkSize(converters.length, params, paramsStdDev);
        for (int i = 0; i < converters.length; i++) {
          add(sb, converters[i].convert(params[i]));
          add(sb, converters[i].convert(paramsStdDev[i]));
        }
      } else {
        checkSize(converters.length, params);
        for (int i = 0; i < converters.length; i++) {
          add(sb, converters[i].convert(params[i]));
          sb.append("\t0");
        }
      }
    } else {
      checkSize(converters.length, params);
      for (int i = 0; i < converters.length; i++) {
        add(sb, converters[i].convert(params[i]));
      }
    }
    if (isShowPrecision()) {
      if (result.hasPrecision()) {
        addPrecision(sb, result.getPrecision(), false);
      } else if (canComputePrecision) {
        addPrecision(sb, calculator.getLSEPrecision(params, result.getNoise()), true);
      } else {
        sb.append("\t0");
      }
    }
    sb.append('\n');
  }

  @Override
  public void addAll(PeakResult[] results) {
    if (fos == null) {
      return;
    }

    int count = 0;

    final StringBuilder sb = new StringBuilder();
    for (final PeakResult result : results) {
      add(sb, result);

      // Flush the output to allow for very large input lists
      if (++count >= 20) {
        writeResult(count, sb.toString());
        if (!isActive()) {
          return;
        }
        sb.setLength(0);
        count = 0;
      }
    }
    writeResult(count, sb.toString());
  }

  /**
   * Output a cluster to the results file. <p> Note: This is not synchronised
   *
   * @param cluster the cluster
   */
  public void addCluster(Cluster cluster) {
    if (fos == null) {
      return;
    }
    if (cluster.size() > 0) {
      final float[] centroid = cluster.getCentroid();
      writeResult(0,
          String.format("#Cluster %f %f (+/-%f) n=%d\n",
              converters[PeakResult.X].convert(centroid[0]),
              converters[PeakResult.X].convert(centroid[1]),
              converters[PeakResult.X].convert(cluster.getStandardDeviation()), cluster.size()));
      addAll(cluster);
    }
  }

  /**
   * Adds all the results from the cluster.
   *
   * @param cluster the cluster
   */
  protected void addAll(Cluster cluster) {
    if (!isShowId() || cluster.getId() == 0) {
      addAll(cluster.getPoints().toArray());
    } else {
      // Store the ID from the trace
      final int id = cluster.getId();
      final ArrayPeakResultStore results2 = new ArrayPeakResultStore(cluster.size());
      cluster.getPoints().forEach(new PeakResultProcedure() {
        @Override
        public void execute(PeakResult result) {
          if (result.getId() == id) {
            results2.add(result);
          } else {
            results2.add(new ExtendedPeakResult(result.getFrame(), result.getOrigX(),
                result.getOrigY(), result.getOrigValue(), result.getError(), result.getNoise(),
                result.getMeanIntensity(), result.getParameters(), result.getParameterDeviations(),
                result.getEndFrame(), id));
          }
        }
      });
      addAll(results2.toArray());
    }
  }

  /**
   * Output a trace to the results file. <p> Note: This is not synchronised
   *
   * @param trace the trace
   */
  public void addTrace(Trace trace) {
    if (fos == null) {
      return;
    }
    if (trace.size() > 0) {
      final float[] centroid = trace.getCentroid();
      writeResult(0,
          String.format("#Trace %f %f (+/-%f) n=%d, b=%d, on=%f, off=%f, signal= %f\n",
              converters[PeakResult.X].convert(centroid[0]),
              converters[PeakResult.X].convert(centroid[1]),
              converters[PeakResult.X].convert(trace.getStandardDeviation()), trace.size(),
              trace.getNBlinks(), trace.getOnTime(), trace.getOffTime(),
              converters[PeakResult.INTENSITY].convert(trace.getSignal())));
      addAll(trace);
    }
  }

  /**
   * Output a comment to the results file. <p> Note: This is not synchronised
   *
   * @param text the text
   */
  public void addComment(String text) {
    if (fos == null) {
      return;
    }
    // Ensure comments are preceded by the comment character
    if (!text.startsWith("#")) {
      text = "#" + text;
    }
    if (text.contains("\n")) {
      text.replace("\n", "\n#");
    }
    if (!text.endsWith("\n")) {
      text += "\n";
    }
    writeResult(0, text);
  }

  @Override
  protected synchronized void writeResult(int count, String result) {
    // In case another thread caused the output to close
    if (fos == null) {
      return;
    }
    size += count;
    try {
      out.write(result);
    } catch (final IOException ioe) {
      closeOutput();
    }
  }

  /** {@inheritDoc} */
  @Override
  protected void sort() throws IOException {
    try (BufferedReader input = new BufferedReader(new FileReader(filename))) {
      final TurboList<Result> results = new TurboList<>(size);
      final StringBuilder header = new StringBuilder();

      String line;
      // Skip the header
      while ((line = input.readLine()) != null) {
        if (line.charAt(0) != '#') {
          // This is the first record
          results.add(new Result(line));
          break;
        }
        header.append(line).append('\n');
      }

      while ((line = input.readLine()) != null) {
        results.add(new Result(line));
      }

      input.close();

      Collections.sort(results);

      try (BufferedWriter output = new BufferedWriter(new FileWriter(filename))) {
        output.write(header.toString());
        for (int i = 0; i < results.size(); i++) {
          output.write(results.getf(i).line);
          output.write("\n");
        }
      }
    }
  }

  private class Result implements Comparable<Result> {
    String line;
    int slice = 0;

    public Result(String line) {
      this.line = line;
      extractSlice();
    }

    private void extractSlice() {
      try (Scanner scanner = new Scanner(line)) {
        scanner.useDelimiter("\t");
        slice = scanner.nextInt();
        if (isShowId()) {
          // The peak is the second column
          slice = scanner.nextInt();
        }
      } catch (final InputMismatchException ex) {
        // Ignore
      } catch (final NoSuchElementException ex) {
        // Ignore
      }
    }

    @Override
    public int compareTo(Result o) {
      // Sort by slice number
      // (Note: peak height is already done in the run(...) method)
      if (slice < o.slice) {
        return -1;
      }
      if (slice > o.slice) {
        return 1;
      }
      return 0;
    }
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
   * @param computePrecision set to true to compute the precision
   */
  public void setComputePrecision(boolean computePrecision) {
    this.computePrecision = computePrecision;
  }
}
