/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.NoSuchElementException;
import java.util.Scanner;
import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.data.utils.Converter;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.ConfigurationException;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

/**
 * Saves the fit results to file.
 */
public class TextFilePeakResults extends SmlmFilePeakResults {
  /** The batch size to write to file in one operation. */
  private static final int BATCH_SIZE = 20;

  private Gaussian2DPeakResultCalculator calculator;

  private PeakResultConversionHelper helper;
  private Converter[] converters;
  private int recordSize;

  private DistanceUnit distanceUnit;
  private IntensityUnit intensityUnit;
  private AngleUnit angleUnit;
  private boolean computePrecision;
  private boolean canComputePrecision;

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

  /**
   * Instantiates a new text file peak results.
   *
   * @param filename the filename
   * @param showDeviations Set to true to show deviations
   * @param showEndFrame Set to true to show the end frame
   * @param showId Set to true to show the id
   * @param showPrecision Set to true to show the precision
   * @param showCategory Set to true to show the category
   */
  public TextFilePeakResults(String filename, boolean showDeviations, boolean showEndFrame,
      boolean showId, boolean showPrecision, boolean showCategory) {
    super(filename, showDeviations, showEndFrame, showId, showPrecision, showCategory);
  }

  @Override
  protected void openOutput() {
    out = new OutputStreamWriter(fos, StandardCharsets.UTF_8);
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

  @Override
  public synchronized void begin() {
    calculator = null;
    canComputePrecision = false;

    if (isShowPrecision() && hasCalibration() && computePrecision) {
      // Determine if we can compute the precision using the current settings
      try {
        calculator = Gaussian2DPeakResultHelper.create(getPsf(), getCalibrationReader(),
            Gaussian2DPeakResultHelper.LSE_PRECISION);
        canComputePrecision = true;
      } catch (final ConfigurationException | ConversionException ex) {
        // Ignore
      }
    }

    // We must correctly convert all the PSF parameter types
    helper = new PeakResultConversionHelper(getCalibration(), getPsf());
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

  @Override
  protected String[] getFieldNames() {
    final String[] unitNames = helper.getUnitNames();

    // Count the field types: int, float, double
    // CHECKSTYLE.OFF: VariableDeclarationUsageDistanceCheck
    int ci = 0;
    int cf = 0;
    int cd = 0;
    // CHECKSTYLE.ON: VariableDeclarationUsageDistanceCheck

    final ArrayList<String> names = new ArrayList<>(20);
    if (isShowId()) {
      names.add("Id");
      ci++;
    }
    if (isShowCategory()) {
      names.add("Category");
      ci++;
    }
    names.add(peakIdColumnName);
    ci++;
    if (isShowEndFrame()) {
      names.add("End " + peakIdColumnName);
      ci++;
    }
    names.add("origX");
    names.add("origY");
    ci += 2;
    names.add("origValue");
    cf++;
    names.add("Error");
    cd++;
    String noiseField = "Noise";
    if (!TextUtils.isNullOrEmpty(unitNames[PeakResult.INTENSITY])) {
      noiseField += " (" + (unitNames[PeakResult.INTENSITY] + ")");
    }
    names.add(noiseField);
    cf++;
    String meanIntensityField = "Mean";
    if (!TextUtils.isNullOrEmpty(unitNames[PeakResult.INTENSITY])) {
      meanIntensityField += " (" + (unitNames[PeakResult.INTENSITY] + ")");
    }
    names.add(meanIntensityField);
    cf++;

    final String[] fields = helper.getNames();

    cf += fields.length * (isShowDeviations() ? 2 : 1);
    for (int i = 0; i < fields.length; i++) {
      String field = fields[i];
      // Add units
      if (!TextUtils.isNullOrEmpty(unitNames[i])) {
        field += " (" + unitNames[i] + ")";
      }
      names.add(field);
      if (isShowDeviations()) {
        names.add("+/-");
      }
    }
    if (isShowPrecision()) {
      names.add("Precision (nm)");
      cf++;
    }

    // Get the estimate record size.
    // Count a tab delimiter for each field and the line separator
    recordSize = ci + cf + cd + System.lineSeparator().length();
    // Integer fields:
    // Integer.toString(Integer.MIN_VALUE).length() == len(-2147483648) == 11
    recordSize += ci * 11;
    // Float fields:
    // Float.toString(-Float.MIN_NORMAL).length() == len(-1.17549435E-38) == 15
    recordSize += cf * 15;
    // Double fields:
    // Double.toString(-Double.MIN_NORMAL).length() == len(-2.2250738585072014E-308) == 24
    recordSize += cd * 24;

    return names.toArray(new String[0]);
  }

  @Override
  public void add(int peak, int origX, int origY, float origValue, double error, float noise,
      float meanIntensity, float[] params, float[] paramsStdDev) {
    if (fos == null) {
      return;
    }

    final StringBuilder sb = new StringBuilder(recordSize);
    addStandardData(sb, 0, 0, peak, peak, origX, origY, origValue, error, noise, meanIntensity);

    // Add the parameters
    if (isShowDeviations()) {
      if (paramsStdDev != null) {
        checkSize(converters.length, params, paramsStdDev);
        for (int i = 0; i < converters.length; i++) {
          addFloat(sb, converters[i].convert(params[i]));
          addFloat(sb, converters[i].convert(paramsStdDev[i]));
        }
      } else {
        checkSize(converters.length, params);
        for (int i = 0; i < converters.length; i++) {
          addFloat(sb, converters[i].convert(params[i]));
          sb.append("\t0");
        }
      }
    } else {
      checkSize(converters.length, params);
      for (int i = 0; i < converters.length; i++) {
        addFloat(sb, converters[i].convert(params[i]));
      }
    }
    if (isShowPrecision()) {
      if (canComputePrecision) {
        addPrecision(sb, calculator.getLsePrecision(params, noise), true);
      } else {
        sb.append("\t0");
      }
    }
    sb.append(System.lineSeparator());
    writeResult(1, sb.toString());
  }

  @Override
  public void add(PeakResult result) {
    if (fos == null) {
      return;
    }

    final StringBuilder sb = new StringBuilder(recordSize);
    add(sb, result);
    writeResult(1, sb.toString());
  }

  private void add(StringBuilder sb, PeakResult result) {
    addStandardData(sb, result.getId(), result.getCategory(), result.getFrame(),
        result.getEndFrame(), result.getOrigX(), result.getOrigY(), result.getOrigValue(),
        result.getError(), result.getNoise(), result.getMeanIntensity());

    // Add the parameters
    final float[] params = result.getParameters();
    if (isShowDeviations()) {
      final float[] paramsStdDev = result.getParameterDeviations();
      if (paramsStdDev != null) {
        checkSize(converters.length, params, paramsStdDev);
        for (int i = 0; i < converters.length; i++) {
          addFloat(sb, converters[i].convert(params[i]));
          addFloat(sb, converters[i].convert(paramsStdDev[i]));
        }
      } else {
        checkSize(converters.length, params);
        for (int i = 0; i < converters.length; i++) {
          addFloat(sb, converters[i].convert(params[i]));
          sb.append("\t0");
        }
      }
    } else {
      checkSize(converters.length, params);
      for (int i = 0; i < converters.length; i++) {
        addFloat(sb, converters[i].convert(params[i]));
      }
    }
    if (isShowPrecision()) {
      if (result.hasPrecision()) {
        addPrecision(sb, result.getPrecision(), false);
      } else if (canComputePrecision) {
        addPrecision(sb, calculator.getLsePrecision(params, result.getNoise()), true);
      } else {
        sb.append("\t0");
      }
    }
    sb.append(System.lineSeparator());
  }

  private void addStandardData(StringBuilder sb, final int id, final int category, final int peak,
      final int endPeak, final int origX, final int origY, final float origValue,
      final double error, final float noise, float meanIntensity) {
    if (isShowId()) {
      sb.append(id).append('\t');
    }
    if (isShowCategory()) {
      sb.append(category).append('\t');
    }
    sb.append(peak);
    if (isShowEndFrame()) {
      sb.append('\t').append(endPeak);
    }
    sb.append('\t').append(origX).append('\t').append(origY);
    addFloat(sb, origValue);
    addDouble(sb, error);
    addFloat(sb, converters[PeakResult.INTENSITY].convert(noise));
    addFloat(sb, converters[PeakResult.INTENSITY].convert(meanIntensity));
  }

  private static void addFloat(StringBuilder sb, float value) {
    sb.append('\t').append(value);
  }

  private static void addDouble(StringBuilder sb, double value) {
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
  public void addAll(PeakResult[] results) {
    if (fos == null) {
      return;
    }

    int count = 0;

    final StringBuilder sb = new StringBuilder(BATCH_SIZE * recordSize);
    for (final PeakResult result : results) {
      add(sb, result);

      // Flush the output to allow for very large input lists
      if (++count >= BATCH_SIZE) {
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
      cluster.getPoints().forEach((PeakResultProcedure) result -> {
        if (result.getId() == id) {
          results2.add(result);
        } else {
          // This will maintain the category but change the ID to the cluster id
          final AttributePeakResult r = new AttributePeakResult(result);
          r.setId(id);
          results2.add(r);
        }
      });
      addAll(results2.toArray());
    }
  }

  /**
   * Output a cluster to the results file.
   *
   * <p>Note: This is not synchronised
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
          String.format("#Cluster %f %f (+/-%f) n=%d%n",
              converters[PeakResult.X].convert(centroid[0]),
              converters[PeakResult.X].convert(centroid[1]),
              converters[PeakResult.X].convert(cluster.getStandardDeviation()), cluster.size()));
      addAll(cluster);
    }
  }

  /**
   * Output a trace to the results file.
   *
   * <p>Note: This is not synchronised
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
          String.format(
              "#Trace %f %f (+/-%f) start=%d, end=%d, n=%d, b=%d, on=%f, off=%f, signal= %f%n",
              converters[PeakResult.X].convert(centroid[0]),
              converters[PeakResult.X].convert(centroid[1]),
              converters[PeakResult.X].convert(trace.getStandardDeviation()),
              trace.getHead().getFrame(), trace.getTail().getEndFrame(), trace.size(),
              trace.getBlinks(), trace.getOnTime(), trace.getOffTime(),
              converters[PeakResult.INTENSITY].convert(trace.getSignal())));
      addAll(trace);
    }
  }

  /**
   * Output a comment to the results file.
   *
   * <p>Note: This is not synchronised
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
    // Remove last line separator
    if (text.endsWith(System.lineSeparator())) {
      text = text.substring(0, text.length() - System.lineSeparator().length());
    }
    // Ensure newline in a comment start with '#'
    if (text.contains(System.lineSeparator())) {
      text = text.replace(System.lineSeparator(), System.lineSeparator() + "#");
    }
    // Write with a new line
    writeResult(0, text + System.lineSeparator());
  }

  @Override
  protected synchronized void writeResult(int count, String result) {
    // In case another thread caused the output to close
    if (fos == null) {
      return;
    }
    try {
      out.write(result);
    } catch (final IOException ex) {
      closeOutput();
    }
    size += count;
  }

  @Override
  protected void sort() throws IOException {
    final LocalList<Result> results = new LocalList<>(size);
    final StringBuilder header = new StringBuilder(2048);

    final Path path = Paths.get(filename);
    try (BufferedReader input = Files.newBufferedReader(path)) {

      // Skip optional columns before the slice
      final int skipCount = (isShowId() ? 1 : 0) + (isShowCategory() ? 1 : 0);

      String line;
      // Skip the header
      while ((line = input.readLine()) != null) {
        if (!line.isEmpty() && line.charAt(0) != '#') {
          // This is the first record
          results.add(new Result(line, skipCount));
          break;
        }
        header.append(line).append(System.lineSeparator());
      }

      while ((line = input.readLine()) != null) {
        results.add(new Result(line, skipCount));
      }
    }

    // Sort by slice number
    Collections.sort(results, (r1, r2) -> Integer.compare(r1.slice, r2.slice));

    try (BufferedWriter output = Files.newBufferedWriter(path)) {
      output.write(header.toString());
      for (int i = 0; i < results.size(); i++) {
        output.write(results.unsafeGet(i).line);
        output.newLine();
      }
    }
  }

  private static class Result {
    String line;
    int slice;

    Result(String line, int skipCount) {
      this.line = line;
      try (Scanner scanner = new Scanner(line)) {
        scanner.useDelimiter("\t");
        // Skip optional columns before the slice
        // CHECKSTYLE.OFF: FallThroughCheck
        switch (skipCount) {
          case 2:
            scanner.nextInt();
            // FALL-THROUGH
          case 1:
            scanner.nextInt();
            // FALL-THROUGH
          default:
            slice = scanner.nextInt();
        }
        // CHECKSTYLE.ON: FallThroughCheck
      } catch (final NoSuchElementException ex) {
        // Ignore
      }
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
