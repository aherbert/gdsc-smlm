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

import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import uk.ac.sussex.gdsc.core.utils.BitFlagUtils;
import uk.ac.sussex.gdsc.core.utils.LocalList;

/**
 * Saves the fit results to a binary file format.
 */
public class BinaryFilePeakResults extends SmlmFilePeakResults {
  /** The constant used to mark the end of the header. */
  public static final String END_HEADER = "END_HEADER";
  /** The batch size to write to file in one operation. */
  private static final int BATCH_SIZE = 20;

  private String[] fieldNames;
  private int fieldCount;
  private int recordSize;

  private BufferedOutputStream out;

  /**
   * Instantiates a new binary file peak results.
   *
   * @param filename the filename
   */
  public BinaryFilePeakResults(String filename) {
    super(filename);
  }

  /**
   * Instantiates a new binary file peak results.
   *
   * @param filename the filename
   * @param showDeviations Set to true to show deviations
   */
  public BinaryFilePeakResults(String filename, boolean showDeviations) {
    super(filename, showDeviations);
  }

  /**
   * Instantiates a new binary file peak results.
   *
   * @param filename the filename
   * @param showDeviations Set to true to show deviations
   * @param showEndFrame Set to true to show the end frame
   */
  public BinaryFilePeakResults(String filename, boolean showDeviations, boolean showEndFrame) {
    super(filename, showDeviations, showEndFrame);
  }

  /**
   * Instantiates a new binary file peak results.
   *
   * @param filename the filename
   * @param showDeviations Set to true to show deviations
   * @param showEndFrame Set to true to show the end frame
   * @param showId Set to true to show the id
   */
  public BinaryFilePeakResults(String filename, boolean showDeviations, boolean showEndFrame,
      boolean showId) {
    super(filename, showDeviations, showEndFrame, showId);
  }

  /**
   * Instantiates a new binary file peak results.
   *
   * @param filename the filename
   * @param showDeviations Set to true to show deviations
   * @param showEndFrame Set to true to show the end frame
   * @param showId Set to true to show the id
   * @param showPrecision Set to true to show the precision
   */
  public BinaryFilePeakResults(String filename, boolean showDeviations, boolean showEndFrame,
      boolean showId, boolean showPrecision) {
    super(filename, showDeviations, showEndFrame, showId, showPrecision);
  }

  /**
   * Instantiates a new binary file peak results.
   *
   * @param filename the filename
   * @param showDeviations Set to true to show deviations
   * @param showEndFrame Set to true to show the end frame
   * @param showId Set to true to show the id
   * @param showPrecision Set to true to show the precision
   * @param showCategory Set to true to show the category
   */
  public BinaryFilePeakResults(String filename, boolean showDeviations, boolean showEndFrame,
      boolean showId, boolean showPrecision, boolean showCategory) {
    super(filename, showDeviations, showEndFrame, showId, showPrecision, showCategory);
  }

  @Override
  protected void openOutput() {
    out = new BufferedOutputStream(fos);
  }

  @Override
  protected void write(String data) {
    try {
      out.write(data.getBytes(StandardCharsets.UTF_8));
    } catch (final IOException ex) {
      closeOutput();
    }
  }

  @Override
  protected String getHeaderEnd() {
    return END_HEADER;
  }

  @Override
  protected String[] getHeaderComments() {
    fieldNames = new PeakResultConversionHelper(null, getPsf()).getNames();
    fieldCount = fieldNames.length;
    int flags = 0;

    final String[] comments = new String[2];
    comments[0] = "Records start after the final comment line";
    final StringBuilder sb = new StringBuilder();
    if (isShowEndFrame()) {
      sb.append('i');
      flags += FLAG_END_FRAME;
    }
    if (isShowId()) {
      sb.append('i');
      flags += FLAG_ID;
    }
    if (isShowCategory()) {
      sb.append('i');
      flags += FLAG_CATEGORY;
    }
    sb.append("iiifdff");
    if (isShowDeviations()) {
      for (int i = 0; i < fieldCount; i++) {
        sb.append("ff");
      }
    } else {
      for (int i = 0; i < fieldCount; i++) {
        sb.append('f');
      }
    }
    if (isShowPrecision()) {
      sb.append("f");
      flags += FLAG_PRECISION;
    }
    recordSize = getDataSize(isShowDeviations(), flags, fieldCount);
    comments[1] = "Binary Format (raw Java bytes) = " + sb.toString();
    return comments;
  }

  @Override
  protected String[] getFieldNames() {
    final ArrayList<String> names = new ArrayList<>(20);
    if (isShowId()) {
      names.add("Id");
    }
    if (isShowCategory()) {
      names.add("Category");
    }
    names.add(peakIdColumnName);
    if (isShowEndFrame()) {
      names.add("End " + peakIdColumnName);
    }
    names.add("origX");
    names.add("origY");
    names.add("origValue");
    names.add("Error");
    names.add("Noise");
    names.add("Signal");
    for (int i = 0; i < fieldCount; i++) {
      names.add(fieldNames[i]);
    }
    if (isShowDeviations()) {
      for (int i = 0; i < fieldCount; i++) {
        names.add(fieldNames[i] + " StdDev");
      }
    }
    if (isShowPrecision()) {
      names.add("Precision (nm)");
    }
    return names.toArray(new String[0]);
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
  public void add(int peak, int origX, int origY, float origValue, double error, float noise,
      float meanIntensity, float[] params, float[] paramsStdDev) {
    if (fos == null) {
      return;
    }

    // Buffer the output for the synchronized write method
    try (ByteArrayOutputStream bytes = new ByteArrayOutputStream(recordSize);
        DataOutputStream buffer = new DataOutputStream(bytes)) {
      addResult(buffer, 0, 0, peak, peak, origX, origY, origValue, error, noise, meanIntensity,
          params, paramsStdDev, 0.0);
      buffer.flush();
      writeResult(1, bytes);
    } catch (final IOException ex) {
      // Do nothing - This result will not be added to the file
    }
  }

  @Override
  public void add(PeakResult result) {
    if (fos == null) {
      return;
    }

    // TODO: Rewrite using a ByteBuffer?

    // Buffer the output for the synchronized write method
    try (ByteArrayOutputStream bytes = new ByteArrayOutputStream(recordSize);
        DataOutputStream buffer = new DataOutputStream(bytes)) {
      addResult(buffer, result.getId(), result.getCategory(), result.getFrame(),
          result.getEndFrame(), result.getOrigX(), result.getOrigY(), result.getOrigValue(),
          result.getError(), result.getNoise(), result.getMeanIntensity(), result.getParameters(),
          result.getParameterDeviations(), result.getPrecision());
      buffer.flush();
      writeResult(1, bytes);
    } catch (final IOException ex) {
      // Do nothing - This result will not be added to the file
    }
  }

  private void addResult(DataOutputStream buffer, final int id, final int category, final int peak,
      final int endPeak, final int origX, final int origY, final float origValue,
      final double error, final float noise, float meanIntensity, final float[] params,
      float[] paramsStdDev, double precision) throws IOException {
    if (isShowId()) {
      buffer.writeInt(id);
    }
    if (isShowCategory()) {
      buffer.writeInt(category);
    }
    buffer.writeInt(peak);
    if (isShowEndFrame()) {
      buffer.writeInt(endPeak);
    }
    buffer.writeInt(origX);
    buffer.writeInt(origY);
    buffer.writeFloat(origValue);
    buffer.writeDouble(error);
    buffer.writeFloat(noise);
    buffer.writeFloat(meanIntensity);

    checkSize(fieldCount, params);
    for (int i = 0; i < fieldCount; i++) {
      buffer.writeFloat(params[i]);
    }
    if (isShowDeviations()) {
      if (paramsStdDev == null) {
        for (int i = 0; i < fieldCount; i++) {
          buffer.writeInt(0); // An empty int is the same size as an empty float
        }
      } else {
        checkSize(fieldCount, paramsStdDev);
        for (int i = 0; i < fieldCount; i++) {
          buffer.writeFloat(paramsStdDev[i]);
        }
      }
    }
    if (isShowPrecision()) {
      buffer.writeFloat((float) precision);
    }
  }

  @Override
  public void addAll(PeakResult[] results) {
    if (fos == null) {
      return;
    }

    int count = 0;

    // Buffer the output for the synchronized write method
    try (ByteArrayOutputStream bytes = new ByteArrayOutputStream(BATCH_SIZE * recordSize);
        DataOutputStream buffer = new DataOutputStream(bytes)) {
      for (final PeakResult result : results) {
        addResult(buffer, result.getId(), result.getCategory(), result.getFrame(),
            result.getEndFrame(), result.getOrigX(), result.getOrigY(), result.getOrigValue(),
            result.getError(), result.getNoise(), result.getMeanIntensity(), result.getParameters(),
            result.getParameterDeviations(), result.getPrecision());

        // Flush the output to allow for very large input lists
        if (++count >= BATCH_SIZE) {
          if (!isActive()) {
            return;
          }
          buffer.flush();
          writeResult(count, bytes);
          bytes.reset();
          count = 0;
        }
      }

      buffer.flush();
      writeResult(count, bytes);
    } catch (final IOException ex) {
      // Do nothing - This result will not be added to the file
    }
  }

  private synchronized void writeResult(int count, ByteArrayOutputStream bytes) {
    // In case another thread caused the output to close
    if (fos == null) {
      return;
    }
    try {
      bytes.writeTo(out);
    } catch (final IOException ioe) {
      closeOutput();
    }
    size += count;
  }

  @Override
  protected void sort() throws IOException {
    final LocalList<Result> results = new LocalList<>(size);
    String header;

    try (DataInputStream input = new DataInputStream(new FileInputStream(filename))) {
      header = readHeader(input);

      int flags = 0;
      if (isShowEndFrame()) {
        flags += FLAG_END_FRAME;
      }
      if (isShowId()) {
        flags += FLAG_ID;
      }
      if (isShowCategory()) {
        flags += FLAG_CATEGORY;
      }
      if (isShowPrecision()) {
        flags += FLAG_PRECISION;
      }
      final byte[] line = new byte[getDataSize(isShowDeviations(), flags, fieldCount)];
      final int offset = (isShowId() ? 4 : 0) + (isShowCategory() ? 4 : 0);
      while (input.read(line) == line.length) {
        results.add(new Result(line, offset));
      }
    }

    // Sort by slice number
    Collections.sort(results, (r1, r2) -> Integer.compare(r1.slice, r2.slice));

    // Must write using the same method as the main code so use a FileOutputStream again
    try (FileOutputStream fos = new FileOutputStream(filename);
        BufferedOutputStream output = new BufferedOutputStream(fos)) {
      output.write(header.getBytes(StandardCharsets.UTF_8));
      for (final Result result : results) {
        output.write(result.line);
      }
    }
  }

  /**
   * Read all lines from the input stream that begin with '#' and collates them into a header. Stops
   * reading if a line contains {@value #END_HEADER}.
   *
   * <p>Lines are defined by the '\n' character. The input stream will read the first non-header
   * line unless the header is terminated by the {@value #END_HEADER} tag.
   *
   * @param input the input
   * @return The header
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static String readHeader(DataInputStream input) throws IOException {
    final StringBuilder header = new StringBuilder(1024);
    final StringBuilder line = new StringBuilder(256);
    do {
      readLine(line, input);
      if (line.charAt(0) == '#') {
        header.append(line);
      } else {
        break;
      }
    } while (!startsWithEndHeader(line));
    return header.toString();
  }

  private static void readLine(StringBuilder sb, DataInputStream input) throws IOException {
    sb.setLength(0);
    byte bi;
    do {
      bi = input.readByte();
      sb.append((char) bi);
    } while (bi != '\n');
  }

  /**
   * Check if the line starts with the end header text. Assumes the first character is the '#'
   * comment character.
   *
   * @param line the line
   * @return true if this is the end header text
   */
  private static boolean startsWithEndHeader(StringBuilder line) {
    if (line.length() - 1 < END_HEADER.length()) {
      return false;
    }
    for (int i = END_HEADER.length(); i > 0; i--) {
      if (line.charAt(i) != END_HEADER.charAt(i - 1)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Gets the data size.
   *
   * @param deviations Set to true to show deviations
   * @param flags the flags
   * @param fieldCount the number of fields
   * @return the data size
   */
  static int getDataSize(boolean deviations, int flags, int fieldCount) {
    // iiifdff + n*f
    // or
    // iiifdff + 2*n*f
    // + Extra i added for the end frame
    // + Extra i added for the id
    // + Extra i added for the category
    // + Extra f added for the precision
    int size = 3 * Integer.BYTES + Float.BYTES + Double.BYTES + Float.BYTES + Float.BYTES
        + fieldCount * Float.BYTES;
    if (deviations) {
      size += fieldCount * Float.BYTES;
    }
    if (BitFlagUtils.areSet(flags, FLAG_END_FRAME)) {
      size += Integer.BYTES;
    }
    if (BitFlagUtils.areSet(flags, FLAG_ID)) {
      size += Integer.BYTES;
    }
    if (BitFlagUtils.areSet(flags, FLAG_CATEGORY)) {
      size += Integer.BYTES;
    }
    if (BitFlagUtils.areSet(flags, FLAG_PRECISION)) {
      size += Float.BYTES;
    }
    return size;
  }

  private static class Result {
    byte[] line;
    int slice;

    Result(byte[] line, int offset) {
      this.line = Arrays.copyOf(line, line.length);
      slice = makeInt(line[offset + 0], line[offset + 1], line[offset + 2], line[offset + 3]);
    }

    private static int makeInt(byte b3, byte b2, byte b1, byte b0) {
      // @formatter:off
      return (( b3         << 24) |
              ((b2 & 0xff) << 16) |
              ((b1 & 0xff) <<  8) |
              ( b0 & 0xff       ));
      // @formatter:on
    }
  }

  @Override
  public boolean isBinary() {
    return true;
  }
}
