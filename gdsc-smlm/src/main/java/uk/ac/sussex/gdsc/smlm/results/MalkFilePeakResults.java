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

package uk.ac.sussex.gdsc.smlm.results;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.NoSuchElementException;
import java.util.Scanner;
import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.data.utils.IdentityTypeConverter;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;

/**
 * Saves the fit results to file using the simple MALK file format (Molecular Accuracy Localisation
 * Keep). This consists of [X,Y,T,Signal] data in a white-space separated format. Comments are
 * allowed with the # character.
 */
public class MalkFilePeakResults extends FilePeakResults {
  /**
   * Converter to change the distances to nm. It is created in {@link #begin()} but may be null.
   */
  protected TypeConverter<DistanceUnit> toNmConverter;
  /**
   * Converter to change the intensity to photons. It is created in {@link #begin()} but may be
   * null.
   */
  protected TypeConverter<IntensityUnit> toPhotonConverter;

  private Writer out;

  /**
   * Instantiates a new MALK file peak results.
   *
   * @param filename the filename
   */
  public MalkFilePeakResults(String filename) {
    super(filename);
  }

  @Override
  protected String getHeaderEnd() {
    return null;
  }

  @Override
  protected String getVersion() {
    return "MALK";
  }

  @Override
  protected void openOutput() {
    out = new BufferedWriter(new OutputStreamWriter(fos, StandardCharsets.UTF_8));
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
    } catch (final Exception ignored) {
      // Ignore exception
    } finally {
      fos = null;
    }
  }

  @Override
  public synchronized void begin() {
    // Ensure we write out in nm and photons if possible.
    if (hasCalibration()) {
      // Copy it so it can be modified
      final CalibrationWriter cw = new CalibrationWriter(getCalibration());

      // Create converters
      try {
        toNmConverter = cw.getDistanceConverter(DistanceUnit.NM);
        cw.setDistanceUnit(DistanceUnit.NM);
      } catch (final ConversionException ignored) {
        // Gracefully fail so ignore this
      }
      try {
        toPhotonConverter = cw.getIntensityConverter(IntensityUnit.PHOTON);
        cw.setIntensityUnit(IntensityUnit.PHOTON);
      } catch (final ConversionException ignored) {
        // Gracefully fail so ignore this
      }

      setCalibration(cw.getCalibration());
    }

    // The data loses PSF information so reset this to a custom type with
    // no additional parameters.
    setPsf(PsfHelper.create(PSFType.CUSTOM));

    super.begin();

    // Create converters to avoid null pointers
    if (toNmConverter == null) {
      toNmConverter = new IdentityTypeConverter<>(null);
    }
    if (toPhotonConverter == null) {
      toPhotonConverter = new IdentityTypeConverter<>(null);
    }
  }

  @Override
  protected String[] getHeaderComments() {
    final String[] comments = new String[3];
    int count = 0;
    if (hasCalibration()) {
      final CalibrationReader cr = getCalibrationReader();
      if (cr.hasNmPerPixel()) {
        comments[count++] =
            String.format("Pixel pitch %s (nm)", MathUtils.rounded(cr.getNmPerPixel()));
      }
      if (cr.hasCountPerPhoton()) {
        comments[count++] =
            String.format("Gain %s (Count/photon)", MathUtils.rounded(cr.getCountPerPhoton()));
      }
      if (cr.hasExposureTime()) {
        comments[count++] = String.format("Exposure time %s (seconds)",
            MathUtils.rounded(cr.getExposureTime() * 1e-3));
      }
    }
    return Arrays.copyOf(comments, count);
  }

  @Override
  protected String[] getFieldNames() {
    final String[] names = {"X", "Y", "Frame", "Signal"};
    if (toNmConverter != null) {
      names[0] += " (nm)";
      names[1] += " (nm)";
    }
    if (toPhotonConverter != null) {
      names[3] += " (photon)";
    }
    return names;
  }

  @Override
  public void add(int peak, int origX, int origY, float origValue, double error, float noise,
      float meanIntensity, float[] params, float[] paramsStdDev) {
    if (fos == null) {
      return;
    }

    final StringBuilder sb = new StringBuilder(100);

    addStandardData(sb, params[PeakResult.X], params[PeakResult.Y], peak,
        params[PeakResult.INTENSITY]);

    writeResult(1, sb.toString());
  }

  @Override
  public void add(PeakResult result) {
    if (fos == null) {
      return;
    }

    final StringBuilder sb = new StringBuilder(100);

    addStandardData(sb, result.getXPosition(), result.getYPosition(), result.getFrame(),
        result.getIntensity());

    writeResult(1, sb.toString());
  }

  private void addStandardData(StringBuilder sb, final float x, final float y, final int frame,
      final float signal) {
    // @formatter:off
    sb.append(toNmConverter.convert(x)).append('\t')
      .append(toNmConverter.convert(y)).append('\t')
      .append(frame).append('\t')
      .append(toPhotonConverter.convert(signal)).append('\n');
    // @formatter:on
  }

  @Override
  public void addAll(PeakResult[] results) {
    if (fos == null) {
      return;
    }

    int count = 0;

    final StringBuilder sb = new StringBuilder(2000);
    for (final PeakResult result : results) {
      // Add the standard data
      addStandardData(sb, result.getXPosition(), result.getYPosition(), result.getFrame(),
          result.getIntensity());

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
   * Adds all the results in the cluster.
   *
   * @param cluster the cluster
   */
  protected void addAll(Cluster cluster) {
    addAll(cluster.getPoints());
  }

  @Override
  protected void sort() throws IOException {
    final LocalList<Result> results = new LocalList<>(size);
    final StringBuilder header = new StringBuilder();

    final Path path = Paths.get(filename);
    try (BufferedReader input = Files.newBufferedReader(path)) {

      String line;
      // Skip the header
      for (line = input.readLine(); line != null; line = input.readLine()) {
        if (line.charAt(0) != '#') {
          // This is the first record
          results.add(new Result(line));
          break;
        }
        header.append(line).append(System.lineSeparator());
      }

      for (line = input.readLine(); line != null; line = input.readLine()) {
        results.add(new Result(line));
      }
    }

    Collections.sort(results, (r1, r2) -> Integer.compare(r1.slice, r2.slice));

    try (BufferedWriter output = Files.newBufferedWriter(path)) {
      output.write(header.toString());
      for (int i = 0; i < results.size(); i++) {
        output.write(results.unsafeGet(i).line);
        output.newLine();
      }
    }
  }

  /**
   * Store the String representation of the result with the frame (slice) number to be used for
   * sorting.
   */
  private static class Result {
    String line;
    int slice;

    Result(String line) {
      this.line = line;
      extractSlice();
    }

    private void extractSlice() {
      try (Scanner scanner = new Scanner(line)) {
        scanner.useDelimiter("\t");
        scanner.nextFloat(); // X
        scanner.nextFloat(); // Y
        slice = scanner.nextInt();
      } catch (final NoSuchElementException ignored) {
        // Ignore
      }
    }
  }
}
