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
import ij.plugin.PlugIn;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Optional;
import java.util.function.Consumer;
import java.util.function.ToDoubleFunction;
import java.util.function.ToIntFunction;
import java.util.stream.Collectors;
import uk.ac.sussex.gdsc.core.data.utils.Converter;
import uk.ac.sussex.gdsc.core.data.utils.IdentityTypeConverter;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.function.ToFloatFunction;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFParameter;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFParameterUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitConverterUtils;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.SaveLocalisationsSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.SaveLocalisationsSettings.Builder;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

/**
 * Saves localisations to a user-specified delimited text file.
 */
public class SaveLocalisations implements PlugIn {
  private static final String TITLE = "Save Localisations";
  private static final char QUOTE = '"';

  /**
   * The type of unit.
   */
  private enum UnitType {
    FRAME, DISTANCE, INTENSITY, ANGLE, NONE;

    static UnitType from(PSFParameterUnit u) {
      if (u == PSFParameterUnit.DISTANCE) {
        return DISTANCE;
      } else if (u == PSFParameterUnit.INTENSITY) {
        return INTENSITY;
      } else if (u == PSFParameterUnit.ANGLE) {
        return ANGLE;
      }
      return NONE;
    }
  }

  /**
   * Specify a field that can be extracted from the results.
   */
  private abstract static class Field {
    private static class DefaultConverter implements Converter {
      static final DefaultConverter INSTANCE = new DefaultConverter();

      @Override
      public double convert(double value) {
        return value;
      }

      @Override
      public float convert(float value) {
        return value;
      }

      @Override
      public double convertBack(double value) {
        return value;
      }

      @Override
      public float convertBack(float value) {
        return value;
      }

      @Override
      public String getFunction() {
        return "f(x) = x";
      }
    }

    final String id;
    final String name;
    final UnitType unit;
    Converter converter = DefaultConverter.INSTANCE;

    Field(String id, String name, UnitType unit) {
      this.id = id;
      this.name = name;
      this.unit = unit;
    }

    @Override
    public String toString() {
      return id + ": " + name;
    }

    /**
     * Sets the converter.
     *
     * @param converter the new converter
     */
    void setConverter(Converter converter) {
      this.converter = converter;
    }

    /**
     * Append to the StringBuilder.
     *
     * @param r the result
     * @param sb the string builder
     * @return the string builder
     */
    abstract StringBuilder appendTo(PeakResult r, StringBuilder sb);
  }

  /**
   * Specify an int field that can be extracted from the results. No convertion is possible.
   */
  private static class IntField extends Field {
    final ToIntFunction<PeakResult> fun;

    IntField(String id, String name, UnitType unit, ToIntFunction<PeakResult> fun) {
      super(id, name, unit);
      this.fun = fun;
    }

    @Override
    StringBuilder appendTo(PeakResult r, StringBuilder sb) {
      return sb.append(fun.applyAsInt(r));
    }
  }

  /**
   * Specify an int field that can be extracted from the results.
   */
  private static class MappedIntField extends IntField {
    MappedIntField(String id, String name, UnitType unit, ToIntFunction<PeakResult> fun) {
      super(id, name, unit, fun);
      converter = null;
    }

    @Override
    StringBuilder appendTo(PeakResult r, StringBuilder sb) {
      final int i = fun.applyAsInt(r);
      if (converter == null) {
        sb.append(i);
      } else {
        sb.append(converter.convert((double) i));
      }
      return sb;
    }
  }

  /**
   * Specify a float field that can be extracted from the results.
   */
  private static class FloatField extends Field {
    final ToFloatFunction<PeakResult> fun;

    FloatField(String id, String name, UnitType unit, ToFloatFunction<PeakResult> fun) {
      super(id, name, unit);
      this.fun = fun;
    }

    @Override
    StringBuilder appendTo(PeakResult r, StringBuilder sb) {
      return sb.append(converter.convert(fun.applyAsFloat(r)));
    }
  }

  /**
   * Specify a double field that can be extracted from the results.
   */
  private static class DoubleField extends Field {
    final ToDoubleFunction<PeakResult> fun;

    DoubleField(String id, String name, UnitType unit, ToDoubleFunction<PeakResult> fun) {
      super(id, name, unit);
      this.fun = fun;
    }

    @Override
    StringBuilder appendTo(PeakResult r, StringBuilder sb) {
      return sb.append(converter.convert(fun.applyAsDouble(r)));
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "There are no fitting results in memory");
      return;
    }

    final SaveLocalisationsSettings.Builder settings =
        SettingsManager.readSaveLocalisationsSettings(0).toBuilder();

    // Get results
    if (!showInputDialog(settings)) {
      return;
    }

    final MemoryPeakResults results =
        ResultsManager.loadInputResults(settings.getInput(), false, null, null);
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(TITLE, "No results could be loaded");
      return;
    }

    final Map<String, Field> fields = createFieldMap(results);

    if (!showDialog(settings, results, fields)) {
      return;
    }

    final Path path = getPath(settings, results);

    if (saveResults(results, settings, fields, path)) {
      final String msg = "Saved " + TextUtils.pleural(results.size(), "localisation");
      IJ.showStatus(msg);
      ImageJUtils.log(msg + " to " + path);
    }
  }

  /**
   * Show a dialog to obtain the name of the input localisations.
   *
   * @param settings the settings
   * @return true, if successful
   */
  private static boolean showInputDialog(Builder settings) {
    final ExtendedGenericDialog gd = createDialog();

    ResultsManager.addInput(gd, settings.getInput(), InputSource.MEMORY);

    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    settings.setInput(ResultsManager.getInputSource(gd));
    SettingsManager.writeSettings(settings.build());
    return true;
  }

  /**
   * Creates a dialog with standard help URL and message.
   *
   * @return the dialog
   */
  private static ExtendedGenericDialog createDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("save-localisations"));
    gd.addMessage("Save localisations in a user-specified text format.");
    return gd;
  }

  /**
   * Creates the field map.
   *
   * @param results the results
   * @return the map
   */
  private static Map<String, Field> createFieldMap(final MemoryPeakResults results) {
    final Map<String, Field> fields = new LinkedHashMap<>(10);
    // Standard fields.
    // Note that even if background and Z are zero they are always stored in a PeakResult
    put(fields, new MappedIntField("T", "Time", UnitType.FRAME, PeakResult::getFrame));
    put(fields, new FloatField("B", "Background", UnitType.INTENSITY, PeakResult::getBackground));
    put(fields, new FloatField("I", "Intensity", UnitType.INTENSITY, PeakResult::getIntensity));
    put(fields, new FloatField("X", "X", UnitType.DISTANCE, PeakResult::getXPosition));
    put(fields, new FloatField("Y", "Y", UnitType.DISTANCE, PeakResult::getYPosition));
    put(fields, new FloatField("Z", "Z", UnitType.DISTANCE, PeakResult::getZPosition));
    // Additional fields. Currently this plugin ignores:
    // hasNoise
    // hasDeviations
    // hasMeanIntensity
    if (results.hasId()) {
      put(fields, new IntField("ID", "Identifier", UnitType.NONE, PeakResult::getId));
    }
    if (results.hasCategory()) {
      put(fields, new IntField("C", "Category", UnitType.NONE, PeakResult::getCategory));
    }
    if (results.hasPrecision()) {
      // Note: precision is in nm. A custom converter is used.
      put(fields, new DoubleField("P", "Precision", UnitType.NONE, PeakResult::getPrecision));
    }
    if (results.hasEndFrame()) {
      put(fields, new MappedIntField("T2", "EndFrame", UnitType.FRAME, PeakResult::getEndFrame));
    }

    final PSF psf = results.getPsf();
    if (psf != null && psf.getParametersList() != null) {
      for (int i = 0; i < psf.getParametersCount(); i++) {
        final PSFParameter param = psf.getParameters(i);
        final int index = i + PeakResult.STANDARD_PARAMETERS;
        put(fields, new FloatField("P" + i, param.getName(), UnitType.from(param.getUnit()),
            r -> r.getParameter(index)));
      }
    }
    return fields;
  }

  /**
   * Put the field into the map using the id as the key.
   *
   * @param map the map
   * @param field the field
   */
  private static void put(Map<String, Field> map, Field field) {
    map.put(field.id, field);
  }

  /**
   * Show a dialog to obtain the output format for the localisations.
   *
   * @param settings the settings
   * @param results the results
   * @param fields the fields
   * @return true, if successful
   */
  private static boolean showDialog(Builder settings, MemoryPeakResults results,
      Map<String, Field> fields) {
    final ExtendedGenericDialog gd = createDialog();

    // Show the available fields
    final StringBuilder sb = new StringBuilder(1024);
    sb.append("Available fields:\n ");
    for (final Field field : fields.values()) {
      sb.append("\n").append(field.toString());
    }
    gd.addMessage(sb.toString());

    // Strip non-existent fields from the current output 'format'
    final String format = Arrays.stream(settings.getFormat().split(" ")).filter(fields::containsKey)
        .collect(Collectors.joining(" "));

    gd.addStringField("Format", format, 20);
    gd.addStringField("Delimiter", escapeDelimiter(settings.getDelimiter()), 3);
    gd.addDirectoryField("Directory", settings.getDirectory());
    gd.addStringField("Suffix", settings.getFileSuffix());
    gd.addCheckbox("Add_header", settings.getAddHeader());

    // Add options to the dialog to choose units if the calibration exists.
    final CalibrationReader cal = results.getCalibrationReader();
    if (cal.getExposureTime() > 0) {
      final String[] items = SettingsManager.getTimeUnitNames();
      gd.addChoice("Time_unit", Arrays.copyOfRange(items, 1, items.length),
          settings.getTimeUnitValue() - 1);
    }
    if (cal.getDistanceUnitValue() > 0) {
      final String[] items = SettingsManager.getDistanceUnitNames();
      gd.addChoice("Distance_unit", Arrays.copyOfRange(items, 1, items.length),
          settings.getDistanceUnitValue() - 1);
    }
    if (cal.getIntensityUnitValue() > 0) {
      final String[] items = SettingsManager.getIntensityUnitNames();
      gd.addChoice("Intensity_unit", Arrays.copyOfRange(items, 1, items.length),
          settings.getIntensityUnitValue() - 1);
    }
    final boolean showAngle = cal.getAngleUnitValue() > 0
        && fields.values().stream().filter(f -> f.unit == UnitType.ANGLE).findAny().isPresent();
    if (showAngle) {
      final String[] items = SettingsManager.getAngleUnitNames();
      gd.addChoice("Angle_unit", Arrays.copyOfRange(items, 1, items.length),
          settings.getAngleUnitValue() - 1);
    }

    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    settings.setFormat(gd.getNextString().toUpperCase(Locale.ROOT));
    settings.setDelimiter(unescapeDelimiter(gd.getNextString()));
    settings.setDirectory(gd.getNextString());
    settings.setFileSuffix(gd.getNextString());
    settings.setAddHeader(gd.getNextBoolean());

    // Read the output units
    if (cal.getExposureTime() > 0) {
      settings.setTimeUnitValue(gd.getNextChoiceIndex() + 1);
    }
    if (cal.getDistanceUnitValue() > 0) {
      settings.setDistanceUnitValue(gd.getNextChoiceIndex() + 1);
    }
    if (cal.getIntensityUnitValue() > 0) {
      settings.setIntensityUnitValue(gd.getNextChoiceIndex() + 1);
    }
    if (showAngle) {
      settings.setAngleUnitValue(gd.getNextChoiceIndex() + 1);
    }

    SettingsManager.writeSettings(settings.build());

    return true;
  }

  /**
   * Escape the delimiter. This converts:
   *
   * <ul>
   *
   * <li>the tab character to {@code \t}
   *
   * </ul>
   *
   * @param delimiter the delimiter
   * @return the string
   */
  private static String escapeDelimiter(String delimiter) {
    return delimiter.replace("\t", "\\t");
  }

  /**
   * Unescape the delimiter. This converts:
   *
   * <ul>
   *
   * <li>{@code \t} to the tab character
   *
   * </ul>
   *
   * @param delimiter the delimiter
   * @return the string
   */
  private static String unescapeDelimiter(String delimiter) {
    return delimiter.replace("\\t", "\t");
  }

  /**
   * Gets the path for the output file.
   *
   * @param settings the settings
   * @param results the results
   * @return the path
   */
  private static Path getPath(final SaveLocalisationsSettings.Builder settings,
      final MemoryPeakResults results) {
    String name = results.getName();
    final String s = settings.getFileSuffix();
    if (s.length() != 0) {
      name += s.charAt(0) == '.' ? s : "." + s;
    }
    return Paths.get(settings.getDirectory(), name);
  }

  /**
   * Save the results.
   *
   * @param results the results
   * @param settings the settings
   * @param fields the fields
   * @param path the path
   * @return true, if successful
   */
  private static boolean saveResults(MemoryPeakResults results, Builder settings,
      Map<String, Field> fields, Path path) {
    // Get the output fields
    final Field[] output =
        Arrays.stream(settings.getFormat().split(" ")).map(fields::get).toArray(Field[]::new);

    if (output.length == 0) {
      IJ.error(TITLE, "No fields to save");
      return false;
    }

    // Create the converters
    final CalibrationReader cal = results.getCalibrationReader();

    if (cal.getExposureTime() > 0) {
      // Do not use cal.getTimeConverter.
      // That converts from the time units set in the TimeCalibration, e.g. for using a custom
      // calibration reader to map time periods to different units.
      // Currently we do not support storing any time values in non-frame units.
      addConverter(output, UnitType.FRAME, UnitConverterUtils.createConverter(TimeUnit.FRAME,
          settings.getTimeUnit(), cal.getExposureTime()));
    }
    if (cal.getDistanceUnitValue() > 0) {
      addConverter(output, UnitType.DISTANCE, cal.getDistanceConverter(settings.getDistanceUnit()));
      // Custom converter for precision (which is always in nm)
      final Optional<Field> result = Arrays.stream(output).filter(x -> "P".equals(x.id)).findAny();
      if (result.isPresent()) {
        result.get().setConverter(UnitConverterUtils.createConverter(DistanceUnit.NM,
            settings.getDistanceUnit(), cal.getNmPerPixel()));
      }
    }
    if (cal.getIntensityUnitValue() > 0) {
      addConverter(output, UnitType.INTENSITY,
          cal.getIntensityConverter(settings.getIntensityUnit()));
    }
    if (cal.getAngleUnitValue() > 0) {
      // If no angle units were present then the angle unit may not be set.
      // Get the converter in a safe method.
      addConverter(output, UnitType.ANGLE, cal.getAngleConverterSafe(settings.getAngleUnit()));
    }

    final StringBuilder sb = new StringBuilder(1024);
    final String delim = settings.getDelimiter();

    // Create a function to add the delimiter
    Consumer<StringBuilder> addDelim;
    if (delim.length() == 1) {
      final char c = delim.charAt(0);
      addDelim = s -> s.append(c);
    } else {
      addDelim = s -> s.append(delim);
    }

    // Create the header. This may raise an exception if the delimiter is in the column names.
    if (settings.getAddHeader()) {
      appendColumnName(sb, output[0].name, delim);
      for (int i = 1; i < output.length; i++) {
        addDelim.accept(sb);
        appendColumnName(sb, output[i].name, delim);
      }
    }

    try (BufferedWriter out = Files.newBufferedWriter(path)) {
      if (settings.getAddHeader()) {
        out.write(sb.toString());
        out.newLine();
      }

      // Create a function to write a record.
      // No use of a ticker here for small data.
      final Ticker ticker = results.size() < (1 << 20) ? Ticker.getDefaultInstance()
          : ImageJUtils.createTicker(results.size(), 1);
      final PeakResultProcedure p = r -> {
        sb.setLength(0);
        output[0].appendTo(r, sb);
        for (int i = 1; i < output.length; i++) {
          addDelim.accept(sb);
          output[i].appendTo(r, sb);
        }
        try {
          out.write(sb.toString());
          out.newLine();
        } catch (final IOException ex) {
          throw new UncheckedIOException(ex);
        }
        ticker.tick();
      };

      results.forEach(p);
      ticker.stop();
    } catch (final IOException ex) {
      throw new UncheckedIOException(ex);
    }

    return true;
  }

  /**
   * Adds the converter for the specified unit type.
   *
   * @param fields the fields
   * @param type the unit type
   * @param c the converter
   */
  private static void addConverter(Field[] fields, UnitType type, Converter c) {
    if (c instanceof IdentityTypeConverter) {
      return;
    }
    for (final Field f : fields) {
      if (f.unit == type) {
        f.setConverter(c);
      }
    }
  }

  private static void appendColumnName(StringBuilder sb, String name, String delim) {
    if (name.contains(delim)) {
      // Try and quote the name
      if (delim.indexOf(QUOTE) != -1) {
        // For now just raise an exception
        throw new IllegalStateException(
            String.format("Column name <%s> contains delimiter <%s> and unable to use quote <%c>",
                name, delim, QUOTE));
      }
      sb.append(QUOTE).append(name).append(QUOTE);
    } else {
      sb.append(name);
    }
  }
}
