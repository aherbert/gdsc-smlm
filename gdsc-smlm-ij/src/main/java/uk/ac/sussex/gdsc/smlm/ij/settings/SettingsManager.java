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

package uk.ac.sussex.gdsc.smlm.ij.settings;

import com.google.protobuf.InvalidProtocolBufferException;
import com.google.protobuf.Message;
import com.google.protobuf.MessageOrBuilder;
import com.google.protobuf.Parser;
import com.google.protobuf.util.JsonFormat;
import com.google.protobuf.util.JsonFormat.Printer;
import ij.Prefs;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Reader;
import java.io.StringReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.EnumSet;
import java.util.Locale;
import java.util.logging.Level;
import java.util.logging.Logger;
import uk.ac.sussex.gdsc.core.utils.BitFlagUtils;
import uk.ac.sussex.gdsc.smlm.data.NamedObject;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraModelSettings;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.DataFilterMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.DataFilterType;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitEngineSettings;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitSolver;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.NoiseEstimatorMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.PrecisionMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.SearchMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.AstigmatismModelSettings;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.CubicSplineSettings;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PsfProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsFileFormat;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageSizeMode;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageType;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsTableFormat;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;
import uk.ac.sussex.gdsc.smlm.engine.FitConfiguration;
import uk.ac.sussex.gdsc.smlm.engine.FitEngineConfiguration;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.AstigmatismModelManagerSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.CameraModelAnalysisSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.CameraModelFisherInformationAnalysisSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.CameraModelManagerSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.ClusteringSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.ConfigurationTemplateSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.CreateDataSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.CropResultsSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.CubicSplineManagerSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.DefaultTemplateSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.FailCountManagerSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.GUIFilterSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.ImageJ3DResultsViewerSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.LoadLocalisationsSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.NucleusMaskSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.OpticsSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.PSFCalculatorSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.PSFCreatorSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.PSFEstimatorSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.SpotFitSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.SummariseResultsSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.TcPalmAnalysisSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.TranslateResultsSettings;

/**
 * Manage the settings for ImageJ plugins.
 */
public final class SettingsManager {

  /** The logger. */
  private static Logger logger = Logger.getLogger(SettingsManager.class.getName());

  /** The Constant DEFAULT_FILENAME. */
  private static final String DEFAULT_FILENAME = System.getProperty("user.home")
      + System.getProperty("file.separator") + "gdsc.smlm.settings.xml";

  /** The Constant DEFAULT_DIRECTORY. */
  private static final String DEFAULT_DIRECTORY =
      System.getProperty("user.home") + File.separatorChar + ".gdsc.smlm";

  /** Use this to suppress warnings. */
  public static final int FLAG_SILENT = 0x01;
  /** Use this flag to suppress returning a default instance. */
  public static final int FLAG_NO_DEFAULT = 0x02;
  /** Use this flag to include insignificant whitespace in JSON output. */
  public static final int FLAG_JSON_WHITESPACE = 0x04;
  /** Use this flag to include default values in JSON output. */
  public static final int FLAG_JSON_DEFAULT_VALUES = 0x08;
  /**
   * Use this to show file-not-found warning when reading settings. The default is to suppress the
   * warning as it is usually usually from a settings file that has not been created.
   */
  public static final int FLAG_SHOW_FILE_NOT_FOUND_ON_READ = 0x00000010;

  /** The settings directory. */
  private static String settingsDirectory;

  static {
    setSettingsDirectory(Prefs.get(Constants.settingsDirectory, DEFAULT_DIRECTORY));
  }

  /**
   * Lazy load a default printer.
   */
  private static class PrinterLoader {
    /** The printer. */
    static final Printer PRINTER = JsonFormat.printer();
  }

  /**
   * Lazy load a default parser.
   */
  private static class ParserLoader {
    /** The parser. */
    static final com.google.protobuf.util.JsonFormat.Parser PARSER = JsonFormat.parser();
  }

  /** No public constructor. */
  private SettingsManager() {}

  // Lazy loading of values and names of enums

  /**
   * Lazy loader for the {@link DistanceUnit} enum.
   */
  private static class DistanceUnitLoader {
    /** The enum VALUES. */
    static final DistanceUnit[] VALUES;
    /** The enum NAMES. */
    static final String[] NAMES;

    static {
      final EnumSet<DistanceUnit> d = EnumSet.allOf(DistanceUnit.class);
      d.remove(DistanceUnit.UNRECOGNIZED);
      VALUES = d.toArray(new DistanceUnit[0]);
      NAMES = new String[VALUES.length];
      for (int i = 0; i < VALUES.length; i++) {
        NAMES[i] = getName(UnitHelper.getName(VALUES[i]), UnitHelper.getShortName(VALUES[i]));
      }
    }
  }

  /**
   * Gets the distance unit values.
   *
   * @return the distance unit values
   */
  public static DistanceUnit[] getDistanceUnitValues() {
    return DistanceUnitLoader.VALUES.clone();
  }

  /**
   * Gets the distance unit names.
   *
   * @return the distance unit names
   */
  public static String[] getDistanceUnitNames() {
    return DistanceUnitLoader.NAMES.clone();
  }

  /**
   * Lazy loader for the {@link IntensityUnit} enum.
   */
  private static class IntensityUnitLoader {
    /** The enum VALUES. */
    static final IntensityUnit[] VALUES;
    /** The enum NAMES. */
    static final String[] NAMES;

    static {
      final EnumSet<IntensityUnit> d = EnumSet.allOf(IntensityUnit.class);
      d.remove(IntensityUnit.UNRECOGNIZED);
      VALUES = d.toArray(new IntensityUnit[0]);
      NAMES = new String[VALUES.length];
      for (int i = 0; i < VALUES.length; i++) {
        NAMES[i] = getName(UnitHelper.getName(VALUES[i]), UnitHelper.getShortName(VALUES[i]));
      }
    }
  }

  /**
   * Gets the intensity unit values.
   *
   * @return the intensity unit values
   */
  public static IntensityUnit[] getIntensityUnitValues() {
    return IntensityUnitLoader.VALUES.clone();
  }

  /**
   * Gets the intensity unit names.
   *
   * @return the intensity unit names
   */
  public static String[] getIntensityUnitNames() {
    return IntensityUnitLoader.NAMES.clone();
  }

  /**
   * Lazy loader for the {@link AngleUnit} enum.
   */
  private static class AngleUnitLoader {
    /** The enum VALUES. */
    static final AngleUnit[] VALUES;
    /** The enum NAMES. */
    static final String[] NAMES;

    static {
      final EnumSet<AngleUnit> d = EnumSet.allOf(AngleUnit.class);
      d.remove(AngleUnit.UNRECOGNIZED);
      VALUES = d.toArray(new AngleUnit[0]);
      NAMES = new String[VALUES.length];
      for (int i = 0; i < VALUES.length; i++) {
        NAMES[i] = getName(UnitHelper.getName(VALUES[i]), UnitHelper.getShortName(VALUES[i]));
      }
    }
  }

  /**
   * Gets the angle unit values.
   *
   * @return the angle unit values
   */
  public static AngleUnit[] getAngleUnitValues() {
    return AngleUnitLoader.VALUES.clone();
  }

  /**
   * Gets the angle unit names.
   *
   * @return the angle unit names
   */
  public static String[] getAngleUnitNames() {
    return AngleUnitLoader.NAMES.clone();
  }

  /**
   * Lazy loader for the {@link TimeUnit} enum.
   */
  private static class TimeUnitLoader {
    /** The enum VALUES. */
    static final TimeUnit[] VALUES;
    /** The enum NAMES. */
    static final String[] NAMES;

    static {
      final EnumSet<TimeUnit> d = EnumSet.allOf(TimeUnit.class);
      d.remove(TimeUnit.UNRECOGNIZED);
      VALUES = d.toArray(new TimeUnit[0]);
      NAMES = new String[VALUES.length];
      for (int i = 0; i < VALUES.length; i++) {
        NAMES[i] = getName(UnitHelper.getName(VALUES[i]), UnitHelper.getShortName(VALUES[i]));
      }
    }
  }

  /**
   * Gets the time unit values.
   *
   * @return the time unit values
   */
  public static TimeUnit[] getTimeUnitValues() {
    return TimeUnitLoader.VALUES.clone();
  }

  /**
   * Gets the time unit names.
   *
   * @return the time unit names
   */
  public static String[] getTimeUnitNames() {
    return TimeUnitLoader.NAMES.clone();
  }

  /**
   * Lazy loader for the {@link ResultsImageType} enum.
   */
  private static class ResultsImageTypeLoader {
    /** The enum VALUES. */
    static final ResultsImageType[] VALUES;
    /** The enum NAMES. */
    static final String[] NAMES;

    static {
      final EnumSet<ResultsImageType> d = EnumSet.allOf(ResultsImageType.class);
      d.remove(ResultsImageType.UNRECOGNIZED);
      VALUES = d.toArray(new ResultsImageType[0]);
      NAMES = new String[VALUES.length];
      for (int i = 0; i < VALUES.length; i++) {
        NAMES[i] = ResultsProtosHelper.getName(VALUES[i]);
      }
    }
  }

  /**
   * Gets the results image type values.
   *
   * @return the results image type values
   */
  public static ResultsImageType[] getResultsImageTypeValues() {
    return ResultsImageTypeLoader.VALUES.clone();
  }

  /**
   * Gets the results image type names.
   *
   * @return the results image type names
   */
  public static String[] getResultsImageTypeNames() {
    return ResultsImageTypeLoader.NAMES.clone();
  }

  /**
   * Lazy loader for the {@link ResultsImageSizeMode} enum.
   */
  private static class ResultsImageSizeModeLoader {
    /** The enum VALUES. */
    static final ResultsImageSizeMode[] VALUES;
    /** The enum NAMES. */
    static final String[] NAMES;

    static {
      final EnumSet<ResultsImageSizeMode> d = EnumSet.allOf(ResultsImageSizeMode.class);
      d.remove(ResultsImageSizeMode.UNRECOGNIZED);
      VALUES = d.toArray(new ResultsImageSizeMode[0]);
      NAMES = new String[VALUES.length];
      for (int i = 0; i < VALUES.length; i++) {
        NAMES[i] = ResultsProtosHelper.getName(VALUES[i]);
      }
    }
  }

  /**
   * Gets the results image type values.
   *
   * @return the results image type values
   */
  public static ResultsImageSizeMode[] getResultsImageSizeModeValues() {
    return ResultsImageSizeModeLoader.VALUES.clone();
  }

  /**
   * Gets the results image type names.
   *
   * @return the results image type names
   */
  public static String[] getResultsImageSizeModeNames() {
    return ResultsImageSizeModeLoader.NAMES.clone();
  }

  /**
   * Lazy loader for the {@link ResultsFileFormat} enum.
   */
  private static class ResultsFileFormatLoader {
    /** The enum VALUES. */
    static final ResultsFileFormat[] VALUES;
    /** The enum NAMES. */
    static final String[] NAMES;

    static {
      final EnumSet<ResultsFileFormat> d = EnumSet.allOf(ResultsFileFormat.class);
      d.remove(ResultsFileFormat.UNRECOGNIZED);
      VALUES = d.toArray(new ResultsFileFormat[0]);
      NAMES = new String[VALUES.length];
      for (int i = 0; i < VALUES.length; i++) {
        NAMES[i] = ResultsProtosHelper.getName(VALUES[i]);
      }
    }
  }

  /**
   * Gets the results file format values.
   *
   * @return the results file format values
   */
  public static ResultsFileFormat[] getResultsFileFormatValues() {
    return ResultsFileFormatLoader.VALUES.clone();
  }

  /**
   * Gets the results file format names.
   *
   * @return the results file format names
   */
  public static String[] getResultsFileFormatNames() {
    return ResultsFileFormatLoader.NAMES.clone();
  }

  /**
   * Lazy loader for the {@link ResultsTableFormat} enum.
   */
  private static class ResultsTableFormatLoader {
    /** The enum VALUES. */
    static final ResultsTableFormat[] VALUES;
    /** The enum NAMES. */
    static final String[] NAMES;

    static {
      final EnumSet<ResultsTableFormat> d = EnumSet.allOf(ResultsTableFormat.class);
      d.remove(ResultsTableFormat.UNRECOGNIZED);
      VALUES = d.toArray(new ResultsTableFormat[0]);
      NAMES = new String[VALUES.length];
      for (int i = 0; i < VALUES.length; i++) {
        NAMES[i] = ResultsProtosHelper.getName(VALUES[i]);
      }
    }
  }

  /**
   * Gets the results table format values.
   *
   * @return the results table format values
   */
  public static ResultsTableFormat[] getResultsTableFormatValues() {
    return ResultsTableFormatLoader.VALUES.clone();
  }

  /**
   * Gets the results table format names.
   *
   * @return the results table format names
   */
  public static String[] getResultsTableFormatNames() {
    return ResultsTableFormatLoader.NAMES.clone();
  }

  /**
   * Lazy loader for the {@link DataFilterType} enum.
   */
  private static class DataFilterTypeLoader {
    /** The enum VALUES. */
    static final DataFilterType[] VALUES;
    /** The enum NAMES. */
    static final String[] NAMES;

    static {
      final EnumSet<DataFilterType> d = EnumSet.allOf(DataFilterType.class);
      d.remove(DataFilterType.UNRECOGNIZED);
      VALUES = d.toArray(new DataFilterType[0]);
      NAMES = new String[VALUES.length];
      for (int i = 0; i < VALUES.length; i++) {
        NAMES[i] = FitProtosHelper.getName(VALUES[i]);
      }
    }
  }

  /**
   * Gets the data filter type values.
   *
   * @return the data filter type values
   */
  public static DataFilterType[] getDataFilterTypeValues() {
    return DataFilterTypeLoader.VALUES.clone();
  }

  /**
   * Gets the data filter type names.
   *
   * @return the data filter type names
   */
  public static String[] getDataFilterTypeNames() {
    return DataFilterTypeLoader.NAMES.clone();
  }

  /**
   * Lazy loader for the {@link DataFilterMethod} enum.
   */
  private static class DataFilterMethodLoader {
    /** The enum VALUES. */
    static final DataFilterMethod[] VALUES;
    /** The enum NAMES. */
    static final String[] NAMES;

    static {
      final EnumSet<DataFilterMethod> d = EnumSet.allOf(DataFilterMethod.class);
      d.remove(DataFilterMethod.UNRECOGNIZED);
      VALUES = d.toArray(new DataFilterMethod[0]);
      NAMES = new String[VALUES.length];
      for (int i = 0; i < VALUES.length; i++) {
        NAMES[i] = FitProtosHelper.getName(VALUES[i]);
      }
    }
  }

  /**
   * Gets the data filter method values.
   *
   * @return the data filter method values
   */
  public static DataFilterMethod[] getDataFilterMethodValues() {
    return DataFilterMethodLoader.VALUES.clone();
  }

  /**
   * Gets the data filter method names.
   *
   * @return the data filter method names
   */
  public static String[] getDataFilterMethodNames() {
    return DataFilterMethodLoader.NAMES.clone();
  }

  /**
   * Lazy loader for the {@link FitSolver} enum.
   */
  private static class FitSolverLoader {
    /** The enum VALUES. */
    static final FitSolver[] VALUES;
    /** The enum NAMES. */
    static final String[] NAMES;

    static {
      final EnumSet<FitSolver> d = EnumSet.allOf(FitSolver.class);
      d.remove(FitSolver.UNRECOGNIZED);
      d.remove(FitSolver.BACKTRACKING_FAST_MLE);
      VALUES = d.toArray(new FitSolver[0]);
      NAMES = new String[VALUES.length];
      for (int i = 0; i < VALUES.length; i++) {
        NAMES[i] = FitProtosHelper.getName(VALUES[i]);
      }
    }
  }

  /**
   * Gets the fit solver values.
   *
   * @return the fit solver values
   */
  public static FitSolver[] getFitSolverValues() {
    return FitSolverLoader.VALUES.clone();
  }

  /**
   * Gets the fit solver names.
   *
   * @return the fit solver names
   */
  public static String[] getFitSolverNames() {
    return FitSolverLoader.NAMES.clone();
  }

  /**
   * Lazy loader for the {@link SearchMethod} enum.
   */
  private static class SearchMethodLoader {
    /** The enum VALUES. */
    static final SearchMethod[] VALUES;
    /** The enum NAMES. */
    static final String[] NAMES;

    static {
      final EnumSet<SearchMethod> d = EnumSet.allOf(SearchMethod.class);
      d.remove(SearchMethod.UNRECOGNIZED);
      d.remove(SearchMethod.BFGS);
      VALUES = d.toArray(new SearchMethod[0]);
      NAMES = new String[VALUES.length];
      for (int i = 0; i < VALUES.length; i++) {
        NAMES[i] = FitProtosHelper.getName(VALUES[i]);
      }
    }
  }

  /**
   * Gets the fit solver values.
   *
   * @return the fit solver values
   */
  public static SearchMethod[] getSearchMethodValues() {
    return SearchMethodLoader.VALUES.clone();
  }

  /**
   * Gets the fit solver names.
   *
   * @return the fit solver names
   */
  public static String[] getSearchMethodNames() {
    return SearchMethodLoader.NAMES.clone();
  }

  /**
   * Lazy loader for the {@link NoiseEstimatorMethod} enum.
   */
  private static class NoiseEstimatorMethodLoader {
    /** The enum VALUES. */
    static final NoiseEstimatorMethod[] VALUES;
    /** The enum NAMES. */
    static final String[] NAMES;

    static {
      final EnumSet<NoiseEstimatorMethod> d = EnumSet.allOf(NoiseEstimatorMethod.class);
      d.remove(NoiseEstimatorMethod.UNRECOGNIZED);
      VALUES = d.toArray(new NoiseEstimatorMethod[0]);
      NAMES = new String[VALUES.length];
      for (int i = 0; i < VALUES.length; i++) {
        NAMES[i] = FitProtosHelper.getName(VALUES[i]);
      }
    }
  }

  /**
   * Gets the noise estimator method values.
   *
   * @return the noise estimator method values
   */
  public static NoiseEstimatorMethod[] getNoiseEstimatorMethodValues() {
    return NoiseEstimatorMethodLoader.VALUES.clone();
  }

  /**
   * Gets the noise estimator method names.
   *
   * @return the noise estimator method names
   */
  public static String[] getNoiseEstimatorMethodNames() {
    return NoiseEstimatorMethodLoader.NAMES.clone();
  }

  /**
   * Lazy loader for the {@link CameraType} enum.
   */
  private static class CameraTypeLoader {
    /** The enum VALUES. */
    static final CameraType[] VALUES;
    /** The enum NAMES. */
    static final String[] NAMES;

    static {
      final EnumSet<CameraType> d = EnumSet.allOf(CameraType.class);
      d.remove(CameraType.UNRECOGNIZED);
      VALUES = d.toArray(new CameraType[0]);
      NAMES = new String[VALUES.length];
      for (int i = 0; i < VALUES.length; i++) {
        NAMES[i] = CalibrationProtosHelper.getName(VALUES[i]);
      }
    }
  }

  /**
   * Gets the camera type values.
   *
   * @return the camera type values
   */
  public static CameraType[] getCameraTypeValues() {
    return CameraTypeLoader.VALUES.clone();
  }

  /**
   * Gets the camera type names.
   *
   * @return the camera type names
   */
  public static String[] getCameraTypeNames() {
    return CameraTypeLoader.NAMES.clone();
  }

  /**
   * Lazy loader for the {@link PrecisionMethod} enum.
   */
  private static class PrecisionMethodLoader {
    /** The enum VALUES. */
    static final PrecisionMethod[] VALUES;
    /** The enum NAMES. */
    static final String[] NAMES;

    static {
      final EnumSet<PrecisionMethod> d = EnumSet.allOf(PrecisionMethod.class);
      d.remove(PrecisionMethod.UNRECOGNIZED);
      VALUES = d.toArray(new PrecisionMethod[0]);
      NAMES = new String[VALUES.length];
      for (int i = 0; i < VALUES.length; i++) {
        NAMES[i] = FitProtosHelper.getName(VALUES[i]);
      }
    }
  }

  /**
   * Gets the precision method values.
   *
   * @return the precision method values
   */
  public static PrecisionMethod[] getPrecisionMethodValues() {
    return PrecisionMethodLoader.VALUES.clone();
  }

  /**
   * Gets the precision method names.
   *
   * @return the precision method names
   */
  public static String[] getPrecisionMethodNames() {
    return PrecisionMethodLoader.NAMES.clone();
  }

  /**
   * Gets the settings directory.
   *
   * @return The settings directory (from the ImageJ preferences or the default under the home
   *         directory)
   */
  public static String getSettingsDirectory() {
    return settingsDirectory;
  }

  /**
   * Sets the settings directory.
   *
   * @param directory the directory
   */
  public static void setSettingsDirectory(String directory) {
    Prefs.set(Constants.settingsDirectory, directory);
    settingsDirectory = directory;
    try {
      Files.createDirectories(Paths.get(settingsDirectory));
    } catch (final IOException ex) {
      logger.log(Level.SEVERE, "Unable create settings directory", ex);
    }
  }

  /**
   * Convert a list of objects into names (e.g. pass in (Object[])enum.getValues()). The first
   * letter is capitalised. The rest of the name is converted to lowercase if it is all uppercase.
   * Remaining mixed case names are left unchanged.
   *
   * <p>Used to convert the settings enumerations into names used with dialogs.
   *
   * @param objects the objects
   * @return the names
   */
  public static String[] getNames(Object... objects) {
    final String[] names = new String[objects.length];
    for (int i = 0; i < names.length; i++) {
      String name;
      if (objects[i] instanceof NamedObject) {
        final NamedObject o = (NamedObject) objects[i];
        name = getName(o.getName(), o.getShortName());
      } else {
        name = objects[i].toString();
      }

      // Capitalise first letter
      if (name.length() > 0 && Character.isLowerCase(name.charAt(0))) {
        name = Character.toUpperCase(name.charAt(0)) + name.substring(1);
      }
      names[i] = name;
    }
    return names;
  }

  /**
   * Gets the name by appending the short name after the name in brackets, only if the short name is
   * different to the name.
   *
   * @param name the name
   * @param shortName the short name
   * @return the name
   */
  public static String getName(String name, String shortName) {
    if (!name.equals(shortName)) {
      name += " (" + shortName + ")";
    }
    return name;
  }

  /**
   * Gets the settings filename (from the ImageJ preferences or the default in the home directory).
   *
   * @return The settings filename
   */
  public static String getSettingsFilename() {
    return Prefs.get(Constants.settingsFilename, DEFAULT_FILENAME);
  }

  /**
   * Save settings filename.
   *
   * @param filename the filename
   */
  public static void saveSettingsFilename(String filename) {
    Prefs.set(Constants.settingsFilename, filename);
  }

  /**
   * Creates the settings file path using the class name to create the filename in the settings
   * directory.
   *
   * @param clazz the clazz
   * @return the path
   */
  private static Path createSettingsPath(Class<?> clazz) {
    return Paths.get(settingsDirectory, clazz.getSimpleName().toLowerCase(Locale.US) + ".settings");
  }

  /**
   * Write a message to a settings file in the settings directory.
   *
   * @param message the message
   * @return true, if successful
   */
  public static boolean writeSettings(Message message) {
    return writeSettings(message, 0);
  }

  /**
   * Write a message to a settings file in the settings directory.
   *
   * @param message the message
   * @param flags the flags
   * @return true, if successful
   */
  public static boolean writeSettings(Message message, int flags) {
    return writeMessage(message, createSettingsPath(message.getClass()), flags);
  }

  /**
   * Write a message to a settings file in the settings directory.
   *
   * @param message the message
   * @return true, if successful
   */
  public static boolean writeSettings(Message.Builder message) {
    return writeSettings(message.build());
  }

  /**
   * Write a message to a settings file in the settings directory.
   *
   * @param message the message
   * @param flags the flags
   * @return true, if successful
   */
  public static boolean writeSettings(Message.Builder message, int flags) {
    return writeSettings(message.build(), flags);
  }

  /**
   * Write to a settings file in the settings directory.
   *
   * @param fitEngineConfiguration the fit engine configuration
   * @param flags the flags
   * @return true, if successful
   */
  public static boolean writeSettings(FitEngineConfiguration fitEngineConfiguration, int flags) {
    final FitConfiguration fitConfig = fitEngineConfiguration.getFitConfiguration();
    // This is fail fast
    boolean result = writeSettings(fitEngineConfiguration.getFitEngineSettings(), flags);
    result &= writeSettings(fitConfig.getCalibration(), flags);
    result &= writeSettings(fitConfig.getPsf(), flags);
    return result;
  }

  /**
   * Clear the settings file for the given class.
   *
   * @param clazz the class
   * @return true, if the settings are cleared
   */
  public static boolean clearSettings(Class<?> clazz) {
    final Path path = createSettingsPath(clazz);
    try {
      if (Files.exists(path)) {
        Files.delete(path);
      }
      return true;
    } catch (final IOException ex) {
      logger.log(Level.WARNING, ex, () -> "Unable to clear the settings file: " + path);
    }
    return false;
  }

  /**
   * Simple class to allow generic message reading.
   *
   * @param <T> the generic message type
   */
  public static class ConfigurationReader<T extends Message> {
    /** the default instance of the message type. */
    private final T mt;

    /**
     * Instantiates a new configuration reader.
     *
     * @param mt the default instance of the message type
     */
    public ConfigurationReader(T mt) {
      this.mt = mt;
    }

    /**
     * Read the message.
     *
     * @return the message
     */
    public T read() {
      return read(0);
    }

    /**
     * Read the message.
     *
     * @param flags the flags
     * @return the message
     */
    @SuppressWarnings("unchecked")
    public T read(int flags) {
      T msg = (T) readMessage(mt.getParserForType(), createSettingsPath(mt.getClass()), flags);
      // If null use the default unless explicitly not requested
      if (msg == null && BitFlagUtils.anyNotSet(flags, FLAG_NO_DEFAULT)) {
        msg = mt;
      }
      return msg;
    }
  }

  /**
   * Read the Calibration from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the Calibration
   */
  public static Calibration readCalibration(int flags) {
    return new ConfigurationReader<>(CalibrationProtosHelper.defaultCalibration).read(flags);
  }

  /**
   * Read the PSF from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the PSF
   */
  public static PSF readPsf(int flags) {
    return new ConfigurationReader<>(PsfProtosHelper.defaultOneAxisGaussian2DPSF).read(flags);
  }

  /**
   * Read the ResultsSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the ResultsSettings
   */
  public static ResultsSettings readResultsSettings(int flags) {
    return new ConfigurationReader<>(ResultsProtosHelper.defaultResultsSettings).read(flags);
  }

  /**
   * Read the FitEngineSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the FitEngineSettings
   */
  private static FitEngineSettings readFitEngineSettings(int flags) {
    return new ConfigurationReader<>(FitProtosHelper.defaultFitEngineSettings).read(flags);
  }

  /**
   * Read the GUIFilterSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the GUIFilterSettings
   */
  public static GUIFilterSettings readGuiFilterSettings(int flags) {
    return new ConfigurationReader<>(GuiSettings.DefaultGUIFilterSettings.INSTANCE).read(flags);
  }

  /**
   * Read the PSFCalculatorSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the PSFCalculatorSettings
   */
  public static PSFCalculatorSettings readPsfCalculatorSettings(int flags) {
    return new ConfigurationReader<>(GuiSettings.DefaultPSFCalculatorSettings.INSTANCE)
        .read(flags);
  }

  /**
   * Read the PSFEstimatorSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the PSFEstimatorSettings
   */
  public static PSFEstimatorSettings readPsfEstimatorSettings(int flags) {
    return new ConfigurationReader<>(GuiSettings.DefaultPSFEstimatorSettings.INSTANCE)
        .read(flags);
  }

  /**
   * Read the CreateDataSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the CreateDataSettings
   */
  public static CreateDataSettings readCreateDataSettings(int flags) {
    return new ConfigurationReader<>(GuiSettings.DefaultCreateDataSettings.INSTANCE)
        .read(flags);
  }

  /**
   * Read the LoadLocalisationsSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the LoadLocalisationsSettings
   */
  public static LoadLocalisationsSettings readLoadLocalisationsSettings(int flags) {
    return new ConfigurationReader<>(GuiSettings.DefaultLoadLocalisationsSettings.INSTANCE)
        .read(flags);
  }

  /**
   * Read the ClusteringSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the ClusteringSettings
   */
  public static ClusteringSettings readClusteringSettings(int flags) {
    return new ConfigurationReader<>(GuiSettings.DefaultClusteringSettings.INSTANCE)
        .read(flags);
  }

  /**
   * Read the OpticsSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the OpticsSettings
   */
  public static OpticsSettings readOpticsSettings(int flags) {
    return new ConfigurationReader<>(GuiSettings.DefaultOpticsSettings.INSTANCE).read(flags);
  }

  /**
   * Read the ConfigurationTemplateSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the ConfigurationTemplateSettings
   */
  public static ConfigurationTemplateSettings readConfigurationTemplateSettings(int flags) {
    return new ConfigurationReader<>(GuiSettings.DefaultConfigurationTemplateSettings.INSTANCE)
        .read(flags);
  }

  /**
   * Read the DefaultTemplateSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the DefaultTemplateSettings
   */
  public static DefaultTemplateSettings readDefaultTemplateSettings(int flags) {
    return new ConfigurationReader<>(DefaultTemplateSettings.getDefaultInstance()).read(flags);
  }

  /**
   * Read the FitEngineConfiguration from the settings file in the settings directory. This loads
   * the current Calibration, PSF and FitEngineSettings.
   *
   * @param flags the flags
   * @return the FitEngineConfiguration
   */
  public static FitEngineConfiguration readFitEngineConfiguration(int flags) {
    final FitEngineSettings fitEngineSettings = readFitEngineSettings(flags);
    final Calibration calibration = readCalibration(flags);
    final PSF psf = readPsf(flags);
    return FitEngineConfiguration.create(fitEngineSettings, calibration, psf);
  }

  /**
   * Read the CameraModelSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the CameraModelSettings
   */
  public static CameraModelSettings readCameraModelSettings(int flags) {
    return new ConfigurationReader<>(CameraModelSettings.getDefaultInstance()).read(flags);
  }

  /**
   * Read the CubicSplineSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the CubicSplineSettings
   */
  public static CubicSplineSettings readCubicSplineSettings(int flags) {
    return new ConfigurationReader<>(CubicSplineSettings.getDefaultInstance()).read(flags);
  }

  /**
   * Read the NucleusMaskSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the NucleusMaskSettings
   */
  public static NucleusMaskSettings readNucleusMaskSettings(int flags) {
    return new ConfigurationReader<>(GuiSettings.DefaultNucleusMaskSettings.INSTANCE)
        .read(flags);
  }

  /**
   * Read the PSFCreatorSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the PSFCreatorSettings
   */
  public static PSFCreatorSettings readPsfCreatorSettings(int flags) {
    return new ConfigurationReader<>(GuiSettings.DefaultPSFCreatorSettings.INSTANCE)
        .read(flags);
  }

  /**
   * Read the CameraModelManagerSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the CameraModelManagerSettings
   */
  public static CameraModelManagerSettings readCameraModelManagerSettings(int flags) {
    return new ConfigurationReader<>(GuiSettings.DefaultCameraModelManagerSettings.INSTANCE)
        .read(flags);
  }

  /**
   * Read the CameraModelAnalysisSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the CameraModelAnalysisSettings
   */
  public static CameraModelAnalysisSettings readCameraModelAnalysisSettings(int flags) {
    return new ConfigurationReader<>(GuiSettings.DefaultCameraModelAnalysisSettings.INSTANCE)
        .read(flags);
  }

  /**
   * Read the CameraModelFisherInformationAnalysisSettings from the settings file in the settings
   * directory.
   *
   * @param flags the flags
   * @return the CameraModelFisherInformationAnalysisSettings
   */
  public static CameraModelFisherInformationAnalysisSettings
      readCameraModelFisherInformationAnalysisSettings(int flags) {
    return new ConfigurationReader<>(
        GuiSettings.DefaultCameraModelFisherInformationAnalysisSettings.INSTANCE).read(flags);
  }

  /**
   * Read the CubicSplineManagerSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the CubicSplineManagerSettings
   */
  public static CubicSplineManagerSettings readCubicSplineManagerSettings(int flags) {
    return new ConfigurationReader<>(GuiSettings.DefaultCubicSplineManagerSettings.INSTANCE)
        .read(flags);
  }

  /**
   * Read the FailCountManagerSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the FailCountManagerSettings
   */
  public static FailCountManagerSettings readFailCountManagerSettings(int flags) {
    return new ConfigurationReader<>(GuiSettings.DefaultFailCountManagerSettings.INSTANCE)
        .read(flags);
  }

  /**
   * Read the AstigmatismModelSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the AstigmatismModelSettings
   */
  public static AstigmatismModelSettings readAstigmatismModelSettings(int flags) {
    return new ConfigurationReader<>(AstigmatismModelSettings.getDefaultInstance()).read(flags);
  }

  /**
   * Read the AstigmatismModelManagerSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the AstigmatismModelManagerSettings
   */
  public static AstigmatismModelManagerSettings readAstigmatismModelManagerSettings(int flags) {
    return new ConfigurationReader<>(
        GuiSettings.DefaultAstigmatismModelManagerSettings.INSTANCE).read(flags);
  }

  /**
   * Read the CropResultsSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the CropResultsSettings
   */
  public static CropResultsSettings readCropResultsSettings(int flags) {
    return new ConfigurationReader<>(CropResultsSettings.getDefaultInstance()).read(flags);
  }

  /**
   * Read the SummariseResultsSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the SummariseResultsSettings
   */
  public static SummariseResultsSettings readSummariseResultsSettings(int flags) {
    return new ConfigurationReader<>(SummariseResultsSettings.getDefaultInstance()).read(flags);
  }

  /**
   * Read the ImageJ3DResultsViewerSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the ImageJ3DResultsViewerSettings
   */
  public static ImageJ3DResultsViewerSettings readImageJ3DResultsViewerSettings(int flags) {
    return new ConfigurationReader<>(GuiSettings.DefaultImageJ3DResultsViewerSettings.INSTANCE)
        .read(flags);
  }

  /**
   * Read the SpotFitSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the SpotFitSettings
   */
  public static SpotFitSettings readSpotFitSettings(int flags) {
    return new ConfigurationReader<>(GuiSettings.DefaultSpotFitSettings.INSTANCE).read(flags);
  }

  /**
   * Read the TranslateResultsSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the TranslateResultsSettings
   */
  public static TranslateResultsSettings readTranslateResultsSettings(int flags) {
    return new ConfigurationReader<>(TranslateResultsSettings.getDefaultInstance()).read(flags);
  }

  /**
   * Read the TcPalmAnalysisSettings from the settings file in the settings directory.
   *
   * @param flags the flags
   * @return the TcPalmAnalysisSettings
   */
  public static TcPalmAnalysisSettings readTcPalmAnalysisSettings(int flags) {
    return new ConfigurationReader<>(GuiSettings.DefaultTcPalmAnalysisSettings.INSTANCE)
        .read(flags);
  }

  /**
   * Write the message to file.
   *
   * <p>If this fails then an error message is written to the ImageJ log
   *
   * @param message the message
   * @param filename the filename
   * @param flags the flags
   * @return True if written
   */
  public static boolean writeMessage(Message message, String filename, int flags) {
    return writeMessage(message, Paths.get(filename), flags);
  }

  /**
   * Write the message to file.
   * 
   * <p>If this fails then an error message is written to the ImageJ log
   *
   * @param message the message
   * @param path the path
   * @param flags the flags
   * @return True if written
   */
  public static boolean writeMessage(Message message, Path path, int flags) {
    try (OutputStream fs = new BufferedOutputStream(Files.newOutputStream(path))) {
      return writeMessage(message, fs, flags);
    } catch (final IOException ex) {
      logWriteError(flags, ex);
    }
    return false;
  }

  /**
   * Write the message.
   *
   * <p>If this fails then an error message is written to the ImageJ log
   *
   * @param message the message
   * @param output the output
   * @param flags the flags
   * @return True if saved
   */
  public static boolean writeMessage(Message message, OutputStream output, int flags) {
    try {
      message.writeDelimitedTo(output);
      return true;
    } catch (final IOException ex) {
      logWriteError(flags, ex);
    }
    return false;
  }

  /**
   * Read the message from file.
   *
   * <p>If this fails then an error message is written to the ImageJ log
   *
   * @param parser the parser
   * @param filename the filename
   * @param flags the flags
   * @return the message
   */
  public static Message readMessage(Parser<? extends Message> parser, String filename, int flags) {
    return readMessage(parser, Paths.get(filename), flags);
  }

  /**
   * Read the message from file.
   * 
   * <p>If this fails then an error message is written to the ImageJ log
   *
   * @param parser the parser
   * @param path the path
   * @param flags the flags
   * @return the message
   */
  public static Message readMessage(Parser<? extends Message> parser, Path path, int flags) {
    try (InputStream is = new BufferedInputStream(Files.newInputStream(path))) {
      return readMessage(parser, is, flags);
    } catch (final IOException ex) {
      // Only print this if the file-not-found flag is present
      // and not silent. This prevents warnings when settings files
      // have yet to be created, i.e. for new users of a settings file.
      if (BitFlagUtils.areSet(flags, FLAG_SHOW_FILE_NOT_FOUND_ON_READ)) {
        logReadError(flags, ex);
      }
    }
    return null;
  }

  /**
   * Read the message from file.
   *
   * <p>If this fails then an error message is written to the ImageJ log
   *
   * @param parser the parser
   * @param input the input
   * @param flags the flags
   * @return the message
   */
  public static Message readMessage(Parser<? extends Message> parser, InputStream input,
      int flags) {
    try {
      return parser.parseDelimitedFrom(input);
    } catch (final InvalidProtocolBufferException ex) {
      logReadError(flags, ex);
    }
    return null;
  }

  /**
   * Log write error.
   *
   * @param flags the flags
   * @param ex the ex
   */
  private static void logWriteError(final int flags, final Exception ex) {
    if (BitFlagUtils.anyNotSet(flags, FLAG_SILENT)) {
      logger.log(Level.SEVERE, "Unable to write message", ex);
    }
  }

  /**
   * Log read error.
   *
   * @param flags the flags
   * @param ex the ex
   */
  private static void logReadError(final int flags, final Exception ex) {
    if (BitFlagUtils.anyNotSet(flags, FLAG_SILENT)) {
      logger.log(Level.SEVERE, "Unable to read message", ex);
    }
  }

  /**
   * Write the message to a JSON string.
   *
   * @param message the message
   * @return the JSON string
   */
  public static String toJson(MessageOrBuilder message) {
    return toJson(message, 0);
  }

  /**
   * Write the message to a JSON string.
   *
   * @param message the message
   * @param flags the flags
   * @return the JSON string
   */
  public static String toJson(MessageOrBuilder message, int flags) {
    final StringBuilder sb = new StringBuilder();
    if (toJson(message, sb, flags)) {
      return sb.toString();
    }
    return null;
  }

  /**
   * Write the message to file in JSON text format.
   *
   * <p>If this fails then an error message is written to the ImageJ log
   *
   * @param message the message
   * @param filename the filename
   * @param flags the flags
   * @return True if written
   */
  public static boolean toJson(MessageOrBuilder message, String filename, int flags) {
    try (BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename))) {
      return toJson(message, writer, flags);
    } catch (final IOException ex) {
      logWriteError(flags, ex);
    }
    return false;
  }

  /**
   * Write the message to file in JSON text format.
   *
   * <p>If this fails then an error message is written to the ImageJ log
   *
   * @param message the message
   * @param file the file
   * @param flags the flags
   * @return True if written
   */
  public static boolean toJson(MessageOrBuilder message, File file, int flags) {
    try (BufferedWriter writer = Files.newBufferedWriter(file.toPath())) {
      return toJson(message, writer, flags);
    } catch (final IOException ex) {
      logWriteError(flags, ex);
    }
    return false;
  }

  /**
   * Write the message to a JSON string.
   *
   * @param message the message
   * @param output the output
   * @param flags the flags
   * @return the JSON string
   */
  public static boolean toJson(MessageOrBuilder message, Appendable output, int flags) {
    try {
      Printer localPrinter = PrinterLoader.PRINTER;
      if (BitFlagUtils.anyNotSet(flags, FLAG_JSON_WHITESPACE)) {
        localPrinter = localPrinter.omittingInsignificantWhitespace();
      }
      if (BitFlagUtils.anySet(flags, FLAG_JSON_DEFAULT_VALUES)) {
        localPrinter = localPrinter.includingDefaultValueFields();
      }
      localPrinter.appendTo(message, output);
      return true;
    } catch (final IOException ex) {
      logWriteError(flags, ex);
    }
    return false;
  }

  /**
   * Read the message from a JSON string.
   *
   * @param json the JSON string
   * @param builder the builder
   * @return true, if successful
   */
  public static boolean fromJson(String json, Message.Builder builder) {
    return fromJson(json, builder, 0);
  }

  /**
   * Read the message from a JSON string.
   *
   * @param json the JSON string
   * @param builder the builder
   * @param flags the flags
   * @return true, if successful
   */
  public static boolean fromJson(String json, Message.Builder builder, int flags) {
    return fromJson(new StringReader(json), builder, flags);
  }

  /**
   * Read the message from a JSON string in a file.
   *
   * @param filename the filename
   * @param builder the builder
   * @param flags the flags
   * @return true, if successful
   */
  public static boolean fromJson(Path filename, Message.Builder builder, int flags) {
    try (Reader reader = Files.newBufferedReader(filename)) {
      return fromJson(reader, builder, flags);
    } catch (final IOException ex) {
      logReadError(flags, ex);
    }
    return false;
  }

  /**
   * Read the message from a JSON string in a file.
   *
   * @param file the file
   * @param builder the builder
   * @param flags the flags
   * @return true, if successful
   */
  public static boolean fromJson(File file, Message.Builder builder, int flags) {
    try (Reader reader = Files.newBufferedReader(file.toPath())) {
      return fromJson(reader, builder, flags);
    } catch (final IOException ex) {
      logReadError(flags, ex);
    }
    return false;
  }

  /**
   * Read the message from a JSON string.
   *
   * @param reader the reader
   * @param builder the builder
   * @param flags the flags
   * @return true, if successful
   */
  public static boolean fromJson(Reader reader, Message.Builder builder, int flags) {
    try {
      ParserLoader.PARSER.merge(reader, builder);
      return true;
    } catch (final IOException ex) {
      logReadError(flags, ex);
    }
    return false;
  }
}
