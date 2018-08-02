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
package uk.ac.sussex.gdsc.smlm.ij.settings;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.Reader;
import java.io.StringReader;
import java.util.EnumSet;

import com.google.protobuf.InvalidProtocolBufferException;
import com.google.protobuf.Message;
import com.google.protobuf.MessageOrBuilder;
import com.google.protobuf.Parser;
import com.google.protobuf.util.JsonFormat;
import com.google.protobuf.util.JsonFormat.Printer;

import ij.IJ;
import ij.Prefs;
import uk.ac.sussex.gdsc.core.utils.BitFlags;
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
import uk.ac.sussex.gdsc.smlm.data.config.FitProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.AstigmatismModelManagerSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.CameraModelAnalysisSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.CameraModelFisherInformationAnalysisSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.CameraModelManagerSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.ClusteringSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.ConfigurationTemplateSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.CreateDataSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.CropResultsSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.CubicSplineManagerSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.DefaultTemplateSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.FailCountManagerSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.GUIFilterSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.ImageJ3DResultsViewerSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.LoadLocalisationsSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.NucleusMaskSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.OpticsSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.PSFCalculatorSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.PSFCreatorSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.PSFEstimatorSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.SpotFitSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.SummariseResultsSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.TranslateResultsSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.AstigmatismModelSettings;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.CubicSplineSettings;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsFileFormat;
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

/**
 * Manage the settings for ImageJ plugins
 */
public class SettingsManager
{
    private static final String DEFAULT_FILENAME = System.getProperty("user.home") +
            System.getProperty("file.separator") + "gdsc.smlm.settings.xml";
    private static final String DEFAULT_DIRECTORY = System.getProperty("user.home") +
            System.getProperty("file.separator") + ".gdsc.smlm";

    /** Use this to suppress warnings. */
    public static final int FLAG_SILENT = 0x00000001;
    /** Use this flag to suppress returning a default instance. */
    public static final int FLAG_NO_DEFAULT = 0x00000002;
    /** Use this flag to include insignificant whitespace in JSON output. */
    public static final int FLAG_JSON_WHITESPACE = 0x00000004;
    /** Use this flag to include default values in JSON output. */
    public static final int FLAG_JSON_DEFAULT_VALUES = 0x00000008;
    /**
     * Use this to show file-not-found warning when reading settings.
     * The default is to suppress the warning as it is usually usually
     * from a settings file that has not been created.
     */
    public static final int FLAG_SHOW_FILE_NOT_FOUND_ON_READ = 0x00000010;

    /** The settings directory. */
    private static File settingsDirectory;
    static
    {
        setSettingsDirectory(Prefs.get(Constants.settingsDirectory, DEFAULT_DIRECTORY));
    }

    /**
     * Gets the settings directory.
     *
     * @return The settings directory (from the ImageJ preferences or the default under the home directory)
     */
    public static String getSettingsDirectory()
    {
        return settingsDirectory.getPath();
    }

    /**
     * Save the settings directory.
     *
     * @param directory
     *            the directory
     */
    public static void setSettingsDirectory(String directory)
    {
        Prefs.set(Constants.settingsDirectory, directory);
        settingsDirectory = new File(directory);
        try
        {
            if (!settingsDirectory.exists())
                settingsDirectory.mkdirs();
        }
        catch (final Exception e)
        {
            IJ.log("Unable create settings directory: " + e.getMessage());
        }
    }

    // Encapsulate the values and names of enums for lazy loading

    private static DistanceUnit[] _DistanceUnitValues;

    /**
     * Gets the distance unit values.
     *
     * @return the distance unit values
     */
    public static DistanceUnit[] getDistanceUnitValues()
    {
        if (_DistanceUnitValues == null)
            initDistanceUnit();
        return _DistanceUnitValues;
    }

    private static String[] _DistanceUnitNames;

    /**
     * Gets the distance unit names.
     *
     * @return the distance unit names
     */
    public static String[] getDistanceUnitNames()
    {
        if (_DistanceUnitNames == null)
            initDistanceUnit();
        return _DistanceUnitNames;
    }

    private static void initDistanceUnit()
    {
        final EnumSet<DistanceUnit> d = EnumSet.allOf(DistanceUnit.class);
        d.remove(DistanceUnit.UNRECOGNIZED);
        //d.remove(DistanceUnit.DISTANCE_UNIT_NA);
        _DistanceUnitValues = d.toArray(new DistanceUnit[d.size()]);
        _DistanceUnitNames = new String[_DistanceUnitValues.length];
        for (int i = 0; i < _DistanceUnitValues.length; i++)
            _DistanceUnitNames[i] = getName(UnitHelper.getName(_DistanceUnitValues[i]),
                    UnitHelper.getShortName(_DistanceUnitValues[i]));
    }

    private static IntensityUnit[] _IntensityUnitValues;

    /**
     * Gets the intensity unit values.
     *
     * @return the intensity unit values
     */
    public static IntensityUnit[] getIntensityUnitValues()
    {
        if (_IntensityUnitValues == null)
            initIntensityUnit();
        return _IntensityUnitValues;
    }

    private static String[] _IntensityUnitNames;

    /**
     * Gets the intensity unit names.
     *
     * @return the intensity unit names
     */
    public static String[] getIntensityUnitNames()
    {
        if (_IntensityUnitNames == null)
            initIntensityUnit();
        return _IntensityUnitNames;
    }

    private static void initIntensityUnit()
    {
        final EnumSet<IntensityUnit> d = EnumSet.allOf(IntensityUnit.class);
        d.remove(IntensityUnit.UNRECOGNIZED);
        //d.remove(IntensityUnit.INTENSITY_UNIT_NA);
        _IntensityUnitValues = d.toArray(new IntensityUnit[d.size()]);
        _IntensityUnitNames = new String[_IntensityUnitValues.length];
        for (int i = 0; i < _IntensityUnitValues.length; i++)
            _IntensityUnitNames[i] = getName(UnitHelper.getName(_IntensityUnitValues[i]),
                    UnitHelper.getShortName(_IntensityUnitValues[i]));
    }

    private static AngleUnit[] _AngleUnitValues;

    /**
     * Gets the angle unit values.
     *
     * @return the angle unit values
     */
    public static AngleUnit[] getAngleUnitValues()
    {
        if (_AngleUnitValues == null)
            initAngleUnit();
        return _AngleUnitValues;
    }

    private static String[] _AngleUnitNames;

    /**
     * Gets the angle unit names.
     *
     * @return the angle unit names
     */
    public static String[] getAngleUnitNames()
    {
        if (_AngleUnitNames == null)
            initAngleUnit();
        return _AngleUnitNames;
    }

    private static void initAngleUnit()
    {
        final EnumSet<AngleUnit> d = EnumSet.allOf(AngleUnit.class);
        d.remove(AngleUnit.UNRECOGNIZED);
        //d.remove(AngleUnit.ANGLE_UNIT_NA);
        _AngleUnitValues = d.toArray(new AngleUnit[d.size()]);
        _AngleUnitNames = new String[_AngleUnitValues.length];
        for (int i = 0; i < _AngleUnitValues.length; i++)
            _AngleUnitNames[i] = getName(UnitHelper.getName(_AngleUnitValues[i]),
                    UnitHelper.getShortName(_AngleUnitValues[i]));
    }

    private static TimeUnit[] _TimeUnitValues;

    /**
     * Gets the time unit values.
     *
     * @return the time unit values
     */
    public static TimeUnit[] getTimeUnitValues()
    {
        if (_TimeUnitValues == null)
            initTimeUnit();
        return _TimeUnitValues;
    }

    private static String[] _TimeUnitNames;

    /**
     * Gets the time unit names.
     *
     * @return the time unit names
     */
    public static String[] getTimeUnitNames()
    {
        if (_TimeUnitNames == null)
            initTimeUnit();
        return _TimeUnitNames;
    }

    private static void initTimeUnit()
    {
        final EnumSet<TimeUnit> d = EnumSet.allOf(TimeUnit.class);
        d.remove(TimeUnit.UNRECOGNIZED);
        //d.remove(TimeUnit.TIME_UNIT_NA);
        _TimeUnitValues = d.toArray(new TimeUnit[d.size()]);
        _TimeUnitNames = new String[_TimeUnitValues.length];
        for (int i = 0; i < _TimeUnitValues.length; i++)
            _TimeUnitNames[i] = getName(UnitHelper.getName(_TimeUnitValues[i]),
                    UnitHelper.getShortName(_TimeUnitValues[i]));
    }

    private static ResultsImageType[] _ResultsImageTypeValues;

    /**
     * Gets the results image type values.
     *
     * @return the results image type values
     */
    public static ResultsImageType[] getResultsImageTypeValues()
    {
        if (_ResultsImageTypeValues == null)
            initResultsImageType();
        return _ResultsImageTypeValues;
    }

    private static String[] _ResultsImageTypeNames;

    /**
     * Gets the results image type names.
     *
     * @return the results image type names
     */
    public static String[] getResultsImageTypeNames()
    {
        if (_ResultsImageTypeNames == null)
            initResultsImageType();
        return _ResultsImageTypeNames;
    }

    private static void initResultsImageType()
    {
        final EnumSet<ResultsImageType> d = EnumSet.allOf(ResultsImageType.class);
        d.remove(ResultsImageType.UNRECOGNIZED);
        _ResultsImageTypeValues = d.toArray(new ResultsImageType[d.size()]);
        _ResultsImageTypeNames = new String[_ResultsImageTypeValues.length];
        for (int i = 0; i < _ResultsImageTypeValues.length; i++)
            _ResultsImageTypeNames[i] = ResultsProtosHelper.getName(_ResultsImageTypeValues[i]);
    }

    private static ResultsFileFormat[] _ResultsFileFormatValues;

    /**
     * Gets the results file format values.
     *
     * @return the results file format values
     */
    public static ResultsFileFormat[] getResultsFileFormatValues()
    {
        if (_ResultsFileFormatValues == null)
            initResultsFileFormat();
        return _ResultsFileFormatValues;
    }

    private static String[] _ResultsFileFormatNames;

    /**
     * Gets the results file format names.
     *
     * @return the results file format names
     */
    public static String[] getResultsFileFormatNames()
    {
        if (_ResultsFileFormatNames == null)
            initResultsFileFormat();
        return _ResultsFileFormatNames;
    }

    private static void initResultsFileFormat()
    {
        final EnumSet<ResultsFileFormat> d = EnumSet.allOf(ResultsFileFormat.class);
        d.remove(ResultsFileFormat.UNRECOGNIZED);
        _ResultsFileFormatValues = d.toArray(new ResultsFileFormat[d.size()]);
        _ResultsFileFormatNames = new String[_ResultsFileFormatValues.length];
        for (int i = 0; i < _ResultsFileFormatValues.length; i++)
            _ResultsFileFormatNames[i] = ResultsProtosHelper.getName(_ResultsFileFormatValues[i]);
    }

    private static ResultsTableFormat[] _ResultsTableFormatValues;

    /**
     * Gets the results table format values.
     *
     * @return the results table format values
     */
    public static ResultsTableFormat[] getResultsTableFormatValues()
    {
        if (_ResultsTableFormatValues == null)
            initResultsTableFormat();
        return _ResultsTableFormatValues;
    }

    private static String[] _ResultsTableFormatNames;

    /**
     * Gets the results table format names.
     *
     * @return the results table format names
     */
    public static String[] getResultsTableFormatNames()
    {
        if (_ResultsTableFormatNames == null)
            initResultsTableFormat();
        return _ResultsTableFormatNames;
    }

    private static void initResultsTableFormat()
    {
        final EnumSet<ResultsTableFormat> d = EnumSet.allOf(ResultsTableFormat.class);
        d.remove(ResultsTableFormat.UNRECOGNIZED);
        _ResultsTableFormatValues = d.toArray(new ResultsTableFormat[d.size()]);
        _ResultsTableFormatNames = new String[_ResultsTableFormatValues.length];
        for (int i = 0; i < _ResultsTableFormatValues.length; i++)
            _ResultsTableFormatNames[i] = ResultsProtosHelper.getName(_ResultsTableFormatValues[i]);
    }

    private static DataFilterType[] _DataFilterTypeValues;

    /**
     * Gets the data filter type values.
     *
     * @return the data filter type values
     */
    public static DataFilterType[] getDataFilterTypeValues()
    {
        if (_DataFilterTypeValues == null)
            initDataFilterType();
        return _DataFilterTypeValues;
    }

    private static String[] _DataFilterTypeNames;

    /**
     * Gets the data filter type names.
     *
     * @return the data filter type names
     */
    public static String[] getDataFilterTypeNames()
    {
        if (_DataFilterTypeNames == null)
            initDataFilterType();
        return _DataFilterTypeNames;
    }

    private static void initDataFilterType()
    {
        final EnumSet<DataFilterType> d = EnumSet.allOf(DataFilterType.class);
        d.remove(DataFilterType.UNRECOGNIZED);
        _DataFilterTypeValues = d.toArray(new DataFilterType[d.size()]);
        _DataFilterTypeNames = new String[_DataFilterTypeValues.length];
        for (int i = 0; i < _DataFilterTypeValues.length; i++)
            _DataFilterTypeNames[i] = FitProtosHelper.getName(_DataFilterTypeValues[i]);
    }

    private static DataFilterMethod[] _DataFilterMethodValues;

    /**
     * Gets the data filter method values.
     *
     * @return the data filter method values
     */
    public static DataFilterMethod[] getDataFilterMethodValues()
    {
        if (_DataFilterMethodValues == null)
            initDataFilterMethod();
        return _DataFilterMethodValues;
    }

    private static String[] _DataFilterMethodNames;

    /**
     * Gets the data filter method names.
     *
     * @return the data filter method names
     */
    public static String[] getDataFilterMethodNames()
    {
        if (_DataFilterMethodNames == null)
            initDataFilterMethod();
        return _DataFilterMethodNames;
    }

    private static void initDataFilterMethod()
    {
        final EnumSet<DataFilterMethod> d = EnumSet.allOf(DataFilterMethod.class);
        d.remove(DataFilterMethod.UNRECOGNIZED);
        _DataFilterMethodValues = d.toArray(new DataFilterMethod[d.size()]);
        _DataFilterMethodNames = new String[_DataFilterMethodValues.length];
        for (int i = 0; i < _DataFilterMethodValues.length; i++)
            _DataFilterMethodNames[i] = FitProtosHelper.getName(_DataFilterMethodValues[i]);
    }

    private static FitSolver[] _FitSolverValues;

    /**
     * Gets the fit solver values.
     *
     * @return the fit solver values
     */
    public static FitSolver[] getFitSolverValues()
    {
        if (_FitSolverValues == null)
            initFitSolver();
        return _FitSolverValues;
    }

    private static String[] _FitSolverNames;

    /**
     * Gets the fit solver names.
     *
     * @return the fit solver names
     */
    public static String[] getFitSolverNames()
    {
        if (_FitSolverNames == null)
            initFitSolver();
        return _FitSolverNames;
    }

    private static void initFitSolver()
    {
        final EnumSet<FitSolver> d = EnumSet.allOf(FitSolver.class);
        d.remove(FitSolver.UNRECOGNIZED);
        _FitSolverValues = d.toArray(new FitSolver[d.size()]);
        _FitSolverNames = new String[_FitSolverValues.length];
        for (int i = 0; i < _FitSolverValues.length; i++)
            _FitSolverNames[i] = FitProtosHelper.getName(_FitSolverValues[i]);
    }

    private static NoiseEstimatorMethod[] _NoiseEstimatorMethodValues;

    /**
     * Gets the noise estimator method values.
     *
     * @return the noise estimator method values
     */
    public static NoiseEstimatorMethod[] getNoiseEstimatorMethodValues()
    {
        if (_NoiseEstimatorMethodValues == null)
            initNoiseEstimatorMethod();
        return _NoiseEstimatorMethodValues;
    }

    private static String[] _NoiseEstimatorMethodNames;

    /**
     * Gets the noise estimator method names.
     *
     * @return the noise estimator method names
     */
    public static String[] getNoiseEstimatorMethodNames()
    {
        if (_NoiseEstimatorMethodNames == null)
            initNoiseEstimatorMethod();
        return _NoiseEstimatorMethodNames;
    }

    private static void initNoiseEstimatorMethod()
    {
        final EnumSet<NoiseEstimatorMethod> d = EnumSet.allOf(NoiseEstimatorMethod.class);
        d.remove(NoiseEstimatorMethod.UNRECOGNIZED);
        _NoiseEstimatorMethodValues = d.toArray(new NoiseEstimatorMethod[d.size()]);
        _NoiseEstimatorMethodNames = new String[_NoiseEstimatorMethodValues.length];
        for (int i = 0; i < _NoiseEstimatorMethodValues.length; i++)
            _NoiseEstimatorMethodNames[i] = FitProtosHelper.getName(_NoiseEstimatorMethodValues[i]);
    }

    private static CameraType[] _CameraTypeValues;

    /**
     * Gets the camera type values.
     *
     * @return the camera type values
     */
    public static CameraType[] getCameraTypeValues()
    {
        if (_CameraTypeValues == null)
            initCameraType();
        return _CameraTypeValues;
    }

    private static String[] _CameraTypeNames;

    /**
     * Gets the camera type names.
     *
     * @return the camera type names
     */
    public static String[] getCameraTypeNames()
    {
        if (_CameraTypeNames == null)
            initCameraType();
        return _CameraTypeNames;
    }

    private static void initCameraType()
    {
        final EnumSet<CameraType> d = EnumSet.allOf(CameraType.class);
        d.remove(CameraType.UNRECOGNIZED);
        _CameraTypeValues = d.toArray(new CameraType[d.size()]);
        _CameraTypeNames = new String[_CameraTypeValues.length];
        for (int i = 0; i < _CameraTypeValues.length; i++)
            _CameraTypeNames[i] = CalibrationProtosHelper.getName(_CameraTypeValues[i]);
    }

    private static PrecisionMethod[] _PrecisionMethodValues;

    /**
     * Gets the precision method values.
     *
     * @return the precision method values
     */
    public static PrecisionMethod[] getPrecisionMethodValues()
    {
        if (_PrecisionMethodValues == null)
            initPrecisionMethod();
        return _PrecisionMethodValues;
    }

    private static String[] _PrecisionMethodNames;

    /**
     * Gets the precision method names.
     *
     * @return the precision method names
     */
    public static String[] getPrecisionMethodNames()
    {
        if (_PrecisionMethodNames == null)
            initPrecisionMethod();
        return _PrecisionMethodNames;
    }

    private static void initPrecisionMethod()
    {
        final EnumSet<PrecisionMethod> d = EnumSet.allOf(PrecisionMethod.class);
        d.remove(PrecisionMethod.UNRECOGNIZED);
        _PrecisionMethodValues = d.toArray(new PrecisionMethod[d.size()]);
        _PrecisionMethodNames = new String[_PrecisionMethodValues.length];
        for (int i = 0; i < _PrecisionMethodValues.length; i++)
            _PrecisionMethodNames[i] = FitProtosHelper.getName(_PrecisionMethodValues[i]);
    }

    /**
     * Convert a list of objects into names (e.g. pass in (Object[])enum.getValues()). The first letter is capitalised.
     * The rest of the name is converted to lowercase if it is all uppercase. Remaining mixed case names are left
     * unchanged.
     * <p>
     * Used to convert the settings enumerations into names used with dialogs.
     *
     * @param objects
     *            the objects
     * @return the names
     */
    public static String[] getNames(Object... objects)
    {
        final String[] names = new String[objects.length];
        for (int i = 0; i < names.length; i++)
        {
            String name;
            if (objects[i] instanceof NamedObject)
            {
                final NamedObject o = (NamedObject) objects[i];
                name = getName(o.getName(), o.getShortName());
            }
            else
                name = objects[i].toString();

            if (name.length() > 0)
                // Capitalise first letter
                if (Character.isLowerCase(name.charAt(0)))
                    name = Character.toUpperCase(name.charAt(0)) + name.substring(1);
            names[i] = name;
        }
        return names;
    }

    /**
     * Gets the name by appending the short name after the name in brackets, only if the short name is different to the
     * name.
     *
     * @param name
     *            the name
     * @param shortName
     *            the short name
     * @return the name
     */
    public static String getName(String name, String shortName)
    {
        if (!name.equals(shortName))
            name += " (" + shortName + ")";
        return name;
    }

    /**
     * @return The settings filename (from the ImageJ preferences or the default in the home directory)
     */
    public static String getSettingsFilename()
    {
        final String filename = Prefs.get(Constants.settingsFilename, DEFAULT_FILENAME);
        return filename;
    }

    /**
     * Save settings filename.
     *
     * @param filename
     *            the filename
     */
    public static void saveSettingsFilename(String filename)
    {
        Prefs.set(Constants.settingsFilename, filename);
    }

    /**
     * Creates the settings file using the class name to create the filename in the settings directory.
     *
     * @param clazz
     *            the clazz
     * @return the file
     */
    private static File createSettingsFile(Class<?> clazz)
    {
        return new File(settingsDirectory, clazz.getSimpleName().toLowerCase() + ".settings");
    }

    /**
     * Write a message to a settings file in the settings directory.
     *
     * @param message
     *            the message
     * @return true, if successful
     */
    public static boolean writeSettings(Message message)
    {
        return writeSettings(message, 0);
    }

    /**
     * Write a message to a settings file in the settings directory.
     *
     * @param message
     *            the message
     * @param flags
     *            the flags
     * @return true, if successful
     */
    public static boolean writeSettings(Message message, int flags)
    {
        return writeMessage(message, createSettingsFile(message.getClass()), flags);
    }

    /**
     * Write a message to a settings file in the settings directory.
     *
     * @param message
     *            the message
     * @return true, if successful
     */
    public static boolean writeSettings(Message.Builder message)
    {
        return writeSettings(message.build());
    }

    /**
     * Write a message to a settings file in the settings directory.
     *
     * @param message
     *            the message
     * @param flags
     *            the flags
     * @return true, if successful
     */
    public static boolean writeSettings(Message.Builder message, int flags)
    {
        return writeSettings(message.build(), flags);
    }

    /**
     * Clear the settings file for the given class.
     *
     * @param clazz
     *            the class
     * @return true, if the settings are cleared
     */
    public static boolean clearSettings(Class<?> clazz)
    {
        final File file = createSettingsFile(clazz);
        try
        {
            if (file.exists())
                return file.delete();
            return true; // Already clear
        }
        catch (final SecurityException e)
        {
            IJ.log("Unable to clear the settings: " + e.getMessage());
        }
        return false;
    }

    /**
     * Simple class to allow generic message reading.
     *
     * @param <T>
     *            the generic message type
     */
    public static class ConfigurationReader<T extends Message>
    {
        /** the default instance of the message type */
        private final T t;

        /**
         * Instantiates a new configuration reader.
         *
         * @param t
         *            the default instance of the message type
         */
        public ConfigurationReader(T t)
        {
            this.t = t;
        }

        /**
         * Read the message.
         *
         * @return the message
         */
        public T read()
        {
            return read(0);
        }

        /**
         * Read the message.
         *
         * @param flags
         *            the flags
         * @return the message
         */
        @SuppressWarnings("unchecked")
        public T read(int flags)
        {
            T c = (T) readMessage(t.getParserForType(), createSettingsFile(t.getClass()), flags);
            if (c == null && BitFlags.anyNotSet(flags, FLAG_NO_DEFAULT))
                c = t;
            return c;
        }
    }

    /**
     * Read the Calibration from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the Calibration
     */
    public static Calibration readCalibration(int flags)
    {
        return new ConfigurationReader<>(CalibrationProtosHelper.defaultCalibration).read(flags);
    }

    /**
     * Read the PSF from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the PSF
     */
    public static PSF readPSF(int flags)
    {
        return new ConfigurationReader<>(PSFProtosHelper.defaultOneAxisGaussian2DPSF).read(flags);
    }

    /**
     * Read the ResultsSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the ResultsSettings
     */
    public static ResultsSettings readResultsSettings(int flags)
    {
        return new ConfigurationReader<>(ResultsProtosHelper.defaultResultsSettings).read(flags);
    }

    /**
     * Read the FitEngineSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the FitEngineSettings
     */
    private static FitEngineSettings readFitEngineSettings(int flags)
    {
        return new ConfigurationReader<>(FitProtosHelper.defaultFitEngineSettings).read(flags);
    }

    /**
     * Read the GUIFilterSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the GUIFilterSettings
     */
    public static GUIFilterSettings readGUIFilterSettings(int flags)
    {
        return new ConfigurationReader<>(GUIProtosHelper.defaultGUIFilterSettings).read(flags);
    }

    /**
     * Read the PSFCalculatorSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the PSFCalculatorSettings
     */
    public static PSFCalculatorSettings readPSFCalculatorSettings(int flags)
    {
        return new ConfigurationReader<>(GUIProtosHelper.defaultPSFCalculatorSettings).read(flags);
    }

    /**
     * Read the PSFEstimatorSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the PSFEstimatorSettings
     */
    public static PSFEstimatorSettings readPSFEstimatorSettings(int flags)
    {
        return new ConfigurationReader<>(GUIProtosHelper.defaultPSFEstimatorSettings).read(flags);
    }

    /**
     * Read the CreateDataSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the CreateDataSettings
     */
    public static CreateDataSettings readCreateDataSettings(int flags)
    {
        return new ConfigurationReader<>(GUIProtosHelper.defaultCreateDataSettings).read(flags);
    }

    /**
     * Read the LoadLocalisationsSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the LoadLocalisationsSettings
     */
    public static LoadLocalisationsSettings readLoadLocalisationsSettings(int flags)
    {
        return new ConfigurationReader<>(GUIProtosHelper.defaultLoadLocalisationsSettings).read(flags);
    }

    /**
     * Read the ClusteringSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the ClusteringSettings
     */
    public static ClusteringSettings readClusteringSettings(int flags)
    {
        return new ConfigurationReader<>(GUIProtosHelper.defaultClusteringSettings).read(flags);
    }

    /**
     * Read the OpticsSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the OpticsSettings
     */
    public static OpticsSettings readOpticsSettings(int flags)
    {
        return new ConfigurationReader<>(GUIProtosHelper.defaultOpticsSettings).read(flags);
    }

    /**
     * Read the ConfigurationTemplateSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the ConfigurationTemplateSettings
     */
    public static ConfigurationTemplateSettings readConfigurationTemplateSettings(int flags)
    {
        return new ConfigurationReader<>(GUIProtosHelper.defaultConfigurationTemplateSettings).read(flags);
    }

    /**
     * Read the DefaultTemplateSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the DefaultTemplateSettings
     */
    public static DefaultTemplateSettings readDefaultTemplateSettings(int flags)
    {
        return new ConfigurationReader<>(DefaultTemplateSettings.getDefaultInstance()).read(flags);
    }

    /**
     * Read the FitEngineConfiguration from the settings file in the settings directory. This loads the current
     * Calibration, PSF and FitEngineSettings.
     *
     * @param flags
     *            the flags
     * @return the FitEngineConfiguration
     */
    public static FitEngineConfiguration readFitEngineConfiguration(int flags)
    {
        final FitEngineSettings fitEngineSettings = readFitEngineSettings(flags);
        final Calibration calibration = readCalibration(flags);
        final PSF psf = readPSF(flags);
        return new FitEngineConfiguration(fitEngineSettings, calibration, psf);
    }

    /**
     * Read the CameraModelSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the CameraModelSettings
     */
    public static CameraModelSettings readCameraModelSettings(int flags)
    {
        return new ConfigurationReader<>(CameraModelSettings.getDefaultInstance()).read(flags);
    }

    /**
     * Read the CubicSplineSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the CubicSplineSettings
     */
    public static CubicSplineSettings readCubicSplineSettings(int flags)
    {
        return new ConfigurationReader<>(CubicSplineSettings.getDefaultInstance()).read(flags);
    }

    /**
     * Read the NucleusMaskSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the NucleusMaskSettings
     */
    public static NucleusMaskSettings readNucleusMaskSettings(int flags)
    {
        return new ConfigurationReader<>(GUIProtosHelper.defaultNucleusMaskSettings).read(flags);
    }

    /**
     * Read the PSFCreatorSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the PSFCreatorSettings
     */
    public static PSFCreatorSettings readPSFCreatorSettings(int flags)
    {
        return new ConfigurationReader<>(GUIProtosHelper.defaultPSFCreatorSettings).read(flags);
    }

    /**
     * Read the CameraModelManagerSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the CameraModelManagerSettings
     */
    public static CameraModelManagerSettings readCameraModelManagerSettings(int flags)
    {
        return new ConfigurationReader<>(GUIProtosHelper.defaultCameraModelManagerSettings).read(flags);
    }

    /**
     * Read the CameraModelAnalysisSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the CameraModelAnalysisSettings
     */
    public static CameraModelAnalysisSettings readCameraModelAnalysisSettings(int flags)
    {
        return new ConfigurationReader<>(GUIProtosHelper.defaultCameraModelAnalysisSettings).read(flags);
    }

    /**
     * Read the CameraModelFisherInformationAnalysisSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the CameraModelFisherInformationAnalysisSettings
     */
    public static CameraModelFisherInformationAnalysisSettings readCameraModelFisherInformationAnalysisSettings(
            int flags)
    {
        return new ConfigurationReader<>(GUIProtosHelper.defaultCameraModelFisherInformationAnalysisSettings)
                .read(flags);
    }

    /**
     * Read the CubicSplineManagerSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the CubicSplineManagerSettings
     */
    public static CubicSplineManagerSettings readCubicSplineManagerSettings(int flags)
    {
        return new ConfigurationReader<>(GUIProtosHelper.defaultCubicSplineManagerSettings).read(flags);
    }

    /**
     * Read the FailCountManagerSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the FailCountManagerSettings
     */
    public static FailCountManagerSettings readFailCountManagerSettings(int flags)
    {
        return new ConfigurationReader<>(GUIProtosHelper.defaultFailCountManagerSettings).read(flags);
    }

    /**
     * Read the AstigmatismModelSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the AstigmatismModelSettings
     */
    public static AstigmatismModelSettings readAstigmatismModelSettings(int flags)
    {
        return new ConfigurationReader<>(AstigmatismModelSettings.getDefaultInstance()).read(flags);
    }

    /**
     * Read the AstigmatismModelManagerSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the AstigmatismModelManagerSettings
     */
    public static AstigmatismModelManagerSettings readAstigmatismModelManagerSettings(int flags)
    {
        return new ConfigurationReader<>(GUIProtosHelper.defaultAstigmatismModelManagerSettings).read(flags);
    }

    /**
     * Read the CropResultsSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the CropResultsSettings
     */
    public static CropResultsSettings readCropResultsSettings(int flags)
    {
        return new ConfigurationReader<>(CropResultsSettings.getDefaultInstance()).read(flags);
    }

    /**
     * Read the SummariseResultsSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the SummariseResultsSettings
     */
    public static SummariseResultsSettings readSummariseResultsSettings(int flags)
    {
        return new ConfigurationReader<>(SummariseResultsSettings.getDefaultInstance()).read(flags);
    }

    /**
     * Read the ImageJ3DResultsViewerSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the ImageJ3DResultsViewerSettings
     */
    public static ImageJ3DResultsViewerSettings readImageJ3DResultsViewerSettings(int flags)
    {
        return new ConfigurationReader<>(GUIProtosHelper.defaultImageJ3DResultsViewerSettings).read(flags);
    }

    /**
     * Read the SpotFitSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the SpotFitSettings
     */
    public static SpotFitSettings readSpotFitSettings(int flags)
    {
        return new ConfigurationReader<>(GUIProtosHelper.defaultSpotFitSettings).read(flags);
    }

    /**
     * Read the TranslateResultsSettings from the settings file in the settings directory.
     *
     * @param flags
     *            the flags
     * @return the TranslateResultsSettings
     */
    public static TranslateResultsSettings readTranslateResultsSettings(int flags)
    {
        return new ConfigurationReader<>(TranslateResultsSettings.getDefaultInstance()).read(flags);
    }

    /**
     * Write to a settings file in the settings directory.
     *
     * @param fitEngineConfiguration
     *            the fit engine configuration
     * @param flags
     *            the flags
     * @return true, if successful
     */
    public static boolean writeSettings(FitEngineConfiguration fitEngineConfiguration, int flags)
    {
        final FitConfiguration fitConfig = fitEngineConfiguration.getFitConfiguration();
        // This is fail fast
        boolean result = writeSettings(fitEngineConfiguration.getFitEngineSettings(), flags);
        result &= writeSettings(fitConfig.getCalibration(), flags);
        result &= writeSettings(fitConfig.getPSF(), flags);
        return result;
    }

    /**
     * Write the message to file.
     * <p>
     * If this fails then an error message is written to the ImageJ log
     *
     * @param message
     *            the message
     * @param filename
     *            the filename
     * @param flags
     *            the flags
     * @return True if written
     */
    public static boolean writeMessage(Message message, String filename, int flags)
    {
        return writeMessage(message, new File(filename), flags);
    }

    /**
     * Write the message to file.
     * <p>
     * If this fails then an error message is written to the ImageJ log
     *
     * @param message
     *            the message
     * @param file
     *            the file
     * @param flags
     *            the flags
     * @return True if written
     */
    public static boolean writeMessage(Message message, File file, int flags)
    {
        try (FileOutputStream fs = new FileOutputStream(file))
        {
            return writeMessage(message, fs, flags);
        }
        catch (final IOException e)
        {
            //e.printStackTrace();
            if (BitFlags.anyNotSet(flags, FLAG_SILENT))
                IJ.log("Unable to write message: " + e.getMessage());
        }
        return false;
    }

    /**
     * Write the message.
     * <p>
     * If this fails then an error message is written to the ImageJ log
     *
     * @param message
     *            the message
     * @param output
     *            the output
     * @param flags
     *            the flags
     * @return True if saved
     */
    public static boolean writeMessage(Message message, OutputStream output, int flags)
    {
        try
        {
            message.writeDelimitedTo(output);
            return true;
        }
        catch (final IOException e)
        {
            if (BitFlags.anyNotSet(flags, FLAG_SILENT))
                IJ.log("Unable to write message: " + e.getMessage());
        }
        return false;
    }

    /**
     * Read the message from file.
     * <p>
     * If this fails then an error message is written to the ImageJ log
     *
     * @param parser
     *            the parser
     * @param filename
     *            the filename
     * @param flags
     *            the flags
     * @return the message
     */
    public static Message readMessage(Parser<? extends Message> parser, String filename, int flags)
    {
        return readMessage(parser, new File(filename), flags);
    }

    /**
     * Read the message from file.
     * <p>
     * If this fails then an error message is written to the ImageJ log
     *
     * @param parser
     *            the parser
     * @param file
     *            the file
     * @param flags
     *            the flags
     * @return the message
     */
    public static Message readMessage(Parser<? extends Message> parser, File file, int flags)
    {
        try (FileInputStream fs = new FileInputStream(file))
        {
            return readMessage(parser, fs, flags);
        }
        catch (final IOException e)
        {
            //e.printStackTrace();
            // Only print this if the file-not-found flag is present
            // and not silent. This prevents warnings when settings files
            // have yet to be created, i.e. for new users of a settings file.
            if (BitFlags.areSet(flags, FLAG_SHOW_FILE_NOT_FOUND_ON_READ) && !BitFlags.anySet(flags, FLAG_SILENT))
                IJ.log("Unable to read message: " + e.getMessage());
        }
        return null;
    }

    /**
     * Read the message from file.
     * <p>
     * If this fails then an error message is written to the ImageJ log
     *
     * @param parser
     *            the parser
     * @param input
     *            the input
     * @param flags
     *            the flags
     * @return the message
     */
    public static Message readMessage(Parser<? extends Message> parser, InputStream input, int flags)
    {
        try
        {
            return parser.parseDelimitedFrom(input);
        }
        catch (final InvalidProtocolBufferException e)
        {
            //e.printStackTrace();
            if (BitFlags.anyNotSet(flags, FLAG_SILENT))
                IJ.log("Unable to read message: " + e.getMessage());
        }
        return null;
    }

    private static Printer printer = null;

    /**
     * Write the message to a JSON string.
     *
     * @param message
     *            the message
     * @return the JSON string
     */
    public static String toJSON(MessageOrBuilder message)
    {
        return toJSON(message, 0);
    }

    /**
     * Write the message to a JSON string.
     *
     * @param message
     *            the message
     * @param flags
     *            the flags
     * @return the JSON string
     */
    public static String toJSON(MessageOrBuilder message, int flags)
    {
        final StringBuilder sb = new StringBuilder();
        if (toJSON(message, sb, flags))
            return sb.toString();
        return null;
    }

    /**
     * Write the message to file in JSON text format.
     * <p>
     * If this fails then an error message is written to the ImageJ log
     *
     * @param message
     *            the message
     * @param filename
     *            the filename
     * @param flags
     *            the flags
     * @return True if written
     */
    public static boolean toJSON(MessageOrBuilder message, String filename, int flags)
    {
        return toJSON(message, new File(filename), flags);
    }

    /**
     * Write the message to file in JSON text format.
     * <p>
     * If this fails then an error message is written to the ImageJ log
     *
     * @param message
     *            the message
     * @param file
     *            the file
     * @param flags
     *            the flags
     * @return True if written
     */
    public static boolean toJSON(MessageOrBuilder message, File file, int flags)
    {
        try (PrintStream fs = new PrintStream(file))
        {
            return toJSON(message, fs, flags);
        }
        catch (final FileNotFoundException e)
        {
            if (BitFlags.anyNotSet(flags, FLAG_SILENT))
                IJ.log("Unable to write message: " + e.getMessage());
        }
        return false;
    }

    /**
     * Write the message to a JSON string.
     *
     * @param message
     *            the message
     * @param output
     *            the output
     * @param flags
     *            the flags
     * @return the JSON string
     */
    public static boolean toJSON(MessageOrBuilder message, Appendable output, int flags)
    {
        try
        {
            if (printer == null)
                printer = JsonFormat.printer();
            Printer p = printer;
            if (BitFlags.anyNotSet(flags, FLAG_JSON_WHITESPACE))
                p = p.omittingInsignificantWhitespace();
            if (BitFlags.anySet(flags, FLAG_JSON_DEFAULT_VALUES))
                p = p.includingDefaultValueFields();
            p.appendTo(message, output);
            return true;
        }
        catch (final IOException e)
        {
            if (BitFlags.anyNotSet(flags, FLAG_SILENT))
                IJ.log("Unable to write message: " + e.getMessage());
        }
        return false;
    }

    private static com.google.protobuf.util.JsonFormat.Parser parser = null;

    /**
     * Read the message from a JSON string.
     *
     * @param json
     *            the JSON string
     * @param builder
     *            the builder
     * @return true, if successful
     */
    public static boolean fromJSON(String json, Message.Builder builder)
    {
        return fromJSON(json, builder, 0);
    }

    /**
     * Read the message from a JSON string.
     *
     * @param json
     *            the JSON string
     * @param builder
     *            the builder
     * @param flags
     *            the flags
     * @return true, if successful
     */
    public static boolean fromJSON(String json, Message.Builder builder, int flags)
    {
        return fromJSON(new StringReader(json), builder, flags);
    }

    /**
     * Read the message from a JSON string in a file.
     *
     * @param file
     *            the file
     * @param builder
     *            the builder
     * @param flags
     *            the flags
     * @return true, if successful
     */
    public static boolean fromJSON(File file, Message.Builder builder, int flags)
    {
        try (Reader reader = new InputStreamReader(new FileInputStream(file)))
        {
            return fromJSON(reader, builder, flags);
        }
        catch (final IOException e)
        {
            if (BitFlags.anyNotSet(flags, FLAG_SILENT))
                IJ.log("Unable to read message: " + e.getMessage());
        }
        return false;
    }

    /**
     * Read the message from a JSON string.
     *
     * @param reader
     *            the reader
     * @param builder
     *            the builder
     * @param flags
     *            the flags
     * @return true, if successful
     */
    public static boolean fromJSON(Reader reader, Message.Builder builder, int flags)
    {
        try
        {
            if (parser == null)
                parser = JsonFormat.parser();
            parser.merge(reader, builder);
            return true;
        }
        catch (final IOException e)
        {
            if (BitFlags.anyNotSet(flags, FLAG_SILENT))
                IJ.log("Unable to read message: " + e.getMessage());
        }
        return false;
    }
}
