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

import uk.ac.sussex.gdsc.core.data.NotImplementedException;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.UnicodeReader;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.PrecisionMethod;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.LoadLocalisationsSettings;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.AttributePeakResult;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.regex.Pattern;

/**
 * Loads generic localisation files into memory.
 */
public class LoadLocalisations implements PlugIn {
  // TODO - Add support for noise and mean signal.
  // Q. Is this required?

  // Time units for the exposure time cannot be in frames as this makes no sense
  private static EnumSet<TimeUnit> set = EnumSet.allOf(TimeUnit.class);
  private static String[] tUnits;
  private static TimeUnit[] tUnitValues;

  static {
    set.remove(TimeUnit.FRAME);
    set.remove(TimeUnit.UNRECOGNIZED);
    tUnits = new String[set.size()];
    tUnitValues = new TimeUnit[set.size()];
    int i = 0;
    for (final TimeUnit t : set) {
      tUnits[i] = SettingsManager.getName(UnitHelper.getName(t), UnitHelper.getShortName(t));
      tUnitValues[i] = t;
      i++;
    }
  }

  /**
   * The localisation.
   */
  public static class Localisation {
    /** The time. */
    int t;
    /** The id. */
    int id;
    /** The x. */
    float x;
    /** The y. */
    float y;
    /** The z. */
    float z;
    /** The intensity. */
    float intensity;
    /** The x standard deviation. */
    float sx = -1;
    /** The y standard deviation. */
    float sy = -1;
    /** The precision. */
    float precision = -1;
  }

  /**
   * A list of {@link Localisation} objects.
   */
  public static class LocalisationList extends ArrayList<Localisation> {
    private static final long serialVersionUID = 6616011992365324247L;

    /** The calibration. */
    public final Calibration calibration;

    /**
     * Instantiates a new localisation list.
     *
     * <p>The time unit for the calibration is expected to be the exposure time in a valid unit of
     * time (e.g. seconds and not frames)
     *
     * @param calibration the calibration
     */
    public LocalisationList(Calibration calibration) {
      this.calibration = calibration;
    }

    /**
     * Convert to peak results.
     *
     * @param name the name
     * @return the memory peak results
     */
    public MemoryPeakResults toPeakResults(String name) {
      final CalibrationWriter calibrationWriter = new CalibrationWriter(this.calibration);

      // Convert exposure time to milliseconds
      final TypeConverter<TimeUnit> timeConverter =
          calibrationWriter.getTimeConverter(TimeUnit.MILLISECOND);
      calibrationWriter.setExposureTime(timeConverter.convert(calibrationWriter.getExposureTime()));
      // This is currently not a method as it is assumed to be milliseconds.
      // calibration.setTimeUnit(TimeUnit.MILLISECOND);

      // Convert precision to nm
      final TypeConverter<DistanceUnit> distanceConverter =
          calibrationWriter.getDistanceConverter(DistanceUnit.NM);

      final MemoryPeakResults results = new MemoryPeakResults();
      results.setName(name);
      if (!hasPrecision()) {
        calibrationWriter.setPrecisionMethod(null);
      }
      results.setCalibration(calibrationWriter.getCalibration());

      if (size() > 0) {
        // Guess the PSF type from the first localisation
        Localisation l = get(0);
        PSFType psfType = PSFType.CUSTOM;
        if (l.sx != -1) {
          psfType = PSFType.ONE_AXIS_GAUSSIAN_2D;
          if (l.sy != -1) {
            psfType = PSFType.TWO_AXIS_GAUSSIAN_2D;
          }
        }
        results.setPsf(PsfHelper.create(psfType));

        for (int i = 0; i < size(); i++) {
          l = get(i);
          final float intensity = (l.intensity <= 0) ? 1 : (float) (l.intensity);
          final float x = (l.x);
          final float y = (l.y);
          final float z = (l.z);

          float[] params;
          switch (psfType) {
            case CUSTOM:
              params = PeakResult.createParams(0, intensity, x, y, z);
              break;
            case ONE_AXIS_GAUSSIAN_2D:
              params = Gaussian2DPeakResultHelper.createOneAxisParams(0, intensity, x, y, z, l.sx);
              break;
            case TWO_AXIS_GAUSSIAN_2D:
              params =
                  Gaussian2DPeakResultHelper.createTwoAxisParams(0, intensity, x, y, z, l.sx, l.sy);
              break;
            default:
              throw new NotImplementedException("Unsupported PSF type: " + psfType);
          }
          final AttributePeakResult peakResult =
              new AttributePeakResult(l.t, (int) x, (int) y, 0, 0, 0, 0, params, null);
          peakResult.setId(l.id);
          if (l.precision > 0) {
            // Convert to nm
            peakResult.setPrecision(distanceConverter.convert(l.precision));
          }
          results.add(peakResult);
        }
      }

      // Convert to preferred units. This can be done even if the results are empty.
      results.convertToPreferredUnits();

      return results;
    }

    private boolean hasPrecision() {
      for (int i = 0; i < size(); i++) {
        if (get(i).precision > 0) {
          return true;
        }
      }
      return false;
    }
  }

  private static final String TITLE = "Load Localisations";

  @Override
  public void run(String arg) {
    SMLMUsageTracker.recordPlugin(this.getClass(), arg);

    final LoadLocalisationsSettings.Builder settings =
        SettingsManager.readLoadLocalisationsSettings(0).toBuilder();

    final String[] path = ImageJUtils.decodePath(settings.getLocalisationsFilename());
    final OpenDialog chooser = new OpenDialog("Localisations_File", path[0], path[1]);
    if (chooser.getFileName() == null) {
      return;
    }

    settings.setLocalisationsFilename(chooser.getDirectory() + chooser.getFileName());

    final LocalisationList localisations = loadLocalisations(settings);

    SettingsManager.writeSettings(settings.build());

    if (localisations == null) {
      // Cancelled
      return;
    }

    if (localisations.isEmpty()) {
      IJ.error(TITLE, "No localisations could be loaded");
      return;
    }

    final MemoryPeakResults results = localisations.toPeakResults(settings.getName());

    // Create the in-memory results
    if (results.size() > 0) {
      MemoryPeakResults.addResults(results);
    }

    final String msg = "Loaded " + TextUtils.pleural(results.size(), "localisation");
    IJ.showStatus(msg);
    ImageJUtils.log(msg);
  }

  /**
   * Load localisations.
   *
   * @param settings the settings
   * @return the localisation list
   */
  static LocalisationList loadLocalisations(LoadLocalisationsSettings.Builder settings) {
    if (!getFields(settings)) {
      return null;
    }

    final LocalisationList localisations = new LocalisationList(settings.getCalibration());

    final String comment = settings.getComment();
    final boolean hasComment = !TextUtils.isNullOrEmpty(comment);
    int errors = 0;
    int count = 0;
    int headerCount = Math.max(0, settings.getHeaderLines());

    try (BufferedReader input = new BufferedReader(
        new UnicodeReader(new FileInputStream(settings.getLocalisationsFilename()), null))) {
      final Pattern p = Pattern.compile(settings.getDelimiter());

      final int it = settings.getFieldT();
      final int iid = settings.getFieldId();
      final int ix = settings.getFieldX();
      final int iy = settings.getFieldY();
      final int iz = settings.getFieldZ();
      final int ii = settings.getFieldI();
      final int isx = settings.getFieldSx();
      final int isy = settings.getFieldSy();
      final int ip = settings.getFieldPrecision();

      String line;
      while ((line = input.readLine()) != null) {
        // Skip header
        if (headerCount-- > 0) {
          continue;
        }
        // Skip empty lines
        if (line.length() == 0) {
          continue;
        }
        // Skip comments
        if (hasComment && line.startsWith(comment)) {
          continue;
        }

        count++;
        final String[] fields = p.split(line);

        final Localisation l = new Localisation();
        try {
          if (it >= 0) {
            l.t = Integer.parseInt(fields[it]);
          }
          if (iid >= 0) {
            l.id = Integer.parseInt(fields[iid]);
          }
          l.x = Float.parseFloat(fields[ix]);
          l.y = Float.parseFloat(fields[iy]);
          if (iz >= 0) {
            l.z = Float.parseFloat(fields[iz]);
          }
          if (ii >= 0) {
            l.intensity = Float.parseFloat(fields[ii]);
          }
          if (isx >= 0) {
            l.sy = l.sx = Integer.parseInt(fields[isx]);
          }
          if (isy >= 0) {
            l.sy = Integer.parseInt(fields[isy]);
          }
          if (ip >= 0) {
            l.precision = Float.parseFloat(fields[ip]);
          }

          localisations.add(l);
        } catch (final NumberFormatException ex) {
          if (errors++ == 0) {
            ImageJUtils.log("%s error on record %d: %s", TITLE, count, ex.getMessage());
          }
        } catch (final IndexOutOfBoundsException ex) {
          if (errors++ == 0) {
            ImageJUtils.log("%s error on record %d: %s", TITLE, count, ex.getMessage());
          }
        }
      }
    } catch (final IOException ex) {
      ImageJUtils.log("%s IO error: %s", TITLE, ex.getMessage());
    }
    if (errors != 0) {
      ImageJUtils.log("%s has %d / %d error lines", TITLE, errors, count);
    }

    return localisations;
  }

  private static boolean getFields(LoadLocalisationsSettings.Builder settings) {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    gd.addMessage("Load delimited localisations");
    gd.addStringField("Dataset_name", settings.getName(), 30);

    gd.addMessage("Calibration:");
    // Allow the full camera type top be captured
    final Calibration.Builder calibrationBuilder = settings.getCalibrationBuilder();
    final CalibrationWriter cw = new CalibrationWriter(calibrationBuilder);
    PeakFit.addCameraOptions(gd, 0, cw);
    // Only primitive support for other calibration
    gd.addNumericField("Pixel_size", cw.getNmPerPixel(), 3, 8, "nm");
    gd.addNumericField("Exposure_time", cw.getExposureTime(), 3, 8, "");

    // This is the unit for the exposure time (used to convert the exposure time to milliseconds).
    // Use the name as the list is a truncated list of the full enum.
    final TimeUnit t = calibrationBuilder.getTimeCalibration().getTimeUnit();
    gd.addChoice("Time_unit", tUnits,
        SettingsManager.getName(UnitHelper.getName(t), UnitHelper.getShortName(t)));

    gd.addMessage("Records:");
    gd.addNumericField("Header_lines", settings.getHeaderLines(), 0);
    gd.addStringField("Comment", settings.getComment());
    gd.addStringField("Delimiter", settings.getDelimiter());
    gd.addChoice("Distance_unit", SettingsManager.getDistanceUnitNames(),
        cw.getDistanceUnitValue());
    gd.addChoice("Intensity_unit", SettingsManager.getIntensityUnitNames(),
        cw.getIntensityUnitValue());

    gd.addMessage("Define the fields:");
    gd.addNumericField("T", settings.getFieldT(), 0);
    gd.addNumericField("ID", settings.getFieldId(), 0);
    gd.addNumericField("X", settings.getFieldX(), 0);
    gd.addNumericField("Y", settings.getFieldY(), 0);
    gd.addNumericField("Z", settings.getFieldZ(), 0);
    gd.addNumericField("Intensity", settings.getFieldI(), 0);
    gd.addNumericField("Sx", settings.getFieldSx(), 0);
    gd.addNumericField("Sy", settings.getFieldSy(), 0);
    gd.addNumericField("Precision", settings.getFieldPrecision(), 0);
    gd.addChoice("Precision_method", SettingsManager.getPrecisionMethodNames(),
        cw.getPrecisionMethodValue());

    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    settings.setName(getNextString(gd, settings.getName()));
    cw.setCameraType(SettingsManager.getCameraTypeValues()[gd.getNextChoiceIndex()]);
    cw.setNmPerPixel(gd.getNextNumber());
    cw.setExposureTime(gd.getNextNumber());

    // The time units used a truncated list so look-up the value from the index
    calibrationBuilder.getTimeCalibrationBuilder()
        .setTimeUnit(tUnitValues[gd.getNextChoiceIndex()]);

    settings.setHeaderLines((int) gd.getNextNumber());
    settings.setComment(gd.getNextString());
    settings.setDelimiter(getNextString(gd, settings.getDelimiter()));

    cw.setDistanceUnit(DistanceUnit.forNumber(gd.getNextChoiceIndex()));
    cw.setIntensityUnit(IntensityUnit.forNumber(gd.getNextChoiceIndex()));

    final int[] columns = new int[9];
    for (int i = 0; i < columns.length; i++) {
      columns[i] = (int) gd.getNextNumber();
    }

    {
      int i = 0;
      settings.setFieldI(columns[i++]);
      settings.setFieldId(columns[i++]);
      settings.setFieldX(columns[i++]);
      settings.setFieldY(columns[i++]);
      settings.setFieldZ(columns[i++]);
      settings.setFieldI(columns[i++]);
      settings.setFieldSx(columns[i++]);
      settings.setFieldSy(columns[i++]);
      settings.setFieldPrecision(columns[i++]);
    }

    cw.setPrecisionMethod(PrecisionMethod.forNumber(gd.getNextChoiceIndex()));

    // Collect the camera calibration
    gd.collectOptions();

    // Validate after reading the dialog (so we store the last entered values)
    if (gd.invalidNumber()) {
      IJ.error(TITLE, "Invalid number in input fields");
      return false;
    }

    for (int i = 0; i < columns.length; i++) {
      if (columns[i] < 0) {
        continue;
      }
      for (int j = i + 1; j < columns.length; j++) {
        if (columns[j] < 0) {
          continue;
        }
        if (columns[i] == columns[j]) {
          IJ.error(TITLE, "Duplicate indicies: " + columns[i]);
          return false;
        }
      }
    }
    if (cw.getNmPerPixel() <= 0) {
      IJ.error(TITLE, "Require positive pixel pitch");
      return false;
    }
    if (cw.isCcdCamera()) {
      if (!cw.hasCountPerPhoton()) {
        IJ.error(TITLE, "Require positive count/photon for CCD camera type");
        return false;
      }
    } else {
      // Q.Validate other camera types?
    }
    if (settings.getFieldX() < 0 || settings.getFieldY() < 0) {
      IJ.error(TITLE, "Require valid X and Y indices");
      return false;
    }
    return true;
  }

  private static String getNextString(GenericDialog gd, String defaultValue) {
    final String value = gd.getNextString();
    if (TextUtils.isNullOrEmpty(value)) {
      return defaultValue;
    }
    return value;
  }
}
