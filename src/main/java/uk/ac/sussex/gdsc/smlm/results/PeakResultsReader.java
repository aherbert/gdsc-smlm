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

import uk.ac.sussex.gdsc.core.logging.TrackProgress;
import uk.ac.sussex.gdsc.core.utils.BitFlagUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.UnicodeReader;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.PSFHelper;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedureX;
import uk.ac.sussex.gdsc.smlm.utils.XStreamUtils;

import com.google.protobuf.InvalidProtocolBufferException;
import com.google.protobuf.util.JsonFormat;

import java.awt.Rectangle;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.util.InputMismatchException;
import java.util.Locale;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Reads the fit results from file.
 */
public class PeakResultsReader {
  // Set up to read two-axis (and theta) Gaussian 2D data into the current format
  private static final int isx;
  private static final int isy;
  private static final int ia;
  private static final int nTwoAxis;
  private static final int nTwoAxisAndTheta;

  static {
    final PSF psf = PSFHelper.create(PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D);
    final int[] indices = PSFHelper.getGaussian2DWxWyIndices(psf);
    isx = indices[0];
    isy = indices[1];
    ia = PSFHelper.getGaussian2DAngleIndex(psf);
    nTwoAxis = PeakResult.STANDARD_PARAMETERS + 2;
    nTwoAxisAndTheta = PeakResult.STANDARD_PARAMETERS + 3;
  }

  /** Index of the background in the parameters array in the legacy GDSC file format. */
  static final int LEGACY_FORMAT_BACKGROUND = 0;
  /** Index of the signal intensity in the parameters array in the legacy GDSC file format. */
  static final int LEGACY_FORMAT_SIGNAL = 1;
  /** Index of the angle in the parameters array in the legacy GDSC file format. */
  static final int LEGACY_FORMAT_ANGLE = 2;
  /** Index of the x-position in the parameters array in the legacy GDSC file format. */
  static final int LEGACY_FORMAT_X_POSITION = 3;
  /** Index of the y-position in the parameters array in the legacy GDSC file format. */
  static final int LEGACY_FORMAT_Y_POSITION = 4;
  /** Index of the x-standard deviation in the parameters array in the legacy GDSC file format. */
  static final int LEGACY_FORMAT_X_SD = 5;
  /** Index of the y-standard deviation in the parameters array in the legacy GDSC file format. */
  static final int LEGACY_FORMAT_Y_SD = 6;

  /**
   * The columns to recognise in the ImageJ table results header for version 1/2. Version 3+ may
   * also have this but we can distinguish because V1/2 had Amplitude/Signal and V3+ does not.
   */
  private static String IMAGEJ_TABLE_RESULTS_HEADER_V1_V2 = "origX\torigY\torigValue\tError\tNoise";

  /** The space pattern. */
  private static Pattern spacePattern = Pattern.compile(" ");
  /** The tab pattern. */
  private static Pattern tabPattern = Pattern.compile("\t");
  /** Simple whitespace pattern for tabs of spaces. */
  private static Pattern whitespacePattern = Pattern.compile("[\t ]");

  private boolean useScanner;
  private boolean rawResults;

  private final String filename;
  private String header;
  private FileFormat format;
  private String version;
  private String name;
  private ImageSource source;
  private Rectangle bounds;
  private CalibrationWriter calibration;
  private PSF psf;
  private String configuration;
  private TrackProgress tracker;
  private ResultOption[] options;

  private boolean deviations;
  private boolean readEndFrame;
  private boolean readId;
  private boolean readPrecision;
  private boolean readSource;
  private int smlmVersion = 4; // Assume the current

  /**
   * Instantiates a new peak results reader.
   *
   * @param filename the filename
   */
  public PeakResultsReader(String filename) {
    this.filename = filename;
  }

  /**
   * Gets the header from the results file.
   *
   * @return The header from the results file
   */
  public String getHeader() {
    if (header == null) {
      try (BufferedReader input =
          new BufferedReader(new UnicodeReader(new FileInputStream(filename), null))) {
        final StringBuilder sb = new StringBuilder();
        String line;
        int count = 0;
        while ((line = input.readLine()) != null) {
          count++;
          if (count == 1) {
            // The NSTORM file format does not have comment characters but does have a single header
            // line
            if (line.startsWith("Channel Name")) {
              sb.append(line).append("\n");
              break;
            }
            // User may try and load the text saved directly from the ImageJ Table Results
            if (line.contains(IMAGEJ_TABLE_RESULTS_HEADER_V1_V2)) {
              sb.append(line).append("\n");
              break;
            }
          }
          if (line.length() == 0) {
            continue;
          }
          if (line.charAt(0) == '#') {
            sb.append(line).append("\n");
          } else {
            break;
          }
        }
        header = sb.toString();

        version = getField("FileVersion");

        format = FileFormat.UNKNOWN;
        if (!guessFormatFromVersion()) {
          guessFormat(line);
        }
      } catch (final IOException ex) {
        // ignore
      }
    }
    return header;
  }

  private boolean guessFormatFromVersion() {
    // Extract information about the file format
    if (version.length() > 0) {
      if (version.startsWith("Binary")) {
        format = FileFormat.SMLM_BINARY;
      } else if (version.startsWith("Text")) {
        format = FileFormat.SMLM_TEXT;
      } else {
        return false;
      }

      if (version.contains(".V4")) {
        smlmVersion = 4;
      } else if (version.contains(".V3")) {
        smlmVersion = 3;
      } else if (version.contains(".V2")) {
        smlmVersion = 2;
      }

      if (smlmVersion == 1) {
        // The original files did not have the version tag
        deviations =
            header.contains((format == FileFormat.SMLM_BINARY) ? "iiiifdfffffffffffffff" : "+/-");
        readEndFrame = readId = readPrecision = false;
      } else {
        deviations = version.contains(".D1");
        // Extended marker has a bit flag within .E[0-9]+.:
        int startIndex = version.indexOf(".E");
        if (startIndex != -1) {
          startIndex += 2;
          final int endIndex = version.indexOf('.', startIndex);
          if (endIndex != -1) {
            // Get the flags
            try {
              final int flags = Integer.parseInt(version.substring(startIndex, endIndex));
              readEndFrame = BitFlagUtils.areSet(flags, SMLMFilePeakResults.FLAG_END_FRAME);
              readId = BitFlagUtils.areSet(flags, SMLMFilePeakResults.FLAG_ID);
              readPrecision = BitFlagUtils.areSet(flags, SMLMFilePeakResults.FLAG_PRECISION);
            } catch (final NumberFormatException ex) {
              // Ignore
            }
          }
        }
      }
      return true;
    }
    return false;
  }

  private void guessFormat(String firstLine) {
    if (header.length() == 0) {
      // No header. Check non-text formats.
      if (TSFPeakResultsReader.isTSF(filename)) {
        format = FileFormat.TSF_BINARY;
        return;
      }

      // We cannot continue guess if there is no non-header data
      if (TextUtils.isNullOrEmpty(firstLine)) {
        return;
      }
    }

    // Support reading old IJ table results.
    // The latest table results have dynamic columns so these must be loaded manually
    // as guessing the column format is not supported.
    if (header.contains(IMAGEJ_TABLE_RESULTS_HEADER_V1_V2)) {
      format = FileFormat.SMLM_TABLE;
      smlmVersion = 2; // V1/V2 doesn't matter
      readId = header.startsWith("#");
      readEndFrame = header.contains("\tEnd ");
      deviations = header.contains("\t+/-\t");
      readSource = header.contains("Source\t");

      // Check for Nikon NSTORM header
    } else if (header.contains("Channel Name")) {
      format = FileFormat.NSTORM;
    } else if (header.contains("<localizations ")) {
      // RapidSTORM can use the MALK format
      if (isMALKFormat(firstLine)) {
        format = FileFormat.MALK;
      } else {
        // There may be many other formats for RapidSTORM. Just support the one we know about.
        format = FileFormat.RAPID_STORM;
      }

      // Check for MALK format: X,Y,T,Signal
    } else if (isMALKFormat(firstLine)) {
      format = FileFormat.MALK;
    } else {
      format = FileFormat.UNKNOWN;
    }
  }

  /**
   * @return True if the results file is binary.
   */
  public boolean isBinary() {
    getFormat();
    return format == FileFormat.SMLM_BINARY;
  }

  /**
   * @return The file format.
   */
  public FileFormat getFormat() {
    if (format == null) {
      getHeader();
    }
    return format;
  }

  /**
   * @return The bounds specified in the results header.
   */
  public Rectangle getBounds() {
    if (bounds == null) {
      getHeader();
      bounds = new Rectangle(0, 0, 0, 0);
      if (header != null) {
        final Pattern pattern = Pattern.compile("x(\\d+) y(\\d+) w(\\d+) h(\\d+)");
        final Matcher match = pattern.matcher(header);
        if (match.find()) {
          bounds.x = Integer.parseInt(match.group(1));
          bounds.y = Integer.parseInt(match.group(2));
          bounds.width = Integer.parseInt(match.group(3));
          bounds.height = Integer.parseInt(match.group(4));
        }
      }
    }
    return bounds;
  }

  /**
   * @return The name specified in the results header.
   */
  public String getName() {
    if (name == null) {
      getHeader();
      name = getField("Name");
    }
    return name;
  }

  /**
   * @return The source specified in the results header.
   */
  public ImageSource getSource() {
    if (source == null) {
      getHeader();
      if (header != null) {
        final String xml = getField("Source");
        if (xml != null && xml.length() > 0 && xml.startsWith("<")) {
          // Convert the XML back
          source = ImageSource.fromXML(xml);
        }
      }
    }
    return source;
  }

  /**
   * @return The calibration specified in the results header.
   */
  @SuppressWarnings("deprecation")
  public Calibration getCalibration() {
    if (calibration == null) {
      getHeader();
      if (header != null && header.length() > 0) {
        if (format == FileFormat.RAPID_STORM) {
          calibration.setDistanceUnit(DistanceUnit.NM);

          // RapidSTORM has a resolution attribute in the header in units of px m^-1
          final Pattern pattern = Pattern.compile("resolution=\"([^ ]+) px m");
          final Matcher match = pattern.matcher(header);
          if (match.find()) {
            try {
              final float resolution = Float.parseFloat(match.group(1));
              if (Double.isFinite(resolution) && resolution > 0) {
                final double nmPerPixel = (float) (1e9 / resolution);
                calibration = new CalibrationWriter();
                calibration.setNmPerPixel(nmPerPixel);
              }
            } catch (final NumberFormatException ex) {
              // Ignore
            }
          }
        } else {
          final String calibrationString = getField("Calibration");
          if (calibrationString != null && calibrationString.length() > 0) {
            // Older formats used XML
            // TODO - Test this still works using older code, e.g. that
            // released to Fiji.
            if (calibrationString.startsWith("<")) {
              // Convert the XML back
              try {
                final uk.ac.sussex.gdsc.smlm.results.Calibration cal =
                    (uk.ac.sussex.gdsc.smlm.results.Calibration) XStreamUtils
                        .fromXML(calibrationString);
                cal.validate();

                // Convert to a calibration helper
                calibration = new CalibrationWriter();
                if (cal.hasNmPerPixel()) {
                  calibration.setNmPerPixel(cal.getNmPerPixel());
                }
                if (cal.hasGain()) {
                  calibration.setCountPerPhoton(cal.getGain());
                }
                if (cal.hasExposureTime()) {
                  calibration.setExposureTime(cal.getExposureTime());
                }
                if (cal.hasReadNoise()) {
                  calibration.setReadNoise(cal.getReadNoise());
                }
                if (cal.hasBias()) {
                  calibration.setBias(cal.getBias());
                }
                if (cal.emCCD) {
                  calibration.setCameraType(CameraType.EMCCD);
                }
                if (cal.hasAmplification()) {
                  if (cal.hasGain()) {
                    calibration.setQuantumEfficiency(cal.getGain() / cal.getAmplification());
                  }
                }

                // Previous version were always in fixed units
                calibration.setDistanceUnit(DistanceUnit.PIXEL);
                calibration.setIntensityUnit(IntensityUnit.COUNT);
                calibration.setAngleUnit(AngleUnit.DEGREE);
              } catch (final ClassCastException ex) {
                ex.printStackTrace();
              } catch (final Exception ex) {
                ex.printStackTrace();
              }
            } else {
              // Assume JSON format
              try {
                final Calibration.Builder calibrationBuilder = Calibration.newBuilder();
                JsonFormat.parser().merge(calibrationString, calibrationBuilder);
                calibration = new CalibrationWriter(calibrationBuilder);
              } catch (final InvalidProtocolBufferException ex) {
                // This should be OK
                System.err.println("Unable to deserialise the Calibration settings");
              }
            }
          }

          if (format == FileFormat.MALK) {
            if (calibration == null) {
              calibration = new CalibrationWriter();
            }
            calibration.setDistanceUnit(DistanceUnit.NM);
            calibration.setIntensityUnit(IntensityUnit.PHOTON);
          }
        }
      }

      // Calibration is a smart object so we can create an empty one
      if (calibration == null) {
        calibration = new CalibrationWriter();
      }
    }
    return calibration.getCalibration();
  }

  /**
   * Gets the PSF specified in the results header.
   *
   * @return the PSF
   */
  public PSF getPSF() {
    if (psf == null) {
      getHeader();
      if (header != null && header.length() > 0) {
        final String json = getField("PSF");
        if (json != null && json.length() > 0
        // && json.startsWith("[")
        ) {
          // Convert the JSON back
          try {
            final PSF.Builder psfBuilder = PSF.newBuilder();
            JsonFormat.parser().merge(json, psfBuilder);
            psf = psfBuilder.build();
          } catch (final InvalidProtocolBufferException ex) {
            // This should be OK
            System.err.println("Unable to deserialise the PSF settings");
          }
        }

        if (psf == null) {
          // Guess base on type
          switch (format) {
            case NSTORM:
            case RAPID_STORM:
              // We currently only support two axis data
              psf = PSFHelper.create(PSFType.TWO_AXIS_GAUSSIAN_2D);
              break;

            // Note: Older GDSC results were all TwoAxisAndTheta
            case SMLM_BINARY:
            case SMLM_TABLE:
            case SMLM_TEXT:
              if (smlmVersion < 3) {
                psf = PSFHelper.create(PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D);
              }
              break;

            case TSF_BINARY:
              // This is read separately so ignore
              break;

            case MALK:
              // Use the custom type with no additional PSF parameters
              psf = PSFHelper.create(PSFType.CUSTOM);
              break;

            case UNKNOWN:
            default:
              break;
          }
        }
      }
    }
    return psf;
  }

  /**
   * Gets the configuration specified in the results header.
   *
   * @return The configuration specified in the results header
   */
  public String getConfiguration() {
    if (configuration == null) {
      getHeader();
      configuration = getField("Configuration");
      if (configuration != null && configuration.length() > 0 && configuration.charAt(0) == '<') {
        // Format the XML back.
        // This is only relevant for configuration serialised as XML.
        // The GDSC SMLM configuration is now a JSON string.
        configuration = uk.ac.sussex.gdsc.core.utils.XmlUtils.formatXml(configuration);
      }
    }
    return configuration;
  }

  private String getField(String name) {
    if (header != null) {
      final Pattern pattern = Pattern.compile(name + " ([^\\n]+)");
      final Matcher match = pattern.matcher(header);
      if (match.find()) {
        return match.group(1);
      }
    }
    return "";
  }

  /**
   * Read the results from the file. The file is read for each invocation of this method, i.e. no
   * caching is done.
   *
   * @return The peak results
   */
  public MemoryPeakResults getResults() {
    getHeader();
    if (header == null || format == null || format == FileFormat.UNKNOWN) {
      return null;
    }

    // Always do this as we can guess the PSF based on the file type
    getPSF();

    // Use a switch statement with no break statements to fall through
    switch (format) {
      case SMLM_BINARY:
      case SMLM_TEXT:
      case MALK:
        // Read SMLM data. We do this for MALK files because we may have written them.
        getName();
        getSource();
        getBounds();
        getConfiguration();
        // RapidSTORM has calibration too
      case RAPID_STORM:
        getCalibration();
      default:
        break;
    }

    MemoryPeakResults results = null;
    switch (format) {
      case SMLM_BINARY:
        results = readBinary();
        break;
      case SMLM_TEXT:
        results = readText();
        break;
      case SMLM_TABLE:
        results = readTable();
        break;
      case RAPID_STORM:
        results = readRapidSTORM();
        break;
      case NSTORM:
        results = readNSTORM();
        break;
      case MALK:
        results = readMALK();
        break;
      case TSF_BINARY:
        results = readTSF();
        break;
      default:
        break;
    }
    if (results != null) {
      results.trimToSize();

      if (!rawResults) {
        if (results.getPSF() != null) {
          // The TSF reader may set the PSF so copy it back
          psf = results.getPSF();
          simplifyPSF(results);

          // Add mean intensity if a Gaussian 2D function
          addMeanIntensity(results);
        }

        // Convert to the preferred units if possible
        results.convertToPreferredUnits();
      }
    }
    return results;
  }

  private void simplifyPSF(MemoryPeakResults results) {
    switch (psf.getPsfType()) {
      case TWO_AXIS_AND_THETA_GAUSSIAN_2D:
        simplifyTwoAxisAndTheta(results);
        break;

      case TWO_AXIS_GAUSSIAN_2D:
        simplifyTwoAxis(results, false);
        break;

      case ASTIGMATIC_GAUSSIAN_2D:
      case ONE_AXIS_GAUSSIAN_2D:
      case UNRECOGNIZED:
      default:
        break;
    }
  }

  /**
   * Simplify two axis and theta Gaussian 2D data.
   *
   * @param results the results
   */
  private void simplifyTwoAxisAndTheta(MemoryPeakResults results) {
    final int[] indices = PSFHelper.getGaussian2DWxWyIndices(psf);
    final int isx = indices[0];
    final int isy = indices[1];
    final int ia = PSFHelper.getGaussian2DAngleIndex(psf);

    // Determine if the angle is non-zero with asymmetric widths
    if (results.forEach(new PeakResultProcedureX() {
      @Override
      public boolean execute(PeakResult peakResult) {
        return (peakResult.getParameter(ia) != 0
            && peakResult.getParameter(isx) != peakResult.getParameter(isy));
      }
    })) {
      // Nothing to simplify
      return;
    }

    simplifyTwoAxis(results, true);
  }

  /**
   * Simplify two axis Gaussian 2D data.
   *
   * @param results the results
   * @param removeTheta the remove theta flag
   */
  private void simplifyTwoAxis(MemoryPeakResults results, boolean removeTheta) {
    final int[] indices = PSFHelper.getGaussian2DWxWyIndices(psf);
    final int isx = indices[0];
    final int isy = indices[1];

    // Columns to remove
    int remove = (removeTheta) ? 1 : 0;
    // New PSF type
    PSFType psfType;

    // Determine if sy is redundant
    if (results.forEach(new PeakResultProcedureX() {
      @Override
      public boolean execute(PeakResult peakResult) {
        return (peakResult.getParameter(isx) != peakResult.getParameter(isy));
      }
    })) {
      if (!removeTheta) {
        // Already a TwoAxis Gaussian
        return;
      }

      // Otherwise this was a TwoAxisAndTheta with 1 column to remove
      // so it should be simplified
      psfType = PSFType.TWO_AXIS_GAUSSIAN_2D;
    } else {
      // sy is redundant so remove another column
      psfType = PSFType.ONE_AXIS_GAUSSIAN_2D;
      remove++;
    }

    // Update the PSF
    final PSF.Builder builder = psf.toBuilder();
    builder.setPsfType(psfType);
    results.setPSF(psf = builder.build());

    // Update the results.
    // We can directly manipulate the params array
    final int newLength = results.getf(0).getNumberOfParameters() - remove;
    if (!deviations) {
      results.forEach(new PeakResultProcedure() {
        @Override
        public void execute(PeakResult peakResult) {
          peakResult.resizeParameters(newLength);
        }
      });
    } else {
      results.forEach(new PeakResultProcedure() {
        @Override
        public void execute(PeakResult peakResult) {
          peakResult.resizeParameters(newLength);
          peakResult.resizeParameterDeviations(newLength);
        }
      });
    }
  }

  /**
   * Add mean intensity if a Gaussian 2D function and no mean intensity has been read.
   *
   * @param results the results
   */
  private void addMeanIntensity(MemoryPeakResults results) {
    // Some formats already read the mean intensity, e.g. SMLM, TSF
    if (PSFHelper.isGaussian2D(psf) && !results.hasMeanIntensity()) {
      final int[] indices = PSFHelper.getGaussian2DWxWyIndices(psf);
      final int isx = indices[0];
      final int isy = indices[1];
      results.forEach(new PeakResultProcedure() {
        @Override
        public void execute(PeakResult peakResult) {
          final float[] p = peakResult.getParameters();
          final float u = (float) Gaussian2DPeakResultHelper
              .getMeanSignalUsingP05(p[PeakResult.INTENSITY], p[isx], p[isy]);
          peakResult.setMeanIntensity(u);
        }
      });
    }
  }

  private int position;

  private MemoryPeakResults readBinary() {
    final MemoryPeakResults results = createResults();

    // Units were added in version 3
    int nFields;
    boolean gaussian2Dformat;
    boolean readMeanIntensity = false;
    if (smlmVersion < 3) {
      nFields = 7;
      gaussian2Dformat = true;
      calibration.setIntensityUnit(IntensityUnit.COUNT);
      calibration.setDistanceUnit(DistanceUnit.PIXEL);
      calibration.setAngleUnit(AngleUnit.DEGREE);
    } else {
      // The mean signal was added in version 4
      readMeanIntensity = smlmVersion > 3;
      gaussian2Dformat = false;
      // The number of fields should be within the PSF object
      nFields = new PeakResultConversionHelper(null, psf).getNames().length;
    }

    try (FileInputStream fis = new FileInputStream(filename)) {
      try (DataInputStream input = new DataInputStream(fis)) {
        // Seek to the start of the binary data by just reading the header again
        BinaryFilePeakResults.readHeader(input);

        @SuppressWarnings("resource")
        final FileChannel channel = fis.getChannel();

        // Format: [i]i[i]iifdf + n*f [+ n*f]
        // where [] are optional and n is the number of fields
        int flags = 0;
        if (readEndFrame) {
          flags += SMLMFilePeakResults.FLAG_END_FRAME;
        }
        if (readId) {
          flags += SMLMFilePeakResults.FLAG_ID;
        }
        if (readPrecision) {
          flags += SMLMFilePeakResults.FLAG_PRECISION;
        }
        int length = BinaryFilePeakResults.getDataSize(deviations, flags, nFields);
        if (!readMeanIntensity) {
          length -= 4; // No float field for mean signal
        }
        final byte[] buffer = new byte[length];

        int c = 0;
        final boolean convert = smlmVersion == 1;
        // Halted by the EOFException
        while (true) {
          // Note: Reading single strips seems fast enough at the moment.
          // This could be modified to read larger blocks of data if necessary.

          final int bytesRead = input.read(buffer);
          if (bytesRead != length) {
            break;
          }

          position = 0;
          final int id = (readId) ? readInt(buffer) : 0;
          final int peak = readInt(buffer);
          final int endPeak = (readEndFrame) ? readInt(buffer) : peak;
          final int origX = readInt(buffer);
          final int origY = readInt(buffer);
          final float origValue = readFloat(buffer);
          final double error = readDouble(buffer);
          final float noise = readFloat(buffer);
          final float meanSignal = (readMeanIntensity) ? readFloat(buffer) : 0;
          float[] params = readData(buffer, new float[nFields]);
          float[] paramsStdDev = (deviations) ? readData(buffer, new float[nFields]) : null;

          // Convert format which had the full Gaussian2D parameters array
          if (gaussian2Dformat) {
            // Convert old binary format with the amplitude to signal
            if (convert) {
              params[LEGACY_FORMAT_SIGNAL] *=
                  2 * Math.PI * params[LEGACY_FORMAT_X_SD] * params[LEGACY_FORMAT_Y_SD];
            }
            params = mapGaussian2DFormatParams(params);
            paramsStdDev = mapGaussian2DFormatDeviations(paramsStdDev);
          }

          if (readPrecision) {
            results.add(createResult(peak, origX, origY, origValue, error, noise, meanSignal,
                params, paramsStdDev, endPeak, id, readFloat(buffer)));
          } else {
            results.add(createResult(peak, origX, origY, origValue, error, noise, meanSignal,
                params, paramsStdDev, endPeak, id));
          }

          if (++c % 512 == 0) {
            showProgress(channel);
          }
        }
      }
    } catch (final EOFException ex) {
      // Ignore. Binary data does not have a size so it is read until the EOF.
    } catch (final IOException ex) {
      // Log this but still return the results that have been read
      ex.printStackTrace();
    }
    return results;
  }

  /**
   * Creates the result. End frame and id will be ignored if the read flags are unset.
   *
   * @return the peak result
   */
  private PeakResult createResult(int startFrame, int origX, int origY, float origValue,
      double error, float noise, float meanIntensity, float[] params, float[] paramsStdDev,
      int endFrame, int id) {
    if (readEndFrame) {
      return new ExtendedPeakResult(startFrame, origX, origY, origValue, error, noise,
          meanIntensity, params, paramsStdDev, endFrame, id);
    }
    if (readId) {
      return new IdPeakResult(startFrame, origX, origY, origValue, error, noise, meanIntensity,
          params, paramsStdDev, id);
    }
    return new PeakResult(startFrame, origX, origY, origValue, error, noise, meanIntensity, params,
        paramsStdDev);
  }

  /**
   * Creates the result. End frame and id will be ignored if the read flags are unset. Only call
   * this when readPrecision is true as it creates an AttributePeakResult and uses the precision.
   *
   * @return the peak result
   */
  private PeakResult createResult(int startFrame, int origX, int origY, float origValue,
      double error, float noise, float meanIntensity, float[] params, float[] paramsStdDev,
      int endFrame, int id, double precision) {
    final AttributePeakResult r = new AttributePeakResult(startFrame, origX, origY, origValue,
        error, noise, meanIntensity, params, paramsStdDev);
    if (readEndFrame) {
      r.setEndFrame(endFrame);
    }
    if (readId) {
      r.setId(id);
    }
    r.setPrecision(precision);
    return r;
  }

  private float[] mapGaussian2DFormatParams(float[] params) {
    final float[] p = mapGaussian2DFormat(params);
    // New format does not have the bias
    p[PeakResult.BACKGROUND] -= calibration.getBias();
    return p;
  }

  private static float[] mapGaussian2DFormatDeviations(float[] params) {
    return (params != null) ? mapGaussian2DFormat(params) : null;
  }

  private static float[] mapGaussian2DFormat(float[] params) {
    final float[] p = new float[nTwoAxisAndTheta];
    p[PeakResult.BACKGROUND] = params[LEGACY_FORMAT_BACKGROUND];
    p[PeakResult.INTENSITY] = params[LEGACY_FORMAT_SIGNAL];
    p[PeakResult.X] = params[LEGACY_FORMAT_X_POSITION];
    p[PeakResult.Y] = params[LEGACY_FORMAT_Y_POSITION];
    p[isx] = params[LEGACY_FORMAT_X_SD];
    p[isy] = params[LEGACY_FORMAT_Y_SD];
    p[ia] = params[LEGACY_FORMAT_ANGLE];
    return p;
  }

  private MemoryPeakResults createResults() {
    final MemoryPeakResults results = new MemoryPeakResults();
    results.setName(name);
    results.setSource(source);
    results.setBounds(bounds);
    if (calibration != null) {
      results.setCalibration(calibration.getCalibration());
    }
    results.setConfiguration(configuration);
    results.setPSF(psf);
    return results;
  }

  private final int readInt(byte[] buffer) {
    final int ch1 = buffer[position] & 0xff;
    final int ch2 = buffer[position + 1] & 0xff;
    final int ch3 = buffer[position + 2] & 0xff;
    final int ch4 = buffer[position + 3] & 0xff;
    position += 4;
    return ((ch1 << 24) + (ch2 << 16) + (ch3 << 8) + (ch4 << 0));
  }

  private final long readLong(byte[] buffer) {
    final long l = (((long) buffer[position + 0] << 56)
        + ((long) (buffer[position + 1] & 255) << 48) + ((long) (buffer[position + 2] & 255) << 40)
        + ((long) (buffer[position + 3] & 255) << 32) + ((long) (buffer[position + 4] & 255) << 24)
        + ((buffer[position + 5] & 255) << 16) + ((buffer[position + 6] & 255) << 8)
        + ((buffer[position + 7] & 255) << 0));
    position += 8;
    return l;
  }

  private final float readFloat(byte[] buffer) {
    return Float.intBitsToFloat(readInt(buffer));
  }

  private final double readDouble(byte[] buffer) {
    return Double.longBitsToDouble(readLong(buffer));
  }

  private void showProgress(FileChannel channel) throws IOException {
    if (tracker != null) {
      tracker.progress(channel.position(), channel.size());
      if (tracker.isEnded()) {
        // Throw an IOException and it will be caught and ignored by all the file reading methods
        throw new IOException("File read was cancelled");
      }
    }
  }

  private float[] readData(byte[] buffer, float[] params) {
    for (int i = 0; i < params.length; i++) {
      params[i] = readFloat(buffer);
    }
    return params;
  }

  private MemoryPeakResults readText() {
    final MemoryPeakResults results = createResults();

    // Units were added in version 3
    int nFields;
    if (smlmVersion < 3) {
      nFields = 7;
      calibration.setIntensityUnit(IntensityUnit.COUNT);
      calibration.setDistanceUnit(DistanceUnit.PIXEL);
      calibration.setAngleUnit(AngleUnit.DEGREE);

      // Note that in older versions the background included the bias.
      // The bias is not included in the table and so the user will have to
      // add this manually. They will also have to add the gain.
    } else {
      // The number of fields should be within the PSF object
      nFields = new PeakResultConversionHelper(null, psf).getNames().length;
    }

    try (FileInputStream fis = new FileInputStream(filename)) {
      try (BufferedReader input = new BufferedReader(new UnicodeReader(fis, null))) {
        @SuppressWarnings("resource")
        final FileChannel channel = fis.getChannel();

        String line;
        int errors = 0;

        final LineReader reader = createLineReader(results, smlmVersion, nFields);

        // Skip the header
        while ((line = input.readLine()) != null) {
          if (line.length() == 0) {
            continue;
          }

          if (line.charAt(0) != '#') {
            // This is the first record
            if (!reader.addPeakResult(line)) {
              errors = 1;
            }
            break;
          }
        }

        int c = 0;
        while ((line = input.readLine()) != null) {
          if (line.length() == 0) {
            continue;
          }
          if (line.charAt(0) == '#') {
            continue;
          }

          if (!reader.addPeakResult(line)) {
            if (++errors >= 10) {
              break;
            }
          }

          if (++c % 512 == 0) {
            showProgress(channel);
          }
        }
      }
    } catch (final IOException ex) {
      // ignore
    }
    return results;
  }

  /**
   * Simple class to call the appropriate method to parse the data.
   */
  private abstract class LineReader {
    final MemoryPeakResults results;

    LineReader(MemoryPeakResults results) {
      this.results = results;
    }

    boolean addPeakResult(String line) {
      final PeakResult result = read(line);
      if (result != null) {
        results.add(result);
        return true;
      }
      return false;
    }

    abstract PeakResult read(String line);
  }

  //@formatter:off
  private class LineReaderV4 extends LineReader {
    int nFields;
    LineReaderV4(MemoryPeakResults results, int nFields) { super(results); this.nFields=nFields; }
    @Override
    PeakResult read(String line) {
      return createPeakResultV4(line, nFields);
    }
  }
  private class LineReaderDV4 extends LineReader {
    int nFields;
    LineReaderDV4(MemoryPeakResults results, int nFields) { super(results); this.nFields=nFields; }
    @Override
    PeakResult read(String line) {
      return createPeakResultDeviationsV4(line, nFields);
    }
  }
  private class LineReaderV3 extends LineReader {
    int nFields;
    LineReaderV3(MemoryPeakResults results, int nFields) { super(results); this.nFields=nFields; }
    @Override
    PeakResult read(String line) {
      return createPeakResultV3(line, nFields);
    }
  }
  private class LineReaderDV3 extends LineReader {
    int nFields;
    LineReaderDV3(MemoryPeakResults results, int nFields) { super(results); this.nFields=nFields; }
    @Override
    PeakResult read(String line) {
      return createPeakResultDeviationsV3(line, nFields);
    }
  }
  private class LineReaderV2 extends LineReader {
    LineReaderV2(MemoryPeakResults results) { super(results); }
    @Override
    PeakResult read(String line) {
      return createPeakResultV2(line);
    }
  }
  private class LineReaderDV2 extends LineReader {
    LineReaderDV2(MemoryPeakResults results) { super(results); }
    @Override
    PeakResult read(String line) {
      return createPeakResultDeviationsV2(line);
    }
  }
  private class LineReaderV1 extends LineReader {
    LineReaderV1(MemoryPeakResults results) { super(results); }
    @Override
    PeakResult read(String line) {
      return createPeakResultV1(line);
    }
  }
  private class LineReaderDV1 extends LineReader {
    LineReaderDV1(MemoryPeakResults results) { super(results); }
    @Override
    PeakResult read(String line) {
      return createPeakResultDeviationsV1(line);
    }
  }
  //@formatter:on

  private LineReader createLineReader(MemoryPeakResults results, int version, int nFields) {
    switch (version) {
      case 4:
        return (deviations) ? new LineReaderDV4(results, nFields)
            : new LineReaderV4(results, nFields);

      case 3:
        return (deviations) ? new LineReaderDV3(results, nFields)
            : new LineReaderV3(results, nFields);

      case 2:
        return (deviations) ? new LineReaderDV2(results) : new LineReaderV2(results);

      case 1:
      default:
        return (deviations) ? new LineReaderDV1(results) : new LineReaderV1(results);
    }
  }

  private PeakResult createPeakResultV1(String line) {
    try {
      float[] params = new float[7];

      if (isUseScanner()) {
        // Code using a Scanner
        try (Scanner scanner = new Scanner(line)) {
          scanner.useDelimiter(tabPattern);
          scanner.useLocale(Locale.US);
          int id = 0;
          int endPeak = 0;
          if (readId) {
            id = scanner.nextInt();
          }
          final int peak = scanner.nextInt();
          if (readEndFrame) {
            endPeak = scanner.nextInt();
          }
          final int origX = scanner.nextInt();
          final int origY = scanner.nextInt();
          final float origValue = scanner.nextFloat();
          final double error = scanner.nextDouble();
          final float noise = scanner.nextFloat();
          final float signal = scanner.nextFloat(); // Ignored but must be read
          for (int i = 0; i < params.length; i++) {
            params[i] = scanner.nextFloat();
          }
          scanner.close();
          params[LEGACY_FORMAT_SIGNAL] = signal;
          params = mapGaussian2DFormatParams(params);
          if (readId || readEndFrame) {
            return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, 0, params,
                null, endPeak, id);
          }
          return new PeakResult(peak, origX, origY, origValue, error, noise, 0, params, null);
        }
      }
      // Code using split and parse
      final String[] fields = tabPattern.split(line);
      int j = 0;
      final int id = (readId) ? Integer.parseInt(fields[j++]) : 0;
      final int peak = Integer.parseInt(fields[j++]);
      final int endPeak = (readEndFrame) ? Integer.parseInt(fields[j++]) : 0;
      final int origX = Integer.parseInt(fields[j++]);
      final int origY = Integer.parseInt(fields[j++]);
      final float origValue = Float.parseFloat(fields[j++]);
      final double error = Double.parseDouble(fields[j++]);
      final float noise = Float.parseFloat(fields[j++]);
      final float signal = Float.parseFloat(fields[j++]);
      for (int i = 0; i < params.length; i++) {
        params[i] = Float.parseFloat(fields[j++]);
      }
      params[LEGACY_FORMAT_SIGNAL] = signal;
      params = mapGaussian2DFormatParams(params);
      if (readId || readEndFrame) {
        return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, 0, params, null,
            endPeak, id);
      }
      return new PeakResult(peak, origX, origY, origValue, error, noise, 0, params, null);
    } catch (final InputMismatchException ex) {
      // Ignore and return null
    } catch (final NoSuchElementException ex) {
      // Ignore and return null
    } catch (final IndexOutOfBoundsException ex) {
      // Ignore and return null
    } catch (final NumberFormatException ex) {
      // Ignore and return null
    }
    return null;
  }

  private PeakResult createPeakResultDeviationsV1(String line) {
    try {
      float[] params = new float[7];
      float[] paramsStdDev = new float[7];

      if (isUseScanner()) {
        // Code using a Scanner
        try (Scanner scanner = new Scanner(line)) {
          scanner.useDelimiter(tabPattern);
          scanner.useLocale(Locale.US);
          int id = 0;
          int endPeak = 0;
          if (readId) {
            id = scanner.nextInt();
          }
          final int peak = scanner.nextInt();
          if (readEndFrame) {
            endPeak = scanner.nextInt();
          }
          final int origX = scanner.nextInt();
          final int origY = scanner.nextInt();
          final float origValue = scanner.nextFloat();
          final double error = scanner.nextDouble();
          final float noise = scanner.nextFloat();
          final float signal = scanner.nextFloat(); // Ignored but must be read
          for (int i = 0; i < params.length; i++) {
            params[i] = scanner.nextFloat();
            paramsStdDev[i] = scanner.nextFloat();
          }
          params[LEGACY_FORMAT_SIGNAL] = signal;

          params = mapGaussian2DFormatParams(params);
          paramsStdDev = mapGaussian2DFormatDeviations(paramsStdDev);
          if (readId || readEndFrame) {
            return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, 0, params,
                paramsStdDev, endPeak, id);
          }
          return new PeakResult(peak, origX, origY, origValue, error, noise, 0, params,
              paramsStdDev);
        }
      }

      // Code using split and parse
      final String[] fields = tabPattern.split(line);
      int j = 0;
      final int id = (readId) ? Integer.parseInt(fields[j++]) : 0;
      final int peak = Integer.parseInt(fields[j++]);
      final int endPeak = (readEndFrame) ? Integer.parseInt(fields[j++]) : 0;
      final int origX = Integer.parseInt(fields[j++]);
      final int origY = Integer.parseInt(fields[j++]);
      final float origValue = Float.parseFloat(fields[j++]);
      final double error = Double.parseDouble(fields[j++]);
      final float noise = Float.parseFloat(fields[j++]);
      final float signal = Float.parseFloat(fields[j++]);
      for (int i = 0; i < params.length; i++) {
        params[i] = Float.parseFloat(fields[j++]);
        paramsStdDev[i] = Float.parseFloat(fields[j++]);
      }
      params[LEGACY_FORMAT_SIGNAL] = signal;
      params = mapGaussian2DFormatParams(params);
      paramsStdDev = mapGaussian2DFormatDeviations(paramsStdDev);
      if (readId || readEndFrame) {
        return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, 0, params,
            paramsStdDev, endPeak, id);
      }
      return new PeakResult(peak, origX, origY, origValue, error, noise, 0, params, paramsStdDev);
    } catch (final InputMismatchException ex) {
      // Ignore and return null
    } catch (final NoSuchElementException ex) {
      // Ignore and return null
    } catch (final IndexOutOfBoundsException ex) {
      // Ignore and return null
    } catch (final NumberFormatException ex) {
      // Ignore and return null
    }
    return null;
  }

  private PeakResult createPeakResultV2(String line) {
    try {
      float[] params = new float[7];

      if (isUseScanner()) {
        // Code using a Scanner
        try (Scanner scanner = new Scanner(line)) {
          scanner.useDelimiter(tabPattern);
          scanner.useLocale(Locale.US);
          int id = 0;
          int endPeak = 0;
          if (readId) {
            id = scanner.nextInt();
          }
          final int peak = scanner.nextInt();
          if (readEndFrame) {
            endPeak = scanner.nextInt();
          }
          final int origX = scanner.nextInt();
          final int origY = scanner.nextInt();
          final float origValue = scanner.nextFloat();
          final double error = scanner.nextDouble();
          final float noise = scanner.nextFloat();
          for (int i = 0; i < params.length; i++) {
            params[i] = scanner.nextFloat();
          }
          params = mapGaussian2DFormatParams(params);
          if (readId || readEndFrame) {
            return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, 0, params,
                null, endPeak, id);
          }
          return new PeakResult(peak, origX, origY, origValue, error, noise, 0, params, null);
        }
      }
      // Code using split and parse
      final String[] fields = tabPattern.split(line);
      int j = 0;
      final int id = (readId) ? Integer.parseInt(fields[j++]) : 0;
      final int peak = Integer.parseInt(fields[j++]);
      final int endPeak = (readEndFrame) ? Integer.parseInt(fields[j++]) : 0;
      final int origX = Integer.parseInt(fields[j++]);
      final int origY = Integer.parseInt(fields[j++]);
      final float origValue = Float.parseFloat(fields[j++]);
      final double error = Double.parseDouble(fields[j++]);
      final float noise = Float.parseFloat(fields[j++]);
      for (int i = 0; i < params.length; i++) {
        params[i] = Float.parseFloat(fields[j++]);
      }
      params = mapGaussian2DFormatParams(params);
      if (readId || readEndFrame) {
        return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, 0, params, null,
            endPeak, id);
      }
      return new PeakResult(peak, origX, origY, origValue, error, noise, 0, params, null);
    } catch (final InputMismatchException ex) {
      // Ignore and return null
    } catch (final NoSuchElementException ex) {
      // Ignore and return null
    } catch (final IndexOutOfBoundsException ex) {
      // Ignore and return null
    } catch (final NumberFormatException ex) {
      // Ignore and return null
    }
    return null;
  }

  private PeakResult createPeakResultDeviationsV2(String line) {
    try {
      float[] params = new float[7];
      float[] paramsStdDev = new float[7];

      if (isUseScanner()) {
        // Code using a Scanner
        try (Scanner scanner = new Scanner(line)) {
          scanner.useDelimiter(tabPattern);
          scanner.useLocale(Locale.US);
          int id = 0;
          int endPeak = 0;
          if (readId) {
            id = scanner.nextInt();
          }
          final int peak = scanner.nextInt();
          if (readEndFrame) {
            endPeak = scanner.nextInt();
          }
          final int origX = scanner.nextInt();
          final int origY = scanner.nextInt();
          final float origValue = scanner.nextFloat();
          final double error = scanner.nextDouble();
          final float noise = scanner.nextFloat();
          for (int i = 0; i < params.length; i++) {
            params[i] = scanner.nextFloat();
            paramsStdDev[i] = scanner.nextFloat();
          }
          params = mapGaussian2DFormatParams(params);
          paramsStdDev = mapGaussian2DFormatDeviations(paramsStdDev);
          if (readId || readEndFrame) {
            return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, 0, params,
                paramsStdDev, endPeak, id);
          }
          return new PeakResult(peak, origX, origY, origValue, error, noise, 0, params,
              paramsStdDev);
        }
      }

      // Code using split and parse
      final String[] fields = tabPattern.split(line);
      int j = 0;
      final int id = (readId) ? Integer.parseInt(fields[j++]) : 0;
      final int peak = Integer.parseInt(fields[j++]);
      final int endPeak = (readEndFrame) ? Integer.parseInt(fields[j++]) : 0;
      final int origX = Integer.parseInt(fields[j++]);
      final int origY = Integer.parseInt(fields[j++]);
      final float origValue = Float.parseFloat(fields[j++]);
      final double error = Double.parseDouble(fields[j++]);
      final float noise = Float.parseFloat(fields[j++]);
      for (int i = 0; i < params.length; i++) {
        params[i] = Float.parseFloat(fields[j++]);
        paramsStdDev[i] = Float.parseFloat(fields[j++]);
      }
      params = mapGaussian2DFormatParams(params);
      paramsStdDev = mapGaussian2DFormatDeviations(paramsStdDev);
      if (readId || readEndFrame) {
        return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, 0, params,
            paramsStdDev, endPeak, id);
      }
      return new PeakResult(peak, origX, origY, origValue, error, noise, 0, params, paramsStdDev);
    } catch (final InputMismatchException ex) {
      // Ignore and return null
    } catch (final NoSuchElementException ex) {
      // Ignore and return null
    } catch (final IndexOutOfBoundsException ex) {
      // Ignore and return null
    } catch (final NumberFormatException ex) {
      // Ignore and return null
    }
    return null;
  }

  private PeakResult createPeakResultV3(String line, int nFields) {
    try {
      final float[] params = new float[nFields];

      if (isUseScanner()) {
        // Code using a Scanner
        try (Scanner scanner = new Scanner(line)) {
          scanner.useDelimiter(tabPattern);
          scanner.useLocale(Locale.US);
          int id = 0;
          int endPeak = 0;
          if (readId) {
            id = scanner.nextInt();
          }
          final int peak = scanner.nextInt();
          if (readEndFrame) {
            endPeak = scanner.nextInt();
          }
          final int origX = scanner.nextInt();
          final int origY = scanner.nextInt();
          final float origValue = scanner.nextFloat();
          final double error = scanner.nextDouble();
          final float noise = scanner.nextFloat();
          for (int i = 0; i < params.length; i++) {
            params[i] = scanner.nextFloat();
          }
          // The format appends a * to computed precision. We ignore these.
          if (readPrecision && !line.endsWith("*")) {
            return createResult(peak, origX, origY, origValue, error, noise, 0, params, null,
                endPeak, id,
                // Read precision here because it is the final field
                scanner.nextFloat());
          }
          return createResult(peak, origX, origY, origValue, error, noise, 0, params, null, endPeak,
              id);
        }
      }
      // Code using split and parse
      final String[] fields = tabPattern.split(line);
      int j = 0;
      final int id = (readId) ? Integer.parseInt(fields[j++]) : 0;
      final int peak = Integer.parseInt(fields[j++]);
      final int endPeak = (readEndFrame) ? Integer.parseInt(fields[j++]) : 0;
      final int origX = Integer.parseInt(fields[j++]);
      final int origY = Integer.parseInt(fields[j++]);
      final float origValue = Float.parseFloat(fields[j++]);
      final double error = Double.parseDouble(fields[j++]);
      final float noise = Float.parseFloat(fields[j++]);
      for (int i = 0; i < params.length; i++) {
        params[i] = Float.parseFloat(fields[j++]);
      }
      // The format appends a * to computed precision. We ignore these.
      if (readPrecision && !line.endsWith("*")) {
        return createResult(peak, origX, origY, origValue, error, noise, 0, params, null, endPeak,
            id,
            // Read precision here because it is the final field
            Float.parseFloat(fields[j]));
      }
      return createResult(peak, origX, origY, origValue, error, noise, 0, params, null, endPeak,
          id);
    } catch (final InputMismatchException ex) {
      // Ignore and return null
    } catch (final NoSuchElementException ex) {
      // Ignore and return null
    } catch (final IndexOutOfBoundsException ex) {
      // Ignore and return null
    } catch (final NumberFormatException ex) {
      // Ignore and return null
    }
    return null;
  }

  private PeakResult createPeakResultDeviationsV3(String line, int nFields) {
    try {
      final float[] params = new float[nFields];
      final float[] paramsStdDev = new float[nFields];

      if (isUseScanner()) {
        // Code using a Scanner
        try (Scanner scanner = new Scanner(line)) {
          scanner.useDelimiter(tabPattern);
          scanner.useLocale(Locale.US);
          int id = 0;
          int endPeak = 0;
          if (readId) {
            id = scanner.nextInt();
          }
          final int peak = scanner.nextInt();
          if (readEndFrame) {
            endPeak = scanner.nextInt();
          }
          final int origX = scanner.nextInt();
          final int origY = scanner.nextInt();
          final float origValue = scanner.nextFloat();
          final double error = scanner.nextDouble();
          final float noise = scanner.nextFloat();
          for (int i = 0; i < params.length; i++) {
            params[i] = scanner.nextFloat();
            paramsStdDev[i] = scanner.nextFloat();
          }
          if (readPrecision && !line.endsWith("*")) {
            return createResult(peak, origX, origY, origValue, error, noise, 0, params,
                paramsStdDev, endPeak, id,
                // Read precision here because it is the final field
                scanner.nextFloat());
          }
          return createResult(peak, origX, origY, origValue, error, noise, 0, params, paramsStdDev,
              endPeak, id);
        }
      }

      // Code using split and parse
      final String[] fields = tabPattern.split(line);
      int j = 0;
      final int id = (readId) ? Integer.parseInt(fields[j++]) : 0;
      final int peak = Integer.parseInt(fields[j++]);
      final int endPeak = (readEndFrame) ? Integer.parseInt(fields[j++]) : 0;
      final int origX = Integer.parseInt(fields[j++]);
      final int origY = Integer.parseInt(fields[j++]);
      final float origValue = Float.parseFloat(fields[j++]);
      final double error = Double.parseDouble(fields[j++]);
      final float noise = Float.parseFloat(fields[j++]);
      for (int i = 0; i < params.length; i++) {
        params[i] = Float.parseFloat(fields[j++]);
        paramsStdDev[i] = Float.parseFloat(fields[j++]);
      }
      if (readPrecision && !line.endsWith("*")) {
        return createResult(peak, origX, origY, origValue, error, noise, 0, params, paramsStdDev,
            endPeak, id,
            // Read precision here because it is the final field
            Float.parseFloat(fields[j]));
      }
      return createResult(peak, origX, origY, origValue, error, noise, 0, params, paramsStdDev,
          endPeak, id);
    } catch (final InputMismatchException ex) {
      // Ignore and return null
    } catch (final NoSuchElementException ex) {
      // Ignore and return null
    } catch (final IndexOutOfBoundsException ex) {
      // Ignore and return null
    } catch (final NumberFormatException ex) {
      // Ignore and return null
    }
    return null;
  }

  private PeakResult createPeakResultV4(String line, int nFields) {
    try {
      final float[] params = new float[nFields];

      if (isUseScanner()) {
        // Code using a Scanner
        try (Scanner scanner = new Scanner(line)) {
          scanner.useDelimiter(tabPattern);
          scanner.useLocale(Locale.US);
          int id = 0;
          int endPeak = 0;
          if (readId) {
            id = scanner.nextInt();
          }
          final int peak = scanner.nextInt();
          if (readEndFrame) {
            endPeak = scanner.nextInt();
          }
          final int origX = scanner.nextInt();
          final int origY = scanner.nextInt();
          final float origValue = scanner.nextFloat();
          final double error = scanner.nextDouble();
          final float noise = scanner.nextFloat();
          final float meanIntensity = scanner.nextFloat();
          for (int i = 0; i < params.length; i++) {
            params[i] = scanner.nextFloat();
          }
          // The format appends a * to computed precision. We ignore these.
          if (readPrecision && !line.endsWith("*")) {
            return createResult(peak, origX, origY, origValue, error, noise, meanIntensity, params,
                null, endPeak, id,
                // Read precision here because it is the final field
                scanner.nextFloat());
          }
          return createResult(peak, origX, origY, origValue, error, noise, meanIntensity, params,
              null, endPeak, id);
        }
      }
      // Code using split and parse
      final String[] fields = tabPattern.split(line);
      int j = 0;
      final int id = (readId) ? Integer.parseInt(fields[j++]) : 0;
      final int peak = Integer.parseInt(fields[j++]);
      final int endPeak = (readEndFrame) ? Integer.parseInt(fields[j++]) : 0;
      final int origX = Integer.parseInt(fields[j++]);
      final int origY = Integer.parseInt(fields[j++]);
      final float origValue = Float.parseFloat(fields[j++]);
      final double error = Double.parseDouble(fields[j++]);
      final float noise = Float.parseFloat(fields[j++]);
      final float meanIntensity = Float.parseFloat(fields[j++]);
      for (int i = 0; i < params.length; i++) {
        params[i] = Float.parseFloat(fields[j++]);
      }
      // The format appends a * to computed precision. We ignore these.
      if (readPrecision && !line.endsWith("*")) {
        return createResult(peak, origX, origY, origValue, error, noise, meanIntensity, params,
            null, endPeak, id,
            // Read precision here because it is the final field
            Float.parseFloat(fields[j]));
      }
      return createResult(peak, origX, origY, origValue, error, noise, meanIntensity, params, null,
          endPeak, id);
    } catch (final InputMismatchException ex) {
      // Ignore and return null
    } catch (final NoSuchElementException ex) {
      // Ignore and return null
    } catch (final IndexOutOfBoundsException ex) {
      // Ignore and return null
    } catch (final NumberFormatException ex) {
      // Ignore and return null
    }
    return null;
  }

  private PeakResult createPeakResultDeviationsV4(String line, int nFields) {
    try {
      final float[] params = new float[nFields];
      final float[] paramsStdDev = new float[nFields];

      if (isUseScanner()) {
        // Code using a Scanner
        try (Scanner scanner = new Scanner(line)) {
          scanner.useDelimiter(tabPattern);
          scanner.useLocale(Locale.US);
          int id = 0;
          int endPeak = 0;
          if (readId) {
            id = scanner.nextInt();
          }
          final int peak = scanner.nextInt();
          if (readEndFrame) {
            endPeak = scanner.nextInt();
          }
          final int origX = scanner.nextInt();
          final int origY = scanner.nextInt();
          final float origValue = scanner.nextFloat();
          final double error = scanner.nextDouble();
          final float noise = scanner.nextFloat();
          final float meanIntensity = scanner.nextFloat();
          for (int i = 0; i < params.length; i++) {
            params[i] = scanner.nextFloat();
            paramsStdDev[i] = scanner.nextFloat();
          }
          if (readPrecision && !line.endsWith("*")) {
            return createResult(peak, origX, origY, origValue, error, noise, meanIntensity, params,
                paramsStdDev, endPeak, id,
                // Read precision here because it is the final field
                scanner.nextFloat());
          }
          return createResult(peak, origX, origY, origValue, error, noise, meanIntensity, params,
              paramsStdDev, endPeak, id);
        }
      }

      // Code using split and parse
      final String[] fields = tabPattern.split(line);
      int j = 0;
      final int id = (readId) ? Integer.parseInt(fields[j++]) : 0;
      final int peak = Integer.parseInt(fields[j++]);
      final int endPeak = (readEndFrame) ? Integer.parseInt(fields[j++]) : 0;
      final int origX = Integer.parseInt(fields[j++]);
      final int origY = Integer.parseInt(fields[j++]);
      final float origValue = Float.parseFloat(fields[j++]);
      final double error = Double.parseDouble(fields[j++]);
      final float noise = Float.parseFloat(fields[j++]);
      final float meanIntensity = Float.parseFloat(fields[j++]);
      for (int i = 0; i < params.length; i++) {
        params[i] = Float.parseFloat(fields[j++]);
        paramsStdDev[i] = Float.parseFloat(fields[j++]);
      }
      if (readPrecision && !line.endsWith("*")) {
        return createResult(peak, origX, origY, origValue, error, noise, meanIntensity, params,
            paramsStdDev, endPeak, id,
            // Read precision here because it is the final field
            Float.parseFloat(fields[j]));
      }
      return createResult(peak, origX, origY, origValue, error, noise, meanIntensity, params,
          paramsStdDev, endPeak, id);
    } catch (final InputMismatchException ex) {
      // Ignore and return null
    } catch (final NoSuchElementException ex) {
      // Ignore and return null
    } catch (final IndexOutOfBoundsException ex) {
      // Ignore and return null
    } catch (final NumberFormatException ex) {
      // Ignore and return null
    }
    return null;
  }

  private MemoryPeakResults readTable() {
    final MemoryPeakResults results = createResults();
    results.setName(new File(filename).getName());

    try (FileInputStream fis = new FileInputStream(filename)) {
      try (BufferedReader input = new BufferedReader(new UnicodeReader(fis, null))) {
        String line;
        int errors = 0;

        // Skip over the single line header
        final String header = input.readLine();

        // V1: had the Signal and Amplitude. Parameters
        // V2: have only the Signal.
        // V3: Has variable columns with units for the PSF parameters. Signal was renamed to
        // Intensity.
        int version;
        final int nFields = 0;
        if (header.contains("Signal")) {
          version = 1;
        } else if (header.contains("Amplitude")) {
          version = 2;
        } else {
          version = 3;

          // We support reading old IJ table results as they had fixed columns.
          // The latest table results have dynamic columns so these must be loaded manually
          // as guessing the column format is not supported.
          return null;

          // This code functioned when the table was not dynamic ...

          // // Get the number of data fields by counting the standard fields
          // String[] columns = header.split("\t");
          // int field = 0;
          // if (readId)
          // field++; // ID #
          // if (readSource)
          // field++; // Source
          // field++; // Frame
          // if (readEndFrame)
          // field++; // End frame
          // field++; // origX
          // field++; // origY
          // field++; // origValue
          // field++; // error
          // field++; // noise
          // field++; // SNR
          //
          // // The remaining fields are PSF parameters with the exception of the final precision
          // field
          //
          // nFields = columns.length - field;
          // if (columns[columns.length - 1].contains("Precision"))
          // nFields--;
          // if (deviations)
          // {
          // nFields /= 2;
          // }
          //
          // // We can guess part of the calibration.
          // if (calibration == null)
          // calibration = new Calibration();
          // int jump = (deviations) ? 2 : 1;
          // // field is currently on Background
          // calibration.setIntensityUnit(UnitHelper.guessIntensityUnitFromShortName(extractUnit(columns[field])));
          // field += jump; // Move to Intensity
          // field += jump; // Move to X
          // calibration.setDistanceUnit(UnitHelper.guessDistanceUnitFromShortName(extractUnit(columns[field])));
          // field += jump; // Move to Y
          // field += jump; // Move to Z
          // // The angle may be used in fields above the standard ones
          // while (field < columns.length)
          // {
          // field += jump;
          // AngleUnit u = UnitHelper.guessAngleUnitFromShortName(extractUnit(columns[field]));
          // if (u != null)
          // {
          // calibration.setAngleUnit(u);
          // break;
          // }
          // }
        }

        @SuppressWarnings("resource")
        final FileChannel channel = fis.getChannel();

        int c = 0;
        while ((line = input.readLine()) != null) {
          if (line.length() == 0) {
            continue;
          }

          if (!addTableResult(results, line, version, nFields)) {
            if (++errors >= 10) {
              break;
            }
          }

          if (++c % 512 == 0) {
            showProgress(channel);
          }
        }
      }
    } catch (final IOException ex) {
      // ignore
    }

    return results;
  }

  /**
   * Extract the unit from a string. The unit is the first string within brackets, e.g. (unit) would
   * return unit.
   *
   * @param string the string
   * @return the unit string (or null)
   */
  public static String extractUnit(String string) {
    if (string != null) {
      final int beginIndex = string.indexOf('(');
      if (beginIndex >= 0) {
        final int endIndex = string.indexOf(')');
        if (endIndex > beginIndex) {
          return string.substring(beginIndex, endIndex);
        }
      }
    }
    return null;
  }

  private boolean addTableResult(MemoryPeakResults results, String line, int version, int nFields) {
    final PeakResult result;
    switch (version) {
      case 3:
        result = createTableResultV3(line, nFields);
        break;

      case 2:
        result = createTableResultV2(line);
        break;

      case 1:
      default:
        result = createTableResultV1(line);
    }
    if (result != null) {
      // Extract the source & bounds from the Source column
      if (results.size() == 0 && readSource) {
        try (Scanner scanner = new Scanner(line)) {
          scanner.useDelimiter(tabPattern);
          scanner.useLocale(Locale.US);
          if (readId) {
            scanner.nextInt();
          }
          final String source = scanner.next();

          if (source.contains(": ")) {
            final String[] fields = source.split(": ");
            results.setName(fields[0]);
            // Bounds is formatted as 'xN yN wN hN'
            final Pattern pattern = Pattern.compile("x(\\d+) y(\\d+) w(\\d+) h(\\d+)");
            final Matcher match = pattern.matcher(fields[1]);
            if (match.find()) {
              final int x = Integer.parseInt(match.group(1));
              final int y = Integer.parseInt(match.group(2));
              final int w = Integer.parseInt(match.group(3));
              final int h = Integer.parseInt(match.group(4));
              results.setBounds(new Rectangle(x, y, w, h));
            }
          } else {
            results.setName(source);
          }
        }
      }

      results.add(result);
      return true;
    }
    return false;
  }

  private PeakResult createTableResultV1(String line) {
    // Text file with fields:
    // [#]
    // [Source]
    // Peak
    // [End Peak]
    // origX
    // origY
    // origValue
    // Error
    // Noise
    // Signal
    // SNR
    // Background
    // [+/-]
    // Amplitude
    // [+/-]
    // Angle
    // [+/-]
    // X
    // [+/-]
    // Y
    // [+/-]
    // X SD
    // [+/-]
    // Y SD
    // [+/-]
    // [Precision]
    try {
      try (Scanner scanner = new Scanner(line)) {
        scanner.useDelimiter(tabPattern);
        scanner.useLocale(Locale.US);
        int id = 0;
        int endPeak = 0;
        if (readId) {
          id = scanner.nextInt();
        }
        if (readSource) {
          scanner.next();
        }
        final int peak = scanner.nextInt();
        if (readEndFrame) {
          endPeak = scanner.nextInt();
        }
        final int origX = scanner.nextInt();
        final int origY = scanner.nextInt();
        final float origValue = scanner.nextFloat();
        final double error = scanner.nextDouble();
        final float noise = scanner.nextFloat();
        @SuppressWarnings("unused")
        final float signal = scanner.nextFloat(); // Ignored but must be read
        @SuppressWarnings("unused")
        final float snr = scanner.nextFloat(); // Ignored but must be read
        float[] params = new float[7];
        float[] paramsStdDev = null;
        if (deviations) {
          paramsStdDev = new float[7];
          for (int i = 0; i < params.length; i++) {
            params[i] = scanner.nextFloat();
            paramsStdDev[i] = scanner.nextFloat();
          }
        } else {
          for (int i = 0; i < params.length; i++) {
            params[i] = scanner.nextFloat();
          }
        }

        // For the new format we store the signal (not the amplitude).
        // Convert the amplitude into a signal
        params[LEGACY_FORMAT_SIGNAL] *=
            2 * Math.PI * params[LEGACY_FORMAT_X_SD] * params[LEGACY_FORMAT_Y_SD];

        params = mapGaussian2DFormatParams(params);
        paramsStdDev = mapGaussian2DFormatDeviations(paramsStdDev);

        if (readId || readEndFrame) {
          return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, 0, params,
              paramsStdDev, endPeak, id);
        }
        return new PeakResult(peak, origX, origY, origValue, error, noise, 0, params, paramsStdDev);
      }
    } catch (final InputMismatchException ex) {
      // Ignore and return null
    } catch (final NoSuchElementException ex) {
      // Ignore and return null
    }
    return null;
  }

  private PeakResult createTableResultV2(String line) {
    // Text file with fields:
    // [#]
    // [Source]
    // Peak
    // [End Peak]
    // origX
    // origY
    // origValue
    // Error
    // Noise
    // SNR
    // Background
    // [+/-]
    // Signal
    // [+/-]
    // Angle
    // [+/-]
    // X
    // [+/-]
    // Y
    // [+/-]
    // X SD
    // [+/-]
    // Y SD
    // [+/-]
    // [Precision]
    try {
      try (Scanner scanner = new Scanner(line)) {
        scanner.useDelimiter(tabPattern);
        scanner.useLocale(Locale.US);
        int id = 0;
        int endPeak = 0;
        if (readId) {
          id = scanner.nextInt();
        }
        if (readSource) {
          scanner.next();
        }
        final int peak = scanner.nextInt();
        if (readEndFrame) {
          endPeak = scanner.nextInt();
        }
        final int origX = scanner.nextInt();
        final int origY = scanner.nextInt();
        final float origValue = scanner.nextFloat();
        final double error = scanner.nextDouble();
        final float noise = scanner.nextFloat();
        @SuppressWarnings("unused")
        final float snr = scanner.nextFloat(); // Ignored but must be read
        float[] params = new float[7];
        float[] paramsStdDev = null;
        if (deviations) {
          paramsStdDev = new float[7];
          for (int i = 0; i < params.length; i++) {
            params[i] = scanner.nextFloat();
            paramsStdDev[i] = scanner.nextFloat();
          }
        } else {
          for (int i = 0; i < params.length; i++) {
            params[i] = scanner.nextFloat();
          }
        }

        params = mapGaussian2DFormatParams(params);
        paramsStdDev = mapGaussian2DFormatDeviations(paramsStdDev);

        if (readId || readEndFrame) {
          return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, 0, params,
              paramsStdDev, endPeak, id);
        }
        return new PeakResult(peak, origX, origY, origValue, error, noise, 0, params, paramsStdDev);
      }
    } catch (final InputMismatchException ex) {
      // Ignore and return null
    } catch (final NoSuchElementException ex) {
      // Ignore and return null
    }
    return null;
  }

  private PeakResult createTableResultV3(String line, int nFields) {
    // Text file with fields:
    // [#]
    // [Source]
    // Peak
    // [End Peak]
    // origX
    // origY
    // origValue
    // Error
    // Noise[ (units)]
    // SNR
    // Background[ (units)]
    // [+/-]
    // Intensity[ (units)]
    // [+/-]
    // X[ (units)]
    // [+/-]
    // Y[ (units)]
    // [+/-]
    // Z[ (units)]
    // [+/-]
    // Repeated:
    // Field[ (units)]
    // [+/-]
    // [Precision]
    try {
      try (Scanner scanner = new Scanner(line)) {
        scanner.useDelimiter(tabPattern);
        scanner.useLocale(Locale.US);
        int id = 0;
        int endPeak = 0;
        if (readId) {
          id = scanner.nextInt();
        }
        if (readSource) {
          scanner.next();
        }
        final int peak = scanner.nextInt();
        if (readEndFrame) {
          endPeak = scanner.nextInt();
        }
        final int origX = scanner.nextInt();
        final int origY = scanner.nextInt();
        final float origValue = scanner.nextFloat();
        final double error = scanner.nextDouble();
        final float noise = scanner.nextFloat();
        @SuppressWarnings("unused")
        final float snr = scanner.nextFloat(); // Ignored but must be read
        final float[] params = new float[nFields];
        float[] paramsStdDev;
        if (deviations) {
          paramsStdDev = new float[params.length];
          for (int i = 0; i < params.length; i++) {
            params[i] = scanner.nextFloat();
            paramsStdDev[i] = scanner.nextFloat();
          }
        } else {
          paramsStdDev = null;
          for (int i = 0; i < params.length; i++) {
            params[i] = scanner.nextFloat();
          }
        }
        if (readId || readEndFrame) {
          return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, 0, params,
              paramsStdDev, endPeak, id);
        }
        return new PeakResult(peak, origX, origY, origValue, error, noise, 0, params, paramsStdDev);
      }
    } catch (final InputMismatchException ex) {
      // Ignore and return null
    } catch (final NoSuchElementException ex) {
      // Ignore and return null
    }
    return null;
  }

  private MemoryPeakResults readRapidSTORM() {
    final MemoryPeakResults results = createResults();
    results.setName(new File(filename).getName());

    try (FileInputStream fis = new FileInputStream(filename)) {
      try (BufferedReader input = new BufferedReader(new UnicodeReader(fis, null))) {
        @SuppressWarnings("resource")
        final FileChannel channel = fis.getChannel();

        String line;
        int errors = 0;

        // Skip the header
        while ((line = input.readLine()) != null) {
          if (line.length() == 0) {
            continue;
          }
          if (line.charAt(0) != '#') {
            // This is the first record
            if (!addRapidSTORMResult(results, line)) {
              errors = 1;
            }
            break;
          }
        }

        int c = 0;
        while ((line = input.readLine()) != null) {
          if (line.length() == 0) {
            continue;
          }
          if (line.charAt(0) == '#') {
            continue;
          }

          if (!addRapidSTORMResult(results, line)) {
            if (++errors >= 10) {
              break;
            }
          }

          if (++c % 512 == 0) {
            showProgress(channel);
          }
        }
      }
    } catch (final IOException ex) {
      // ignore
    }
    return results;
  }

  private static boolean addRapidSTORMResult(MemoryPeakResults results, String line) {
    final PeakResult result = createRapidSTORMResult(line);
    if (result != null) {
      results.add(result);
      return true;
    }
    return false;
  }

  private static PeakResult createRapidSTORMResult(String line) {
    // Text file with fields:
    // X (nm)
    // Y (nm)
    // Frame
    // Amplitude*
    // sx^2 (pm^2)
    // sy^2 (pm^2)
    // 2 kernel improvement
    // Fit residues chi square
    // *Note that the RapidSTORM Amplitude is the signal.
    try {
      try (Scanner scanner = new Scanner(line)) {
        scanner.useDelimiter(spacePattern);
        scanner.useLocale(Locale.US);
        final float x = scanner.nextFloat();
        final float y = scanner.nextFloat();
        final int peak = scanner.nextInt();
        final float signal = scanner.nextFloat();
        final float sx2 = scanner.nextFloat();
        final float sy2 = scanner.nextFloat();
        @SuppressWarnings("unused")
        final float kernelImprovement = scanner.nextFloat();
        final double error = scanner.nextDouble();

        // Convert from pm^2 to nm
        final float sx = (float) (Math.sqrt(sx2) * 1000);
        final float sy = (float) (Math.sqrt(sy2) * 1000);

        final float[] params = new float[nTwoAxis];
        params[PeakResult.INTENSITY] = signal;
        params[PeakResult.X] = x;
        params[PeakResult.Y] = y;
        params[isx] = sx;
        params[isy] = sy;

        // Store the signal as the original value
        return new PeakResult(peak, (int) x, (int) y, signal, error, 0.0f, 0, params, null);
      }
    } catch (final InputMismatchException ex) {
      // Ignore and return null
    } catch (final NoSuchElementException ex) {
      // Ignore and return null
    }
    return null;
  }

  private MemoryPeakResults readNSTORM() {
    final MemoryPeakResults results = createResults();
    results.setName(new File(filename).getName());

    try (FileInputStream fis = new FileInputStream(filename)) {
      try (BufferedReader input = new BufferedReader(new UnicodeReader(fis, null))) {
        @SuppressWarnings("resource")
        final FileChannel channel = fis.getChannel();

        String line;
        int errors = 0;

        // Skip the single line header
        input.readLine();

        int c = 0;
        while ((line = input.readLine()) != null) {
          if (line.length() == 0) {
            continue;
          }

          if (!addNSTORMResult(results, line)) {
            if (++errors >= 10) {
              break;
            }
          }

          if (++c % 512 == 0) {
            showProgress(channel);
          }
        }
      }
    } catch (final IOException ex) {
      // ignore
    }

    // The following relationship holds when length == 1:
    // intensity = height * 2 * pi * sd0 * sd1 / pixel_pitch^2
    // => Pixel_pitch = sqrt(height * 2 * pi * sd0 * sd1 / intensity)
    // Try and create a calibration
    final Statistics pixelPitch = new Statistics();
    results.forEach(new PeakResultProcedureX() {
      final double twoPi = 2 * Math.PI;

      @Override
      public boolean execute(PeakResult p) {
        if (p.getFrame() == p.getEndFrame()) {
          final float height = p.getOrigValue();
          final float intensity = p.getParameter(PeakResult.INTENSITY);
          final float sd0 = p.getParameter(isx);
          final float sd1 = p.getParameter(isy);
          pixelPitch.add(Math.sqrt(height * twoPi * sd0 * sd1 / intensity));
          // Stop when we have enough for a good guess
          return (pixelPitch.getN() > 100);
        }
        return false;
      }
    });

    // TODO - Support all the NSTORM formats: one-axis, two-axis, rotated, 3D.
    // Is this information in the header?
    // We could support setting the PSF as a Gaussian2D with one/two axis SD.
    // This would mean updating all the result params if it is a one axis PSF.
    // For now just record it as a 2 axis PSF.

    // Create a calibration
    calibration = new CalibrationWriter();

    // Q. Is NSTORM in photons?
    calibration.setIntensityUnit(IntensityUnit.COUNT);
    calibration.setDistanceUnit(DistanceUnit.NM);

    if (pixelPitch.getN() > 0) {
      final double nmPerPixel = pixelPitch.getMean();
      calibration.setNmPerPixel(nmPerPixel);
    }

    results.setCalibration(getCalibration());

    return results;
  }

  private static boolean addNSTORMResult(MemoryPeakResults results, String line) {
    final PeakResult result = createNSTORMResult(line);
    if (result != null) {
      results.add(result);
      return true;
    }
    return false;
  }

  // So that the fields can be named
  @SuppressWarnings("unused")
  private static PeakResult createNSTORMResult(String line) {
    // Note that the NSTORM file contains traced molecules hence the Frame
    // and Length fields.

    // Text file with tab delimited fields:
    // Channel Name: This is the name of the Channel where the molecule was
    // detected.
    // X: The X position of the centroid of the molecule in
    // nanometers. Similar to the conventional image, molecules
    // positions in the image are relative to the upper left corner of
    // the image.
    // Y: The Y position of the centroid of the molecule in
    // nanometers. Similar to the conventional image, molecules
    // positions in the image are relative to the upper left corner of
    // the image.
    // Xc: The X position of the centroid of the molecule (in
    // nanometers) with drift correction applied. If no drift
    // correction was applied to this data then Xc= X.
    // Yc: The Y position of the centroid of the molecule (in
    // nanometers) with drift correction applied. If no drift
    // correction was applied to this data then Yc= Y.
    // Height: Intensity of the peak height in the detection frame (after the
    // detection process)
    // Area: Volume under the peak. Units are intensity * pixel^2
    // Width: Square Root(Wx*Wy)
    // Phi: The Angle of the molecule. This is the axial angle for 2D
    // and distance from Z calibration curve in nm in Wx,Wy
    // space for 3D)
    // Ax: Axial ratio of Wy/Wx
    // BG: The local background for the molecule
    // I: Accumulated intensity
    // Frame: The sequential frame number where the molecule was first
    // detected.
    // Length: The number of consecutive frames the molecule was
    // detected.
    // Link: For Internal Use only
    // Valid: For Internal Use only
    // Z: The Z coordinate of the molecule in Nanometers (origin is
    // the cover glass).
    // Zc: The Z position of the molecule (in nanometers) with drift
    // correction applied. If no drift correction was applied to this
    // data then Zc= Z.
    try {
      try (Scanner scanner = new Scanner(line)) {
        scanner.useDelimiter(tabPattern);
        scanner.useLocale(Locale.US);
        final String channelName = scanner.next();
        final float x = scanner.nextFloat();
        final float y = scanner.nextFloat();
        final float xc = scanner.nextFloat();
        final float yc = scanner.nextFloat();
        final float height = scanner.nextFloat();
        final float area = scanner.nextFloat();
        float width = scanner.nextFloat();
        final float phi = scanner.nextFloat();
        float ax = scanner.nextFloat();
        final float bg = scanner.nextFloat();
        final float i = scanner.nextFloat();
        final int frame = scanner.nextInt();
        final int length = scanner.nextInt();
        // These are not needed
        // float link = scanner.nextFloat();
        // float valid = scanner.nextFloat();
        // float z = scanner.nextFloat();
        // float zc = scanner.nextFloat();

        // The coordinates are in nm
        // The values are in ADUs. The area value is the signal.

        // The following relationship holds when length == 1:
        // Area = Height * 2 * pi * (Width / (pixel_pitch*2) )^2
        // => Pixel_pitch = 0.5 * Width / sqrt(Area / (Height * 2 * pi))

        final float[] params = new float[nTwoAxis];
        params[PeakResult.BACKGROUND] = bg;
        // params[ia] = ax;
        params[PeakResult.INTENSITY] = area;
        params[PeakResult.X] = xc;
        params[PeakResult.Y] = yc;

        // Convert width (2*SD) to SD
        width /= 2f;

        // Convert to separate XY widths using the axial ratio
        if (ax == 1) {
          params[isx] = width;
          params[isy] = width;
        } else {
          // Ensure the axial ratio is long/short
          if (ax < 1) {
            ax = 1.0f / ax;
          }
          final double a = Math.sqrt(ax);

          params[isx] = (float) (width * a);
          params[isy] = (float) (width / a);
        }

        // Store the signal as the original value
        return new ExtendedPeakResult(frame, (int) xc, (int) yc, height, 0.0, 0.0f, 0, params, null,
            frame + length - 1, 0);
      }
    } catch (final InputMismatchException ex) {
      // Ignore
    } catch (final NoSuchElementException ex) {
      // Ignore
    }
    return null;
  }

  private MemoryPeakResults readMALK() {
    final MemoryPeakResults results = createResults();
    if (TextUtils.isNullOrEmpty(name)) {
      results.setName(new File(filename).getName());
    }

    try (FileInputStream fis = new FileInputStream(filename)) {
      try (BufferedReader input = new BufferedReader(new UnicodeReader(fis, null))) {
        @SuppressWarnings("resource")
        final FileChannel channel = fis.getChannel();

        String line;
        int errors = 0;

        // Skip the header
        while ((line = input.readLine()) != null) {
          if (line.length() == 0) {
            continue;
          }

          if (line.charAt(0) != '#') {
            // This is the first record
            if (!addMALKResult(results, line)) {
              errors = 1;
            }
            break;
          }
        }

        int c = 0;
        while ((line = input.readLine()) != null) {
          if (line.length() == 0) {
            continue;
          }
          if (line.charAt(0) == '#') {
            continue;
          }

          if (!addMALKResult(results, line)) {
            if (++errors >= 10) {
              break;
            }
          }

          if (++c % 512 == 0) {
            showProgress(channel);
          }
        }
      }
    } catch (final IOException ex) {
      // ignore
    }

    // Set default calibration for MALK format.
    // The calibration may not be null if this was a GDSC MALK file since that has a header.
    if (calibration == null) {
      calibration = new CalibrationWriter();
      // Default assumption is nm
      calibration.setDistanceUnit(DistanceUnit.NM);
      // MALK uses photons
      calibration.setIntensityUnit(IntensityUnit.PHOTON);

      results.setCalibration(getCalibration());
    }

    return results;
  }

  private boolean addMALKResult(MemoryPeakResults results, String line) {
    final PeakResult result = createMALKResult(line);
    if (result != null) {
      results.add(result);
      return true;
    }
    return false;
  }

  private boolean isMALKFormat(String firstLine) {
    // The MALK file format is very simple: X,Y,T,Signal
    final String[] fields = whitespacePattern.split(firstLine);
    if (fields.length != 4) {
      return false;
    }
    if (createMALKResult(firstLine) == null) {
      return false;
    }
    return true;
  }

  private PeakResult createMALKResult(String line) {
    try {
      final float[] params = new float[PeakResult.STANDARD_PARAMETERS];

      if (isUseScanner()) {
        // Code using a Scanner
        try (Scanner scanner = new Scanner(line)) {
          scanner.useDelimiter(whitespacePattern);
          scanner.useLocale(Locale.US);
          params[PeakResult.X] = scanner.nextFloat();
          params[PeakResult.Y] = scanner.nextFloat();
          final int peak = scanner.nextInt();
          params[PeakResult.INTENSITY] = scanner.nextFloat();

          return new PeakResult(peak, 0, 0, 0, 0, 0, 0, params, null);
        }
      }
      // Code using split and parse
      final String[] fields = whitespacePattern.split(line);

      params[PeakResult.X] = Float.parseFloat(fields[0]);
      params[PeakResult.Y] = Float.parseFloat(fields[1]);
      final int peak = Integer.parseInt(fields[2]);
      params[PeakResult.INTENSITY] = Float.parseFloat(fields[3]);

      return new PeakResult(peak, 0, 0, 0, 0, 0, 0, params, null);
    } catch (final InputMismatchException ex) {
      // Ignore and return null
    } catch (final NoSuchElementException ex) {
      // Ignore and return null
    } catch (final IndexOutOfBoundsException ex) {
      // Ignore and return null
    } catch (final NumberFormatException ex) {
      // Ignore and return null
    }
    return null;
  }

  private MemoryPeakResults readTSF() {
    final TSFPeakResultsReader reader = new TSFPeakResultsReader(filename);
    reader.setOptions(options);
    return reader.read();
  }

  /**
   * Gets the tracker.
   *
   * @return the tracker
   */
  public TrackProgress getTracker() {
    return tracker;
  }

  /**
   * Sets the tracker.
   *
   * @param tracker the tracker to set
   */
  public void setTracker(TrackProgress tracker) {
    this.tracker = tracker;
  }

  /**
   * Checks if using a {@link Scanner} to read the text. The default is
   * {@link Pattern#split(CharSequence)}.
   *
   * @return true, if using a scanner
   */
  public boolean isUseScanner() {
    return useScanner;
  }

  /**
   * Set to true to use a {@link Scanner} to read the text. The default is
   * {@link Pattern#split(CharSequence)}.
   *
   * @param useScanner Set to true to use a {@link Scanner} to read the text
   */
  public void setUseScanner(boolean useScanner) {
    this.useScanner = useScanner;
  }

  /**
   * Checks if returning raw results (i.e. not converted to the preferred units).
   *
   * @return true, if returning raw results
   */
  public boolean isRawResults() {
    return rawResults;
  }

  /**
   * Set to true to return raw results. The default is to convert to the preferred units.
   *
   * @param rawResults Set to true to return raw results
   */
  public void setRawResults(boolean rawResults) {
    this.rawResults = rawResults;
  }

  /**
   * Gets the options for reading the results. Allows specific file formats to provide options for
   * how to read the data.
   *
   * @return the options
   */
  public ResultOption[] getOptions() {
    getHeader();
    if (header == null || format == null) {
      return null;
    }
    if (format == FileFormat.TSF_BINARY) {
      final TSFPeakResultsReader reader = new TSFPeakResultsReader(filename);
      return reader.getOptions();
    }
    return null;
  }

  /**
   * Sets the options for reading the results.
   *
   * @param options the new options for reading the results
   */
  public void setOptions(ResultOption[] options) {
    this.options = options;
  }
}
