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

package uk.ac.sussex.gdsc.smlm.results;

import uk.ac.sussex.gdsc.core.data.NotImplementedException;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.ConfigurationException;
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.CameraType;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.FitMode;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.FluorophoreType;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.IntensityUnits;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.LocationUnits;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.ROI;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.Spot;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.Spot.Builder;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.SpotList;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.ThetaUnits;

import com.google.protobuf.InvalidProtocolBufferException;
import com.google.protobuf.util.JsonFormat;
import com.google.protobuf.util.JsonFormat.Printer;

import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.RandomAccessFile;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Saves the fit results to file using the Tagged Spot File (TSF) format.
 *
 * <p>Write out a TSF file assuming the results are in the standard GSDC SMLM format (intensity in
 * counts, angles in degrees).
 *
 * <p>To satisfy the format the calibration must be set including the amplification
 * (electrons/count) and camera bias. The bias is removed from the background. If amplification is
 * not strictly positive then the calibration gain will be written to the TSF 'electron conversion
 * factor' field.
 */
public class TsfPeakResultsWriter extends AbstractPeakResults {
  private static Logger logger = Logger.getLogger(TsfPeakResultsWriter.class.getName());

  private static final CameraType[] cameraTypeMap;
  private static final ThetaUnits[] thetaUnitsMap;
  private static final LocationUnits[] locationUnitsMap;
  private static final IntensityUnits[] intensityUnitsMap;

  static {
    // These should have 1:1 mapping. We can extends the TSF proto if necessary.
    cameraTypeMap = new CameraType[uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType
        .values().length];
    cameraTypeMap[uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType.CCD.ordinal()] =
        CameraType.CCD;
    cameraTypeMap[uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType.EMCCD.ordinal()] =
        CameraType.EMCCD;
    cameraTypeMap[uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType.SCMOS.ordinal()] =
        CameraType.SCMOS;
    thetaUnitsMap =
        new ThetaUnits[uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit.values().length];
    thetaUnitsMap[uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit.RADIAN.ordinal()] =
        ThetaUnits.RADIANS;
    thetaUnitsMap[uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit.DEGREE.ordinal()] =
        ThetaUnits.DEGREES;
    locationUnitsMap = new LocationUnits[uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit
        .values().length];
    locationUnitsMap[uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit.NM.ordinal()] =
        LocationUnits.NM;
    locationUnitsMap[uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit.UM.ordinal()] =
        LocationUnits.UM;
    locationUnitsMap[uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit.PIXEL.ordinal()] =
        LocationUnits.PIXELS;
    intensityUnitsMap =
        new IntensityUnits[uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit
            .values().length];
    intensityUnitsMap[uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit.COUNT.ordinal()] =
        IntensityUnits.COUNTS;
    intensityUnitsMap[uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit.PHOTON
        .ordinal()] = IntensityUnits.PHOTONS;
  }

  /**
   * Application ID assigned to GDSC SMLM ImageJ plugins.
   */
  public static final int APPLICATION_ID = 4;

  private OutputStream out;

  private final String filename;

  private int size;
  private AtomicInteger id;

  private int isx;
  private int isy;
  private int ia;
  private FitMode fitMode;

  private int boxSize;

  /**
   * Instantiates a new TSF peak results writer.
   *
   * @param filename the filename
   */
  public TsfPeakResultsWriter(String filename) {
    // Note: it would be good to be able to use the ability to write to any output stream. However
    // the TSF format requires a seek at the end of writing to record the offset. seek() is not
    // supported by OutputStream. It is supported by: RandomAccessFile, RandomAccessStream (for
    // input). So for now just support write to file.
    this.filename = filename;
  }

  @Override
  public void begin() {
    out = null;
    size = 0;

    // Only support Gaussian 2D data
    if (getPsf() == null || !PsfHelper.isGaussian2D(getPsf())) {
      logger.log(Level.SEVERE, "TSF format requires a Gaussian 2D PSF");
      closeOutput();
      return;
    }
    final int[] indices = PsfHelper.getGaussian2DWxWyIndices(getPsf());
    isx = indices[0];
    isy = indices[1];
    try {
      ia = PsfHelper.getGaussian2DAngleIndex(getPsf());
      fitMode = FitMode.TWOAXISANDTHETA;
    } catch (final ConfigurationException ex) {
      // This is not an angled PSF. Revert to 1/2 axis:
      fitMode = (isx == isy) ? FitMode.ONEAXIS : FitMode.TWOAXIS;
    }

    id = new AtomicInteger();
    try {
      out = new BufferedOutputStream(new FileOutputStream(filename));
    } catch (final Exception ex) {
      logger.log(Level.SEVERE, ex, () -> "Failed to write open TSF file: " + filename);
      closeOutput();
      return;
    }

    // Write the offsets used in the TSF format
    try {
      final DataOutputStream dos = new DataOutputStream(out);
      dos.writeInt(0);
      dos.writeLong(0);
    } catch (final IOException ex) {
      logger.log(Level.SEVERE, "Failed to write TSF offset fields", ex);
      closeOutput();
    }
  }

  private void closeOutput() {
    if (out == null) {
      return;
    }

    try {
      out.close();
    } catch (final Exception ex) {
      logger.log(Level.SEVERE, "Failed to close output file", ex);
    } finally {
      out = null;
    }
  }

  @Override
  public boolean isActive() {
    return out != null;
  }

  @Override
  public void add(int peak, int origX, int origY, float origValue, double error, float noise,
      float meanIntensity, float[] params, float[] paramsStdDev) {
    if (out == null) {
      return;
    }

    final Spot.Builder builder = Spot.newBuilder();
    builder.setMolecule(id.incrementAndGet());
    builder.setChannel(1);
    builder.setFluorophoreType(1);
    builder.setFrame(peak);
    builder.setXPosition(origX);
    builder.setYPosition(origY);

    setParam(params, builder);

    builder.setError(error);
    builder.setNoise(noise);
    builder.setMeanIntensity(meanIntensity);
    builder.setOriginalValue(origValue);
    if (paramsStdDev != null) {
      addNewParamStdDevs(builder, paramsStdDev);
    }

    final Spot spot = builder.build();

    writeResult(1, spot);
  }

  @Override
  public void add(PeakResult result) {
    final float[] params = result.getParameters();

    final Spot.Builder builder = Spot.newBuilder();
    builder.setMolecule(id.incrementAndGet());
    builder.setChannel(1);
    builder.setFluorophoreType(1);
    builder.setFrame(result.getFrame());
    builder.setXPosition(result.getOrigX());
    builder.setYPosition(result.getOrigY());

    setParam(params, builder);

    if (result.hasPrecision()) {
      // Use the actual precision
      final float precision = (float) result.getPrecision();
      builder.setXPrecision(precision);
      builder.setYPrecision(precision);
    }

    if (result.hasId()) {
      builder.setCluster(result.getId());
    }

    builder.setError(result.getError());
    builder.setNoise(result.getNoise());
    builder.setMeanIntensity(result.getMeanIntensity());
    if (result.hasEndFrame()) {
      builder.setEndFrame(result.getEndFrame());
    }
    builder.setOriginalValue(result.getOrigValue());
    if (result.hasParameterDeviations()) {
      addNewParamStdDevs(builder, result.getParameterDeviations());
    }

    final Spot spot = builder.build();

    writeResult(1, spot);
  }

  /**
   * Sets the width. Convert the X/Y widths used in GDSC SMLM to the single width and shape
   * parameters used in TSF.
   *
   * @param params the params
   * @param builder the builder
   */
  private void setParam(float[] params, Spot.Builder builder) {
    builder.setBackground(params[PeakResult.BACKGROUND]);
    builder.setIntensity(params[PeakResult.INTENSITY]);
    builder.setX(params[PeakResult.X]);
    builder.setY(params[PeakResult.Y]);
    builder.setZ(params[PeakResult.Z]);

    switch (fitMode) {
      case ONEAXIS:
        builder.setWidth((float) (Gaussian2DFunction.SD_TO_FWHM_FACTOR * params[isx]));
        break;

      case TWOAXIS:
        builder.setWidth((float) (Gaussian2DFunction.SD_TO_FWHM_FACTOR
            * Gaussian2DPeakResultHelper.getStandardDeviation(params[isx], params[isy])));
        builder.setA(params[isx] / params[isy]);
        break;

      default:
        // Fit mode is validated as a 2D gaussian in the begin() method
        // TWOAXISANDTHETA
        builder.setWidth((float) (Gaussian2DFunction.SD_TO_FWHM_FACTOR
            * Gaussian2DPeakResultHelper.getStandardDeviation(params[isx], params[isy])));
        builder.setA(params[isx] / params[isy]);
        builder.setTheta(params[ia]);
        break;
    }
  }

  /**
   * Adds the params std dev assuming a new builder.
   *
   * @param builder the builder
   * @param paramStdDev the params std dev
   */
  private static void addNewParamStdDevs(Builder builder, float[] paramStdDev) {
    // Note: paramsStdDev for X/Y could be set into the X/Y Precision field.
    for (int i = 0; i < paramStdDev.length; i++) {
      builder.addParamStdDevs(paramStdDev[i]);
    }
  }

  @Override
  public void addAll(PeakResult[] results) {
    if (out == null) {
      return;
    }

    final Spot[] spots = new Spot[20];
    int count = 0;
    final Spot.Builder builder = Spot.newBuilder();
    for (final PeakResult result : results) {
      final float[] params = result.getParameters();

      builder.setMolecule(id.incrementAndGet());
      builder.setChannel(1);
      builder.setFluorophoreType(1);
      builder.setFrame(result.getFrame());
      builder.setXPosition(result.getOrigX());
      builder.setYPosition(result.getOrigY());

      setParam(params, builder);

      if (result.hasPrecision()) {
        // Use the actual precision
        final float precision = (float) result.getPrecision();
        builder.setXPrecision(precision);
        builder.setYPrecision(precision);
      }

      if (result.hasId()) {
        builder.setCluster(result.getId());
      } else {
        builder.clearCluster();
      }

      builder.setError(result.getError());
      builder.setNoise(result.getNoise());
      builder.setMeanIntensity(result.getMeanIntensity());
      if (result.hasEndFrame()) {
        builder.setEndFrame(result.getEndFrame());
      } else {
        builder.clearEndFrame();
      }
      builder.setOriginalValue(result.getOrigValue());
      addParamStdDevs(builder, result.getParameterDeviations());

      spots[count++] = builder.build();

      // Flush the output to allow for very large input lists
      if (count >= spots.length) {
        writeResult(count, spots);
        if (!isActive()) {
          return;
        }
        count = 0;
      }
    }
    writeResult(count, spots);
  }

  /**
   * Adds the params std dev assuming an existing builder (allowing re-use of the space).
   *
   * @param builder the builder
   * @param paramStdDev the params std dev
   */
  private static void addParamStdDevs(Builder builder, float[] paramStdDev) {
    // Note: paramsStdDev for X/Y could be set into the X/Y Precision field.
    if (paramStdDev == null) {
      if (builder.getParamStdDevsCount() != 0) {
        builder.clearParamStdDevs();
      }
      return;
    }

    // Reuse the space
    if (builder.getParamStdDevsCount() == paramStdDev.length) {
      for (int i = 0; i < paramStdDev.length; i++) {
        builder.setParamStdDevs(i, paramStdDev[i]);
      }
    } else {
      builder.clearParamStdDevs();
      addNewParamStdDevs(builder, paramStdDev);
    }
  }

  private synchronized void writeResult(int count, Spot... spots) {
    // In case another thread caused the output to close
    if (out == null) {
      return;
    }
    size += count;
    try {
      for (int i = 0; i < count; i++) {
        spots[i].writeDelimitedTo(out);
      }
    } catch (final IOException ex) {
      logger.log(Level.SEVERE, "Failed to write Spot message", ex);
      closeOutput();
    }
  }

  @Override
  public int size() {
    return size;
  }

  @Override
  public void end() {
    // Close the buffered output
    closeOutput();

    // Write the spot list and the offset to the SpotList message into the offset position.
    // Re-open the file for random access.
    try (RandomAccessFile f = new RandomAccessFile(new File(filename), "rw")) {
      // The offset is the amount to skip forward after reading the int
      // magic number (4 bytes) and long offset (8 bytes).
      final long offset = f.length() - 12;

      // Write the spotlist message at the end
      f.seek(f.length());
      final ByteArrayOutputStream baos = new ByteArrayOutputStream();
      createSpotList().writeDelimitedTo(baos);
      f.write(baos.toByteArray());

      // Write the offset
      f.seek(4);
      f.writeLong(offset);
    } catch (final Exception ex) {
      logger.log(Level.SEVERE, "Failed to record offset for SpotList message", ex);
    }
  }

  private SpotList createSpotList() {
    final SpotList.Builder builder = SpotList.newBuilder();

    builder.setApplicationId(APPLICATION_ID);

    builder.setNrSpots(size);

    // Add the standard details the TSF supports. We use extensions to add GDSC SMLM data.
    if (!TextUtils.isNullOrEmpty(getName())) {
      builder.setName(getName());
    }
    if (getSource() != null) {
      builder.setNrPixelsX(getSource().width);
      builder.setNrPixelsY(getSource().height);
      builder.setNrFrames(getSource().frames);

      builder.setSource(singleLine(getSource().toXml()));
    }
    if (getBounds() != null) {
      final ROI.Builder roiBuilder = builder.getRoiBuilder();
      roiBuilder.setX(getBounds().x);
      roiBuilder.setY(getBounds().y);
      roiBuilder.setXWidth(getBounds().width);
      roiBuilder.setYWidth(getBounds().height);
      builder.setRoi(roiBuilder.build());
    }
    if (hasCalibration()) {
      final CalibrationReader cr = getCalibrationReader();
      if (cr.hasNmPerPixel()) {
        builder.setPixelSize((float) cr.getNmPerPixel());
      }

      if (cr.hasExposureTime()) {
        builder.setExposureTime(cr.getExposureTime());
      }
      if (cr.hasReadNoise()) {
        builder.setReadNoise(cr.getReadNoise());
      }
      if (cr.hasBias()) {
        builder.setBias(cr.getBias());
      }
      if (cr.hasCameraType()) {
        builder.setCameraType(cameraTypeMap[cr.getCameraType().ordinal()]);
      }

      if (cr.hasDistanceUnit()) {
        builder.setLocationUnits(locationUnitsMap[cr.getDistanceUnit().ordinal()]);
      }
      if (cr.hasIntensityUnit()) {
        builder.setIntensityUnits(intensityUnitsMap[cr.getIntensityUnit().ordinal()]);
      }
      if (cr.hasAngleUnit()) {
        builder.setThetaUnits(thetaUnitsMap[cr.getAngleUnit().ordinal()]);
      }

      // We can use some logic here to get the QE
      if (cr.hasCountPerPhoton()) {
        builder.setGain(cr.getCountPerPhoton());

        final double qe = (cr.hasQuantumEfficiency()) ? cr.getQuantumEfficiency() : 1;
        // e-/photon / count/photon => e-/count
        final double ecf = qe / cr.getCountPerPhoton();

        builder.addEcf(ecf);
        builder.addQe(qe);
      }
    }
    if (!TextUtils.isNullOrEmpty(getConfiguration())) {
      builder.setConfiguration(singleLine(getConfiguration()));
    }
    if (getPsf() != null) {
      try {
        final Printer printer = JsonFormat.printer().omittingInsignificantWhitespace();
        builder.setPSF(printer.print(getPsf()));
      } catch (final InvalidProtocolBufferException ex) {
        // This shouldn't happen so throw it
        throw new NotImplementedException("Unable to serialise the PSF settings", ex);
      }
    }

    // Have a property so the boxSize can be set
    if (boxSize > 0) {
      builder.setBoxSize(boxSize);
    }

    builder.setFitMode(fitMode);

    final FluorophoreType.Builder typeBuilder = FluorophoreType.newBuilder();
    typeBuilder.setId(1);
    typeBuilder.setDescription("Default fluorophore");
    typeBuilder.setIsFiducial(false);
    builder.addFluorophoreTypes(typeBuilder.build());

    return builder.build();
  }

  private static String singleLine(String text) {
    return text.replaceAll("\n *", "");
  }

  /**
   * Gets the box size.
   *
   * @return the box size (in pixels) of rectangular box used in Gaussian fitting
   */
  public int getBoxSize() {
    return boxSize;
  }

  /**
   * Sets the box size.
   *
   * @param boxSize the box size (in pixels) of rectangular box used in Gaussian fitting
   */
  public void setBoxSize(int boxSize) {
    this.boxSize = boxSize;
  }
}
