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

import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;

import java.awt.Rectangle;
import java.util.Collection;

/**
 * Abstract base class for peak results.
 */
public abstract class AbstractPeakResults implements PeakResults {
  /** The default for nm./pixel */
  public static final double DEFAULT_NM_PER_PIXEL = 0;
  /** The default for gain. */
  public static final double DEFAULT_GAIN = 0;
  /** The default for emCCD. */
  public static final boolean DEFAULT_EMCCD = true;

  private ImageSource source;
  private Rectangle bounds;
  private Calibration calibration;
  private PSF psf;
  private String configuration = "";
  private String name = "";

  /** The calibration reader. This is encapsulated */
  private CalibrationReader calibrationReader;

  /** {@inheritDoc} */
  @Override
  public void addAll(Collection<PeakResult> results) {
    // Utility function
    addAll(results.toArray(new PeakResult[results.size()]));
  }

  /** {@inheritDoc} */
  @Override
  public void addAll(PeakResultStore results) {
    // Utility function
    addAll(results.toArray());
  }

  /** {@inheritDoc} */
  @Override
  public void setSource(ImageSource source) {
    this.source = source;
  }

  /** {@inheritDoc} */
  @Override
  public ImageSource getSource() {
    return source;
  }

  /** {@inheritDoc} */
  @Override
  public void setBounds(Rectangle bounds) {
    this.bounds = bounds;
  }

  /** {@inheritDoc} */
  @Override
  public Rectangle getBounds() {
    return bounds;
  }

  /**
   * Gets the bounds as a string.
   *
   * @return the bounds string
   */
  public String getBoundsString() {
    if (bounds != null) {
      return String.format("x%d y%d w%d h%d", bounds.x, bounds.y, bounds.width, bounds.height);
    }
    return "";
  }

  /** {@inheritDoc} */
  @Override
  public void setCalibration(Calibration calibration) {
    this.calibration = calibration;
    calibrationReader = (calibration != null) ? new CalibrationReader(calibration) : null;
  }

  /** {@inheritDoc} */
  @Override
  public Calibration getCalibration() {
    return calibration;
  }

  /**
   * Checks for calibration.
   *
   * @return true, if successful
   */
  public boolean hasCalibration() {
    return (calibration != null);
  }

  /**
   * Gets the calibration reader with the current calibration.
   *
   * @return the calibration reader (or null)
   */
  public CalibrationReader getCalibrationReader() {
    return calibrationReader;
  }

  /**
   * Gets the calibration writer with the current calibration (must not be null). The writer can be
   * used to update the calibration but changes are not saved until
   * {@link uk.ac.sussex.gdsc.smlm.results.PeakResults#setCalibration(uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration)}
   * is called with the new calibration.
   *
   * @return the calibration writer
   * @throws IllegalArgumentException if the calibration is null
   */
  public CalibrationWriter getCalibrationWriter() {
    return new CalibrationWriter(calibration);
  }

  /**
   * Gets the calibration writer with the current calibration, or a default calibration. The writer
   * can be used to update the calibration but changes are not saved until
   * {@link uk.ac.sussex.gdsc.smlm.results.PeakResults#setCalibration(uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration)}
   * is called with the new calibration.
   *
   * @return the calibration writer
   */
  public CalibrationWriter getCalibrationWriterSafe() {
    return (calibration != null) ? new CalibrationWriter(calibration) : new CalibrationWriter();
  }

  /** {@inheritDoc} */
  @Override
  public PSF getPSF() {
    return psf;
  }

  /** {@inheritDoc} */
  @Override
  public void setPSF(PSF psf) {
    this.psf = psf;
  }

  /** {@inheritDoc} */
  @Override
  public void setConfiguration(String configuration) {
    this.configuration = configuration;
  }

  /** {@inheritDoc} */
  @Override
  public String getConfiguration() {
    return configuration;
  }

  /**
   * @return The name of the results set (or the source if empty).
   */
  @Override
  public String getName() {
    if (name.length() > 0) {
      return name;
    }
    return (getSource() != null) ? getSource().getName() : "";
  }

  /**
   * @param name The name of the results set
   */
  @Override
  public void setName(String name) {
    if (name == null) {
      this.name = "";
    } else {
      this.name = name;
    }
  }

  /**
   * Get the nm-per-pixel from the calibration, or if not available, return the
   * {@link #DEFAULT_NM_PER_PIXEL}
   *
   * @return the nmPerPixel
   */
  public double getNmPerPixel() {
    return (calibration != null) ? calibrationReader.getNmPerPixel() : DEFAULT_NM_PER_PIXEL;
  }

  /**
   * Get the gain from the calibration, or if not available, return the {@link #DEFAULT_GAIN}
   *
   * @return the gain
   */
  public double getGain() {
    return (calibration != null) ? calibrationReader.getCountPerPhoton() : DEFAULT_GAIN;
  }

  /**
   * Checks for a CCD camera.
   *
   * @return true, if successful
   */
  public boolean isCCDCamera() {
    return (calibration != null) ? calibrationReader.isCCDCamera() : false;
  }

  /**
   * Get the EMCCD flag from the calibration, or if not available, return the {@link #DEFAULT_EMCCD}
   *
   * @return the EMCCD flag
   * @deprecated Replaced by the camera type
   */
  @Deprecated
  public boolean isEMCCD() {
    return (calibration != null && calibrationReader.isCCDCamera()) ? calibrationReader.isEMCCD()
        : DEFAULT_EMCCD;
  }

  /**
   * Checks if the results have a valid calibration to compute the localisation precision. This
   * requires the pixel size and camera gain, or alternatively the units to be in nm and photons,
   * and camera CCD type.
   *
   * @return true, if is calibrated for precision
   */
  public boolean isCalibratedForPrecision() {
    if (calibration != null) {
      if (!calibrationReader.isCCDCamera()) {
        return false;
      }
      final DistanceUnit du = calibrationReader.getDistanceUnit();
      final IntensityUnit iu = calibrationReader.getIntensityUnit();
      if (du == DistanceUnit.NM && iu == IntensityUnit.PHOTON) {
        return true;
      }
      return isCalibrated();
    }
    return false;
  }

  /**
   * Checks if the results have a basic calibration. This requires the pixel size and camera gain
   * with the distance and intensity units.
   *
   * @return true, if is calibrated
   */
  public boolean isCalibrated() {
    if (calibration != null) {
      final DistanceUnit du = calibrationReader.getDistanceUnit();
      final IntensityUnit iu = calibrationReader.getIntensityUnit();
      //@formatter:off
      return (du != null && calibrationReader.getNmPerPixel() > 0) &&
           (iu != null && calibrationReader.getCountPerPhoton() > 0);
      //@formatter:on
    }
    return false;
  }

  /** {@inheritDoc} */
  @Override
  public void copySettings(PeakResults peakResults) {
    this.setSource(peakResults.getSource());
    this.setBounds(peakResults.getBounds());
    this.setCalibration(peakResults.getCalibration());
    this.setPSF(peakResults.getPSF());
    this.setConfiguration(peakResults.getConfiguration());
    this.setName(peakResults.getName());
  }
}
