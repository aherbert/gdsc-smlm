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

package uk.ac.sussex.gdsc.smlm.data.config;

import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType;

/**
 * Contains helper functions for the CalibrationProtos class.
 */
public final class CalibrationProtosHelper {
  /** The default Calibration. */
  public static final Calibration defaultCalibration;

  static {
    final Calibration.Builder builder = Calibration.newBuilder();
    // Note: Ideally we would set QE to be 1 but this will involve creating a
    // camera calibration and it is more useful to have the default as null.
    defaultCalibration = builder.build();
  }

  /**
   * Gets the name.
   *
   * @param value the value
   * @return the name
   */
  public static String getName(CameraType value) {
    switch (value) {
      case CAMERA_TYPE_NA:
        return "NA";
      case CCD:
        return "CCD";
      case EMCCD:
        return "EMCCD";
      case SCMOS:
        return "sCMOS";
      case UNRECOGNIZED:
        return ProtosHelperUtils.UNKNOWN;
      default:
        throw new IllegalArgumentException("Unknown name: " + value);
    }
  }

  /**
   * Checks if is CCD camera type.
   *
   * @param cameraType the camera type
   * @return true, if is CCD camera type
   */
  public static boolean isCcdCameraType(CameraType cameraType) {
    return cameraType == CameraType.EMCCD || cameraType == CameraType.CCD;
  }

  /** No public constructor. */
  private CalibrationProtosHelper() {}
}
