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

import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;

import java.util.ArrayList;

/**
 * Contain the configuration file settings for the batch fitting plugin.
 */
public class BatchSettings {
  /** The images. */
  public ArrayList<String> images = new ArrayList<>();

  /** The parameters. */
  public ArrayList<ParameterSettings> parameters = new ArrayList<>();

  /** The results directory. */
  public String resultsDirectory = null;

  /** Set to true to run the peak fit plugin. */
  public boolean runPeakFit = true;

  private Calibration calibration = null;

  /**
   * Gets the calibration.
   *
   * @return the calibration
   */
  public Calibration getCalibration() {
    return (calibration != null) ? calibration : Calibration.getDefaultInstance();
  }

  /**
   * Sets the calibration.
   *
   * @param calibration the calibration to set
   */
  public void setCalibration(Calibration calibration) {
    this.calibration = calibration;
  }
}
