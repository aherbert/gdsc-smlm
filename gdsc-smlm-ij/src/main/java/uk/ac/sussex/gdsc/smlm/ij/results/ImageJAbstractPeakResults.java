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

package uk.ac.sussex.gdsc.smlm.ij.results;

import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.results.AbstractPeakResults;
import uk.ac.sussex.gdsc.smlm.results.ThreadSafePeakResults;

/**
 * Wrap the fit results to provide convenience methods.
 */
public abstract class ImageJAbstractPeakResults extends AbstractPeakResults
    implements ThreadSafePeakResults {
  /**
   * Sets the calibration.
   *
   * <p>The calibration distance unit is set to pixels and the intensity unit to photons.
   *
   * @param nmPerPixel the nm per pixel
   * @param gain the gain
   */
  public void setCalibration(double nmPerPixel, double gain) {
    final CalibrationWriter cw = getCalibrationWriterSafe();
    cw.setNmPerPixel(nmPerPixel);
    cw.setCountPerPhoton(gain);
    cw.setDistanceUnit(DistanceUnit.PIXEL);
    cw.setIntensityUnit(IntensityUnit.PHOTON);
    setCalibration(cw.getCalibration());
  }
}
