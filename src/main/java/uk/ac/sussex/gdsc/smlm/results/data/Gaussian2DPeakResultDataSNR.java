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

package uk.ac.sussex.gdsc.smlm.results.data;

import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Gets the SNR from the PeakResult assuming a Gaussian 2D PSF. This is ratio of the average signal
 * value to the standard deviation of the background. The result must have the standard deviation
 * for each dimension in the first two additional parameters of the PeakResult parameter array.
 *
 * <p>Assumes that the mean signal is half the total signal divided by the ellipsoid area of the
 * Gaussian that contains half the signal.
 */
public class Gaussian2DPeakResultDataSNR extends PeakResultDataFloat {
  /** The index of the x width. */
  static final int i = PeakResult.STANDARD_PARAMETERS;
  /** The index of the y width. */
  static final int j = i + 1;

  /** {@inheritDoc} */
  @Override
  public Float getValue(PeakResult result) {
    return new Float(Gaussian2DPeakResultHelper.getMeanSignalUsingP05(result.getIntensity(),
        result.getParameter(i), result.getParameter(j)) / result.getNoise());
  }

  /** {@inheritDoc} */
  @Override
  public String getValueName() {
    return "Gaussian2D SNR";
  }
}
