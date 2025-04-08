/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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
 * Gets the mean signal from the PeakResult assuming a Gaussian 2D PSF. The result must have the
 * standard deviation for each dimension in the first two additional parameters of the PeakResult
 * parameter array.
 *
 * <p>Assumes that the mean signal is the total signal within 1 standard deviation of the centre
 * divided by the elliptical area of the Gaussian.
 */
public class Gaussian2DPeakResultDataMeanSignal implements PeakResultDataFloat {
  /** The index of the x width. */
  private static final int I = PeakResult.STANDARD_PARAMETERS;
  /** The index of the y width. */
  private static final int J = I + 1;

  @Override
  public Float getValue(PeakResult result) {
    return Float.valueOf(
        (float) Gaussian2DPeakResultHelper.getMeanSignalUsingR1(result.getIntensity(),
            result.getParameter(I), result.getParameter(J)));
  }

  @Override
  public String getValueName() {
    return "Gaussian2D mean signal";
  }
}
