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

package uk.ac.sussex.gdsc.smlm.results.filter;

import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Filter results using a width range. Assumes width is different on the X and Y axis and they are
 * combined using s = sqrt(s0*s1).
 */
public class XyWidthFilter2 extends WidthFilter2 implements IMultiFilter {

  /**
   * Instantiates a new XY width filter 2.
   *
   * @param minWidth the min width
   * @param maxWidth the max width
   */
  public XyWidthFilter2(double minWidth, double maxWidth) {
    super(minWidth, maxWidth);
  }

  @Override
  protected String generateName() {
    return "XYWidth " + minWidth + "-" + maxWidth;
  }

  @Override
  public void setup(MemoryPeakResults peakResults) {
    calculator =
        Gaussian2DPeakResultHelper.create(peakResults.getPsf(), peakResults.getCalibration(), 0);

    // Set the width limit
    final double[] s = PsfHelper.getGaussian2DWxWy(peakResults.getPsf());
    final double s2 = s[0] * s[1];
    lowerSigmaThreshold = (float) (s2 * minWidth * minWidth);
    upperSigmaThreshold = Filter.getUpperLimit(s2 * maxWidth * maxWidth);
  }

  @Override
  protected void setup(final double minWidth, double maxWidth) {
    widthEnabled = false;
    if (maxWidth > 1 && maxWidth != Double.POSITIVE_INFINITY) {
      upperSigmaThreshold = Filter.getUpperLimit(maxWidth * maxWidth);
      widthEnabled = upperSigmaThreshold != Float.POSITIVE_INFINITY;
    } else {
      upperSigmaThreshold = Float.POSITIVE_INFINITY;
    }
    if (minWidth < 1) {
      widthEnabled = true;
      lowerSigmaThreshold = (float) (minWidth * minWidth);
    } else {
      lowerSigmaThreshold = 0f;
    }
  }

  @Override
  public boolean accept(PeakResult peak) {
    final float sd2 = calculator.getStandardDeviation2(peak.getParameters());
    return sd2 <= upperSigmaThreshold && sd2 >= lowerSigmaThreshold;
  }

  @Override
  public int getValidationFlags() {
    return V_X_SD_FACTOR | V_Y_SD_FACTOR;
  }

  @Override
  public int validate(final PreprocessedPeakResult peak) {
    if (widthEnabled) {
      final float s2 = peak.getXSdFactor() * peak.getYSdFactor();
      if (s2 > upperSigmaThreshold || s2 < lowerSigmaThreshold) {
        return V_X_SD_FACTOR | V_Y_SD_FACTOR;
      }
    }
    return 0;
  }

  @Override
  public String getDescription() {
    return "Filter results using an XY width range. (Width is relative to initial peak width.)";
  }

  @Override
  public Filter adjustParameter(int index, double delta) {
    checkIndex(index);
    if (index == 0) {
      return new XyWidthFilter2(updateParameter(minWidth, delta, DEFAULT_MIN_RANGE), maxWidth);
    }
    return new XyWidthFilter2(minWidth,
        updateParameter(maxWidth, delta, WidthFilter.DEFAULT_RANGE));
  }

  @Override
  public Filter create(double... parameters) {
    return new XyWidthFilter2(parameters[0], parameters[1]);
  }
}
