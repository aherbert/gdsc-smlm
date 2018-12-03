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
package uk.ac.sussex.gdsc.smlm.results.filter;

import uk.ac.sussex.gdsc.smlm.data.config.PSFHelper;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Filter results using an upper width factor. Assumes width is different on the X and Y axis and
 * they are combined using s = sqrt(s0*s1).
 */
public class XYWidthFilter extends WidthFilter implements IMultiFilter {
  /**
   * Instantiates a new XY width filter.
   *
   * @param width the width
   */
  public XYWidthFilter(double width) {
    super(width);
  }

  @Override
  public void setup(MemoryPeakResults peakResults) {
    calculator =
        Gaussian2DPeakResultHelper.create(peakResults.getPSF(), peakResults.getCalibration(), 0);

    // Set the width limit
    final double[] s = PSFHelper.getGaussian2DWxWy(peakResults.getPSF());
    upperSigmaThreshold = Filter.getUpperLimit(s[0] * s[1] * width * width);
  }

  @Override
  protected void setup(final double width) {
    upperSigmaThreshold = Filter.getUpperLimit(width * width);
    widthEnabled = (width != Float.POSITIVE_INFINITY);
  }

  @Override
  public boolean accept(PeakResult peak) {
    return calculator.getStandardDeviation2(peak.getParameters()) <= upperSigmaThreshold;
  }

  @Override
  public int getValidationFlags() {
    return V_X_SD_FACTOR | V_Y_SD_FACTOR;
  }

  @Override
  public int validate(final PreprocessedPeakResult peak) {
    if (widthEnabled) {
      if (peak.getXSDFactor() * peak.getYSDFactor() > upperSigmaThreshold) {
        return V_X_SD_FACTOR | V_Y_SD_FACTOR;
      }
    }
    return 0;
  }

  /** {@inheritDoc} */
  @Override
  public String getDescription() {
    return "Filter results using an upper XY width factor. (Width is relative to initial peak width.)";
  }

  /** {@inheritDoc} */
  @Override
  public Filter adjustParameter(int index, double delta) {
    checkIndex(index);
    return new XYWidthFilter(updateParameter(width, delta, DEFAULT_RANGE));
  }

  /** {@inheritDoc} */
  @Override
  public Filter create(double... parameters) {
    return new XYWidthFilter(parameters[0]);
  }
}
