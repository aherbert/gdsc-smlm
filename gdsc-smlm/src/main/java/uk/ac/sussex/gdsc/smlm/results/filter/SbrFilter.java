/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Filter results using a signal-to-background ratio (SBR) threshold.
 *
 * <p>Requires the bias to be configured at or above zero. If the background is below the configured
 * bias or there is no bias then the filter resorts to a signal-to-noise filter. If there is a
 * background level above the bias then this is assumed to be the variance of the photon shot noise
 * and the noise is taken at the square root of the background level.
 */
public class SbrFilter extends DirectFilter {
  /**
   * The signal-to-background ratio (SBR).
   */
  @XStreamAsAttribute
  final float sbr;

  /**
   * Instantiates a new signal-to-background ratio (SBR) filter.
   *
   * @param sbr the signal-to-background ratio (SBR)
   */
  public SbrFilter(float sbr) {
    this.sbr = Math.max(0, sbr);
  }

  @Override
  public void setup(MemoryPeakResults peakResults) {
    // Do nothing
  }

  @Override
  public boolean accept(PeakResult peak) {
    final double background = peak.getBackground();
    if (background > 0) {
      return peak.getMeanIntensity() / Math.sqrt(background) >= this.sbr;
    }
    return peak.getSnr() >= this.sbr;
  }

  @Override
  public int getValidationFlags() {
    return FilterValidationFlag.PHOTONS | FilterValidationFlag.BACKGROUND
        | FilterValidationFlag.SNR;
  }

  @Override
  public int validate(final PreprocessedPeakResult peak) {
    final double background = peak.getBackground();
    if (background > 0) {
      // Get the mean signal assuming the integral / area of 1 SD of the Gaussian
      if (Gaussian2DPeakResultHelper.getMeanSignalUsingR1(peak.getSignal(), peak.getXSd(),
          peak.getYSd()) / Math.sqrt(background) < this.sbr) {
        return FilterValidationFlag.PHOTONS | FilterValidationFlag.BACKGROUND;
      }
      return 0;
    }
    if (peak.getSnr() < this.sbr) {
      return FilterValidationFlag.SNR;
    }
    return 0;
  }

  @Override
  public String getDescription() {
    return "Filter results using a lower SBR threshold.";
  }

  @Override
  public int getNumberOfParameters() {
    return 1;
  }

  @Override
  protected double getParameterValueInternal(int index) {
    return sbr;
  }

  @Override
  public double getParameterIncrement(int index) {
    checkIndex(index);
    return SnrFilter.DEFAULT_INCREMENT;
  }

  @Override
  public ParameterType getParameterType(int index) {
    checkIndex(index);
    return ParameterType.SBR;
  }

  @Override
  public Filter adjustParameter(int index, double delta) {
    checkIndex(index);
    return new SbrFilter(updateParameter(sbr, delta, SnrFilter.DEFAULT_RANGE));
  }

  @Override
  public Filter create(double... parameters) {
    return new SbrFilter((float) parameters[0]);
  }

  @Override
  public void weakestParameters(double[] parameters) {
    setMin(parameters, 0, sbr);
  }

  @Override
  public double[] mutationStepRange() {
    return new double[] {SnrFilter.DEFAULT_RANGE};
  }
}
