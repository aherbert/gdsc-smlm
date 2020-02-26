/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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
import com.thoughtworks.xstream.annotations.XStreamOmitField;
import uk.ac.sussex.gdsc.smlm.data.config.ConfigurationException;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Filter results using a precision threshold. Calculates the precision using the true fitted
 * background if a bias is provided.
 */
public class PrecisionFilter2 extends DirectFilter implements IMultiFilter {
  @XStreamAsAttribute
  private final double precision;
  @XStreamOmitField
  private double variance;
  @XStreamOmitField
  private boolean useBackground;
  @XStreamOmitField
  private Gaussian2DPeakResultCalculator calculator;

  /**
   * Instantiates a new precision filter 2.
   *
   * @param precision the precision
   */
  public PrecisionFilter2(double precision) {
    this.precision = Math.max(0, precision);
  }

  @Override
  public void setup(MemoryPeakResults peakResults) {
    try {
      calculator = Gaussian2DPeakResultHelper.create(peakResults.getPsf(),
          peakResults.getCalibration(), Gaussian2DPeakResultHelper.LSE_PRECISION_X);
      useBackground = true;
    } catch (final ConfigurationException ex) {
      calculator = Gaussian2DPeakResultHelper.create(peakResults.getPsf(),
          peakResults.getCalibration(), Gaussian2DPeakResultHelper.LSE_PRECISION);
      useBackground = false;
    }
    variance = Filter.getDUpperSquaredLimit(precision);
  }

  @Override
  public boolean accept(PeakResult peak) {
    if (useBackground) {
      // Use the estimated background for the peak
      return calculator.getLsePrecision(peak.getParameters()) <= variance;
    }
    // Use the background noise to estimate precision
    return calculator.getLsePrecision(peak.getParameters(), peak.getNoise()) <= variance;
  }

  @Override
  public int getValidationFlags() {
    return FilterValidationFlag.LOCATION_VARIANCE2;
  }

  @Override
  public int validate(final PreprocessedPeakResult peak) {
    if (peak.getLocationVariance2() > variance) {
      return getValidationFlags();
    }
    return 0;
  }

  @Override
  public String getDescription() {
    return "Filter results using an upper precision threshold (uses fitted background"
        + " to set noise).";
  }

  @Override
  public int getNumberOfParameters() {
    return 1;
  }

  @Override
  protected double getParameterValueInternal(int index) {
    return precision;
  }

  @Override
  public double getParameterIncrement(int index) {
    checkIndex(index);
    return PrecisionFilter.DEFAULT_INCREMENT;
  }

  @Override
  public ParameterType getParameterType(int index) {
    checkIndex(index);
    return ParameterType.PRECISION2;
  }

  @Override
  public Filter adjustParameter(int index, double delta) {
    checkIndex(index);
    return new PrecisionFilter2(updateParameter(precision, delta, PrecisionFilter.DEFAULT_RANGE));
  }

  @Override
  public Filter create(double... parameters) {
    return new PrecisionFilter2(parameters[0]);
  }

  @Override
  public void weakestParameters(double[] parameters) {
    setMax(parameters, 0, precision);
  }

  @Override
  public int lowerBoundOrientation(int index) {
    return 1;
  }

  @Override
  public double[] upperLimit() {
    return new double[] {PrecisionFilter.UPPER_LIMIT};
  }

  @Override
  public double[] mutationStepRange() {
    return new double[] {PrecisionFilter.DEFAULT_RANGE};
  }

  @Override
  public double getSignal() {
    return 0;
  }

  @Override
  public double getSnr() {
    return 0;
  }

  @Override
  public double getMinWidth() {
    return 0;
  }

  @Override
  public double getMaxWidth() {
    return 0;
  }

  @Override
  public double getShift() {
    return 0;
  }

  @Override
  public double getEShift() {
    return 0;
  }

  @Override
  public double getPrecision() {
    return precision;
  }

  @Override
  public PrecisionType getPrecisionType() {
    return PrecisionType.ESTIMATE_USING_LOCAL_BACKGROUND;
  }

  @Override
  public double getMinZ() {
    return 0;
  }

  @Override
  public double getMaxZ() {
    return 0;
  }
}
