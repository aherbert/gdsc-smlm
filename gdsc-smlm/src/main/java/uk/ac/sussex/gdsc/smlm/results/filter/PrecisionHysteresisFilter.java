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

package uk.ac.sussex.gdsc.smlm.results.filter;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Filter results using a precision threshold. Any results below the lower precision limit are
 * included. Any results above the upper precision limit are excluded. Any results between the
 * limits are included only if they can be traced through time, optionally via other candidates, to
 * a valid result.
 */
public class PrecisionHysteresisFilter extends HysteresisFilter {
  /** The default range. */
  static final double[] DEFAULT_RANGE =
      {0, 0, 0, 0, PrecisionFilter.DEFAULT_RANGE, PrecisionFilter.DEFAULT_RANGE};

  @XStreamAsAttribute
  private final double strictPrecision;
  @XStreamAsAttribute
  private final double range;
  @XStreamOmitField
  private double lowerVariance;
  @XStreamOmitField
  private double upperVariance;
  @XStreamOmitField
  private Gaussian2DPeakResultCalculator calculator;

  /**
   * Instantiates a new precision hysteresis filter.
   *
   * @param searchDistance the search distance
   * @param searchDistanceMode 0 = relative to the precision of the candidates; 1 = Absolute (in nm)
   * @param timeThreshold the time threshold
   * @param timeThresholdMode 0 = frames; 1 = seconds
   * @param strictPrecision the strict precision
   * @param range the range
   */
  public PrecisionHysteresisFilter(double searchDistance, int searchDistanceMode,
      double timeThreshold, int timeThresholdMode, double strictPrecision, double range) {
    super(searchDistance, searchDistanceMode, timeThreshold, timeThresholdMode);
    this.strictPrecision = Math.max(0, strictPrecision);
    this.range = Math.max(0, range);
  }

  @Override
  protected String generateName() {
    return String.format("Precision Hysteresis %.2f +%.2f (%s)", strictPrecision, range,
        getTraceParameters());
  }

  @Override
  public void setup(MemoryPeakResults peakResults) {
    calculator = Gaussian2DPeakResultHelper.create(peakResults.getPsf(),
        peakResults.getCalibration(), Gaussian2DPeakResultHelper.LSE_PRECISION);
    lowerVariance = Filter.getDUpperSquaredLimit(strictPrecision);
    upperVariance = Filter.getDUpperSquaredLimit(strictPrecision + range);
    super.setup(peakResults);
  }

  @Override
  protected PeakStatus getStatus(PeakResult result) {
    // Use the background noise to estimate precision
    final double variance = calculator.getLseVariance(result.getParameters(), result.getNoise());
    if (variance <= lowerVariance) {
      return PeakStatus.OK;
    } else if (variance <= upperVariance) {
      return PeakStatus.CANDIDATE;
    }
    return PeakStatus.REJECT;
  }

  @Override
  public double getNumericalValue() {
    return strictPrecision;
  }

  @Override
  public String getNumericalValueName() {
    return ParameterType.PRECISION.toString() + " +" + range;
  }

  @Override
  public String getDescription() {
    return "Filter results using a precision threshold. Any results below the lower precision "
        + "limit are included. Any results above the upper precision limit are excluded. "
        + super.getDescription();
  }

  @Override
  public int getNumberOfParameters() {
    return 2 + super.getNumberOfParameters();
  }

  @Override
  protected double getParameterValueInternal(int index) {
    if (index < super.getNumberOfParameters()) {
      return super.getParameterValueInternal(index);
    }
    index -= super.getNumberOfParameters();
    return (index == 0) ? strictPrecision : range;
  }

  @Override
  public double getParameterIncrement(int index) {
    checkIndex(index);
    if (index < super.getNumberOfParameters()) {
      return super.getParameterIncrement(index);
    }
    return PrecisionFilter.DEFAULT_INCREMENT;
  }

  @Override
  public ParameterType getParameterType(int index) {
    checkIndex(index);
    if (index < super.getNumberOfParameters()) {
      return super.getParameterType(index);
    }
    index -= super.getNumberOfParameters();
    return (index == 0) ? ParameterType.PRECISION : ParameterType.PRECISION_RANGE;
  }

  @Override
  public Filter adjustParameter(int index, double delta) {
    checkIndex(index);
    // No adjustment of the mode parameters
    if (index == 1 || index == 3) {
      return this;
    }
    final double[] parameters = {searchDistance, searchDistanceMode, timeThreshold,
        timeThresholdMode, strictPrecision, range};
    if (index == 0) {
      parameters[0] = updateParameter(parameters[0], delta, getDefaultSearchRange());
    } else if (index == 2) {
      parameters[2] = updateParameter(parameters[2], delta, getDefaultTimeRange());
    } else {
      parameters[index] = updateParameter(parameters[index], delta, DEFAULT_RANGE[index]);
    }
    return create(parameters);
  }

  @Override
  public Filter create(double... parameters) {
    return new PrecisionHysteresisFilter(parameters[0], (int) parameters[1], parameters[2],
        (int) parameters[3], parameters[4], parameters[5]);
  }

  @Override
  public void weakestParameters(double[] parameters) {
    super.weakestParameters(parameters);

    // Hysteresis filters require all the potential candidates, so disable hysteresis above the
    // candidate threshold
    setMax(parameters, 4, strictPrecision + range);
    parameters[5] = 0;
  }

  @Override
  public double[] upperLimit() {
    return new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY,
        PrecisionFilter.UPPER_LIMIT, PrecisionFilter.UPPER_LIMIT};
  }

  @Override
  public double[] mutationStepRange() {
    return new double[] {getDefaultSearchRange(), getDefaultTimeRange(),
        PrecisionFilter.DEFAULT_RANGE, PrecisionFilter.DEFAULT_RANGE};
  }
}
