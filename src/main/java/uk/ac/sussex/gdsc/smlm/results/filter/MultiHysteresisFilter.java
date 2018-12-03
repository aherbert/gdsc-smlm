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
import uk.ac.sussex.gdsc.smlm.ga.Chromosome;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

/**
 * Filter results using multiple thresholds: Signal, SNR, width, coordinate shift and precision. <p>
 * Any results with the strict limits are included. Any results outside the weak limits are
 * excluded. Any results between the strict and weak limits are included only if they can be traced
 * through time, optionally via other candidates, to a valid result.
 */
public class MultiHysteresisFilter extends HysteresisFilter {
  /** The strict signal. */
  @XStreamAsAttribute
  final double strictSignal;

  /** The strict snr. */
  @XStreamAsAttribute
  final float strictSnr;

  /** The strict min width. */
  @XStreamAsAttribute
  final double strictMinWidth;

  /** The strict max width. */
  @XStreamAsAttribute
  final double strictMaxWidth;

  /** The strict shift. */
  @XStreamAsAttribute
  final double strictShift;

  /** The strict precision. */
  @XStreamAsAttribute
  final double strictPrecision;

  /** The range signal. */
  @XStreamAsAttribute
  final double rangeSignal;

  /** The range snr. */
  @XStreamAsAttribute
  final float rangeSnr;

  /** The range min width. */
  @XStreamAsAttribute
  final double rangeMinWidth;

  /** The range max width. */
  @XStreamAsAttribute
  final double rangeMaxWidth;

  /** The range shift. */
  @XStreamAsAttribute
  final double rangeShift;

  /** The range precision. */
  @XStreamAsAttribute
  final double rangePrecision;

  /** The strict signal threshold. */
  @XStreamOmitField
  float strictSignalThreshold;

  /** The weak signal threshold. */
  @XStreamOmitField
  float weakSignalThreshold;

  /** The weak snr. */
  @XStreamOmitField
  float weakSnr;

  /** The strict min sigma threshold. */
  @XStreamOmitField
  float strictMinSigmaThreshold;

  /** The weak min sigma threshold. */
  @XStreamOmitField
  float weakMinSigmaThreshold;

  /** The strict max sigma threshold. */
  @XStreamOmitField
  float strictMaxSigmaThreshold;

  /** The weak max sigma threshold. */
  @XStreamOmitField
  float weakMaxSigmaThreshold;

  /** The strict offset. */
  @XStreamOmitField
  float strictOffset;

  /** The weak offset. */
  @XStreamOmitField
  float weakOffset;

  /** The strict variance. */
  @XStreamOmitField
  double strictVariance;

  /** The weak variance. */
  @XStreamOmitField
  double weakVariance;

  /** The calculator. */
  @XStreamOmitField
  Gaussian2DPeakResultCalculator calculator;

  /**
   * Instantiates a new multi hysteresis filter.
   *
   * @param searchDistance the search distance
   * @param searchDistanceMode 0 = relative to the precision of the candidates; 1 = Absolute (in nm)
   * @param timeThreshold the time threshold
   * @param timeThresholdMode 0 = frames; 1 = seconds
   * @param strictSignal the strict signal
   * @param rangeSignal the range signal
   * @param strictSnr the strict snr
   * @param rangeSnr the range snr
   * @param strictMinWidth the strict min width
   * @param rangeMinWidth the range min width
   * @param strictMaxWidth the strict max width
   * @param rangeMaxWidth the range max width
   * @param strictShift the strict shift
   * @param rangeShift the range shift
   * @param strictPrecision the strict precision
   * @param rangePrecision the range precision
   */
  public MultiHysteresisFilter(double searchDistance, int searchDistanceMode, double timeThreshold,
      int timeThresholdMode, double strictSignal, double rangeSignal, float strictSnr,
      float rangeSnr, double strictMinWidth, double rangeMinWidth, double strictMaxWidth,
      double rangeMaxWidth, double strictShift, double rangeShift, double strictPrecision,
      double rangePrecision) {
    super(searchDistance, searchDistanceMode, timeThreshold, timeThresholdMode);
    this.strictSignal = Math.max(0, strictSignal);
    this.rangeSignal = Math.max(0, rangeSignal);
    this.strictSnr = Math.max(0, strictSnr);
    this.rangeSnr = Math.max(0, rangeSnr);
    this.strictMinWidth = Math.max(0, strictMinWidth);
    this.rangeMinWidth = Math.max(0, rangeMinWidth);
    this.strictMaxWidth = Math.max(0, strictMaxWidth);
    this.rangeMaxWidth = Math.max(0, rangeMaxWidth);
    this.strictShift = Math.max(0, strictShift);
    this.rangeShift = Math.max(0, rangeShift);
    this.strictPrecision = Math.max(0, strictPrecision);
    this.rangePrecision = Math.max(0, rangePrecision);
  }

  @Override
  protected String generateName() {
    return String.format(
        "Multi Hysteresis: Signal=%.1f-%.1f, SNR=%.1f-%.1f, MinWidth=%.2f-%.2f, MaxWidth=%.2f+%.2f, Shift=%.2f+%.2f, Precision=%.1f+%.1f (%s)",
        strictSignal, rangeSignal, strictSnr, rangeSnr, strictMinWidth, rangeMinWidth,
        strictMaxWidth, rangeMaxWidth, strictShift, rangeShift, strictPrecision, rangePrecision,
        getTraceParameters());
  }

  @Override
  public void setup(MemoryPeakResults peakResults) {
    setupCalculator(peakResults);

    // Set the signal limit using the gain
    strictSignalThreshold = (float) (strictSignal);
    weakSignalThreshold = (float) (strictSignal - rangeSignal);

    weakSnr = strictSnr - rangeSnr;

    // Set the width limit
    strictMinSigmaThreshold = weakMinSigmaThreshold = 0;
    strictMaxSigmaThreshold = weakMaxSigmaThreshold = Float.POSITIVE_INFINITY;
    // Set the shift limit
    strictOffset = weakOffset = Float.POSITIVE_INFINITY;

    final double s = PSFHelper.getGaussian2DWx(peakResults.getPSF());
    strictMinSigmaThreshold = (float) (s * strictMinWidth);
    strictMaxSigmaThreshold = Filter.getUpperLimit(s * strictMaxWidth);
    weakMinSigmaThreshold = (float) (s * (strictMinWidth - rangeMinWidth));
    weakMaxSigmaThreshold = Filter.getUpperLimit(s * (strictMaxWidth + rangeMaxWidth));

    strictOffset = Filter.getUpperLimit(s * strictShift);
    weakOffset = Filter.getUpperLimit(s * (strictShift + rangeShift));

    // Configure the precision limit
    strictVariance = Filter.getDUpperSquaredLimit(strictPrecision);
    weakVariance = Filter.getDUpperSquaredLimit(strictPrecision + rangePrecision);

    super.setup(peakResults);
  }

  /**
   * Set up the calculator.
   *
   * @param peakResults the results
   */
  protected void setupCalculator(MemoryPeakResults peakResults) {
    calculator = Gaussian2DPeakResultHelper.create(peakResults.getPSF(),
        peakResults.getCalibration(), Gaussian2DPeakResultHelper.LSE_PRECISION);
  }

  @Override
  protected PeakStatus getStatus(PeakResult result) {
    // Check weak thresholds
    if (result.getIntensity() < weakSignalThreshold) {
      return PeakStatus.REJECT;
    }
    final float snr = result.getSNR();
    if (snr < weakSnr) {
      return PeakStatus.REJECT;
    }
    final float sd = calculator.getStandardDeviation(result.getParameters());
    if (sd < weakMinSigmaThreshold || sd > weakMaxSigmaThreshold) {
      return PeakStatus.REJECT;
    }
    if (Math.abs(result.getXPosition()) > weakOffset
        || Math.abs(result.getYPosition()) > weakOffset) {
      return PeakStatus.REJECT;
    }
    final double variance = getVariance(result);
    if (variance > weakVariance) {
      return PeakStatus.REJECT;
    }

    // Check the strict thresholds
    if (result.getIntensity() < strictSignalThreshold) {
      return PeakStatus.CANDIDATE;
    }
    if (snr < strictSnr) {
      return PeakStatus.CANDIDATE;
    }
    if (sd < strictMinSigmaThreshold || sd > strictMaxSigmaThreshold) {
      return PeakStatus.CANDIDATE;
    }
    if (Math.abs(result.getXPosition()) > strictOffset
        || Math.abs(result.getYPosition()) > strictOffset) {
      return PeakStatus.CANDIDATE;
    }
    if (variance > strictVariance) {
      return PeakStatus.CANDIDATE;
    }

    return PeakStatus.OK;
  }

  /**
   * Gets the variance.
   *
   * @param result the result
   * @return the variance
   */
  protected double getVariance(PeakResult result) {
    return calculator.getLSEVariance(result.getParameters(), result.getNoise());
  }

  @Override
  public double getNumericalValue() {
    return strictSnr;
  }

  @Override
  public String getNumericalValueName() {
    return ParameterType.SNR.toString() + " +" + rangeSnr;
  }

  @Override
  public String getDescription() {
    return "Filter results using a multiple thresholds: Signal, SNR, width, shift, precision. Any results within the "
        + "strict limits are included. Any results outside the weak limits are excluded. "
        + super.getDescription();
  }

  /** {@inheritDoc} */
  @Override
  public int getNumberOfParameters() {
    return 12 + super.getNumberOfParameters();
  }

  /** {@inheritDoc} */
  @Override
  protected double getParameterValueInternal(int index) {
    if (index < super.getNumberOfParameters()) {
      return super.getParameterValueInternal(index);
    }
    index -= super.getNumberOfParameters();
    switch (index) {
      case 0:
        return strictSignal;
      case 1:
        return rangeSignal;
      case 2:
        return strictSnr;
      case 3:
        return rangeSnr;
      case 4:
        return strictMinWidth;
      case 5:
        return rangeMinWidth;
      case 6:
        return strictMaxWidth;
      case 7:
        return rangeMaxWidth;
      case 8:
        return strictShift;
      case 9:
        return rangeShift;
      case 10:
        return strictPrecision;
      default:
        return rangePrecision;
    }
  }

  @Override
  public double getParameterIncrement(int index) {
    if (index < super.getNumberOfParameters()) {
      return super.getParameterValueInternal(index);
    }
    index -= super.getNumberOfParameters();
    switch (index) {
      case 0:
        return SignalFilter.DEFAULT_INCREMENT;
      case 1:
        return SignalFilter.DEFAULT_INCREMENT;
      case 2:
        return SNRFilter.DEFAULT_INCREMENT;
      case 3:
        return SNRFilter.DEFAULT_INCREMENT;
      case 4:
        return WidthFilter2.DEFAULT_MIN_INCREMENT;
      case 5:
        return WidthFilter2.DEFAULT_MIN_INCREMENT;
      case 6:
        return WidthFilter.DEFAULT_INCREMENT;
      case 7:
        return WidthFilter.DEFAULT_INCREMENT;
      case 8:
        return ShiftFilter.DEFAULT_INCREMENT;
      case 9:
        return ShiftFilter.DEFAULT_INCREMENT;
      case 10:
        return PrecisionFilter.DEFAULT_INCREMENT;
      default:
        return PrecisionFilter.DEFAULT_INCREMENT;
    }
  }

  /** {@inheritDoc} */
  @Override
  public ParameterType getParameterType(int index) {
    checkIndex(index);
    if (index < super.getNumberOfParameters()) {
      return super.getParameterType(index);
    }
    index -= super.getNumberOfParameters();
    switch (index) {
      case 0:
        return ParameterType.SIGNAL;
      case 1:
        return ParameterType.SIGNAL_RANGE;
      case 2:
        return ParameterType.SNR;
      case 3:
        return ParameterType.SNR_RANGE;
      case 4:
        return ParameterType.MIN_WIDTH;
      case 5:
        return ParameterType.MIN_WIDTH_RANGE;
      case 6:
        return ParameterType.MAX_WIDTH;
      case 7:
        return ParameterType.MAX_WIDTH_RANGE;
      case 8:
        return ParameterType.SHIFT;
      case 9:
        return ParameterType.SHIFT_RANGE;
      case 10:
        return getPrecisionParamaterType();
      default:
        return getPrecisionRangeParamaterType();
    }
  }

  /**
   * Gets the precision paramater type.
   *
   * @return the precision paramater type
   */
  protected ParameterType getPrecisionParamaterType() {
    return ParameterType.PRECISION;
  }

  /**
   * Gets the precision range paramater type.
   *
   * @return the precision range paramater type
   */
  protected ParameterType getPrecisionRangeParamaterType() {
    return ParameterType.PRECISION_RANGE;
  }

  /** The default range. */
  static double[] defaultRange = new double[] {0, 0, 0, 0, SignalFilter.DEFAULT_RANGE,
      SignalFilter.DEFAULT_RANGE, SNRFilter.DEFAULT_RANGE, SNRFilter.DEFAULT_RANGE,
      WidthFilter2.DEFAULT_MIN_RANGE, WidthFilter2.DEFAULT_MIN_RANGE, WidthFilter.DEFAULT_RANGE,
      WidthFilter.DEFAULT_RANGE, ShiftFilter.DEFAULT_RANGE, ShiftFilter.DEFAULT_RANGE,
      PrecisionFilter.DEFAULT_RANGE, PrecisionFilter.DEFAULT_RANGE};

  /** {@inheritDoc} */
  @Override
  public Filter adjustParameter(int index, double delta) {
    checkIndex(index);
    // No adjustment of the mode parameters
    if (index == 1 || index == 3) {
      return this;
    }
    final double[] parameters = new double[] {searchDistance, searchDistanceMode, timeThreshold,
        timeThresholdMode, strictSignal, rangeSignal, strictSnr, rangeSnr, strictMinWidth,
        rangeMinWidth, strictMaxWidth, rangeMaxWidth, strictShift, rangeShift, strictPrecision,
        rangePrecision};
    if (index == 0) {
      parameters[0] = updateParameter(parameters[0], delta, getDefaultSearchRange());
    } else if (index == 2) {
      parameters[2] = updateParameter(parameters[2], delta, getDefaultTimeRange());
    } else {
      parameters[index] = updateParameter(parameters[index], delta, defaultRange[index]);
    }
    return create(parameters);
  }

  /** {@inheritDoc} */
  @Override
  public Filter create(double... parameters) {
    return new MultiHysteresisFilter(parameters[0], (int) parameters[1], parameters[2],
        (int) parameters[3], parameters[4], parameters[5], (float) parameters[6],
        (float) parameters[7], parameters[8], parameters[9], parameters[10], parameters[11],
        parameters[12], parameters[13], parameters[14], parameters[15]);
  }

  /** {@inheritDoc} */
  @Override
  public void weakestParameters(double[] parameters) {
    super.weakestParameters(parameters);

    // Hysteresis filters require all the potential candidates, so disable hysteresis above the
    // candidate threshold
    setMin(parameters, 4, strictSignal - rangeSignal);
    parameters[5] = 0;
    setMin(parameters, 6, strictSnr - rangeSnr);
    parameters[7] = 0;
    setMin(parameters, 8, strictMinWidth - rangeMinWidth);
    parameters[9] = 0;
    setMax(parameters, 10, strictMaxWidth + rangeMaxWidth);
    parameters[11] = 0;
    setMax(parameters, 12, strictShift + rangeShift);
    parameters[13] = 0;
    setMax(parameters, 14, strictPrecision + rangePrecision);
    parameters[15] = 0;
  }

  /** {@inheritDoc} */
  @Override
  public Chromosome<FilterScore> newChromosome(double[] sequence) {
    // Override the default Hysteresis filter implementation for speed since this is the filter we
    // will most likely optimise using the genetic algorithm
    return new MultiHysteresisFilter(sequence[0], searchDistanceMode, sequence[1],
        timeThresholdMode, sequence[2], sequence[3], (float) sequence[4], (float) sequence[5],
        sequence[6], sequence[7], sequence[8], sequence[9], sequence[10], sequence[11],
        sequence[12], sequence[13]);
  }

  /** {@inheritDoc} */
  @Override
  public double[] upperLimit() {
    return new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY,
        Double.POSITIVE_INFINITY, WidthFilter.UPPER_LIMIT, WidthFilter.UPPER_LIMIT,
        WidthFilter.UPPER_LIMIT, WidthFilter.UPPER_LIMIT, ShiftFilter.UPPER_LIMIT,
        ShiftFilter.UPPER_LIMIT, PrecisionFilter.UPPER_LIMIT, PrecisionFilter.UPPER_LIMIT};
  }

  /** {@inheritDoc} */
  @Override
  public double[] mutationStepRange() {
    return new double[] {getDefaultSearchRange(), getDefaultTimeRange(), SignalFilter.DEFAULT_RANGE,
        SignalFilter.DEFAULT_RANGE, SNRFilter.DEFAULT_RANGE, SNRFilter.DEFAULT_RANGE,
        WidthFilter2.DEFAULT_MIN_RANGE, WidthFilter2.DEFAULT_MIN_RANGE, WidthFilter.DEFAULT_RANGE,
        WidthFilter.DEFAULT_RANGE, ShiftFilter.DEFAULT_RANGE, ShiftFilter.DEFAULT_RANGE,
        PrecisionFilter.DEFAULT_RANGE, PrecisionFilter.DEFAULT_RANGE};
  }
}
