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

import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

/**
 * Filter results using a signal-to-noise (SNR) threshold. Any results above the upper SNR limit are
 * included. Any results below the lower SNR limit are excluded. Any results between the limits are
 * included only if they can be traced through time, optionally via other candidates, to a valid
 * result.
 */
public class SNRHysteresisFilter extends HysteresisFilter {
  @XStreamAsAttribute
  private final float strictSnr;
  @XStreamAsAttribute
  private final float range;
  @XStreamOmitField
  private float weakSnr;

  /**
   * Instantiates a new SNR hysteresis filter.
   *
   * @param searchDistance the search distance
   * @param searchDistanceMode 0 = relative to the precision of the candidates; 1 = Absolute (in nm)
   * @param timeThreshold the time threshold
   * @param timeThresholdMode 0 = frames; 1 = seconds
   * @param strictSnr the strict snr
   * @param range the range
   */
  public SNRHysteresisFilter(double searchDistance, int searchDistanceMode, double timeThreshold,
      int timeThresholdMode, float strictSnr, float range) {
    super(searchDistance, searchDistanceMode, timeThreshold, timeThresholdMode);
    this.strictSnr = Math.max(0, strictSnr);
    this.range = Math.max(0, range);
  }

  @Override
  protected String generateName() {
    return String.format("SNR Hysteresis %.2f -%.2f (%s)", strictSnr, range, getTraceParameters());
  }

  @Override
  public void setup(MemoryPeakResults peakResults) {
    weakSnr = strictSnr - range;
    super.setup(peakResults);
  }

  @Override
  protected PeakStatus getStatus(PeakResult result) {
    final float snr = result.getSNR();
    if (snr >= strictSnr) {
      return PeakStatus.OK;
    } else if (snr >= weakSnr) {
      return PeakStatus.CANDIDATE;
    }
    return PeakStatus.REJECT;
  }

  @Override
  public double getNumericalValue() {
    return strictSnr;
  }

  @Override
  public String getNumericalValueName() {
    return ParameterType.SNR.toString() + " +" + range;
  }

  @Override
  public String getDescription() {
    return "Filter results using a signal-to-noise (SNR) threshold. Any results above the upper SNR "
        + "limit are included. Any results below the lower SNR limit are excluded. "
        + super.getDescription();
  }

  /** {@inheritDoc} */
  @Override
  public int getNumberOfParameters() {
    return 2 + super.getNumberOfParameters();
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
        return strictSnr;
      default:
        return range;
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
        return ParameterType.SNR;
      default:
        return ParameterType.SNR_RANGE;
    }
  }

  private static double[] defaultRange =
      new double[] {0, 0, 0, 0, SNRFilter.DEFAULT_RANGE, SNRFilter.DEFAULT_RANGE};

  /** {@inheritDoc} */
  @Override
  public Filter adjustParameter(int index, double delta) {
    checkIndex(index);
    // No adjustment of the mode parameters
    if (index == 1 || index == 3) {
      return this;
    }
    final double[] parameters = new double[] {searchDistance, searchDistanceMode, timeThreshold,
        timeThresholdMode, strictSnr, range};
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
    return new SNRHysteresisFilter(parameters[0], (int) parameters[1], parameters[2],
        (int) parameters[3], (float) parameters[4], (float) parameters[5]);
  }

  /** {@inheritDoc} */
  @Override
  public void weakestParameters(double[] parameters) {
    super.weakestParameters(parameters);

    // Hysteresis filters require all the potential candidates, so disable hysteresis above the
    // candidate threshold
    setMin(parameters, 4, strictSnr);
    parameters[5] = 0;
  }

  /** {@inheritDoc} */
  @Override
  public double[] mutationStepRange() {
    return new double[] {getDefaultSearchRange(), getDefaultTimeRange(), SNRFilter.DEFAULT_RANGE,
        SNRFilter.DEFAULT_RANGE};
  }
}
