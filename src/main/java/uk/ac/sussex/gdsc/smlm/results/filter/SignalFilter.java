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

import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

/**
 * Filter results using a signal threshold.
 */
public class SignalFilter extends DirectFilter implements IMultiFilter {
  /**
   * The default increment. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome} interface.
   */
  public static final double DEFAULT_INCREMENT = 5;
  /**
   * The default range. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome} interface.
   */
  public static final double DEFAULT_RANGE = 30;

  @XStreamAsAttribute
  private final double signal;
  @XStreamOmitField
  private float signalThreshold;

  /**
   * Instantiates a new signal filter.
   *
   * @param signal the signal
   */
  public SignalFilter(double signal) {
    this.signal = Math.max(0, signal);
  }

  @Override
  public void setup(MemoryPeakResults peakResults) {
    // Set the signal limit using the gain
    signalThreshold = (float) (signal * peakResults.getGain());
  }

  @Override
  public boolean accept(PeakResult peak) {
    return peak.getIntensity() >= signalThreshold;
  }

  @Override
  public int getValidationFlags() {
    return V_PHOTONS;
  }

  @Override
  public int validate(final PreprocessedPeakResult peak) {
    if (peak.getSignal() < signal) {
      return V_PHOTONS;
    }
    return 0;
  }

  @Override
  public String getDescription() {
    return "Filter results using a lower signal threshold. The threshold is applied in "
        + "photons (i.e. the signal is divided by the calibrated gain).";
  }

  @Override
  public int getNumberOfParameters() {
    return 1;
  }

  @Override
  protected double getParameterValueInternal(int index) {
    return signal;
  }

  @Override
  public double getParameterIncrement(int index) {
    checkIndex(index);
    return SignalFilter.DEFAULT_INCREMENT;
  }

  @Override
  public ParameterType getParameterType(int index) {
    checkIndex(index);
    return ParameterType.SIGNAL;
  }

  @Override
  public Filter adjustParameter(int index, double delta) {
    checkIndex(index);
    return new SignalFilter(updateParameter(signal, delta, DEFAULT_RANGE));
  }

  @Override
  public Filter create(double... parameters) {
    return new SignalFilter(parameters[0]);
  }

  @Override
  public void weakestParameters(double[] parameters) {
    setMin(parameters, 0, signal);
  }

  @Override
  public double[] mutationStepRange() {
    return new double[] {DEFAULT_RANGE};
  }

  @Override
  public double getSignal() {
    return signal;
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
    return 0;
  }

  @Override
  public PrecisionType getPrecisionType() {
    return PrecisionType.NONE;
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
