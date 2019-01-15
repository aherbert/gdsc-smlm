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

import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

/**
 * Filter results using an amplitude-to-noise ratio (ANR) threshold.
 */
public class ANRFilter extends DirectFilter {
  /** The amplitude-to-noise ratio (ANR). */
  @XStreamAsAttribute
  private final float anr;

  @XStreamOmitField
  private Gaussian2DPeakResultCalculator calculator;

  /**
   * Instantiates a new amplitude-to-noise ratio (ANR) filter.
   *
   * @param anr the amplitude-to-noise ratio (ANR) threshold
   */
  public ANRFilter(float anr) {
    this.anr = Math.max(0, anr);
  }

  @Override
  public void setup(MemoryPeakResults peakResults) {
    calculator = Gaussian2DPeakResultHelper.create(peakResults.getPSF(),
        peakResults.getCalibration(), Gaussian2DPeakResultHelper.AMPLITUDE);
  }

  @Override
  public boolean accept(PeakResult peak) {
    return getANR(calculator, peak) >= this.anr;
  }

  /**
   * Gets the amplitude-to-noise ratio (ANR).
   *
   * @param calculator the calculator
   * @param peak the peak
   * @return the amplitude-to-noise ratio (ANR)
   */
  static float getANR(Gaussian2DPeakResultCalculator calculator, PeakResult peak) {
    return (peak.getNoise() > 0) ? calculator.getAmplitude(peak.getParameters()) / peak.getNoise()
        : Float.POSITIVE_INFINITY;
  }

  /**
   * Gets the amplitude-to-noise ratio (ANR).
   *
   * @param peak the peak
   * @return the amplitude-to-noise ratio (ANR)
   */
  static float getANR(PreprocessedPeakResult peak) {
    return (peak.getNoise() > 0) ? peak.getAmplitude() / peak.getNoise() : Float.POSITIVE_INFINITY;
  }

  @Override
  public int getValidationFlags() {
    return V_AMPLITUDE | V_NOISE;
  }

  @Override
  public int validate(final PreprocessedPeakResult peak) {
    if (getANR(peak) < this.anr) {
      return V_AMPLITUDE | V_NOISE;
    }
    return 0;
  }

  @Override
  public String getDescription() {
    return "Filter results using a lower ANR threshold.";
  }

  @Override
  public int getNumberOfParameters() {
    return 1;
  }

  @Override
  protected double getParameterValueInternal(int index) {
    return anr;
  }

  @Override
  public double getParameterIncrement(int index) {
    checkIndex(index);
    return SNRFilter.DEFAULT_INCREMENT;
  }

  @Override
  public ParameterType getParameterType(int index) {
    checkIndex(index);
    return ParameterType.ANR;
  }

  @Override
  public Filter adjustParameter(int index, double delta) {
    checkIndex(index);
    return new ANRFilter(updateParameter(anr, delta, SNRFilter.DEFAULT_RANGE));
  }

  @Override
  public Filter create(double... parameters) {
    return new ANRFilter((float) parameters[0]);
  }

  @Override
  public void weakestParameters(double[] parameters) {
    setMin(parameters, 0, anr);
  }

  @Override
  public double[] mutationStepRange() {
    return new double[] {SNRFilter.DEFAULT_RANGE};
  }
}
