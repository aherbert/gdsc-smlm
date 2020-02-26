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

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Filter results using an upper width factor. Assumes width is identical on the X and Y axis.
 */
public class WidthFilter extends DirectFilter implements IMultiFilter {
  /**
   * The default increment. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome} interface.
   */
  public static final double DEFAULT_INCREMENT = 0.05;
  /**
   * The default range. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome} interface.
   */
  public static final double DEFAULT_RANGE = 1;
  /**
   * The default limit. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome} interface.
   */
  public static final double UPPER_LIMIT = 5;

  /** The width. */
  @XStreamAsAttribute
  protected final double width;

  /** The upper sigma threshold. */
  @XStreamOmitField
  protected float upperSigmaThreshold;

  /** The width enabled. */
  @XStreamOmitField
  protected boolean widthEnabled;

  /** The calculator. */
  @XStreamOmitField
  protected Gaussian2DPeakResultCalculator calculator;

  /**
   * Instantiates a new width filter.
   *
   * @param width the width
   */
  public WidthFilter(double width) {
    this.width = Math.max(0, width);
  }

  @Override
  public void setup(MemoryPeakResults peakResults) {
    calculator =
        Gaussian2DPeakResultHelper.create(peakResults.getPsf(), peakResults.getCalibration(), 0);

    // Set the width limit
    final double s = PsfHelper.getGaussian2DWx(peakResults.getPsf());
    upperSigmaThreshold = Filter.getUpperLimit(s * width);
  }

  @Override
  public void setup() {
    setup(width);
  }

  @Override
  public void setup(int flags) {
    if (areSet(flags, FilterValidationOption.NO_WIDTH)) {
      widthEnabled = false;
    } else {
      setup(width);
    }
  }

  @Override
  public void setup(int flags, FilterSetupData... filterSetupData) {
    setup(flags);
  }

  /**
   * Sets up the filter.
   *
   * @param width the new up
   */
  protected void setup(final double width) {
    upperSigmaThreshold = Filter.getUpperLimit(width);
    widthEnabled = (width != Float.POSITIVE_INFINITY);
  }

  @Override
  public int getFilterSetupFlags() {
    return (widthEnabled) ? 0 : FilterValidationOption.NO_WIDTH;
  }

  @Override
  public boolean accept(PeakResult peak) {
    return calculator.getStandardDeviation(peak.getParameters()) <= upperSigmaThreshold;
  }

  @Override
  public int getValidationFlags() {
    return FilterValidationFlag.X_SD_FACTOR;
  }

  @Override
  public int validate(final PreprocessedPeakResult peak) {
    if (widthEnabled && peak.getXSdFactor() > upperSigmaThreshold) {
      return getValidationFlags();
    }
    return 0;
  }

  @Override
  public String getDescription() {
    return "Filter results using an upper width factor. (Width is relative to initial peak width.)";
  }

  @Override
  public int getNumberOfParameters() {
    return 1;
  }

  @Override
  protected double getParameterValueInternal(int index) {
    return width;
  }

  @Override
  public double getParameterIncrement(int index) {
    checkIndex(index);
    return DEFAULT_INCREMENT;
  }

  @Override
  public ParameterType getParameterType(int index) {
    checkIndex(index);
    return ParameterType.MAX_WIDTH;
  }

  @Override
  public Filter adjustParameter(int index, double delta) {
    checkIndex(index);
    return new WidthFilter(updateParameter(width, delta, DEFAULT_RANGE));
  }

  @Override
  public Filter create(double... parameters) {
    return new WidthFilter(parameters[0]);
  }

  @Override
  public void weakestParameters(double[] parameters) {
    setMax(parameters, 0, width);
  }

  @Override
  public int lowerBoundOrientation(int index) {
    return 1;
  }

  @Override
  public double[] upperLimit() {
    return new double[] {UPPER_LIMIT};
  }

  @Override
  public double[] mutationStepRange() {
    return new double[] {DEFAULT_RANGE};
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
    return width;
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
