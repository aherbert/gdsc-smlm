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
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Filter results using a width range. Assumes width is identical on the X and Y axis.
 */
public class WidthFilter2 extends DirectFilter implements IMultiFilter {
  /**
   * The default increment for the min width. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome}
   * interface.
   */
  public static final double DEFAULT_MIN_INCREMENT = 0.02;
  /**
   * The default range for the min width. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome}
   * interface.
   */
  public static final double DEFAULT_MIN_RANGE = 1;

  /** The min width. */
  @XStreamAsAttribute
  protected final double minWidth;

  /** The max width. */
  @XStreamAsAttribute
  protected final double maxWidth;

  /** The lower sigma threshold. */
  @XStreamOmitField
  protected float lowerSigmaThreshold;

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
   * Instantiates a new width filter 2.
   *
   * @param minWidth the min width
   * @param maxWidth the max width
   */
  public WidthFilter2(double minWidth, double maxWidth) {
    // Only swap if max width is enabled
    if (maxWidth != 0 && maxWidth < minWidth) {
      final double f = maxWidth;
      maxWidth = minWidth;
      minWidth = f;
    }
    this.minWidth = Math.max(0, minWidth);
    this.maxWidth = Math.max(0, maxWidth);
  }

  @Override
  protected String generateName() {
    return "Width " + minWidth + "-" + maxWidth;
  }

  @Override
  public void setup(MemoryPeakResults peakResults) {
    calculator =
        Gaussian2DPeakResultHelper.create(peakResults.getPsf(), peakResults.getCalibration(), 0);

    // Set the width limit
    lowerSigmaThreshold = 0;
    upperSigmaThreshold = Float.POSITIVE_INFINITY;
    final double s = PsfHelper.getGaussian2DWx(peakResults.getPsf());
    lowerSigmaThreshold = (float) (s * minWidth);
    upperSigmaThreshold = Filter.getUpperLimit(s * maxWidth);
  }

  @Override
  public void setup() {
    setup(minWidth, maxWidth);
  }

  @Override
  public void setup(int flags) {
    if (areSet(flags, FilterValidationOption.NO_WIDTH)) {
      widthEnabled = false;
    } else {
      setup(minWidth, maxWidth);
    }
  }

  @Override
  public void setup(int flags, FilterSetupData... filterSetupData) {
    setup(flags);
  }

  /**
   * Setup the filter.
   *
   * @param minWidth the min width
   * @param maxWidth the max width
   */
  protected void setup(final double minWidth, double maxWidth) {
    widthEnabled = false;
    if (maxWidth > 1 && maxWidth != Double.POSITIVE_INFINITY) {
      upperSigmaThreshold = Filter.getUpperLimit(maxWidth);
      widthEnabled = upperSigmaThreshold != Float.POSITIVE_INFINITY;
    } else {
      upperSigmaThreshold = Float.POSITIVE_INFINITY;
    }
    if (minWidth < 1) {
      widthEnabled = true;
      lowerSigmaThreshold = (float) minWidth;
    } else {
      lowerSigmaThreshold = 0f;
    }
  }

  @Override
  public int getFilterSetupFlags() {
    return (widthEnabled) ? 0 : FilterValidationOption.NO_WIDTH;
  }

  @Override
  public boolean accept(PeakResult peak) {
    final float sd = calculator.getStandardDeviation(peak.getParameters());
    return sd <= upperSigmaThreshold && sd >= lowerSigmaThreshold;
  }

  @Override
  public int getValidationFlags() {
    return FilterValidationFlag.X_SD_FACTOR;
  }

  @Override
  public int validate(final PreprocessedPeakResult peak) {
    if (widthEnabled && (peak.getXSdFactor() > upperSigmaThreshold
        || peak.getXSdFactor() < lowerSigmaThreshold)) {
      return getValidationFlags();
    }
    return 0;
  }

  @Override
  public String getDescription() {
    return "Filter results using a width range. (Width is relative to initial peak width.)";
  }

  @Override
  public int getNumberOfParameters() {
    return 2;
  }

  @Override
  protected double getParameterValueInternal(int index) {
    return (index == 0) ? minWidth : maxWidth;
  }

  @Override
  public double getParameterIncrement(int index) {
    checkIndex(index);
    return (index == 0) ? WidthFilter2.DEFAULT_MIN_INCREMENT : WidthFilter.DEFAULT_INCREMENT;
  }

  @Override
  public ParameterType getParameterType(int index) {
    checkIndex(index);
    return (index == 0) ? ParameterType.MIN_WIDTH : ParameterType.MAX_WIDTH;
  }

  @Override
  public Filter adjustParameter(int index, double delta) {
    checkIndex(index);
    return (index == 0)
        ? new WidthFilter2(updateParameter(minWidth, delta, DEFAULT_MIN_RANGE), maxWidth)
        : new WidthFilter2(minWidth, updateParameter(maxWidth, delta, WidthFilter.DEFAULT_RANGE));
  }

  @Override
  public Filter create(double... parameters) {
    return new WidthFilter2(parameters[0], parameters[1]);
  }

  @Override
  public void weakestParameters(double[] parameters) {
    setMin(parameters, 0, minWidth);
    setMax(parameters, 1, maxWidth);
  }

  @Override
  public int lowerBoundOrientation(int index) {
    return (index == 1) ? 1 : -1;
  }

  @Override
  public double[] upperLimit() {
    return new double[] {WidthFilter.UPPER_LIMIT, WidthFilter.UPPER_LIMIT};
  }

  @Override
  public double[] mutationStepRange() {
    return new double[] {WidthFilter2.DEFAULT_MIN_RANGE, WidthFilter.DEFAULT_RANGE};
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
    return minWidth;
  }

  @Override
  public double getMaxWidth() {
    return maxWidth;
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
