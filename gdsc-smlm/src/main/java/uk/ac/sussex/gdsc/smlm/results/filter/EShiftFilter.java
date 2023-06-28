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
import com.thoughtworks.xstream.annotations.XStreamOmitField;
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Filter results using a X/Y coordinate shift as a Euclidian distance. This filter requires that
 * the result X and Y coordinates are reported relative to their initial positions.
 */
public class EShiftFilter extends DirectFilter implements IMultiFilter {
  /**
   * The default increment. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome} interface.
   */
  public static final double DEFAULT_INCREMENT = 0.05;
  /**
   * The default range. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome} interface.
   */
  public static final double DEFAULT_RANGE = 10;
  /**
   * The default limit. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome} interface.
   */
  public static final double UPPER_LIMIT = 5;

  @XStreamAsAttribute
  private final double eshift;
  @XStreamOmitField
  private float eoffset;
  @XStreamOmitField
  private float eshift2;
  @XStreamOmitField
  private boolean shiftEnabled;

  /**
   * Instantiates a new e shift filter.
   *
   * @param eshift the eshift
   */
  public EShiftFilter(double eshift) {
    this.eshift = Math.max(0, eshift);
  }

  @Override
  public void setup(MemoryPeakResults peakResults) {
    // Set the shift limit
    final double[] s = PsfHelper.getGaussian2DWxWy(peakResults.getPsf());
    eoffset = getUpperLimit(s[0] * s[1] * eshift * eshift);
  }

  @Override
  public void setup() {
    setup(eshift);
  }

  @Override
  public void setup(int flags) {
    if (areSet(flags, FilterValidationOption.NO_SHIFT)) {
      shiftEnabled = false;
    } else {
      setup(eshift);
    }
  }

  @Override
  public void setup(int flags, FilterSetupData... filterSetupData) {
    if (areSet(flags, FilterValidationOption.NO_SHIFT)) {
      shiftEnabled = false;
      return;
    }

    for (int i = filterSetupData.length; i-- > 0;) {
      if (filterSetupData[i] instanceof ShiftFilterSetupData) {
        // Convert standard shift to Euclidian for a 2D?
        // Leaving it creates a circle with radius at the box edge.
        // Updating it creates a circle with radius at the box corner.
        final double shift = ((ShiftFilterSetupData) filterSetupData[i]).shift;
        // Leave for now
        // shift = Math.sqrt(shift * shift * 2)
        setup(shift);
        return;
      }
    }
    // Default
    setup(eshift);
  }

  private void setup(final double eshift) {
    eshift2 = getUpperSquaredLimit(eshift);
    shiftEnabled = (eshift2 != Float.POSITIVE_INFINITY);
  }

  @Override
  public FilterSetupData[] getFilterSetupData() {
    if (shiftEnabled && eshift2 != Float.POSITIVE_INFINITY) {
      if (eshift2 == getUpperSquaredLimit(eshift)) {
        // This is the default so ignore
        return null;
      }
      return getFilterSetupData(new ShiftFilterSetupData(Math.sqrt(eshift2)));
    }
    return null;
  }

  @Override
  public int getFilterSetupFlags() {
    return (shiftEnabled) ? 0 : FilterValidationOption.NO_SHIFT;
  }

  @Override
  public boolean accept(PeakResult peak) {
    final float dx = peak.getXPosition();
    final float dy = peak.getYPosition();
    return dx * dx + dy * dy <= eoffset;
  }

  @Override
  public int getValidationFlags() {
    return FilterValidationFlag.X_RELATIVE_SHIFT | FilterValidationFlag.Y_RELATIVE_SHIFT;
  }

  @Override
  public int validate(final PreprocessedPeakResult peak) {
    if (shiftEnabled && (peak.getXRelativeShift2() + peak.getYRelativeShift2()) > eshift2) {
      return getValidationFlags();
    }
    return 0;
  }

  @Override
  public String getDescription() {
    return "Filter results using a Euclidian shift factor. (Euclidian shift is relative to "
        + "initial peak width.)";
  }

  @Override
  public int getNumberOfParameters() {
    return 1;
  }

  @Override
  protected double getParameterValueInternal(int index) {
    return eshift;
  }

  @Override
  public double getParameterIncrement(int index) {
    checkIndex(index);
    return DEFAULT_INCREMENT;
  }

  @Override
  public ParameterType getParameterType(int index) {
    checkIndex(index);
    return ParameterType.ESHIFT;
  }

  @Override
  public Filter adjustParameter(int index, double delta) {
    checkIndex(index);
    return new EShiftFilter(updateParameter(eshift, delta, DEFAULT_RANGE));
  }

  @Override
  public Filter create(double... parameters) {
    return new EShiftFilter(parameters[0]);
  }

  @Override
  public void weakestParameters(double[] parameters) {
    setMax(parameters, 0, eshift);
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
    return 0;
  }

  @Override
  public double getShift() {
    return 0;
  }

  @Override
  public double getEShift() {
    return eshift;
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
