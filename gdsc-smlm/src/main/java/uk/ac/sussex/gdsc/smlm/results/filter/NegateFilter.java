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

import com.thoughtworks.xstream.annotations.XStreamOmitField;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Filter results using the negation of a filter.
 *
 * <p>Note that if the filter is not a {@link DirectFilter} then the result of filtering a
 * PreprocessedPeakResult is always false using that filter.
 * 
 * <p>Warning: This filter requires a {@link DirectFilter} to implement the weakest filter
 * functionality. This is performed using the {@link DirectFilter#lowerBoundOrientation(int)} to
 * identify the strongest parameters. Otherwise it will raise an exception if comparison of
 * parameters are made to identify the weakest.
 */
public class NegateFilter extends DirectFilter {
  private static final String NEGATE_PREFIX = "Negate ";

  /** The first filter. */
  protected Filter filter;

  /**
   * The first filter, if it is a {@link DirectFilter}.
   */
  @XStreamOmitField
  protected DirectFilter dfilter;

  /**
   * Create an instance.
   *
   * @param filter1 the first filter
   */
  public NegateFilter(Filter filter1) {
    this.filter = filter1;
    init();
  }

  private void init() {
    if (filter instanceof DirectFilter) {
      dfilter = (DirectFilter) filter;
    }
  }

  @Override
  protected void initialiseState() {
    init();
  }

  @Override
  public NegateFilter clone() {
    final NegateFilter filter = (NegateFilter) super.clone();
    filter.filter = this.filter.clone();
    filter.initialiseState();
    return filter;
  }

  @Override
  protected String generateName() {
    return NEGATE_PREFIX + filter.getName();
  }

  @Override
  protected String generateType() {
    return NEGATE_PREFIX + filter.getType();
  }

  @Override
  public String getDescription() {
    return NEGATE_PREFIX + filter.getDescription();
  }

  @Override
  public boolean requiresParameterDeviations() {
    return filter.requiresParameterDeviations();
  }

  @Override
  public int getValidationFlags() {
    if (dfilter != null) {
      return dfilter.getValidationFlags();
    }
    return 0;
  }

  @Override
  public boolean accept(PeakResult peak) {
    return !filter.accept(peak);
  }

  @Override
  public int validate(final PreprocessedPeakResult peak) {
    if (dfilter != null) {
      if (dfilter.validate(peak) == 0) {
        return dfilter.getValidationFlags();
      }
      return 0;
    }
    // We have to negate the standard result of 0.
    // Use a non-zero result not corresponding to any FilterValidationFlag.
    return Integer.MIN_VALUE;
  }

  @Override
  public void setup(MemoryPeakResults peakResults) {
    filter.setup(peakResults);
  }

  @Override
  public void setup() {
    if (dfilter != null) {
      dfilter.setup();
    }
  }

  @Override
  public void setup(int flags) {
    if (dfilter != null) {
      dfilter.setup(flags);
    }
  }

  @Override
  public void setup(int flags, FilterSetupData... filterSetupData) {
    if (dfilter != null) {
      dfilter.setup(flags, filterSetupData);
    }
  }

  @Override
  public int getFilterSetupFlags() {
    if (dfilter != null) {
      return dfilter.getFilterSetupFlags();
    }
    return 0;
  }

  @Override
  public FilterSetupData[] getFilterSetupData() {
    if (dfilter != null) {
      return dfilter.getFilterSetupData();
    }
    return null;
  }

  @Override
  public void end() {
    filter.end();
  }

  @Override
  public double getNumericalValue() {
    return filter.getNumericalValue();
  }

  @Override
  public String getNumericalValueName() {
    return NEGATE_PREFIX + filter.getNumericalValueName();
  }

  @Override
  public int getNumberOfParameters() {
    return filter.getNumberOfParameters();
  }

  @Override
  protected double getParameterValueInternal(int index) {
    return filter.getParameterValueInternal(index);
  }

  @Override
  public double getParameterIncrement(int index) {
    return filter.getParameterIncrement(index);
  }

  @Override
  public double getDisabledParameterValue(int index) {
    return filter.getDisabledParameterValue(index);
  }

  @Override
  public ParameterType getParameterType(int index) {
    return filter.getParameterType(index);
  }

  @Override
  public Filter adjustParameter(int index, double delta) {
    return new NegateFilter(filter.adjustParameter(index, delta));
  }

  @Override
  public int lowerBoundOrientation(int index) {
    if (dfilter != null) {
      return -dfilter.lowerBoundOrientation(index);
    }
    return 0;
  }

  @Override
  public Filter create(double... parameters) {
    return new NegateFilter(filter.create(parameters));
  }

  /**
   * {@inheritDoc}.
   *
   * <p>Note: The negation of the weakest parameters requires the
   * {@link DirectFilter#lowerBoundOrientation(int)} to determine the strongest parameters of the
   * underlying filter.
   *
   * @throws UnsupportedOperationException if the underlying filter is not a {@link DirectFilter}
   */
  @Override
  public void weakestParameters(double[] parameters) {
    // To perform the opposite of weakest requires all filters to implement
    // the negation (i.e. the strongest parameters).
    // We can do this for direct filters using the opposite of the lower bound orientation.
    if (dfilter != null) {
      // Assume parameters is the correct length and access using getParameterValueInternal
      for (int i = 0; i < parameters.length; i++) {
        final int orientation = dfilter.lowerBoundOrientation(i);
        if (orientation < 0) {
          setMax(parameters, i, filter.getParameterValueInternal(i));
        } else if (orientation > 0) {
          setMin(parameters, i, filter.getParameterValueInternal(i));
        }
      }
    }
    throw new UnsupportedOperationException(
        "Unable to identify strongest parameters in filter " + filter.getName());
  }

  @Override
  public boolean subsetWithFailCount() {
    return filter.subsetWithFailCount();
  }

  @Override
  public int length() {
    return filter.length();
  }

  @Override
  public double[] lowerLimit() {
    return filter.lowerLimit();
  }

  @Override
  public double[] upperLimit() {
    return filter.upperLimit();
  }

  @Override
  public double[] sequence() {
    return filter.sequence();
  }

  @Override
  public double[] mutationStepRange() {
    return filter.mutationStepRange();
  }

  @Override
  public int[] getChromosomeParameters() {
    return filter.getChromosomeParameters();
  }

  @Override
  public int compareTo(Filter o) {
    return -filter.compareTo(o);
  }
}
