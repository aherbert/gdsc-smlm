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

import com.thoughtworks.xstream.annotations.XStreamOmitField;

import org.apache.commons.lang3.NotImplementedException;

import java.util.Arrays;

/**
 * Filter results using the combination of two filters.
 *
 * <p>Note that is the filter is not a DirectFilter then the result of filtering a
 * PreprocessedPeakResult is always true using that filter.
 */
public abstract class CombinedFilter extends DirectFilter {
  /** The first filter. */
  protected Filter filter1;

  /** The second filter. */
  protected Filter filter2;

  /**
   * The first filter, if it is a {@link DirectFilter}.
   */
  @XStreamOmitField
  protected DirectFilter dfilter1;

  /**
   * The second filter, if it is a {@link DirectFilter}.
   */
  @XStreamOmitField
  protected DirectFilter dfilter2;

  /** The result from the first filter. */
  @XStreamOmitField
  protected int result1;

  /** The result from the second filter. */
  @XStreamOmitField
  protected int result2;

  /**
   * Instantiates a new combined filter.
   *
   * @param filter1 the first filter
   * @param filter2 the second filter
   */
  public CombinedFilter(Filter filter1, Filter filter2) {
    this.filter1 = filter1;
    this.filter2 = filter2;
    initialiseState();
  }

  @Override
  protected void initialiseState() {
    if (filter1 instanceof DirectFilter) {
      dfilter1 = (DirectFilter) filter1;
    }
    if (filter2 instanceof DirectFilter) {
      dfilter2 = (DirectFilter) filter2;
    }
  }

  /** {@inheritDoc} */
  @Override
  public Filter clone() {
    // Add a reminder to implement clone
    throw new NotImplementedException("Derived classes must clone filter1 and filter2");
  }

  @Override
  protected String generateName() {
    final StringBuilder sb = new StringBuilder();
    addText(sb, filter1, filter1.getName());
    sb.append(" ").append(getOperator()).append(" ");
    addText(sb, filter2, filter2.getName());
    return sb.toString();
  }

  private static void addText(StringBuilder sb, Filter filter, String text) {
    if (filter instanceof CombinedFilter) {
      sb.append("(");
    }
    sb.append(text);
    if (filter instanceof CombinedFilter) {
      sb.append(")");
    }
  }

  @Override
  protected String generateType() {
    final StringBuilder sb = new StringBuilder();
    addText(sb, filter1, filter1.getType());
    sb.append(" ").append(getOperator()).append(" ");
    addText(sb, filter2, filter2.getType());
    return sb.toString();
  }

  @Override
  public String getDescription() {
    final StringBuilder sb = new StringBuilder();
    addText(sb, filter1, filter1.getDescription());
    sb.append(" ").append(getOperator()).append(" ");
    addText(sb, filter2, filter2.getDescription());
    return sb.toString();
  }

  /**
   * Get the string representation of the operator used to combine the two filters. This is used in
   * the filter name.
   *
   * @return The operator
   */
  protected abstract String getOperator();

  @Override
  public boolean requiresParameterDeviations() {
    return filter1.requiresParameterDeviations() || filter2.requiresParameterDeviations();
  }

  @Override
  public int getValidationFlags() {
    int flags = 0;
    if (dfilter1 != null) {
      flags |= dfilter1.getValidationFlags();
    }
    if (dfilter2 != null) {
      flags |= dfilter2.getValidationFlags();
    }
    return flags;
  }

  /**
   * Filter the result using filter1.
   *
   * @param peak The result
   * @return The filter result
   */
  public boolean accept1(PeakResult peak) {
    return filter1.accept(peak);
  }

  /**
   * Filter the result using filter1 if it is a DirectFilter, otherwise return true.
   *
   * @param peak The result
   * @return The filter result, or true
   */
  public boolean accept1(PreprocessedPeakResult peak) {
    result1 = (dfilter1 != null) ? dfilter1.validate(peak) : 0;
    return result1 == 0;
  }

  /**
   * Filter the result using filter2.
   *
   * @param peak The result
   * @return The filter result
   */
  public boolean accept2(PeakResult peak) {
    return filter2.accept(peak);
  }

  /**
   * Filter the result using filter2 if it is a DirectFilter, otherwise return true.
   *
   * @param peak The result
   * @return The filter result, or true
   */
  public boolean accept2(PreprocessedPeakResult peak) {
    result2 = (dfilter2 != null) ? dfilter2.validate(peak) : 0;
    return result2 == 0;
  }

  /** {@inheritDoc} */
  @Override
  public void setup(MemoryPeakResults peakResults) {
    filter1.setup(peakResults);
    filter2.setup(peakResults);
  }

  @Override
  public void setup() {
    if (dfilter1 != null) {
      dfilter1.setup();
    }
    if (dfilter2 != null) {
      dfilter2.setup();
    }
  }

  @Override
  public void setup(int flags) {
    if (dfilter1 != null) {
      dfilter1.setup(flags);
    }
    if (dfilter2 != null) {
      dfilter2.setup(flags);
    }
  }

  @Override
  public void setup(int flags, FilterSetupData... filterSetupData) {
    if (dfilter1 != null) {
      dfilter1.setup(flags, filterSetupData);
    }
    if (dfilter2 != null) {
      dfilter2.setup(flags, filterSetupData);
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>The method from the combined filter doesn't throw, leaving that to the underlying filters.
   * This does mean that a combined filter can be created from two already initialised filters and
   * the flags returned may not exactly recreate the state, since they are joined.
   */
  @Override
  public int getFilterSetupFlags() {
    int flags = 0;
    if (dfilter1 != null) {
      flags |= dfilter1.getFilterSetupFlags();
    }
    if (dfilter2 != null) {
      flags |= dfilter2.getFilterSetupFlags();
    }
    return flags;
  }

  /**
   * {@inheritDoc}
   *
   * <p>The method from the combined filter doesn't throw, leaving that to the underlying filters.
   * This does mean that a combined filter can be created from two already initialised filters and
   * the data returned may not exactly recreate the state, since they are joined.
   */
  @Override
  public FilterSetupData[] getFilterSetupData() {
    final FilterSetupData[] data1 = (dfilter1 == null) ? null : dfilter1.getFilterSetupData();
    if (data1 != null) {
      if (dfilter2 != null) {
        final FilterSetupData[] data2 = dfilter1.getFilterSetupData();
        if (data2 != null) {
          // Merge
          final int size = data1.length + data2.length;
          final FilterSetupData[] merge = Arrays.copyOf(data1, size);
          System.arraycopy(data2, 0, merge, data1.length, data2.length);
          return merge;
        }
      }
      return data1;
    }
    if (dfilter2 != null) {
      return dfilter2.getFilterSetupData();
    }
    return null;
  }

  /** {@inheritDoc} */
  @Override
  public void end() {
    filter1.end();
    filter2.end();
  }

  /** {@inheritDoc} */
  @Override
  public double getNumericalValue() {
    return filter1.getNumericalValue();
  }

  /** {@inheritDoc} */
  @Override
  public String getNumericalValueName() {
    return filter1.getNumericalValueName();
  }

  @Override
  public int getNumberOfParameters() {
    return filter1.getNumberOfParameters() + filter2.getNumberOfParameters();
  }

  @Override
  protected double getParameterValueInternal(int index) {
    if (index < filter1.getNumberOfParameters()) {
      return filter1.getParameterValueInternal(index);
    }
    return filter2.getParameterValueInternal(index - filter1.getNumberOfParameters());
  }

  @Override
  public double getParameterIncrement(int index) {
    checkIndex(index);
    if (index < filter1.getNumberOfParameters()) {
      return filter1.getParameterIncrement(index);
    }
    return filter2.getParameterIncrement(index - filter1.getNumberOfParameters());
  }

  @Override
  public double getDisabledParameterValue(int index) {
    checkIndex(index);
    if (index < filter1.getNumberOfParameters()) {
      return filter1.getDisabledParameterValue(index);
    }
    return filter2.getDisabledParameterValue(index - filter1.getNumberOfParameters());
  }

  @Override
  public ParameterType getParameterType(int index) {
    checkIndex(index);
    if (index < filter1.getNumberOfParameters()) {
      return filter1.getParameterType(index);
    }
    return filter2.getParameterType(index - filter1.getNumberOfParameters());
  }

  @Override
  public Filter adjustParameter(int index, double delta) {
    checkIndex(index);
    Filter f1 = filter1;
    Filter f2 = filter2;
    if (index < filter1.getNumberOfParameters()) {
      f1 = filter1.adjustParameter(index, delta);
    } else {
      f2 = filter2.adjustParameter(index - filter1.getNumberOfParameters(), delta);
    }
    return createFilter(f1, f2);
  }

  @Override
  public int lowerBoundOrientation(int index) {
    if (dfilter1 == null) {
      return 0;
    }

    if (index < dfilter1.getNumberOfParameters()) {
      return dfilter1.lowerBoundOrientation(index);
    }

    if (dfilter2 == null) {
      return 0;
    }

    return dfilter2.lowerBoundOrientation(index - dfilter1.getNumberOfParameters());
  }

  /**
   * Create a new combined filter from the two input filters.
   *
   * @param f1 the first filter
   * @param f2 the second filter
   * @return the filter
   */
  protected abstract Filter createFilter(Filter f1, Filter f2);

  /** {@inheritDoc} */
  @Override
  public Filter create(double... parameters) {
    final double[] p1 = Arrays.copyOf(parameters, filter1.getNumberOfParameters());
    final double[] p2 =
        Arrays.copyOfRange(parameters, filter1.getNumberOfParameters(), parameters.length);
    return createFilter(filter1.create(p1), filter2.create(p2));
  }

  /** {@inheritDoc} */
  @Override
  public void weakestParameters(double[] parameters) {
    final double[] p1 = Arrays.copyOf(parameters, filter1.getNumberOfParameters());
    final double[] p2 =
        Arrays.copyOfRange(parameters, filter1.getNumberOfParameters(), parameters.length);
    filter1.weakestParameters(p1);
    filter2.weakestParameters(p2);
    System.arraycopy(p1, 0, parameters, 0, p1.length);
    System.arraycopy(p2, 0, parameters, p1.length, p2.length);
  }

  /** {@inheritDoc} */
  @Override
  public boolean subsetWithFailCount() {
    return filter1.subsetWithFailCount() && filter2.subsetWithFailCount();
  }

  /** {@inheritDoc} */
  @Override
  public int length() {
    return filter1.length() + filter2.length();
  }

  /** {@inheritDoc} */
  @Override
  public double[] lowerLimit() {
    final double[] l1 = filter1.lowerLimit();
    final double[] l2 = filter2.lowerLimit();
    if (l1 == null && l2 == null) {
      return null;
    }
    return combine(getLowerLimit(filter1, l1), getLowerLimit(filter2, l2));
  }

  private static double[] getLowerLimit(Filter filter, double[] lower) {
    if (lower == null) {
      // Default to zero on the lower so no need to fill
      lower = new double[filter.getNumberOfParameters()];
    }
    return lower;
  }

  /** {@inheritDoc} */
  @Override
  public double[] upperLimit() {
    final double[] u1 = filter1.upperLimit();
    final double[] u2 = filter2.upperLimit();
    if (u1 == null && u2 == null) {
      return null;
    }
    return combine(getUpperLimit(filter1, u1), getUpperLimit(filter2, u2));
  }

  private static double[] getUpperLimit(Filter filter, double[] upper) {
    if (upper == null) {
      upper = new double[filter.getNumberOfParameters()];
      Arrays.fill(upper, Double.POSITIVE_INFINITY);
    }
    return upper;
  }

  /** {@inheritDoc} */
  @Override
  public double[] sequence() {
    return combine(filter1.sequence(), filter2.sequence());
  }

  /** {@inheritDoc} */
  @Override
  public double[] mutationStepRange() {
    return combine(filter1.mutationStepRange(), filter2.mutationStepRange());
  }

  private static double[] combine(double[] s1, double[] s2) {
    final double[] s = new double[s1.length + s2.length];
    System.arraycopy(s1, 0, s, 0, s1.length);
    System.arraycopy(s2, 0, s, s1.length, s2.length);
    return s;
  }

  /** {@inheritDoc} */
  @Override
  public int[] getChromosomeParameters() {
    final int[] s1 = filter1.getChromosomeParameters();
    final int[] s2 = filter2.getChromosomeParameters();
    final int[] s = new int[s1.length + s2.length];
    System.arraycopy(s1, 0, s, 0, s1.length);
    // Copy the next array but offset the index by the number of parameters in first filter
    // so that getParameterName(int) works OK
    final int n1 = filter1.getNumberOfParameters();
    for (int i = 0, j = s1.length; i < s2.length; i++, j++) {
      s[j] = s2[i] + n1;
    }
    return s;
  }
}
