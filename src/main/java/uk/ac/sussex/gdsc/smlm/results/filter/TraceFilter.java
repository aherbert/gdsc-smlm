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
import uk.ac.sussex.gdsc.smlm.results.Trace;
import uk.ac.sussex.gdsc.smlm.results.TraceManager;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

import java.util.HashSet;
import java.util.Set;

/**
 * Filter results that can be traced over time frames.
 */
public class TraceFilter extends Filter {
  /**
   * The default distance increment. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome}
   * interface.
   */
  private static final double DEFAULT_DISTANCE_INCREMENT = 0.05;
  /**
   * The default time increment. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome} interface.
   */
  private static final int DEFAULT_TIME_INCREMENT = 1;
  /**
   * The default distance range. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome} interface.
   */
  private static final double DEFAULT_DISTANCE_RANGE = 2;
  /**
   * The default time range. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome} interface.
   */
  private static final int DEFAULT_TIME_RANGE = 10;

  // The old names must be preserved for XStream serialisation
  // @CHECKSTYLE.OFF: MemberName

  /**
   * The distance.
   */
  @XStreamAsAttribute
  private final double d;

  /** The time. */
  @XStreamAsAttribute
  private final int t;

  @XStreamOmitField
  private Set<PeakResult> ok;

  /**
   * Instantiates a new trace filter.
   *
   * @param distance the distance
   * @param time the time
   */
  public TraceFilter(double distance, int time) {
    this.d = Math.max(0, distance);
    this.t = Math.max(0, time);
  }

  @Override
  protected String generateName() {
    return String.format("Trace distance=%.2f, time=%d", d, t);
  }

  @Override
  public void setup(MemoryPeakResults peakResults) {
    ok = new HashSet<>();

    // Trace molecules. Anything that is part of a trace is OK
    final TraceManager tm = new TraceManager(peakResults);
    tm.traceMolecules(d, t);
    final Trace[] traces = tm.getTraces();
    for (final Trace trace : traces) {
      if (trace.size() > 1) {
        for (int i = 0; i < trace.size(); i++) {
          ok.add(trace.get(i));
        }
      }
    }
  }

  /**
   * {@inheritDoc}
   *
   * @throws NullPointerException if not first initialised with a call to
   *         {@link #setup(MemoryPeakResults)}
   */
  @Override
  public boolean accept(PeakResult peak) {
    return ok.contains(peak);
  }

  @Override
  public double getNumericalValue() {
    return t;
  }

  @Override
  public String getNumericalValueName() {
    return ParameterType.TIME_THRESHOLD.toString();
  }

  @Override
  public String getDescription() {
    return "Filter results that can be traced over time frames.";
  }

  /** {@inheritDoc} */
  @Override
  public int getNumberOfParameters() {
    return 2;
  }

  /** {@inheritDoc} */
  @Override
  protected double getParameterValueInternal(int index) {
    switch (index) {
      case 0:
        return d;
      default:
        return t;
    }
  }

  /** {@inheritDoc} */
  @Override
  public double getParameterIncrement(int index) {
    checkIndex(index);
    switch (index) {
      case 0:
        return DEFAULT_DISTANCE_INCREMENT;
      default:
        return DEFAULT_TIME_INCREMENT;
    }
  }

  /** {@inheritDoc} */
  @Override
  public ParameterType getParameterType(int index) {
    checkIndex(index);
    switch (index) {
      case 0:
        return ParameterType.DISTANCE_THRESHOLD;
      default:
        return ParameterType.TIME_THRESHOLD;
    }
  }

  /** {@inheritDoc} */
  @Override
  public Filter adjustParameter(int index, double delta) {
    checkIndex(index);
    switch (index) {
      case 0:
        return new TraceFilter(updateParameter(d, delta, DEFAULT_DISTANCE_RANGE), t);
      default:
        return new TraceFilter(d, updateParameter(t, delta, DEFAULT_TIME_RANGE));
    }
  }

  /** {@inheritDoc} */
  @Override
  public Filter create(double... parameters) {
    return new TraceFilter(parameters[0], (int) parameters[1]);
  }

  /** {@inheritDoc} */
  @Override
  public void weakestParameters(double[] parameters) {
    setMax(parameters, 0, d);
    setMax(parameters, 1, t);
  }

  /** {@inheritDoc} */
  @Override
  public double[] mutationStepRange() {
    return new double[] {DEFAULT_DISTANCE_RANGE, DEFAULT_TIME_RANGE};
  }
}
